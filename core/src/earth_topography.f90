!> Earth's land fraction and surface geopotential from the ETOPO 2022 intermediate
!> file, smoothed with a Gaussian kernel whose width follows the truncation
!> (docs/dynamics/earth-topography.md).  Same outputs as the analytic topography.
module earth_topography
  use iso_fortran_env, only: real64
  use harmonics, only: harmonic_transform
  use planet_parameters, only: earth_gravity
  implicit none
  private

  integer, parameter, public :: earth_data_longitudes = 720, earth_data_latitudes = 360
  real(real64), parameter :: earth_data_cell_degrees = 0.5_real64
  real(real64), parameter :: pi = acos(-1.0_real64)
  real(real64), parameter :: degrees_to_radians = pi/180.0_real64

  type, public :: earth_topography_config
    !> Intermediate file written by scripts/prepare_earth_topography.py, relative to
    !> the working directory (core/) with core/ as the fallback prefix.
    character(len=128) :: data_path = 'data/earth_topography_0p5deg.bin'
    !> Kernel half width s = kernel_scale_factor * 180 deg / T.
    real(real64) :: kernel_scale_factor = 1.25_real64
    !> The kernel is cut off at window_factor * s.
    real(real64) :: window_factor = 3.0_real64
    !> Residual spectral filter exp(-kappa (n(n+1)/(T(T+1)))^2) on Phi_s; 0 disables it.
    real(real64) :: residual_filter_strength = 0.0_real64
    !> Grid points with a smaller land fraction count as open ocean in the diagnostics.
    real(real64) :: open_ocean_land_fraction = 0.01_real64
  end type earth_topography_config

  type, public :: earth_topography_diagnostics
    character(len=256) :: resolved_data_path = ''
    real(real64) :: kernel_half_width_degrees = 0.0_real64
    real(real64) :: source_land_fraction = 0.0_real64
    real(real64) :: source_maximum_height_metres = 0.0_real64
    real(real64) :: global_land_fraction = 0.0_real64
    real(real64) :: minimum_truncated_height_metres = 0.0_real64
    real(real64) :: minimum_longitude_degrees = 0.0_real64
    real(real64) :: minimum_latitude_degrees = 0.0_real64
    real(real64) :: maximum_truncated_height_metres = 0.0_real64
    real(real64) :: maximum_longitude_degrees = 0.0_real64
    real(real64) :: maximum_latitude_degrees = 0.0_real64
    !> RMS of truncated minus smoothed target height.
    real(real64) :: truncation_rms_metres = 0.0_real64
    !> RMS of truncated minus raw cell heights box-averaged to the grid.
    real(real64) :: total_rms_metres = 0.0_real64
    real(real64) :: open_ocean_rms_metres = 0.0_real64
    real(real64) :: open_ocean_minimum_metres = 0.0_real64
  end type earth_topography_diagnostics

  public :: generate_earth_topography, read_earth_topography_data, earth_data_cell_centre

contains

  !> Longitude (0..360) and latitude of the centre of intermediate cell (i, j), degrees.
  pure subroutine earth_data_cell_centre(i, j, longitude, latitude)
    integer, intent(in) :: i, j
    real(real64), intent(out) :: longitude, latitude
    longitude = (real(i, real64) - 0.5_real64)*earth_data_cell_degrees
    latitude = -90.0_real64 + (real(j, real64) - 0.5_real64)*earth_data_cell_degrees
  end subroutine earth_data_cell_centre

  !> Reads the land fraction and land height blocks of the intermediate file.
  subroutine read_earth_topography_data(data_path, land_fraction, land_height, resolved_path)
    character(*), intent(in) :: data_path
    real(real64), allocatable, intent(out) :: land_fraction(:, :), land_height(:, :)
    character(len=:), allocatable, intent(out) :: resolved_path
    integer :: unit, i, j, status
    integer(kind=8) :: file_size
    logical :: exists

    resolved_path = trim(data_path)
    inquire (file=resolved_path, exist=exists)
    if (.not. exists) then
      resolved_path = 'core/'//trim(data_path)
      inquire (file=resolved_path, exist=exists)
    end if
    if (.not. exists) error stop 'earth topography data file not found: '//trim(data_path)
    inquire (file=resolved_path, size=file_size)
    if (file_size /= 2_8*8_8*int(earth_data_longitudes, 8)*int(earth_data_latitudes, 8)) then
      error stop 'earth topography data file has an unexpected size: '//resolved_path
    end if
    allocate (land_fraction(earth_data_longitudes, earth_data_latitudes))
    allocate (land_height(earth_data_longitudes, earth_data_latitudes))
    open (newunit=unit, file=resolved_path, status='old', access='stream', form='unformatted', &
          action='read', convert='little_endian', iostat=status)
    if (status /= 0) error stop 'earth topography data file could not be opened: '//resolved_path
    do j = 1, earth_data_latitudes
      read (unit, iostat=status) (land_fraction(i, j), i=1, earth_data_longitudes)
      if (status /= 0) error stop 'earth topography land fraction block could not be read'
    end do
    do j = 1, earth_data_latitudes
      read (unit, iostat=status) (land_height(i, j), i=1, earth_data_longitudes)
      if (status /= 0) error stop 'earth topography land height block could not be read'
    end do
    close (unit)
    if (any(land_fraction < 0.0_real64) .or. any(land_fraction > 1.0_real64)) then
      error stop 'earth topography land fraction is outside [0,1]'
    end if
    if (any(land_height < 0.0_real64) .or. any(land_height > 9000.0_real64)) then
      error stop 'earth topography land height is outside [0, 9000] m'
    end if
  end subroutine read_earth_topography_data

  !> Smooth the intermediate fields onto the grid, truncate g z_s, and diagnose the result.
  subroutine generate_earth_topography(transform, config, land_fraction, target_height, &
                                       surface_geopotential, truncated_height, diagnostics)
    type(harmonic_transform), intent(in) :: transform
    type(earth_topography_config), intent(in) :: config
    real(real64), allocatable, intent(out) :: land_fraction(:, :), target_height(:, :)
    complex(real64), allocatable, intent(out) :: surface_geopotential(:, :)
    real(real64), allocatable, intent(out) :: truncated_height(:, :)
    type(earth_topography_diagnostics), intent(out) :: diagnostics
    real(real64), allocatable :: data_land(:, :), data_height(:, :), raw_height(:, :), geopotential_grid(:, :)
    real(real64), allocatable :: weights(:), latitude(:)
    character(len=:), allocatable :: resolved_path
    integer, allocatable :: nlon(:)
    integer :: truncation, i, j, n, m
    real(real64) :: half_width, area_weight, truncation_error, total_error, ocean_error, ocean_area
    real(real64) :: laplacian_ratio, longitude

    if (config%kernel_scale_factor <= 0.0_real64) error stop 'earth topography kernel scale must be positive'
    if (config%window_factor <= 0.0_real64) error stop 'earth topography window factor must be positive'
    if (config%residual_filter_strength < 0.0_real64) error stop 'earth topography filter strength must be >= 0'

    call read_earth_topography_data(config%data_path, data_land, data_height, resolved_path)
    nlon = transform%get_nlon()
    truncation = size(nlon)/2 - 1
    weights = transform%get_gaussian_weights()
    allocate (latitude(size(nlon)))
    do j = 1, size(nlon)
      latitude(j) = asin(transform%mu(j))/degrees_to_radians
    end do
    half_width = config%kernel_scale_factor*180.0_real64/real(truncation, real64)

    call transform%allocate_field(land_fraction)
    allocate (target_height, mold=land_fraction)
    allocate (geopotential_grid, mold=land_fraction)
    land_fraction = 0.0_real64
    target_height = 0.0_real64
    geopotential_grid = 0.0_real64
    call smooth_onto_grid(data_land, data_height, nlon, latitude, half_width, config%window_factor, &
                          land_fraction, target_height)
    do j = 1, size(nlon)
      geopotential_grid(1:nlon(j), j) = earth_gravity*target_height(1:nlon(j), j)
    end do
    call transform%grid_to_spectral(geopotential_grid, surface_geopotential)
    if (config%residual_filter_strength > 0.0_real64) then
      do m = 0, truncation
        do n = m, truncation
          laplacian_ratio = real(n*(n + 1), real64)/real(truncation*(truncation + 1), real64)
          surface_geopotential(n, m) = surface_geopotential(n, m)* &
            exp(-config%residual_filter_strength*laplacian_ratio**2)
        end do
      end do
    end if
    call transform%spectral_to_grid(surface_geopotential, truncated_height)
    truncated_height = truncated_height/earth_gravity
    call box_average_onto_grid(data_height, nlon, latitude, raw_height)

    diagnostics = earth_topography_diagnostics()
    diagnostics%resolved_data_path = resolved_path
    diagnostics%kernel_half_width_degrees = half_width
    diagnostics%source_land_fraction = data_area_mean(data_land)
    diagnostics%source_maximum_height_metres = maxval(data_height)
    diagnostics%minimum_truncated_height_metres = huge(0.0_real64)
    diagnostics%maximum_truncated_height_metres = -huge(0.0_real64)
    diagnostics%open_ocean_minimum_metres = huge(0.0_real64)
    truncation_error = 0.0_real64
    total_error = 0.0_real64
    ocean_error = 0.0_real64
    ocean_area = 0.0_real64
    do j = 1, size(nlon)
      area_weight = 0.5_real64*weights(j)/real(nlon(j), real64)
      do i = 1, nlon(j)
        longitude = 360.0_real64*real(i - 1, real64)/real(nlon(j), real64)
        diagnostics%global_land_fraction = diagnostics%global_land_fraction + area_weight*land_fraction(i, j)
        if (truncated_height(i, j) < diagnostics%minimum_truncated_height_metres) then
          diagnostics%minimum_truncated_height_metres = truncated_height(i, j)
          diagnostics%minimum_longitude_degrees = longitude
          diagnostics%minimum_latitude_degrees = latitude(j)
        end if
        if (truncated_height(i, j) > diagnostics%maximum_truncated_height_metres) then
          diagnostics%maximum_truncated_height_metres = truncated_height(i, j)
          diagnostics%maximum_longitude_degrees = longitude
          diagnostics%maximum_latitude_degrees = latitude(j)
        end if
        truncation_error = truncation_error + area_weight*(truncated_height(i, j) - target_height(i, j))**2
        total_error = total_error + area_weight*(truncated_height(i, j) - raw_height(i, j))**2
        if (land_fraction(i, j) < config%open_ocean_land_fraction) then
          ocean_area = ocean_area + area_weight
          ocean_error = ocean_error + area_weight*truncated_height(i, j)**2
          diagnostics%open_ocean_minimum_metres = min(diagnostics%open_ocean_minimum_metres, truncated_height(i, j))
        end if
      end do
    end do
    diagnostics%truncation_rms_metres = sqrt(max(truncation_error, 0.0_real64))
    diagnostics%total_rms_metres = sqrt(max(total_error, 0.0_real64))
    if (ocean_area > 0.0_real64) then
      diagnostics%open_ocean_rms_metres = sqrt(max(ocean_error/ocean_area, 0.0_real64))
    else
      diagnostics%open_ocean_minimum_metres = 0.0_real64
    end if
  end subroutine generate_earth_topography

  !> Area-weighted mean of an intermediate field over the sphere.
  real(real64) function data_area_mean(field) result(mean)
    real(real64), intent(in) :: field(:, :)
    real(real64) :: cell_longitude, cell_latitude, weight, weight_sum
    integer :: j

    mean = 0.0_real64
    weight_sum = 0.0_real64
    do j = 1, earth_data_latitudes
      call earth_data_cell_centre(1, j, cell_longitude, cell_latitude)
      weight = cos(cell_latitude*degrees_to_radians)
      mean = mean + weight*sum(field(:, j))
      weight_sum = weight_sum + weight*real(earth_data_longitudes, real64)
    end do
    mean = mean/weight_sum
  end function data_area_mean

  !> Gaussian-kernel average of the intermediate fields at every grid point.  The
  !> kernel exp(-(theta/s)^2) in great-circle angle theta is cut off at
  !> window_factor * s and normalized over the included cells, so a land fraction
  !> in [0,1] and a non-negative height are preserved.
  subroutine smooth_onto_grid(data_land, data_height, nlon, latitude, half_width, window_factor, &
                              land_fraction, target_height)
    real(real64), intent(in) :: data_land(:, :), data_height(:, :), latitude(:), half_width, window_factor
    integer, intent(in) :: nlon(:)
    real(real64), intent(inout) :: land_fraction(:, :), target_height(:, :)
    real(real64) :: cell_longitude(earth_data_longitudes), cell_latitude(earth_data_latitudes)
    real(real64) :: cell_cos(earth_data_latitudes), cell_sin(earth_data_latitudes)
    real(real64) :: cell_cos_lon(earth_data_longitudes), cell_sin_lon(earth_data_longitudes)
    real(real64) :: s, cutoff, sin_half_cutoff, band_cos_minimum, sin_half_dlon_limit, point_cos, point_sin
    real(real64) :: point_lon, point_cos_lon, point_sin_lon, cos_dlon, sin_half_dlon, cos_theta, theta
    real(real64) :: weight, weight_sum, land_sum, height_sum, cos_lat_j, band_edge
    integer :: i, j, ic, jc, jc_low, jc_high
    real(real64) :: dummy

    do ic = 1, earth_data_longitudes
      call earth_data_cell_centre(ic, 1, cell_longitude(ic), dummy)
      cell_cos_lon(ic) = cos(cell_longitude(ic)*degrees_to_radians)
      cell_sin_lon(ic) = sin(cell_longitude(ic)*degrees_to_radians)
    end do
    do jc = 1, earth_data_latitudes
      call earth_data_cell_centre(1, jc, dummy, cell_latitude(jc))
      cell_cos(jc) = cos(cell_latitude(jc)*degrees_to_radians)
      cell_sin(jc) = sin(cell_latitude(jc)*degrees_to_radians)
    end do
    s = half_width*degrees_to_radians
    cutoff = window_factor*s
    sin_half_cutoff = sin(0.5_real64*min(cutoff, pi))

    !$omp parallel do default(shared) schedule(dynamic, 1) &
    !$omp private(j, i, ic, jc, jc_low, jc_high, cos_lat_j, band_edge, band_cos_minimum, sin_half_dlon_limit) &
    !$omp private(point_cos, point_sin, point_lon, point_cos_lon, point_sin_lon, cos_dlon, sin_half_dlon) &
    !$omp private(cos_theta, theta, weight, weight_sum, land_sum, height_sum)
    do j = 1, size(nlon)
      cos_lat_j = cos(latitude(j)*degrees_to_radians)
      point_sin = sin(latitude(j)*degrees_to_radians)
      point_cos = cos_lat_j
      ! Latitude band of cells that can lie within the cutoff.
      jc_low = max(1, floor((latitude(j) - window_factor*half_width + 90.0_real64)/earth_data_cell_degrees) + 1)
      jc_high = min(earth_data_latitudes, &
                    floor((latitude(j) + window_factor*half_width + 90.0_real64)/earth_data_cell_degrees) + 1)
      ! Longitude prefilter: sin(theta/2) >= cos(phi_max) sin(dlon/2), so cells with
      ! cos(phi_max) sin(dlon/2) > sin(cutoff/2) are outside the cutoff.
      band_edge = min(abs(latitude(j)) + window_factor*half_width, 90.0_real64)
      band_cos_minimum = cos(band_edge*degrees_to_radians)
      if (band_cos_minimum > 0.0_real64) then
        sin_half_dlon_limit = sin_half_cutoff/band_cos_minimum
      else
        sin_half_dlon_limit = 2.0_real64
      end if
      do i = 1, nlon(j)
        point_lon = 2.0_real64*pi*real(i - 1, real64)/real(nlon(j), real64)
        point_cos_lon = cos(point_lon)
        point_sin_lon = sin(point_lon)
        weight_sum = 0.0_real64
        land_sum = 0.0_real64
        height_sum = 0.0_real64
        do jc = jc_low, jc_high
          do ic = 1, earth_data_longitudes
            cos_dlon = cell_cos_lon(ic)*point_cos_lon + cell_sin_lon(ic)*point_sin_lon
            if (sin_half_dlon_limit < 1.0_real64) then
              sin_half_dlon = sqrt(max(0.0_real64, 0.5_real64*(1.0_real64 - cos_dlon)))
              if (sin_half_dlon > sin_half_dlon_limit) cycle
            end if
            cos_theta = point_sin*cell_sin(jc) + point_cos*cell_cos(jc)*cos_dlon
            theta = acos(max(-1.0_real64, min(1.0_real64, cos_theta)))
            if (theta > cutoff) cycle
            weight = cell_cos(jc)*exp(-(theta/s)**2)
            weight_sum = weight_sum + weight
            land_sum = land_sum + weight*data_land(ic, jc)
            height_sum = height_sum + weight*data_height(ic, jc)
          end do
        end do
        if (weight_sum <= 0.0_real64) error stop 'earth topography kernel window contains no cells'
        land_fraction(i, j) = max(0.0_real64, min(1.0_real64, land_sum/weight_sum))
        target_height(i, j) = max(0.0_real64, height_sum/weight_sum)
      end do
    end do
    !$omp end parallel do
  end subroutine smooth_onto_grid

  !> Raw cell heights averaged over each grid point's latitude band (between ring
  !> midpoints) and nearest-longitude sector, for the total-error diagnostic.
  subroutine box_average_onto_grid(data_height, nlon, latitude, raw_height)
    real(real64), intent(in) :: data_height(:, :), latitude(:)
    integer, intent(in) :: nlon(:)
    real(real64), allocatable, intent(out) :: raw_height(:, :)
    real(real64), allocatable :: weight_sum(:, :), edges(:)
    real(real64) :: cell_longitude, cell_latitude, area
    integer :: i, j, ic, jc, rings

    rings = size(nlon)
    allocate (raw_height(maxval(nlon), rings), weight_sum(maxval(nlon), rings), edges(0:rings))
    raw_height = 0.0_real64
    weight_sum = 0.0_real64
    edges(0) = -90.0_real64
    edges(rings) = 90.0_real64
    do j = 1, rings - 1
      edges(j) = 0.5_real64*(latitude(j) + latitude(j + 1))
    end do
    j = 1
    do jc = 1, earth_data_latitudes
      call earth_data_cell_centre(1, jc, cell_longitude, cell_latitude)
      do while (j < rings .and. cell_latitude > edges(j))
        j = j + 1
      end do
      area = cos(cell_latitude*degrees_to_radians)
      do ic = 1, earth_data_longitudes
        call earth_data_cell_centre(ic, jc, cell_longitude, cell_latitude)
        i = modulo(nint(cell_longitude*real(nlon(j), real64)/360.0_real64), nlon(j)) + 1
        raw_height(i, j) = raw_height(i, j) + area*data_height(ic, jc)
        weight_sum(i, j) = weight_sum(i, j) + area
      end do
    end do
    do j = 1, rings
      do i = 1, nlon(j)
        if (weight_sum(i, j) <= 0.0_real64) error stop 'earth topography box average has an empty grid cell'
        raw_height(i, j) = raw_height(i, j)/weight_sum(i, j)
      end do
    end do
  end subroutine box_average_onto_grid

end module earth_topography
