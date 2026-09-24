!> Starts a land--sea case from a yearly snapshot of the same planet at a lower
!> truncation (docs/cases/land-sea.md#t63-の初期値).  The spectral prognostic
!> fields are zero-padded, log surface pressure is moved hydrostatically onto the
!> higher-resolution terrain, and the grid surface fields are interpolated between
!> the two octahedral grids with land and ocean masks.
module land_sea_restart
  use iso_fortran_env, only: real64
  use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
  use harmonics, only: harmonic_transform
  use dry_atmosphere, only: dry_atmosphere_solver
  use dry_physics_config, only: dry_model_physics_config
  use dry_vertical_coordinate, only: hybrid_sigma_coordinate, dry_air_gas_constant
  use moist_thermodynamics, only: virtual_temperature
  use planet_parameters, only: planet_config
  use sea_ice, only: equilibrate_sea_ice
  implicit none
  private

  public :: set_land_sea_state_from_snapshot
  public :: read_ring_field, read_spectral_snapshot_field, read_snapshot_time
  public :: pad_spectral_field, regrid_ring_field, terrain_log_surface_pressure_change

  !> Tolerance for the land fraction stored with the snapshot against the one
  !> regenerated at the snapshot truncation.
  real(real64), parameter :: land_fraction_tolerance = 1.0e-12_real64

  type, public :: land_sea_restart_summary
    integer :: source_truncation = -1
    real(real64) :: source_time_seconds = 0.0_real64
    !> Largest |Delta ln p_s| of the terrain adjustment on the target grid.
    real(real64) :: maximum_log_surface_pressure_change = 0.0_real64
    !> Target points with land (ocean) whose four interpolation neighbours hold
    !> no land (ocean); they take the nearest source point that does.
    integer :: land_fallback_points = 0
    integer :: ocean_fallback_points = 0
    !> Target ocean points that start with sea ice.
    integer :: ice_points = 0
  end type land_sea_restart_summary

contains

  !> Installs the state of year `year` from `source_directory` on `solver`, whose
  !> truncation must not be lower than that of `source_transform`.  The source
  !> terrain (`source_surface_geopotential`, `source_land_fraction`) is regenerated
  !> by the caller at the source truncation; the stored land fraction must match it.
  !> The snapshot time must be a whole number of `calendar_period`s, so that the
  !> target run can start at time zero with the same season and time of day.
  subroutine set_land_sea_state_from_snapshot(solver, transform, physics, planet, surface_geopotential, &
                                              land_fraction, source_transform, source_surface_geopotential, &
                                              source_land_fraction, source_directory, year, calendar_period, &
                                              summary)
    type(dry_atmosphere_solver), intent(inout) :: solver
    type(harmonic_transform), intent(inout) :: transform
    type(dry_model_physics_config), intent(in) :: physics
    type(planet_config), intent(in) :: planet
    complex(real64), intent(in) :: surface_geopotential(0:, 0:)
    real(real64), intent(in) :: land_fraction(:, :)
    type(harmonic_transform), intent(in) :: source_transform
    complex(real64), intent(in) :: source_surface_geopotential(0:, 0:)
    real(real64), intent(in) :: source_land_fraction(:, :)
    character(*), intent(in) :: source_directory
    integer, intent(in) :: year
    real(real64), intent(in) :: calendar_period
    type(land_sea_restart_summary), intent(out) :: summary
    type(hybrid_sigma_coordinate) :: coordinate
    complex(real64), allocatable :: zeta(:, :, :), delta(:, :, :), temperature(:, :, :), humidity(:, :, :)
    complex(real64), allocatable :: log_ps(:, :), source_spectral(:, :), log_ps_change(:, :)
    real(real64), allocatable :: stored_land(:, :), source(:, :), target_template(:, :)
    real(real64), allocatable :: surface_temperature(:, :), deep_temperature(:, :), surface_water(:, :)
    real(real64), allocatable :: land_temperature(:, :), ocean_temperature(:, :), snow_water(:, :)
    real(real64), allocatable :: sea_ice_fraction(:, :), sea_ice_volume(:, :)
    logical, allocatable :: source_land(:, :), source_ocean(:, :)
    integer, allocatable :: source_nlon(:), target_nlon(:)
    character(len=4) :: year_text
    character(len=2) :: level_text
    character(len=:), allocatable :: prefix
    integer :: source_truncation, target_truncation, levels, k, i, j, fallback
    real(real64) :: ocean_heat_capacity, periods

    source_nlon = source_transform%get_nlon()
    target_nlon = transform%get_nlon()
    source_truncation = size(source_nlon)/2 - 1
    target_truncation = solver%get_truncation()
    if (target_truncation < source_truncation) error stop 'restart target truncation is lower than the source'
    coordinate = solver%get_coordinate()
    levels = coordinate%number_of_levels
    summary%source_truncation = source_truncation
    write (year_text, '(i4.4)') year

    summary%source_time_seconds = read_snapshot_time(source_directory//'/yearly_time_y'//year_text//'.json')
    periods = summary%source_time_seconds/calendar_period
    if (abs(periods - real(nint(periods), real64)) > 1.0e-9_real64) then
      error stop 'restart snapshot time is not a whole number of calendar years'
    end if

    ! The stored land fraction must be the one this code regenerates for the source,
    ! otherwise the snapshot belongs to a different terrain.
    call source_transform%allocate_field(stored_land)
    call read_ring_field(source_directory//'/land_fraction.bin', source_nlon, stored_land)
    do j = 1, size(source_nlon)
      if (maxval(abs(stored_land(1:source_nlon(j), j) - source_land_fraction(1:source_nlon(j), j))) > &
          land_fraction_tolerance) then
        error stop 'restart source land fraction differs from the regenerated source terrain; rerun the source case'
      end if
    end do

    ! Spectral state: the source coefficients are the target coefficients with n <= source T.
    allocate (zeta(0:target_truncation + 1, 0:target_truncation, levels))
    allocate (delta, temperature, humidity, mold=zeta)
    do k = 1, levels
      write (level_text, '(i2.2)') k
      prefix = source_directory//'/yearly_'
      call read_spectral_snapshot_field(prefix//'zeta_spectral_y'//year_text//'_l'//level_text//'.bin', &
                                        source_truncation, source_spectral)
      call pad_spectral_field(source_spectral, target_truncation, log_ps)
      zeta(:, :, k) = log_ps
      call read_spectral_snapshot_field(prefix//'delta_spectral_y'//year_text//'_l'//level_text//'.bin', &
                                        source_truncation, source_spectral)
      call pad_spectral_field(source_spectral, target_truncation, log_ps)
      delta(:, :, k) = log_ps
      call read_spectral_snapshot_field(prefix//'temperature_spectral_y'//year_text//'_l'//level_text//'.bin', &
                                        source_truncation, source_spectral)
      call pad_spectral_field(source_spectral, target_truncation, log_ps)
      temperature(:, :, k) = log_ps
      call read_spectral_snapshot_field(prefix//'specific_humidity_spectral_y'//year_text//'_l'//level_text// &
                                        '.bin', source_truncation, source_spectral)
      call pad_spectral_field(source_spectral, target_truncation, log_ps)
      humidity(:, :, k) = log_ps
    end do
    call read_spectral_snapshot_field(source_directory//'/yearly_log_surface_pressure_spectral_y'//year_text// &
                                      '.bin', source_truncation, source_spectral)
    call pad_spectral_field(source_spectral, target_truncation, log_ps)
    call terrain_log_surface_pressure_change(transform, target_nlon, source_surface_geopotential, &
      surface_geopotential, temperature(:, :, levels), humidity(:, :, levels), log_ps_change, &
      summary%maximum_log_surface_pressure_change)
    log_ps = log_ps + log_ps_change

    ! Grid surface fields.  Land quantities are interpolated from source points with
    ! land, ocean quantities from source points with ocean.
    allocate (source_land(size(source_land_fraction, 1), size(source_land_fraction, 2)))
    allocate (source_ocean, mold=source_land)
    source_land = source_land_fraction > 0.0_real64
    source_ocean = source_land_fraction < 1.0_real64
    call transform%allocate_field(target_template)
    call source_transform%allocate_field(source)

    call regrid_snapshot('surface_temperature', surface_temperature)
    call regrid_snapshot('deep_temperature', deep_temperature, source_land, land_fraction > 0.0_real64, fallback)
    summary%land_fallback_points = fallback
    call regrid_snapshot('land_temperature', land_temperature, source_land, land_fraction > 0.0_real64)
    call regrid_snapshot('surface_water', surface_water, source_land, land_fraction > 0.0_real64)
    call regrid_snapshot('snow_water', snow_water, source_land, land_fraction > 0.0_real64)
    call regrid_snapshot('ocean_temperature', ocean_temperature, source_ocean, land_fraction < 1.0_real64, fallback)
    summary%ocean_fallback_points = fallback
    call regrid_snapshot('sea_ice_fraction', sea_ice_fraction, source_ocean, land_fraction < 1.0_real64)
    call regrid_snapshot('sea_ice_volume', sea_ice_volume, source_ocean, land_fraction < 1.0_real64)

    where (land_fraction > 0.0_real64)
      surface_water = min(physics%bucket%capacity, max(0.0_real64, surface_water))
      snow_water = max(0.0_real64, snow_water)
    elsewhere
      surface_water = 0.0_real64
      snow_water = 0.0_real64
    end where
    ! Interpolation is linear in T_o and V, hence in the ocean--ice energy
    ! C_o (T_o - T_f) - rho_i L_f V; the phase equilibrium is restored from it.
    ocean_heat_capacity = physics%radiation%seawater_density*physics%radiation%seawater_specific_heat* &
                          physics%radiation%slab_ocean_depth
    do j = 1, size(target_nlon)
      do i = 1, target_nlon(j)
        if (land_fraction(i, j) >= 1.0_real64) then
          sea_ice_fraction(i, j) = 0.0_real64
          sea_ice_volume(i, j) = 0.0_real64
          cycle
        end if
        sea_ice_fraction(i, j) = min(1.0_real64, max(0.0_real64, sea_ice_fraction(i, j)))
        sea_ice_volume(i, j) = max(0.0_real64, sea_ice_volume(i, j))
        if (physics%sea_ice%enabled) then
          call equilibrate_sea_ice(physics%sea_ice, ocean_heat_capacity, ocean_temperature(i, j), &
                                   sea_ice_fraction(i, j), sea_ice_volume(i, j))
        else
          sea_ice_fraction(i, j) = 0.0_real64
          sea_ice_volume(i, j) = 0.0_real64
        end if
        if (sea_ice_fraction(i, j) > 0.0_real64) summary%ice_points = summary%ice_points + 1
      end do
    end do

    call solver%set_planet(planet)
    call solver%set_initial_state(zeta, delta, temperature, log_ps, surface_geopotential, &
                                  specific_humidity=humidity, land_fraction=land_fraction)
    ! Sea ice is validated against the active physics, so the physics precede the surface.
    call solver%set_physics(physics)
    call solver%set_surface_state(surface_temperature, deep_temperature, surface_water, &
                                  land_temperature=land_temperature, ocean_temperature=ocean_temperature, &
                                  sea_ice_fraction=sea_ice_fraction, sea_ice_volume=sea_ice_volume, &
                                  snow_water=snow_water)

  contains

    subroutine regrid_snapshot(name, field, valid, needed, fallback_points)
      character(*), intent(in) :: name
      real(real64), allocatable, intent(out) :: field(:, :)
      logical, intent(in), optional :: valid(:, :), needed(:, :)
      integer, intent(out), optional :: fallback_points

      call read_ring_field(source_directory//'/yearly_'//name//'_y'//year_text//'.bin', source_nlon, source)
      field = target_template
      call regrid_ring_field(source_transform%mu, source_nlon, source, transform%mu, target_nlon, field, &
                             valid, needed, fallback_points)
    end subroutine regrid_snapshot

  end subroutine set_land_sea_state_from_snapshot

  !> Change of ln p_s that keeps the pressure at the source surface when the
  !> surface moves from Phi_s^src to Phi_s: Delta ln p_s = -(Phi_s - Phi_s^src)/(R_d T_v,N),
  !> with the lowest-level virtual temperature standing for the thin layer between
  !> the two surfaces.  The change is evaluated on the target grid and analysed on
  !> its own, so the zero-padded source ln p_s is not passed through a transform.
  subroutine terrain_log_surface_pressure_change(transform, nlon, source_surface_geopotential, &
                                                 surface_geopotential, lowest_temperature, lowest_humidity, &
                                                 change, maximum_change)
    type(harmonic_transform), intent(inout) :: transform
    integer, intent(in) :: nlon(:)
    complex(real64), intent(in) :: source_surface_geopotential(0:, 0:), surface_geopotential(0:, 0:)
    complex(real64), intent(in) :: lowest_temperature(0:, 0:), lowest_humidity(0:, 0:)
    complex(real64), allocatable, intent(out) :: change(:, :)
    real(real64), intent(out) :: maximum_change
    complex(real64), allocatable :: padded(:, :)
    real(real64), allocatable :: source_phi(:, :), target_phi(:, :), temperature(:, :), humidity(:, :)
    real(real64), allocatable :: grid_change(:, :)
    integer :: truncation, i, j

    truncation = ubound(surface_geopotential, 2)
    call pad_spectral_field(source_surface_geopotential, truncation, padded)
    call transform%spectral_to_grid(padded, source_phi)
    call transform%spectral_to_grid(surface_geopotential, target_phi)
    call transform%spectral_to_grid(lowest_temperature, temperature)
    call transform%spectral_to_grid(lowest_humidity, humidity)
    allocate (grid_change, mold=target_phi)
    grid_change = 0.0_real64
    maximum_change = 0.0_real64
    do j = 1, size(nlon)
      do i = 1, nlon(j)
        grid_change(i, j) = -(target_phi(i, j) - source_phi(i, j))/ &
                            (dry_air_gas_constant*virtual_temperature(temperature(i, j), humidity(i, j)))
        maximum_change = max(maximum_change, abs(grid_change(i, j)))
      end do
    end do
    call transform%grid_to_spectral(grid_change, change)
  end subroutine terrain_log_surface_pressure_change

  !> Copies the coefficients n, m <= source T into a zeroed (0:T+1, 0:T) array.
  subroutine pad_spectral_field(source, truncation, padded)
    complex(real64), intent(in) :: source(0:, 0:)
    integer, intent(in) :: truncation
    complex(real64), allocatable, intent(out) :: padded(:, :)
    integer :: source_truncation

    source_truncation = ubound(source, 2)
    if (source_truncation > truncation) error stop 'cannot pad a spectral field to a lower truncation'
    allocate (padded(0:truncation + 1, 0:truncation))
    padded = cmplx(0.0_real64, 0.0_real64, kind=real64)
    padded(0:source_truncation, 0:source_truncation) = source(0:source_truncation, 0:source_truncation)
  end subroutine pad_spectral_field

  !> Bilinear interpolation between octahedral Gaussian grids: linear in longitude
  !> on the two source rings that bracket the target latitude and linear in latitude
  !> between them.  Beyond the outermost source ring only that ring is used.
  !> With `valid`, invalid source points are dropped and the remaining weights
  !> renormalized; a target point in `needed` (all points by default) left without a
  !> valid neighbour takes the nearest valid source point on the sphere, and one
  !> outside `needed` keeps the unmasked interpolation.  `destination` must be allocated
  !> with the padded target shape; the padding is left untouched.
  subroutine regrid_ring_field(source_mu, source_nlon, source, target_mu, target_nlon, destination, valid, needed, &
                               fallback_points)
    real(real64), intent(in) :: source_mu(:), source(:, :), target_mu(:)
    integer, intent(in) :: source_nlon(:), target_nlon(:)
    real(real64), intent(inout) :: destination(:, :)
    logical, intent(in), optional :: valid(:, :), needed(:, :)
    integer, intent(out), optional :: fallback_points
    real(real64), parameter :: pi = acos(-1.0_real64)
    real(real64), allocatable :: source_latitude(:)
    real(real64) :: latitude, longitude, ring_weight(2), weight, value_sum, weight_sum, plain_sum, x, f
    integer :: ring(2), jt, it, r, js, n, i0, points(2), p, fallback
    logical :: point_needed

    source_latitude = asin(source_mu)
    fallback = 0
    do jt = 1, size(target_nlon)
      latitude = asin(target_mu(jt))
      if (latitude <= source_latitude(1)) then
        ring = [1, 1]
        ring_weight = [1.0_real64, 0.0_real64]
      else if (latitude >= source_latitude(size(source_nlon))) then
        ring = size(source_nlon)
        ring_weight = [1.0_real64, 0.0_real64]
      else
        js = 1
        do while (source_latitude(js + 1) <= latitude)
          js = js + 1
        end do
        ring = [js, js + 1]
        ring_weight(2) = (latitude - source_latitude(js))/(source_latitude(js + 1) - source_latitude(js))
        ring_weight(1) = 1.0_real64 - ring_weight(2)
      end if
      do it = 1, target_nlon(jt)
        longitude = 2.0_real64*pi*real(it - 1, real64)/real(target_nlon(jt), real64)
        value_sum = 0.0_real64
        weight_sum = 0.0_real64
        plain_sum = 0.0_real64
        do r = 1, 2
          if (ring_weight(r) <= 0.0_real64) cycle
          js = ring(r)
          n = source_nlon(js)
          x = real(it - 1, real64)/real(target_nlon(jt), real64)*real(n, real64)
          i0 = floor(x)
          f = x - real(i0, real64)
          points = [modulo(i0, n) + 1, modulo(i0 + 1, n) + 1]
          do p = 1, 2
            weight = ring_weight(r)*merge(1.0_real64 - f, f, p == 1)
            if (weight <= 0.0_real64) cycle
            plain_sum = plain_sum + weight*source(points(p), js)
            if (present(valid)) then
              if (.not. valid(points(p), js)) cycle
            end if
            value_sum = value_sum + weight*source(points(p), js)
            weight_sum = weight_sum + weight
          end do
        end do
        point_needed = .true.
        if (present(needed)) point_needed = needed(it, jt)
        if (weight_sum > 0.0_real64) then
          destination(it, jt) = value_sum/weight_sum
        else if (point_needed) then
          destination(it, jt) = nearest_valid_value(latitude, longitude)
          fallback = fallback + 1
        else
          destination(it, jt) = plain_sum
        end if
      end do
    end do
    if (present(fallback_points)) fallback_points = fallback

  contains

    real(real64) function nearest_valid_value(latitude, longitude) result(value)
      real(real64), intent(in) :: latitude, longitude
      real(real64) :: best, closeness, source_longitude
      integer :: i, j
      logical :: found

      found = .false.
      best = -2.0_real64
      value = 0.0_real64
      do j = 1, size(source_nlon)
        do i = 1, source_nlon(j)
          if (.not. valid(i, j)) cycle
          source_longitude = 2.0_real64*pi*real(i - 1, real64)/real(source_nlon(j), real64)
          ! Cosine of the great-circle angle; larger is closer.
          closeness = sin(latitude)*sin(source_latitude(j)) + &
                      cos(latitude)*cos(source_latitude(j))*cos(longitude - source_longitude)
          if (closeness > best) then
            best = closeness
            value = source(i, j)
            found = .true.
          end if
        end do
      end do
      if (.not. found) error stop 'regridding found no valid source point'
    end function nearest_valid_value

  end subroutine regrid_ring_field

  !> Reads a ring-ordered field (field_binary_writer%write_field) into the padded
  !> array `field`; the padding is zeroed.
  subroutine read_ring_field(path, nlon, field)
    character(*), intent(in) :: path
    integer, intent(in) :: nlon(:)
    real(real64), intent(inout) :: field(:, :)
    integer :: unit, j, status
    integer(8) :: bytes
    logical :: exists

    inquire (file=path, exist=exists, size=bytes)
    if (.not. exists) then
      write (*, '(2a)') 'missing restart file: ', path
      error stop 'restart source file is missing; run the lower-resolution case first'
    end if
    if (bytes /= 8_8*int(sum(nlon), 8)) then
      write (*, '(2a)') 'unexpected size of restart file: ', path
      error stop 'restart grid field has the wrong size for the source truncation'
    end if
    field = 0.0_real64
    open (newunit=unit, file=path, status='old', access='stream', form='unformatted', action='read', &
          convert='little_endian', iostat=status)
    if (status /= 0) error stop 'cannot open restart grid field'
    do j = 1, size(nlon)
      read (unit) field(1:nlon(j), j)
      if (.not. all(ieee_is_finite(field(1:nlon(j), j)))) error stop 'nonfinite value in restart grid field'
    end do
    close (unit)
  end subroutine read_ring_field

  !> Reads a spectral field (field_binary_writer%write_spectral_field) of truncation
  !> T into a (0:T+1, 0:T) array whose row n = T+1 is zero.
  subroutine read_spectral_snapshot_field(path, truncation, field)
    character(*), intent(in) :: path
    integer, intent(in) :: truncation
    complex(real64), allocatable, intent(out) :: field(:, :)
    real(real64) :: parts(2)
    integer :: unit, n, m, status
    integer(8) :: bytes
    logical :: exists

    inquire (file=path, exist=exists, size=bytes)
    if (.not. exists) then
      write (*, '(2a)') 'missing restart file: ', path
      error stop 'restart source file is missing; run the lower-resolution case first'
    end if
    if (bytes /= 16_8*int(truncation + 1, 8)**2) then
      write (*, '(2a)') 'unexpected size of restart file: ', path
      error stop 'restart spectral field has the wrong size for the source truncation'
    end if
    allocate (field(0:truncation + 1, 0:truncation))
    field = cmplx(0.0_real64, 0.0_real64, kind=real64)
    open (newunit=unit, file=path, status='old', access='stream', form='unformatted', action='read', &
          convert='little_endian', iostat=status)
    if (status /= 0) error stop 'cannot open restart spectral field'
    do m = 0, truncation
      do n = 0, truncation
        read (unit) parts
        field(n, m) = cmplx(parts(1), parts(2), kind=real64)
      end do
    end do
    close (unit)
    if (.not. all(ieee_is_finite(real(field, real64))) .or. .not. all(ieee_is_finite(aimag(field)))) then
      error stop 'nonfinite value in restart spectral field'
    end if
  end subroutine read_spectral_snapshot_field

  !> The "time_seconds" entry of a yearly_time_y*.json file.
  real(real64) function read_snapshot_time(path) result(time_seconds)
    character(*), intent(in) :: path
    character(len=256) :: line
    integer :: unit, status, colon, comma
    logical :: exists, found

    inquire (file=path, exist=exists)
    if (.not. exists) then
      write (*, '(2a)') 'missing restart file: ', path
      error stop 'restart source file is missing; run the lower-resolution case first'
    end if
    found = .false.
    time_seconds = 0.0_real64
    open (newunit=unit, file=path, status='old', action='read')
    do
      read (unit, '(a)', iostat=status) line
      if (status /= 0) exit
      if (index(line, '"time_seconds"') == 0) cycle
      colon = index(line, ':')
      comma = index(line, ',', back=.true.)
      if (comma <= colon) comma = len_trim(line) + 1
      read (line(colon + 1:comma - 1), *, iostat=status) time_seconds
      if (status /= 0) error stop 'cannot parse restart snapshot time'
      found = .true.
      exit
    end do
    close (unit)
    if (.not. found) error stop 'restart snapshot time is missing'
  end function read_snapshot_time

end module land_sea_restart
