!> Analytic continents and mountain ranges used by the land--sea case.
module topography
  use iso_fortran_env, only: real64
  use harmonics, only: harmonic_transform
  use planet_parameters, only: earth_gravity
  implicit none
  private

  type, public :: ellipse_continent_config
    real(real64) :: longitude_degrees = 0.0_real64
    real(real64) :: latitude_degrees = 0.0_real64
    real(real64) :: semi_axis_east_degrees = 1.0_real64
    real(real64) :: semi_axis_north_degrees = 1.0_real64
    real(real64) :: orientation_degrees = 0.0_real64
  end type ellipse_continent_config

  type, public :: polar_cap_config
    real(real64) :: edge_latitude_degrees = 70.0_real64
    integer :: hemisphere_sign = -1
  end type polar_cap_config

  type, public :: mountain_range_config
    real(real64) :: longitude_degrees = 0.0_real64
    real(real64) :: latitude_degrees = 0.0_real64
    real(real64) :: length_degrees = 20.0_real64
    real(real64) :: half_width_degrees = 6.0_real64
    real(real64) :: orientation_degrees = 0.0_real64
    real(real64) :: height_metres = 1000.0_real64
  end type mountain_range_config

  type, public :: topography_config
    type(ellipse_continent_config) :: continent_a = &
      ellipse_continent_config(60.0_real64, 45.0_real64, 55.0_real64, 28.0_real64, 0.0_real64)
    type(ellipse_continent_config) :: continent_b = &
      ellipse_continent_config(300.0_real64, -10.0_real64, 28.0_real64, 40.0_real64, 0.0_real64)
    type(ellipse_continent_config) :: continent_c = &
      ellipse_continent_config(150.0_real64, -30.0_real64, 22.0_real64, 16.0_real64, 0.0_real64)
    type(polar_cap_config) :: south_polar_cap = polar_cap_config(70.0_real64, -1)
    type(mountain_range_config) :: mountain_a = &
      mountain_range_config(75.0_real64, 35.0_real64, 50.0_real64, 6.0_real64, 0.0_real64, 2500.0_real64)
    type(mountain_range_config) :: mountain_b = &
      mountain_range_config(295.0_real64, -10.0_real64, 60.0_real64, 6.0_real64, 90.0_real64, 2500.0_real64)
    real(real64) :: coast_width_degrees = 6.0_real64
    real(real64) :: base_height_metres = 300.0_real64
  end type topography_config

  type, public :: topography_diagnostics
    real(real64) :: global_land_fraction = 0.0_real64
    real(real64) :: minimum_truncated_height_metres = 0.0_real64
    real(real64) :: maximum_truncated_height_metres = 0.0_real64
    real(real64) :: height_rms_error_metres = 0.0_real64
  end type topography_diagnostics

  public :: generate_topography

contains

  !> Generate the untruncated land fraction and analytic height, then project the
  !> geopotential to the model truncation and diagnose the projection error.
  subroutine generate_topography(transform, config, land_fraction, analytic_height, &
                                 surface_geopotential, truncated_height, diagnostics)
    type(harmonic_transform), intent(in) :: transform
    type(topography_config), intent(in) :: config
    real(real64), allocatable, intent(out) :: land_fraction(:, :), analytic_height(:, :)
    complex(real64), allocatable, intent(out) :: surface_geopotential(:, :)
    real(real64), allocatable, intent(out) :: truncated_height(:, :)
    type(topography_diagnostics), intent(out) :: diagnostics
    real(real64), allocatable :: shape_a(:, :), shape_b(:, :), shape_c(:, :), shape_cap(:, :)
    real(real64), allocatable :: ridge_a(:, :), ridge_b(:, :), geopotential_grid(:, :)
    real(real64), allocatable :: weights(:)
    integer, allocatable :: nlon(:)
    integer :: i, j
    real(real64) :: area_weight, squared_error

    if (config%coast_width_degrees <= 0.0_real64) error stop 'topography coast width must be positive'
    call transform%allocate_field(land_fraction)
    allocate (analytic_height, mold=land_fraction)
    allocate (shape_a, mold=land_fraction)
    allocate (shape_b, mold=land_fraction)
    allocate (shape_c, mold=land_fraction)
    allocate (shape_cap, mold=land_fraction)
    allocate (ridge_a, mold=land_fraction)
    allocate (ridge_b, mold=land_fraction)
    allocate (geopotential_grid, mold=land_fraction)
    nlon = transform%get_nlon()
    weights = transform%get_gaussian_weights()
    land_fraction = 0.0_real64
    analytic_height = 0.0_real64
    shape_a = 0.0_real64
    shape_b = 0.0_real64
    shape_c = 0.0_real64
    shape_cap = 0.0_real64
    ridge_a = 0.0_real64
    ridge_b = 0.0_real64

    call fill_ellipse(transform, nlon, config%continent_a, config%coast_width_degrees, shape_a)
    call fill_ellipse(transform, nlon, config%continent_b, config%coast_width_degrees, shape_b)
    call fill_ellipse(transform, nlon, config%continent_c, config%coast_width_degrees, shape_c)
    call fill_polar_cap(transform, nlon, config%south_polar_cap, config%coast_width_degrees, shape_cap)
    call fill_mountain(transform, nlon, config%mountain_a, ridge_a)
    call fill_mountain(transform, nlon, config%mountain_b, ridge_b)

    do j = 1, size(nlon)
      do i = 1, nlon(j)
        land_fraction(i, j) = 1.0_real64 - (1.0_real64 - shape_a(i, j))* &
          (1.0_real64 - shape_b(i, j))*(1.0_real64 - shape_c(i, j))*(1.0_real64 - shape_cap(i, j))
        land_fraction(i, j) = max(0.0_real64, min(1.0_real64, land_fraction(i, j)))
        analytic_height(i, j) = land_fraction(i, j)*(config%base_height_metres + ridge_a(i, j) + ridge_b(i, j))
        geopotential_grid(i, j) = earth_gravity*analytic_height(i, j)
      end do
    end do
    call transform%grid_to_spectral(geopotential_grid, surface_geopotential)
    call transform%spectral_to_grid(surface_geopotential, truncated_height)
    truncated_height = truncated_height/earth_gravity

    diagnostics = topography_diagnostics()
    diagnostics%minimum_truncated_height_metres = huge(0.0_real64)
    diagnostics%maximum_truncated_height_metres = -huge(0.0_real64)
    squared_error = 0.0_real64
    do j = 1, size(nlon)
      area_weight = 0.5_real64*weights(j)/real(nlon(j), real64)
      do i = 1, nlon(j)
        diagnostics%global_land_fraction = diagnostics%global_land_fraction + area_weight*land_fraction(i, j)
        diagnostics%minimum_truncated_height_metres = min(diagnostics%minimum_truncated_height_metres, &
                                                           truncated_height(i, j))
        diagnostics%maximum_truncated_height_metres = max(diagnostics%maximum_truncated_height_metres, &
                                                           truncated_height(i, j))
        squared_error = squared_error + area_weight*(truncated_height(i, j) - analytic_height(i, j))**2
      end do
    end do
    diagnostics%height_rms_error_metres = sqrt(max(squared_error, 0.0_real64))
  end subroutine generate_topography

  subroutine local_coordinates(longitude_degrees, latitude_degrees, centre_longitude, centre_latitude, &
                               orientation, x_rotated, y_rotated)
    real(real64), intent(in) :: longitude_degrees, latitude_degrees, centre_longitude, centre_latitude, orientation
    real(real64), intent(out) :: x_rotated, y_rotated
    real(real64) :: delta_longitude, x, y, theta

    delta_longitude = modulo(longitude_degrees - centre_longitude + 180.0_real64, 360.0_real64) - 180.0_real64
    x = cos(latitude_degrees*acos(-1.0_real64)/180.0_real64)*delta_longitude
    y = latitude_degrees - centre_latitude
    theta = orientation*acos(-1.0_real64)/180.0_real64
    x_rotated = x*cos(theta) + y*sin(theta)
    y_rotated = -x*sin(theta) + y*cos(theta)
  end subroutine local_coordinates

  subroutine fill_ellipse(transform, nlon, config, coast_width, field)
    type(harmonic_transform), intent(in) :: transform
    integer, intent(in) :: nlon(:)
    type(ellipse_continent_config), intent(in) :: config
    real(real64), intent(in) :: coast_width
    real(real64), intent(inout) :: field(:, :)
    integer :: i, j
    real(real64) :: longitude, latitude, x, y, distance

    if (min(config%semi_axis_east_degrees, config%semi_axis_north_degrees) <= 0.0_real64) then
      error stop 'topography ellipse axes must be positive'
    end if
    do j = 1, size(nlon)
      latitude = asin(transform%mu(j))*180.0_real64/acos(-1.0_real64)
      do i = 1, nlon(j)
        longitude = 360.0_real64*real(i - 1, real64)/real(nlon(j), real64)
        call local_coordinates(longitude, latitude, config%longitude_degrees, config%latitude_degrees, &
                               config%orientation_degrees, x, y)
        distance = sqrt((x/config%semi_axis_east_degrees)**2 + (y/config%semi_axis_north_degrees)**2)
        field(i, j) = 0.5_real64*(1.0_real64 - tanh((distance - 1.0_real64)* &
          min(config%semi_axis_east_degrees, config%semi_axis_north_degrees)/coast_width))
      end do
    end do
  end subroutine fill_ellipse

  subroutine fill_polar_cap(transform, nlon, config, coast_width, field)
    type(harmonic_transform), intent(in) :: transform
    integer, intent(in) :: nlon(:)
    type(polar_cap_config), intent(in) :: config
    real(real64), intent(in) :: coast_width
    real(real64), intent(inout) :: field(:, :)
    integer :: j
    real(real64) :: latitude, value

    if (abs(config%hemisphere_sign) /= 1) error stop 'polar cap hemisphere sign must be +1 or -1'
    do j = 1, size(nlon)
      latitude = asin(transform%mu(j))*180.0_real64/acos(-1.0_real64)
      value = 0.5_real64*(1.0_real64 + tanh((real(config%hemisphere_sign, real64)*latitude - &
        config%edge_latitude_degrees)/coast_width))
      field(1:nlon(j), j) = value
    end do
  end subroutine fill_polar_cap

  subroutine fill_mountain(transform, nlon, config, field)
    type(harmonic_transform), intent(in) :: transform
    integer, intent(in) :: nlon(:)
    type(mountain_range_config), intent(in) :: config
    real(real64), intent(inout) :: field(:, :)
    integer :: i, j
    real(real64) :: longitude, latitude, x, y

    if (config%length_degrees <= 0.0_real64 .or. config%half_width_degrees <= 0.0_real64) then
      error stop 'mountain length and half width must be positive'
    end if
    do j = 1, size(nlon)
      latitude = asin(transform%mu(j))*180.0_real64/acos(-1.0_real64)
      do i = 1, nlon(j)
        longitude = 360.0_real64*real(i - 1, real64)/real(nlon(j), real64)
        call local_coordinates(longitude, latitude, config%longitude_degrees, config%latitude_degrees, &
                               config%orientation_degrees, x, y)
        field(i, j) = config%height_metres*exp(-(y/config%half_width_degrees)**2)* &
          0.5_real64*(1.0_real64 - tanh((abs(x) - 0.5_real64*config%length_degrees)/ &
                                        config%half_width_degrees))
      end do
    end do
  end subroutine fill_mountain

end module topography
