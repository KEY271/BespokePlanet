!> Unit checks for docs/dynamics/earth-topography.md and docs/cases/land-sea-earth.md.
program check_earth_topography
  use iso_fortran_env, only: real64
  use harmonics, only: harmonic_transform
  use earth_topography, only: earth_topography_config, earth_topography_diagnostics, generate_earth_topography, &
                              read_earth_topography_data, earth_data_longitudes, earth_data_latitudes
  use dry_vertical_coordinate, only: hybrid_sigma_coordinate
  use dry_initial_conditions, only: jablonowski_williamson_topographic_initial_state, topographic_surface_pressure
  use dry_physics_config, only: dry_model_physics_config
  use dry_case_initial_conditions, only: land_sea_case_physics, radiation_case_planet, set_land_sea_case_state
  use planet_parameters, only: planet_config, earth_gravity
  use dry_atmosphere, only: dry_atmosphere_solver
  use topography, only: topography_config, topography_diagnostics, generate_topography
  use numerics_config, only: model_numerics_config
  use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
  implicit none

  type :: generated_terrain
    type(harmonic_transform) :: transform
    real(real64), allocatable :: land_fraction(:, :), target_height(:, :), truncated_height(:, :)
    complex(real64), allocatable :: surface_geopotential(:, :)
    type(earth_topography_diagnostics) :: diagnostics
  end type generated_terrain

  type(generated_terrain) :: t31, t63

  call check_intermediate_file()
  call generate(t31, 31)
  call generate(t63, 63)
  call check_generation(t31)
  call check_generation(t63)
  call check_resolution_dependence()
  call check_surface_pressure(t31)
  call check_surface_pressure(t63)
  call check_dry_integration(t31, 31)
  call check_dry_integration(t63, 63)
  call check_land_integration(t31, 31)
  write (*, '(a)') 'check_earth_topography: all checks passed'

contains

  subroutine check_intermediate_file()
    real(real64), allocatable :: land(:, :), height(:, :)
    character(len=:), allocatable :: path

    call read_earth_topography_data('data/earth_topography_0p5deg.bin', land, height, path)
    if (size(land, 1) /= earth_data_longitudes .or. size(land, 2) /= earth_data_latitudes) then
      error stop 'intermediate file has the wrong shape'
    end if
    if (minval(land) < 0.0_real64 .or. maxval(land) > 1.0_real64 .or. minval(height) < 0.0_real64) then
      error stop 'intermediate file values are out of range'
    end if
    if (maxval(height) < 5000.0_real64 .or. maxval(height) > 6500.0_real64) then
      error stop 'intermediate file maximum height is not the Karakoram cell mean'
    end if
    write (*, '(a,a)') 'intermediate file: ', path
  end subroutine check_intermediate_file

  subroutine generate(terrain, truncation)
    type(generated_terrain), intent(inout) :: terrain
    integer, intent(in) :: truncation
    type(earth_topography_config) :: config

    call terrain%transform%init(truncation)
    call generate_earth_topography(terrain%transform, config, terrain%land_fraction, terrain%target_height, &
                                   terrain%surface_geopotential, terrain%truncated_height, terrain%diagnostics)
    write (*, '(a,i0,a,f8.4,a,f8.2,a,f8.2,a,f8.2,a,f8.2,a,f8.2)') 'T', truncation, &
      ': <f_L> = ', terrain%diagnostics%global_land_fraction, &
      '  min z_s = ', terrain%diagnostics%minimum_truncated_height_metres, &
      '  max z_s = ', terrain%diagnostics%maximum_truncated_height_metres, &
      '  open-ocean RMS = ', terrain%diagnostics%open_ocean_rms_metres, &
      '  truncation RMS = ', terrain%diagnostics%truncation_rms_metres, &
      '  total RMS = ', terrain%diagnostics%total_rms_metres
    write (*, '(a,f8.2,a,f8.2,a,f8.2,a)') '  ocean (f_L < 0.01) z_s: RMS ', terrain%diagnostics%ocean_height_rms_metres, &
      ' m, max ', terrain%diagnostics%ocean_maximum_height_metres, ' m; kernel half width ', &
      terrain%diagnostics%kernel_half_width_degrees, ' deg'
  end subroutine generate

  real(real64) function nearest_value(terrain, field, longitude, latitude) result(value)
    type(generated_terrain), intent(in) :: terrain
    real(real64), intent(in) :: field(:, :), longitude, latitude
    integer, allocatable :: nlon(:)
    integer :: j, best_j, i
    real(real64) :: best_distance, ring_latitude

    nlon = terrain%transform%get_nlon()
    best_j = 1
    best_distance = huge(0.0_real64)
    do j = 1, size(nlon)
      ring_latitude = asin(terrain%transform%mu(j))*180.0_real64/acos(-1.0_real64)
      if (abs(ring_latitude - latitude) < best_distance) then
        best_distance = abs(ring_latitude - latitude)
        best_j = j
      end if
    end do
    i = modulo(nint(longitude*real(nlon(best_j), real64)/360.0_real64), nlon(best_j)) + 1
    value = field(i, best_j)
  end function nearest_value

  subroutine check_generation(terrain)
    type(generated_terrain), intent(in) :: terrain
    integer, allocatable :: nlon(:)
    integer :: j
    real(real64) :: tibet, antarctica, pacific_height, pacific_land

    nlon = terrain%transform%get_nlon()
    do j = 1, size(nlon)
      if (minval(terrain%land_fraction(1:nlon(j), j)) < 0.0_real64 .or. &
          maxval(terrain%land_fraction(1:nlon(j), j)) > 1.0_real64) then
        error stop 'smoothed land fraction is outside [0,1]'
      end if
      if (minval(terrain%target_height(1:nlon(j), j)) < 0.0_real64) error stop 'smoothed target height is negative'
    end do
    if (terrain%diagnostics%minimum_truncated_height_metres < -50.0_real64) then
      error stop 'truncated Earth terrain has excessive negative overshoot'
    end if
    if (terrain%diagnostics%open_ocean_rms_metres > 10.0_real64) error stop 'open-ocean ripple RMS exceeds 10 m'
    if (abs(terrain%diagnostics%global_land_fraction - terrain%diagnostics%source_land_fraction) > 0.01_real64) then
      error stop 'grid land fraction differs from the intermediate file by more than 0.01'
    end if
    if (terrain%diagnostics%truncation_rms_metres > 10.0_real64) error stop 'truncation RMS error exceeds 10 m'
    tibet = nearest_value(terrain, terrain%truncated_height, 90.0_real64, 33.0_real64)
    antarctica = nearest_value(terrain, terrain%truncated_height, 90.0_real64, -80.0_real64)
    pacific_height = nearest_value(terrain, terrain%truncated_height, 200.0_real64, 0.0_real64)
    pacific_land = nearest_value(terrain, terrain%land_fraction, 200.0_real64, 0.0_real64)
    write (*, '(a,f8.1,a,f8.1,a,f8.2,a,es10.2)') '  Tibet ', tibet, ' m, Antarctica ', antarctica, &
      ' m, central Pacific ', pacific_height, ' m, f_L ', pacific_land
    if (tibet < 2000.0_real64) error stop 'Tibet is lower than 2000 m'
    if (antarctica < 2500.0_real64) error stop 'the Antarctic plateau is lower than 2500 m'
    if (abs(pacific_height) > 20.0_real64 .or. pacific_land > 0.01_real64) error stop 'central Pacific is not flat ocean'
  end subroutine check_generation

  subroutine check_resolution_dependence()
    if (t63%diagnostics%total_rms_metres >= t31%diagnostics%total_rms_metres) then
      error stop 'T63 total RMS error is not smaller than at T31'
    end if
    if (abs(t63%diagnostics%global_land_fraction - t31%diagnostics%global_land_fraction) > &
        0.005_real64*t31%diagnostics%global_land_fraction) then
      error stop 'T31 and T63 global land fractions differ by more than 0.5 percent'
    end if
    if (t63%diagnostics%maximum_truncated_height_metres <= t31%diagnostics%maximum_truncated_height_metres) then
      error stop 'T63 does not retain more of the highest terrain than T31'
    end if
  end subroutine check_resolution_dependence

  subroutine check_surface_pressure(terrain)
    type(generated_terrain), intent(in) :: terrain
    type(dry_model_physics_config) :: physics
    type(planet_config) :: planet
    integer, allocatable :: nlon(:)
    integer :: i, j
    real(real64) :: pressure, minimum_pressure

    physics = land_sea_case_physics()
    planet = radiation_case_planet(physics)
    nlon = terrain%transform%get_nlon()
    minimum_pressure = huge(0.0_real64)
    do j = 1, size(nlon)
      do i = 1, nlon(j)
        pressure = topographic_surface_pressure(terrain%transform%mu(j), earth_gravity*terrain%truncated_height(i, j), &
                                                planet%rotation_rate)
        if (.not. ieee_is_finite(pressure)) error stop 'topographic surface pressure is not finite'
        minimum_pressure = min(minimum_pressure, pressure)
      end do
    end do
    write (*, '(a,f10.1,a)') '  minimum surface pressure ', minimum_pressure/100.0_real64, ' hPa'
    if (minimum_pressure < 5.0e4_real64) error stop 'surface pressure falls below 500 hPa'
  end subroutine check_surface_pressure

  !> One day of the dry, physics-free atmosphere over Earth terrain must stay within an
  !> order of magnitude of the same integration over the analytic land--sea terrain
  !> (docs/cases/land-sea.md).  The flat surface is printed for reference only: the
  !> Jablonowski-Williamson state has no meridional wind, so a ratio to it is not meaningful.
  subroutine check_dry_integration(terrain, truncation)
    type(generated_terrain), intent(inout) :: terrain
    integer, intent(in) :: truncation
    type(topography_config) :: analytic_config
    type(topography_diagnostics) :: analytic_diagnostics
    real(real64), allocatable :: analytic_land(:, :), analytic_height(:, :), analytic_truncated(:, :)
    complex(real64), allocatable :: analytic_geopotential(:, :), flat_geopotential(:, :)
    real(real64) :: earth_divergence, earth_wind, analytic_divergence, analytic_wind, flat_divergence, flat_wind

    allocate (flat_geopotential(0:truncation + 1, 0:truncation))
    flat_geopotential = cmplx(0.0_real64, 0.0_real64, kind=real64)
    call generate_topography(terrain%transform, analytic_config, analytic_land, analytic_height, &
                             analytic_geopotential, analytic_truncated, analytic_diagnostics)
    call integrate_dry_day(terrain, truncation, flat_geopotential, flat_divergence, flat_wind)
    call integrate_dry_day(terrain, truncation, analytic_geopotential, analytic_divergence, analytic_wind)
    call integrate_dry_day(terrain, truncation, terrain%surface_geopotential, earth_divergence, earth_wind)
    write (*, '(a,i0,a,3es10.2,a,3f7.2)') 'T', truncation, ' one dry day (flat, analytic, Earth): divergence RMS', &
      flat_divergence, analytic_divergence, earth_divergence, '; max |v|', flat_wind, analytic_wind, earth_wind
    if (earth_divergence > 10.0_real64*analytic_divergence) then
      error stop 'divergence over Earth terrain is an order of magnitude larger than over the analytic terrain'
    end if
    if (earth_wind > 10.0_real64*analytic_wind) then
      error stop 'meridional wind over Earth terrain is an order of magnitude larger than over the analytic terrain'
    end if
  end subroutine check_dry_integration

  subroutine integrate_dry_day(terrain, truncation, surface_geopotential, divergence_rms, maximum_v)
    type(generated_terrain), intent(inout) :: terrain
    integer, intent(in) :: truncation
    complex(real64), intent(in) :: surface_geopotential(0:, 0:)
    real(real64), intent(out) :: divergence_rms, maximum_v
    type(dry_atmosphere_solver) :: solver
    type(model_numerics_config) :: numerics
    type(hybrid_sigma_coordinate) :: coordinate
    type(dry_model_physics_config) :: physics
    type(planet_config) :: planet
    complex(real64), allocatable :: zeta(:, :, :), delta(:, :, :), temperature(:, :, :), log_ps(:, :)
    real(real64), allocatable :: zeta_grid(:, :, :), delta_grid(:, :, :), temperature_grid(:, :, :)
    real(real64), allocatable :: surface_pressure(:, :), u(:, :, :), v(:, :, :), weights(:)
    integer, allocatable :: nlon(:)
    integer :: step, j, k
    real(real64) :: total

    physics = land_sea_case_physics()
    planet = radiation_case_planet(physics)
    numerics = model_numerics_config()
    numerics%truncation = truncation
    numerics%time_step = 1200.0_real64
    call solver%init_with_config(numerics)
    coordinate = solver%get_coordinate()
    call jablonowski_williamson_topographic_initial_state(terrain%transform, truncation, coordinate, &
      surface_geopotential, planet%rotation_rate, zeta, delta, temperature, log_ps)
    call solver%set_planet(planet)
    call solver%set_initial_state(zeta, delta, temperature, log_ps, surface_geopotential)
    do step = 1, 72
      call solver%advance()
    end do
    call solver%get_fields(zeta_grid, delta_grid, temperature_grid, surface_pressure, u, v)
    if (.not. all(ieee_is_finite(delta_grid)) .or. .not. all(ieee_is_finite(v))) then
      error stop 'dry integration over terrain produced a non-finite field'
    end if
    nlon = terrain%transform%get_nlon()
    weights = terrain%transform%get_gaussian_weights()
    total = 0.0_real64
    maximum_v = 0.0_real64
    do k = 1, size(delta_grid, 3)
      do j = 1, size(nlon)
        total = total + 0.5_real64*weights(j)*sum(delta_grid(1:nlon(j), j, k)**2)/real(nlon(j), real64)
        maximum_v = max(maximum_v, maxval(abs(v(1:nlon(j), j, k))))
      end do
    end do
    divergence_rms = sqrt(total/real(size(delta_grid, 3), real64))
  end subroutine integrate_dry_day

  subroutine check_land_integration(terrain, truncation)
    type(generated_terrain), intent(inout) :: terrain
    integer, intent(in) :: truncation
    type(dry_model_physics_config) :: physics
    type(planet_config) :: planet
    type(model_numerics_config) :: numerics
    type(dry_atmosphere_solver) :: solver
    real(real64), allocatable :: zeta(:, :, :), delta(:, :, :), temperature(:, :, :)
    real(real64), allocatable :: surface_pressure(:, :), u(:, :, :), v(:, :, :)
    real(real64), allocatable :: surface_temperature(:, :), deep_temperature(:, :), humidity(:, :, :)
    integer :: step

    physics = land_sea_case_physics()
    planet = radiation_case_planet(physics)
    numerics = model_numerics_config()
    numerics%truncation = truncation
    numerics%time_step = 1200.0_real64
    call solver%init_with_config(numerics)
    call set_land_sea_case_state(solver, terrain%transform, physics, planet, terrain%surface_geopotential, &
                                 terrain%land_fraction)
    do step = 1, 3
      call solver%advance()
    end do
    call solver%get_fields(zeta, delta, temperature, surface_pressure, u, v, &
      surface_temperature=surface_temperature, deep_temperature=deep_temperature, specific_humidity=humidity)
    if (.not. all(ieee_is_finite(temperature)) .or. .not. all(ieee_is_finite(surface_pressure)) .or. &
        .not. all(ieee_is_finite(surface_temperature)) .or. .not. all(ieee_is_finite(deep_temperature)) .or. &
        .not. all(ieee_is_finite(humidity))) then
      error stop 'short Earth land-sea integration produced a non-finite field'
    end if
  end subroutine check_land_integration

end program check_earth_topography
