!> Unit checks for docs/cases/land-sea.md.
program check_land_sea
  use iso_fortran_env, only: real64
  use harmonics, only: harmonic_transform
  use topography, only: topography_config, topography_diagnostics, generate_topography
  use dry_vertical_coordinate, only: hybrid_sigma_coordinate, reference_surface_pressure
  use dry_initial_conditions, only: jablonowski_williamson_initial_state, &
                                    jablonowski_williamson_topographic_initial_state, &
                                    topographic_surface_pressure
  use dry_physics_config, only: dry_model_physics_config, mixed_surface_properties
  use dry_case_initial_conditions, only: land_sea_case_physics, radiation_case_planet
  use dry_radiation, only: radiation_tendency
  use planet_parameters, only: earth_gravity, planet_config
  use dry_atmosphere, only: dry_atmosphere_solver
  use numerics_config, only: model_numerics_config
  use dry_case_initial_conditions, only: set_land_sea_case_state
  use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
  implicit none

  call check_default_topography()
  call check_flat_initial_state_identity()
  call check_surface_pressure()
  call check_mixed_surface_budget()
  call check_land_integration()
  write (*, '(a)') 'check_land_sea: all checks passed'

contains

  subroutine check_default_topography()
    type(harmonic_transform) :: transform31, transform63
    type(topography_config) :: config
    type(topography_diagnostics) :: d31, d63
    real(real64), allocatable :: land31(:, :), analytic31(:, :), height31(:, :)
    real(real64), allocatable :: land63(:, :), analytic63(:, :), height63(:, :)
    complex(real64), allocatable :: phi31(:, :), phi63(:, :)

    call transform31%init(31)
    call transform63%init(63)
    call generate_topography(transform31, config, land31, analytic31, phi31, height31, d31)
    call generate_topography(transform63, config, land63, analytic63, phi63, height63, d63)
    if (minval(land31) < 0.0_real64 .or. maxval(land31) > 1.0_real64 .or. &
        minval(land63) < 0.0_real64 .or. maxval(land63) > 1.0_real64) then
      error stop 'generated land fraction is outside [0,1]'
    end if
    if (minval(analytic31) < 0.0_real64 .or. minval(analytic63) < 0.0_real64) then
      error stop 'analytic terrain contains negative height'
    end if
    if (d31%minimum_truncated_height_metres < -50.0_real64 .or. &
        d63%minimum_truncated_height_metres < -50.0_real64) then
      error stop 'truncated terrain has excessive negative overshoot'
    end if
    if (d63%height_rms_error_metres >= d31%height_rms_error_metres) then
      error stop 'T63 terrain RMS error is not smaller than T31 error'
    end if
    if (abs(d63%global_land_fraction - d31%global_land_fraction) > &
        0.005_real64*d31%global_land_fraction) then
      error stop 'T31 and T63 global land fractions differ by more than 0.5 percent'
    end if
    if (config%coast_width_degrees < 180.0_real64/31.0_real64 .or. &
        config%mountain_a%half_width_degrees < 180.0_real64/31.0_real64 .or. &
        config%mountain_b%half_width_degrees < 180.0_real64/31.0_real64) then
      error stop 'default terrain does not meet the T31 width condition'
    end if
    write (*, '(a,2f10.5,a,2f10.3)') 'land fractions T31/T63 = ', &
      d31%global_land_fraction, d63%global_land_fraction, '; RMS heights = ', &
      d31%height_rms_error_metres, d63%height_rms_error_metres
  end subroutine check_default_topography

  subroutine check_flat_initial_state_identity()
    type(harmonic_transform) :: transform
    type(hybrid_sigma_coordinate) :: coordinate
    type(dry_model_physics_config) :: physics
    type(planet_config) :: planet
    complex(real64), allocatable :: zeta_old(:, :, :), delta_old(:, :, :), temperature_old(:, :, :)
    complex(real64), allocatable :: log_ps_old(:, :), phi_old(:, :)
    complex(real64), allocatable :: zeta_new(:, :, :), delta_new(:, :, :), temperature_new(:, :, :)
    complex(real64), allocatable :: log_ps_new(:, :), zero_phi(:, :)

    call transform%init(31)
    call coordinate%init_default()
    physics = land_sea_case_physics()
    planet = radiation_case_planet(physics)
    call jablonowski_williamson_initial_state(transform, 31, coordinate, .false., zeta_old, delta_old, &
      temperature_old, log_ps_old, phi_old, planet%rotation_rate, flat_terrain=.true.)
    allocate (zero_phi(0:32, 0:31))
    zero_phi = cmplx(0.0_real64, 0.0_real64, kind=real64)
    call jablonowski_williamson_topographic_initial_state(transform, 31, coordinate, zero_phi, &
      planet%rotation_rate, zeta_new, delta_new, temperature_new, log_ps_new)
    if (any(zeta_new /= zeta_old) .or. any(delta_new /= delta_old) .or. &
        any(temperature_new /= temperature_old) .or. any(log_ps_new /= log_ps_old)) then
      error stop 'zero-topography initial state differs from the flat radiation state'
    end if
  end subroutine check_flat_initial_state_identity

  subroutine check_surface_pressure()
    type(dry_model_physics_config) :: physics
    type(planet_config) :: planet
    real(real64) :: pressure

    physics = land_sea_case_physics()
    planet = radiation_case_planet(physics)
    pressure = topographic_surface_pressure(0.3_real64, 0.0_real64, planet%rotation_rate)
    if (pressure /= reference_surface_pressure) error stop 'zero terrain does not give exactly p0'
    pressure = topographic_surface_pressure(0.3_real64, 2000.0_real64*earth_gravity, planet%rotation_rate)
    if (pressure < 7.5e4_real64 .or. pressure > 8.2e4_real64) then
      error stop '2000 m topographic surface pressure is outside the expected range'
    end if
  end subroutine check_surface_pressure

  subroutine check_mixed_surface_budget()
    type(dry_model_physics_config) :: physics
    real(real64), parameter :: pressure_half(0:4) = [100.0_real64, 10000.0_real64, 40000.0_real64, &
                                                     70000.0_real64, 100000.0_real64]
    real(real64), parameter :: temperature(4) = [220.0_real64, 245.0_real64, 270.0_real64, 285.0_real64]
    real(real64) :: atmospheric_tendency(4), surface_tendency, deep_tendency
    real(real64) :: incoming, reflected, outgoing, capacity, albedo, wetness, exchange
    real(real64) :: capacity0, albedo0, wetness0, exchange0
    real(real64) :: capacity1, albedo1, wetness1, exchange1, total_energy_tendency
    real(real64), parameter :: land = 0.4_real64, latent_heat_flux = 100.0_real64

    physics = land_sea_case_physics()
    call mixed_surface_properties(physics%radiation, physics%evaporation, 0.0_real64, &
                                  capacity0, albedo0, wetness0, exchange0)
    call mixed_surface_properties(physics%radiation, physics%evaporation, 1.0_real64, &
                                  capacity1, albedo1, wetness1, exchange1)
    call mixed_surface_properties(physics%radiation, physics%evaporation, land, &
                                  capacity, albedo, wetness, exchange)
    if (capacity /= land*capacity1 + (1.0_real64 - land)*capacity0 .or. &
        albedo /= land*albedo1 + (1.0_real64 - land)*albedo0 .or. &
        wetness /= land*wetness1 + (1.0_real64 - land)*wetness0 .or. &
        exchange /= land*exchange1 + (1.0_real64 - land)*exchange0) then
      error stop 'mixed surface coefficients are not linear in land fraction'
    end if
    call radiation_tendency(physics%radiation, pressure_half, temperature, 288.0_real64, 286.0_real64, &
      3.0_real64, 4.0_real64, 0.0_real64, acos(-1.0_real64), 0.0_real64, atmospheric_tendency, &
      surface_tendency, deep_tendency, incoming, reflected, outgoing, latent_heat_flux=latent_heat_flux, &
      land_fraction=land)
    total_energy_tendency = sum(physics%radiation%dry_air_specific_heat* &
      (pressure_half(1:4) - pressure_half(0:3))/earth_gravity*atmospheric_tendency) + &
      capacity*surface_tendency + physics%radiation%deep_ground_heat_capacity*deep_tendency + latent_heat_flux
    if (abs(total_energy_tendency - (incoming - reflected - outgoing)) > 1.0e-9_real64) then
      error stop 'mixed surface column energy budget does not close'
    end if
  end subroutine check_mixed_surface_budget

  subroutine check_land_integration()
    type(harmonic_transform) :: transform
    type(topography_config) :: terrain
    type(topography_diagnostics) :: terrain_diagnostics
    type(dry_model_physics_config) :: physics
    type(planet_config) :: planet
    type(model_numerics_config) :: numerics
    type(dry_atmosphere_solver) :: solver
    real(real64), allocatable :: land_fraction(:, :), analytic_height(:, :), height(:, :)
    real(real64), allocatable :: zeta(:, :, :), delta(:, :, :), temperature(:, :, :)
    real(real64), allocatable :: surface_pressure(:, :), u(:, :, :), v(:, :, :)
    real(real64), allocatable :: surface_temperature(:, :), deep_temperature(:, :), humidity(:, :, :)
    complex(real64), allocatable :: surface_geopotential(:, :)
    integer :: step

    call transform%init(31)
    call generate_topography(transform, terrain, land_fraction, analytic_height, surface_geopotential, &
                             height, terrain_diagnostics)
    physics = land_sea_case_physics()
    planet = radiation_case_planet(physics)
    numerics = model_numerics_config()
    numerics%truncation = 31
    numerics%time_step = 1200.0_real64
    call solver%init_with_config(numerics)
    call set_land_sea_case_state(solver, transform, physics, planet, surface_geopotential, land_fraction)
    do step = 1, 3
      call solver%advance()
    end do
    call solver%get_fields(zeta, delta, temperature, surface_pressure, u, v, &
      surface_temperature=surface_temperature, deep_temperature=deep_temperature, specific_humidity=humidity)
    if (.not. all(ieee_is_finite(temperature)) .or. .not. all(ieee_is_finite(surface_pressure)) .or. &
        .not. all(ieee_is_finite(surface_temperature)) .or. .not. all(ieee_is_finite(deep_temperature)) .or. &
        .not. all(ieee_is_finite(humidity))) then
      error stop 'short land-sea integration produced a non-finite field'
    end if
  end subroutine check_land_integration

end program check_land_sea
