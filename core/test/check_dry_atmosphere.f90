program check_dry_atmosphere
  use iso_fortran_env, only: real64
  use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
  use harmonics, only: harmonic_transform
  use dry_vertical_coordinate, only: hybrid_sigma_coordinate, reference_surface_pressure, &
                                     dry_air_gas_constant, dry_air_kappa
  use dry_gravity_wave, only: dry_gravity_wave_solver, dry_gravity_wave_implicitness
  use dry_atmosphere, only: dry_atmosphere_solver
  use dry_initial_conditions, only: jablonowski_williamson_initial_state
  use dry_convection, only: dry_convective_adjustment_tendency, dry_convective_adjustment_time
  use dry_held_suarez, only: held_suarez_forcing, held_suarez_initial_temperature, &
                              held_suarez_temperature_perturbation, held_suarez_sigma_boundary, &
                              held_suarez_upper_thermal_rate, held_suarez_lower_thermal_rate, &
                              held_suarez_friction_rate, held_suarez_minimum_equilibrium_temperature, &
                              held_suarez_equatorial_temperature, held_suarez_equator_to_pole_difference, &
                              held_suarez_vertical_difference
  use dry_radiation, only: radiation_tendency, shortwave_downward_flux, stefan_boltzmann_constant, &
                           longwave_surface_optical_depth, dry_air_specific_heat, &
                           dry_gravity_acceleration, solar_constant, surface_shortwave_albedo, &
                           ultraviolet_shortwave_fraction, ozone_shortwave_optical_depth, &
                           ozone_longwave_optical_depth, ozone_pressure_lower_bound, &
                           ozone_pressure_upper_bound, ozone_peak_pressure, ozone_log_pressure_width, &
                           ozone_layer_optical_depth, ozone_longwave_layer_optical_depth, &
                           surface_exchange_coefficient, gustiness_speed, &
                           surface_heat_capacity, deep_ground_heat_capacity, ground_exchange_coefficient, &
                           axial_tilt, orbital_period, &
                           solar_day, days_per_month, months_per_year, days_per_year, planetary_rotation_rate, &
                           radiation_calendar_date, radiation_diagnostics, radiation_daily_accumulator, &
                           radiation_top_rayleigh_levels, radiation_top_rayleigh_rate, radiation_rayleigh_rate
  implicit none

  integer, parameter :: truncation = 5
  real(real64), parameter :: time_step = 1200.0_real64

  call check_reference_atmosphere()
  call check_precomputed_gravity_wave_inverse()
  call check_resting_atmosphere()
  call check_held_suarez_forcing()
  call check_held_suarez_state()
  call check_dry_convective_adjustment()
  call check_radiation_top_rayleigh_friction()
  call check_solar_geometry()
  call check_ozone_absorption()
  call check_radiation_column()
  call check_radiation_daily_accumulator()
  call check_radiation_state()
  call check_jablonowski_state()
  call check_jablonowski_steady_state()
  call check_flat_terrain_balanced_state()

contains

  subroutine check_radiation_top_rayleigh_friction()
    if (radiation_top_rayleigh_levels /= 2 .or. &
        abs(radiation_top_rayleigh_rate - 1.0_real64/solar_day) > 1.0e-15_real64 .or. &
        abs(radiation_rayleigh_rate(1) - radiation_top_rayleigh_rate) > 1.0e-15_real64 .or. &
        abs(radiation_rayleigh_rate(2) - radiation_top_rayleigh_rate) > 1.0e-15_real64 .or. &
        abs(radiation_rayleigh_rate(3)) > 1.0e-15_real64) then
      error stop 'radiation top-level Rayleigh friction is incorrect'
    end if
  end subroutine check_radiation_top_rayleigh_friction

  subroutine check_dry_convective_adjustment()
    real(real64), parameter :: pressure_half(0:4) = [ &
      1000.0_real64, 20000.0_real64, 50000.0_real64, 80000.0_real64, 100000.0_real64]
    real(real64), parameter :: unstable_potential_temperature(4) = [ &
      340.0_real64, 280.0_real64, 300.0_real64, 250.0_real64]
    real(real64), parameter :: stable_potential_temperature(4) = [ &
      340.0_real64, 300.0_real64, 280.0_real64, 250.0_real64]
    real(real64) :: exner(4), delta_p(4), temperature(4), tendency(4)
    real(real64) :: reference_temperature(4), reference_potential_temperature(4)
    real(real64) :: layer_log_pressure, alpha
    integer :: k

    do k = 1, 4
      delta_p(k) = pressure_half(k) - pressure_half(k - 1)
      layer_log_pressure = log(pressure_half(k)/pressure_half(k - 1))
      alpha = 1.0_real64 - pressure_half(k - 1)*layer_log_pressure/delta_p(k)
      exner(k) = (pressure_half(k)*exp(-alpha)/reference_surface_pressure)**dry_air_kappa
    end do

    temperature = exner*stable_potential_temperature
    call dry_convective_adjustment_tendency(pressure_half, temperature, tendency)
    if (maxval(abs(tendency)) > 1.0e-15_real64) then
      error stop 'dry convective adjustment changed a stable column'
    end if

    temperature = exner*unstable_potential_temperature
    call dry_convective_adjustment_tendency(pressure_half, temperature, tendency)
    reference_temperature = temperature + dry_convective_adjustment_time*tendency
    reference_potential_temperature = reference_temperature/exner
    if (any(reference_potential_temperature(1:3) < &
            reference_potential_temperature(2:4) - 1.0e-12_real64)) then
      error stop 'dry convective adjustment reference profile is unstable'
    end if
    if (abs(reference_potential_temperature(2) - reference_potential_temperature(3)) > 1.0e-12_real64 .or. &
        abs(tendency(1)) > 1.0e-15_real64 .or. abs(tendency(4)) > 1.0e-15_real64) then
      error stop 'dry convective adjustment pooled the wrong layers'
    end if
    if (abs(sum(tendency*delta_p)) > 1.0e-10_real64) then
      error stop 'dry convective adjustment does not conserve column enthalpy'
    end if
  end subroutine check_dry_convective_adjustment

  subroutine check_solar_geometry()
    type(harmonic_transform) :: transform
    integer, allocatable :: nlon(:)
    real(real64), allocatable :: weights(:)
    real(real64) :: global_mean, longitude, seconds_of_day
    integer :: i, j, calendar_year, calendar_month, calendar_day

    if (abs(solar_day - 86400.0_real64) > 1.0e-12_real64 .or. days_per_month /= 30 .or. &
        months_per_year /= 12 .or. days_per_year /= 360 .or. &
        abs(orbital_period - real(days_per_year, real64)*solar_day) > 1.0e-12_real64 .or. &
        abs(planetary_rotation_rate - 2.0_real64*acos(-1.0_real64)* &
          (1.0_real64/solar_day + 1.0_real64/(360.0_real64*solar_day))) > 1.0e-18_real64) then
      error stop 'radiation 360-day calendar constants are incorrect'
    end if
    call radiation_calendar_date(0.0_real64, calendar_year, calendar_month, calendar_day, seconds_of_day)
    if (calendar_year /= 1 .or. calendar_month /= 4 .or. calendar_day /= 1 .or. &
        abs(seconds_of_day) > 1.0e-12_real64) error stop 'radiation calendar does not start on year 1 April 1'
    call radiation_calendar_date(30.0_real64*solar_day, calendar_year, calendar_month, &
                                 calendar_day, seconds_of_day)
    if (calendar_year /= 1 .or. calendar_month /= 5 .or. calendar_day /= 1 .or. &
        abs(seconds_of_day) > 1.0e-12_real64) error stop 'radiation 30-day month boundary is incorrect'
    call radiation_calendar_date(270.0_real64*solar_day, calendar_year, calendar_month, &
                                 calendar_day, seconds_of_day)
    if (calendar_year /= 2 .or. calendar_month /= 1 .or. calendar_day /= 1) then
      error stop 'radiation 30-day month rollover is incorrect'
    end if
    call radiation_calendar_date(360.0_real64*solar_day, calendar_year, calendar_month, &
                                 calendar_day, seconds_of_day)
    if (calendar_year /= 2 .or. calendar_month /= 4 .or. calendar_day /= 1) then
      error stop 'radiation 360-day year rollover is incorrect'
    end if

    if (abs(shortwave_downward_flux(0.0_real64, 0.0_real64, 0.0_real64) - solar_constant) > 1.0e-12_real64 .or. &
        abs(shortwave_downward_flux(0.0_real64, acos(-1.0_real64), 0.0_real64)) > 1.0e-12_real64 .or. &
        abs(shortwave_downward_flux(0.0_real64, 0.0_real64, orbital_period) - solar_constant) > 1.0e-12_real64) then
      error stop 'radiation equinox day/night geometry is incorrect'
    end if
    if (abs(shortwave_downward_flux(1.0_real64, 0.0_real64, 0.25_real64*orbital_period) - &
            solar_constant*sin(axial_tilt)) > 1.0e-10_real64) then
      error stop 'radiation solstice declination is incorrect'
    end if

    call transform%init(31)
    nlon = transform%get_nlon()
    weights = transform%get_gaussian_weights()
    global_mean = 0.0_real64
    do j = 1, size(nlon)
      do i = 1, nlon(j)
        longitude = 2.0_real64*acos(-1.0_real64)*real(i - 1, real64)/real(nlon(j), real64)
        global_mean = global_mean + 0.5_real64*weights(j)/real(nlon(j), real64)* &
          shortwave_downward_flux(transform%mu(j), longitude, 0.137_real64*orbital_period)
      end do
    end do
    if (abs(global_mean - solar_constant/4.0_real64) > 5.0_real64) then
      error stop 'discrete global shortwave input is inconsistent with S0/4'
    end if
  end subroutine check_solar_geometry

  subroutine check_ozone_absorption()
    real(real64) :: upper_optical_depth, lower_optical_depth, expected_optical_depth
    real(real64) :: upper_longwave_optical_depth, lower_longwave_optical_depth
    real(real64) :: pressure_lower, pressure_upper, x_lower, x_upper, x_profile_lower, x_profile_upper

    upper_optical_depth = ozone_layer_optical_depth(ozone_pressure_lower_bound, ozone_peak_pressure)
    lower_optical_depth = ozone_layer_optical_depth(ozone_peak_pressure, ozone_pressure_upper_bound)
    upper_longwave_optical_depth = &
      ozone_longwave_layer_optical_depth(ozone_pressure_lower_bound, ozone_peak_pressure)
    lower_longwave_optical_depth = &
      ozone_longwave_layer_optical_depth(ozone_peak_pressure, ozone_pressure_upper_bound)
    if (abs(upper_optical_depth - 0.5_real64*ozone_shortwave_optical_depth) > 1.0e-15_real64 .or. &
        abs(lower_optical_depth - 0.5_real64*ozone_shortwave_optical_depth) > 1.0e-15_real64 .or. &
        abs(upper_longwave_optical_depth - 0.5_real64*ozone_longwave_optical_depth) > 1.0e-15_real64 .or. &
        abs(lower_longwave_optical_depth - 0.5_real64*ozone_longwave_optical_depth) > 1.0e-15_real64 .or. &
        abs(ozone_layer_optical_depth(1.0_real64, ozone_pressure_lower_bound)) > 1.0e-15_real64 .or. &
        abs(ozone_layer_optical_depth(ozone_pressure_upper_bound, 1.0e5_real64)) > 1.0e-15_real64 .or. &
        abs(ozone_longwave_layer_optical_depth(1.0_real64, ozone_pressure_lower_bound)) > 1.0e-15_real64 .or. &
        abs(ozone_longwave_layer_optical_depth(ozone_pressure_upper_bound, 1.0e5_real64)) > 1.0e-15_real64) then
      error stop 'ozone optical-depth profile has incorrect bounds or normalization'
    end if

    pressure_lower = 300.0_real64
    pressure_upper = 1000.0_real64
    x_lower = log(pressure_lower/ozone_peak_pressure)/(sqrt(2.0_real64)*ozone_log_pressure_width)
    x_upper = log(pressure_upper/ozone_peak_pressure)/(sqrt(2.0_real64)*ozone_log_pressure_width)
    x_profile_lower = log(ozone_pressure_lower_bound/ozone_peak_pressure)/ &
      (sqrt(2.0_real64)*ozone_log_pressure_width)
    x_profile_upper = log(ozone_pressure_upper_bound/ozone_peak_pressure)/ &
      (sqrt(2.0_real64)*ozone_log_pressure_width)
    expected_optical_depth = ozone_shortwave_optical_depth*(erf(x_upper) - erf(x_lower))/ &
      (erf(x_profile_upper) - erf(x_profile_lower))
    if (abs(ozone_layer_optical_depth(pressure_lower, pressure_upper) - expected_optical_depth) > 1.0e-15_real64) then
      error stop 'ozone optical depth does not use the documented error-function integral'
    end if
  end subroutine check_ozone_absorption

  subroutine check_radiation_column()
    real(real64), parameter :: pressure_half(0:2) = [1000.0_real64, 40000.0_real64, 100000.0_real64]
    ! Every state-dependent term receives the same RAW-filtered previous-time column.
    real(real64), parameter :: temperature(2) = [252.0_real64, 279.0_real64]
    real(real64), parameter :: surface_temperature = 291.0_real64
    real(real64), parameter :: deep_temperature = 285.0_real64
    real(real64) :: temperature_tendency(2), expected_temperature_tendency(2), surface_tendency, deep_tendency
    real(real64) :: transmission(2), emission(2), upward_longwave(0:2), downward_longwave(0:2), net_longwave(0:2)
    real(real64) :: surface_shortwave, absorbed_shortwave, sensible_heat, pressure_thickness
    real(real64) :: total_energy_tendency
    real(real64) :: incoming_shortwave, reflected_shortwave, outgoing_longwave
    integer :: k

    call radiation_tendency(pressure_half, temperature, surface_temperature, deep_temperature, &
                            3.0_real64, 4.0_real64, 0.0_real64, 0.0_real64, 0.0_real64, &
                            temperature_tendency, surface_tendency, deep_tendency, &
                            incoming_shortwave, reflected_shortwave, outgoing_longwave)
    absorbed_shortwave = ultraviolet_shortwave_fraction*solar_constant*(1.0_real64 - &
      exp(-ozone_layer_optical_depth(pressure_half(0), pressure_half(1))))
    surface_shortwave = solar_constant - absorbed_shortwave
    do k = 1, 2
      pressure_thickness = pressure_half(k) - pressure_half(k - 1)
      transmission(k) = exp(-(longwave_surface_optical_depth*pressure_thickness/ &
        (pressure_half(2) - pressure_half(0)) + &
        ozone_longwave_layer_optical_depth(pressure_half(k - 1), pressure_half(k))))
      emission(k) = (1.0_real64 - transmission(k))*stefan_boltzmann_constant*temperature(k)**4
    end do
    downward_longwave(0) = 0.0_real64
    do k = 1, 2
      downward_longwave(k) = transmission(k)*downward_longwave(k - 1) + emission(k)
    end do
    upward_longwave(2) = stefan_boltzmann_constant*surface_temperature**4
    do k = 2, 1, -1
      upward_longwave(k - 1) = transmission(k)*upward_longwave(k) + emission(k)
    end do
    net_longwave = upward_longwave - downward_longwave
    do k = 1, 2
      pressure_thickness = pressure_half(k) - pressure_half(k - 1)
      expected_temperature_tendency(k) = dry_gravity_acceleration/ &
        (dry_air_specific_heat*pressure_thickness)*(net_longwave(k) - net_longwave(k - 1))
    end do
    expected_temperature_tendency(1) = expected_temperature_tendency(1) + &
      dry_gravity_acceleration/(dry_air_specific_heat*(pressure_half(1) - pressure_half(0)))*absorbed_shortwave
    sensible_heat = pressure_half(2)/(dry_air_gas_constant*temperature(2))* &
      dry_air_specific_heat*surface_exchange_coefficient*sqrt(3.0_real64**2 + 4.0_real64**2 + &
      gustiness_speed**2)*(surface_temperature - temperature(2))
    expected_temperature_tendency(2) = expected_temperature_tendency(2) + &
      dry_gravity_acceleration/(dry_air_specific_heat*(pressure_half(2) - pressure_half(1)))*sensible_heat
    if (maxval(abs(temperature_tendency - expected_temperature_tendency)) > 1.0e-16_real64 .or. &
        abs(deep_tendency - ground_exchange_coefficient*(surface_temperature - deep_temperature)/ &
          deep_ground_heat_capacity) > 1.0e-18_real64) then
      error stop 'radiation column tendencies are incorrect'
    end if
    if (abs(incoming_shortwave - solar_constant) > 1.0e-12_real64 .or. &
        abs(reflected_shortwave - surface_shortwave_albedo*surface_shortwave) > 1.0e-12_real64) then
      error stop 'radiation column shortwave fluxes are incorrect'
    end if
    total_energy_tendency = sum(dry_air_specific_heat* &
      (pressure_half(1:2) - pressure_half(0:1))/dry_gravity_acceleration*temperature_tendency) + &
      surface_heat_capacity*surface_tendency + deep_ground_heat_capacity*deep_tendency
    if (abs(total_energy_tendency - &
        (incoming_shortwave - reflected_shortwave - upward_longwave(0))) > 1.0e-10_real64 .or. &
        abs(outgoing_longwave - upward_longwave(0)) > 1.0e-12_real64) then
      error stop 'radiation column does not conserve energy'
    end if

    call radiation_tendency(pressure_half, temperature, surface_temperature, deep_temperature, &
                            3.0_real64, 4.0_real64, 0.0_real64, acos(0.5_real64), 0.0_real64, &
                            temperature_tendency, surface_tendency, deep_tendency, &
                            incoming_shortwave, reflected_shortwave, outgoing_longwave)
    surface_shortwave = 0.5_real64*(solar_constant - absorbed_shortwave)
    if (abs(incoming_shortwave - 0.5_real64*solar_constant) > 1.0e-12_real64 .or. &
        abs(reflected_shortwave - surface_shortwave_albedo*surface_shortwave) > 1.0e-12_real64) then
      error stop 'radiation column does not use the prescribed vertical ozone path'
    end if
  end subroutine check_radiation_column

  subroutine check_radiation_daily_accumulator()
    type(radiation_daily_accumulator) :: accumulator
    type(radiation_diagnostics) :: sample, means

    sample%time_seconds = 86400.0_real64
    sample%mean_atmospheric_temperature = 250.0_real64
    sample%mean_surface_temperature = 280.0_real64
    sample%mean_deep_temperature = 281.0_real64
    sample%mean_kinetic_energy = 100.0_real64
    sample%mean_surface_pressure = 1.0e5_real64
    sample%mean_incoming_shortwave = 300.0_real64
    sample%mean_reflected_shortwave = 90.0_real64
    sample%mean_outgoing_longwave = 240.0_real64
    call accumulator%add(sample)
    sample%time_seconds = 87600.0_real64
    sample%mean_atmospheric_temperature = 254.0_real64
    sample%mean_surface_temperature = 290.0_real64
    sample%mean_deep_temperature = 283.0_real64
    sample%mean_kinetic_energy = 120.0_real64
    sample%mean_surface_pressure = 1.0002e5_real64
    sample%mean_incoming_shortwave = 400.0_real64
    sample%mean_reflected_shortwave = 120.0_real64
    sample%mean_outgoing_longwave = 250.0_real64
    call accumulator%add(sample)
    call accumulator%take(means)
    if (abs(means%time_seconds - 86400.0_real64) > 0.0_real64 .or. &
        abs(means%mean_atmospheric_temperature - 252.0_real64) > 1.0e-12_real64 .or. &
        abs(means%mean_surface_temperature - 285.0_real64) > 1.0e-12_real64 .or. &
        abs(means%mean_deep_temperature - 282.0_real64) > 1.0e-12_real64 .or. &
        abs(means%mean_kinetic_energy - 110.0_real64) > 1.0e-12_real64 .or. &
        abs(means%mean_surface_pressure - 1.0001e5_real64) > 1.0e-9_real64 .or. &
        abs(means%mean_incoming_shortwave - 350.0_real64) > 1.0e-12_real64 .or. &
        abs(means%mean_reflected_shortwave - 105.0_real64) > 1.0e-12_real64 .or. &
        abs(means%mean_outgoing_longwave - 245.0_real64) > 1.0e-12_real64) then
      error stop 'radiation daily accumulator does not return equal-weight interval means'
    end if
    if (accumulator%count() /= 0) error stop 'radiation daily accumulator was not reset after take'
  end subroutine check_radiation_daily_accumulator

  subroutine check_radiation_state()
    type(dry_atmosphere_solver) :: solver
    real(real64), allocatable :: zeta(:, :, :), delta(:, :, :), temperature(:, :, :)
    real(real64), allocatable :: surface_pressure(:, :), u(:, :, :), v(:, :, :)
    real(real64), allocatable :: surface_temperature(:, :), deep_temperature(:, :)
    real(real64), allocatable :: monthly_surface_temperature(:, :), monthly_surface_pressure(:, :)
    real(real64), allocatable :: zonal_temperature(:, :), zonal_u(:, :), zonal_v(:, :), eddy_uv(:, :), eddy_vt(:, :)
    complex(real64), allocatable :: snapshot_zeta_spectral(:, :, :), snapshot_delta_spectral(:, :, :)
    complex(real64), allocatable :: snapshot_temperature_spectral(:, :, :), snapshot_log_ps_spectral(:, :)
    integer, allocatable :: nlon(:)
    type(harmonic_transform) :: transform
    integer :: j, step
    real(real64) :: maximum_initial_difference
    real(real64) :: diagnostic_time, mean_atmospheric_temperature, mean_surface_temperature
    real(real64) :: mean_deep_temperature, mean_kinetic_energy, mean_surface_pressure
    real(real64) :: mean_incoming_shortwave, mean_reflected_shortwave, mean_outgoing_longwave

    call transform%init(truncation)
    nlon = transform%get_nlon()
    call solver%init(truncation, time_step)
    call solver%set_radiation_state()
    call solver%get_fields(zeta, delta, temperature, surface_pressure, u, v, &
                           surface_temperature=surface_temperature, deep_temperature=deep_temperature)
    maximum_initial_difference = 0.0_real64
    do j = 1, size(nlon)
      maximum_initial_difference = max(maximum_initial_difference, &
        maxval(abs(surface_temperature(1:nlon(j), j) - temperature(1:nlon(j), j, size(temperature, 3)))), &
        maxval(abs(deep_temperature(1:nlon(j), j) - temperature(1:nlon(j), j, size(temperature, 3)))))
    end do
    if (maximum_initial_difference > 1.0e-10_real64) then
      error stop 'radiation ground temperatures do not match the initial lowest atmospheric level'
    end if

    do step = 1, 4
      call solver%advance()
    end do
    ! Four advances collect samples at t = 0, 1200, 2400 and 3600 s; the mean is stamped with t = 0.
    call solver%take_radiation_daily_means(diagnostic_time, mean_atmospheric_temperature, &
      mean_surface_temperature, mean_deep_temperature, mean_kinetic_energy, mean_surface_pressure, &
      mean_incoming_shortwave, mean_reflected_shortwave, mean_outgoing_longwave)
    if (abs(diagnostic_time) > 1.0e-12_real64 .or. &
        abs(mean_incoming_shortwave - solar_constant/4.0_real64) > 5.0_real64 .or. &
        abs(mean_reflected_shortwave - surface_shortwave_albedo* &
          (1.0_real64 - ultraviolet_shortwave_fraction* &
          (1.0_real64 - exp(-ozone_shortwave_optical_depth)))* &
          mean_incoming_shortwave) > 1.0e-10_real64 .or. &
        min(mean_atmospheric_temperature, mean_surface_temperature, mean_deep_temperature, &
            mean_surface_pressure, mean_outgoing_longwave) <= 0.0_real64 .or. mean_kinetic_energy < 0.0_real64) then
      error stop 'radiation daily diagnostics are incorrect'
    end if
    call solver%take_radiation_monthly_means(monthly_surface_temperature, monthly_surface_pressure, &
                                             zonal_temperature, zonal_u, zonal_v, eddy_uv, eddy_vt)
    call solver%get_spectral_state(snapshot_zeta_spectral, snapshot_delta_spectral, &
                                   snapshot_temperature_spectral, snapshot_log_ps_spectral)
    if (.not. all(ieee_is_finite(real(snapshot_temperature_spectral, real64))) .or. &
        .not. all(ieee_is_finite(real(snapshot_log_ps_spectral, real64)))) then
      error stop 'radiation instantaneous snapshot produced a non-finite value'
    end if
    call solver%get_fields(zeta, delta, temperature, surface_pressure, u, v, &
                           surface_temperature=surface_temperature, deep_temperature=deep_temperature)
    if (.not. all(ieee_is_finite(temperature)) .or. .not. all(ieee_is_finite(surface_temperature)) .or. &
        .not. all(ieee_is_finite(deep_temperature))) then
      error stop 'radiation state produced a non-finite temperature'
    end if
    if (.not. all(ieee_is_finite(monthly_surface_temperature)) .or. &
        .not. all(ieee_is_finite(monthly_surface_pressure)) .or. &
        .not. all(ieee_is_finite(zonal_temperature)) .or. .not. all(ieee_is_finite(eddy_uv)) .or. &
        .not. all(ieee_is_finite(eddy_vt))) then
      error stop 'radiation monthly diagnostics produced a non-finite value'
    end if
  end subroutine check_radiation_state

  subroutine check_reference_atmosphere()
    type(hybrid_sigma_coordinate) :: coordinate
    real(real64), parameter :: expected_a(0:12) = [ &
      100.0_real64, 300.0_real64, 1000.0_real64, 5000.0_real64, &
      10000.0_real64, 8000.0_real64, 8000.0_real64, 10000.0_real64, &
      12000.0_real64, 10000.0_real64, 7000.0_real64, 3000.0_real64, &
      0.0_real64]
    real(real64), parameter :: expected_b(0:12) = [ &
      0.0_real64, 0.0_real64, 0.0_real64, 0.0_real64, 0.0_real64, &
      0.1_real64, 0.2_real64, 0.3_real64, 0.4_real64, 0.55_real64, &
      0.7_real64, 0.85_real64, 1.0_real64]
    real(real64), parameter :: expected_eta(12) = [ &
      0.002_real64, 0.0065_real64, 0.03_real64, 0.075_real64, &
      0.14_real64, 0.23_real64, 0.34_real64, 0.46_real64, &
      0.585_real64, 0.71_real64, 0.825_real64, 0.94_real64]

    call coordinate%init_default()
    if (coordinate%number_of_levels /= 12) error stop 'default dry atmosphere does not have twelve levels'
    if (maxval(abs(coordinate%a_half - expected_a)) > 1.0e-12_real64 .or. &
        maxval(abs(coordinate%b_half - expected_b)) > 1.0e-14_real64) then
      error stop 'dry atmosphere has the wrong default hybrid-sigma coefficients'
    end if
    if (abs(coordinate%reference_p_half(0) - 100.0_real64) > 1.0e-12_real64) then
      error stop 'dry atmosphere has the wrong top pressure'
    end if
    if (abs(coordinate%reference_p_half(12) - reference_surface_pressure) > 1.0e-12_real64) then
      error stop 'dry atmosphere has the wrong reference surface pressure'
    end if
    if (any(coordinate%reference_delta_p <= 0.0_real64) .or. &
        any(coordinate%reference_temperature <= 0.0_real64)) then
      error stop 'dry reference atmosphere contains a nonphysical value'
    end if
    if (maxval(abs(coordinate%full_level_eta - expected_eta)) > 1.0e-14_real64) then
      error stop 'dry full-level eta values do not match the hybrid-sigma pressures'
    end if
  end subroutine check_reference_atmosphere

  subroutine check_precomputed_gravity_wave_inverse()
    type(hybrid_sigma_coordinate) :: coordinate
    type(dry_gravity_wave_solver) :: gravity_wave
    complex(real64), allocatable :: surface_pressure(:, :), delta(:, :, :), temperature(:, :, :)
    complex(real64), allocatable :: rhs_surface_pressure(:, :), rhs_delta(:, :, :), rhs_temperature(:, :, :)
    complex(real64), allocatable :: next_surface_pressure(:, :), next_delta(:, :, :), next_temperature(:, :, :)
    complex(real64), allocatable :: g_next_surface_pressure(:, :), g_next_delta(:, :, :), g_next_temperature(:, :, :)
    complex(real64), allocatable :: residual_surface_pressure(:, :), residual_delta(:, :, :)
    complex(real64), allocatable :: residual_temperature(:, :, :)
    real(real64) :: centered_interval

    call coordinate%init_default()
    call gravity_wave%init(truncation, time_step, coordinate)
    call allocate_zero_state(coordinate%number_of_levels, surface_pressure, delta, temperature)
    delta(3, 2, 4) = cmplx(2.0e-5_real64, -1.0e-5_real64, real64)
    temperature(3, 2, 7) = cmplx(0.7_real64, 0.2_real64, real64)
    surface_pressure(3, 2) = cmplx(1.0e-4_real64, -2.0e-4_real64, real64)
    call gravity_wave%apply(surface_pressure, delta, temperature, &
                            rhs_surface_pressure, rhs_delta, rhs_temperature)
    centered_interval = 2.0_real64*time_step
    call gravity_wave%solve(centered_interval, surface_pressure, delta, temperature, &
                            surface_pressure, delta, temperature, &
                            rhs_surface_pressure, rhs_delta, rhs_temperature, &
                            next_surface_pressure, next_delta, next_temperature)
    call gravity_wave%apply(next_surface_pressure, next_delta, next_temperature, &
                            g_next_surface_pressure, g_next_delta, g_next_temperature)
    residual_surface_pressure = next_surface_pressure - &
      centered_interval*dry_gravity_wave_implicitness*g_next_surface_pressure - &
      (surface_pressure + centered_interval*(1.0_real64 - dry_gravity_wave_implicitness)*rhs_surface_pressure)
    residual_delta = next_delta - centered_interval*dry_gravity_wave_implicitness*g_next_delta - &
      (delta + centered_interval*(1.0_real64 - dry_gravity_wave_implicitness)*rhs_delta)
    residual_temperature = next_temperature - centered_interval*dry_gravity_wave_implicitness*g_next_temperature - &
      (temperature + centered_interval*(1.0_real64 - dry_gravity_wave_implicitness)*rhs_temperature)
    if (maxval(abs(residual_surface_pressure)) > 1.0e-10_real64 .or. &
        maxval(abs(residual_delta)) > 1.0e-10_real64 .or. &
        maxval(abs(residual_temperature)) > 1.0e-10_real64) then
      error stop 'precomputed dry gravity-wave inverse does not solve its implicit system'
    end if
  end subroutine check_precomputed_gravity_wave_inverse

  subroutine check_resting_atmosphere()
    type(hybrid_sigma_coordinate) :: coordinate
    type(harmonic_transform) :: transform
    type(dry_atmosphere_solver) :: solver
    complex(real64), allocatable :: zeta(:, :, :), delta(:, :, :), temperature(:, :, :), log_ps(:, :)
    complex(real64), allocatable :: initial_temperature(:, :, :), initial_log_ps(:, :)
    complex(real64), allocatable :: level_spectral(:, :)
    real(real64), allocatable :: grid(:, :)
    integer :: k

    call coordinate%init_default()
    call transform%init(truncation)
    call allocate_zero_state(coordinate%number_of_levels, log_ps, delta, temperature)
    allocate (zeta, mold=delta)
    zeta = 0.0_real64
    call transform%allocate_field(grid)
    do k = 1, coordinate%number_of_levels
      grid = coordinate%reference_temperature(k)
      call transform%grid_to_spectral(grid, level_spectral)
      temperature(:, :, k) = level_spectral
    end do
    grid = log(reference_surface_pressure)
    call transform%grid_to_spectral(grid, log_ps)
    initial_temperature = temperature
    initial_log_ps = log_ps
    call solver%init(truncation, time_step)
    call solver%set_initial_state(zeta, delta, temperature, log_ps)
    call solver%advance()
    call solver%advance()
    call solver%get_spectral_state(zeta, delta, temperature, log_ps)
    if (maxval(abs(zeta)) > 1.0e-11_real64 .or. maxval(abs(delta)) > 1.0e-11_real64) then
      error stop 'resting dry reference atmosphere developed wind'
    end if
    if (maxval(abs(temperature - initial_temperature)) > 1.0e-8_real64 .or. &
        maxval(abs(log_ps - initial_log_ps)) > 1.0e-10_real64) then
      error stop 'resting dry reference atmosphere changed'
    end if
  end subroutine check_resting_atmosphere

  subroutine check_held_suarez_forcing()
    real(real64) :: forcing_u, forcing_v, temperature_tendency
    real(real64) :: sinphi, pressure, surface_pressure, temperature, u, v
    real(real64) :: pressure_ratio, cosphi_squared, sigma_weight, equilibrium_temperature
    real(real64) :: thermal_rate

    call held_suarez_forcing(0.0_real64, reference_surface_pressure, reference_surface_pressure, &
                             323.0_real64, 12.0_real64, -6.0_real64, &
                             forcing_u, forcing_v, temperature_tendency)
    if (abs(forcing_u + 12.0_real64*held_suarez_friction_rate) > 1.0e-15_real64 .or. &
        abs(forcing_v - 6.0_real64*held_suarez_friction_rate) > 1.0e-15_real64 .or. &
        abs(temperature_tendency + 8.0_real64*held_suarez_lower_thermal_rate) > 1.0e-15_real64) then
      error stop 'Held-Suarez lower-boundary forcing is incorrect'
    end if

    sinphi = 0.5_real64
    pressure = 0.5_real64*reference_surface_pressure
    surface_pressure = reference_surface_pressure
    temperature = 250.0_real64
    u = 12.0_real64
    v = -6.0_real64
    call held_suarez_forcing(sinphi, pressure, surface_pressure, temperature, u, v, &
                             forcing_u, forcing_v, temperature_tendency)
    pressure_ratio = pressure/reference_surface_pressure
    cosphi_squared = 1.0_real64 - sinphi**2
    sigma_weight = max(0.0_real64, (pressure/surface_pressure - held_suarez_sigma_boundary)/ &
                       (1.0_real64 - held_suarez_sigma_boundary))
    equilibrium_temperature = max(held_suarez_minimum_equilibrium_temperature, &
      (held_suarez_equatorial_temperature - held_suarez_equator_to_pole_difference*sinphi**2 - &
       held_suarez_vertical_difference*log(pressure_ratio)*cosphi_squared)*pressure_ratio**dry_air_kappa)
    thermal_rate = held_suarez_upper_thermal_rate + &
      (held_suarez_lower_thermal_rate - held_suarez_upper_thermal_rate)*sigma_weight*cosphi_squared**2
    if (abs(forcing_u) > 1.0e-15_real64 .or. abs(forcing_v) > 1.0e-15_real64 .or. &
        abs(temperature_tendency + thermal_rate*(temperature - equilibrium_temperature)) > 1.0e-15_real64) then
      error stop 'Held-Suarez free-atmosphere forcing is incorrect'
    end if
  end subroutine check_held_suarez_forcing

  subroutine check_held_suarez_state()
    type(harmonic_transform) :: transform
    type(dry_atmosphere_solver) :: forced_solver, unforced_solver
    complex(real64), allocatable :: zeta_spectral(:, :, :), delta_spectral(:, :, :)
    complex(real64), allocatable :: temperature_spectral(:, :, :), log_ps_spectral(:, :)
    complex(real64), allocatable :: unforced_temperature(:, :, :)
    real(real64), allocatable :: zeta(:, :, :), delta(:, :, :), temperature(:, :, :)
    real(real64), allocatable :: surface_pressure(:, :), u(:, :, :), v(:, :, :)
    integer, allocatable :: nlon(:)
    integer :: i, j, k
    real(real64) :: longitude, cosphi, expected_temperature

    call transform%init(truncation)
    nlon = transform%get_nlon()
    call forced_solver%init(truncation, time_step)
    call forced_solver%set_held_suarez_state()
    call forced_solver%get_fields(zeta, delta, temperature, surface_pressure, u, v)
    do k = 1, size(temperature, 3)
      do j = 1, size(nlon)
        cosphi = sqrt(max(0.0_real64, 1.0_real64 - transform%mu(j)**2))
        do i = 1, nlon(j)
          longitude = 2.0_real64*acos(-1.0_real64)*real(i - 1, real64)/real(nlon(j), real64)
          expected_temperature = held_suarez_initial_temperature + &
            held_suarez_temperature_perturbation*cosphi*sin(longitude)
          if (abs(temperature(i, j, k) - expected_temperature) > 1.0e-10_real64) then
            error stop 'Held-Suarez initial temperature is incorrect'
          end if
          if (abs(zeta(i, j, k)) > 1.0e-12_real64 .or. abs(delta(i, j, k)) > 1.0e-12_real64 .or. &
              abs(u(i, j, k)) > 1.0e-10_real64 .or. abs(v(i, j, k)) > 1.0e-10_real64 .or. &
              abs(surface_pressure(i, j) - reference_surface_pressure) > 1.0e-8_real64) then
            error stop 'Held-Suarez atmosphere does not start at rest with uniform surface pressure'
          end if
        end do
      end do
    end do

    call forced_solver%get_spectral_state(zeta_spectral, delta_spectral, temperature_spectral, log_ps_spectral)
    call unforced_solver%init(truncation, time_step)
    call unforced_solver%set_initial_state(zeta_spectral, delta_spectral, temperature_spectral, log_ps_spectral)
    call forced_solver%advance()
    call unforced_solver%advance()
    call forced_solver%get_spectral_state(zeta_spectral, delta_spectral, temperature_spectral, log_ps_spectral)
    call unforced_solver%get_spectral_state(zeta_spectral, delta_spectral, unforced_temperature, log_ps_spectral)
    if (maxval(abs(temperature_spectral - unforced_temperature)) <= 1.0e-8_real64) then
      error stop 'Held-Suarez thermal forcing was not enabled by its initial-state setup'
    end if
    call forced_solver%get_fields(zeta, delta, temperature, surface_pressure, u, v)
    if (.not. all(ieee_is_finite(temperature)) .or. .not. all(ieee_is_finite(surface_pressure)) .or. &
        .not. all(ieee_is_finite(u)) .or. .not. all(ieee_is_finite(v))) then
      error stop 'Held-Suarez state produced a non-finite value'
    end if
  end subroutine check_held_suarez_state

  subroutine check_jablonowski_state()
    type(dry_atmosphere_solver) :: solver
    real(real64), allocatable :: zeta(:, :, :), delta(:, :, :), temperature(:, :, :)
    real(real64), allocatable :: surface_pressure(:, :), u(:, :, :), v(:, :, :)
    real(real64) :: cfl
    integer :: step

    call solver%init(truncation, time_step)
    call solver%set_jablonowski_williamson_state(.true.)
    do step = 1, 4
      call solver%advance()
    end do
    call solver%get_fields(zeta, delta, temperature, surface_pressure, u, v, cfl)
    if (.not. all(ieee_is_finite(zeta)) .or. .not. all(ieee_is_finite(delta)) .or. &
        .not. all(ieee_is_finite(temperature)) .or. .not. all(ieee_is_finite(surface_pressure)) .or. &
        .not. all(ieee_is_finite(u)) .or. .not. all(ieee_is_finite(v)) .or. &
        .not. ieee_is_finite(cfl)) then
      error stop 'Jablonowski-Williamson dry state produced a non-finite value'
    end if
    if (maxval(abs(u)) <= 1.0_real64) error stop 'Jablonowski-Williamson jet is absent'
  end subroutine check_jablonowski_state

  subroutine check_jablonowski_steady_state()
    ! Without the wind perturbation the Jablonowski-Williamson base state is a steady
    ! solution.  It stays balanced only when eta is consistent with the hybrid-sigma
    ! pressures and the surface geopotential is supplied; otherwise |v| grows to
    ! several m/s and ps oscillates by several hPa within a day.
    integer, parameter :: steady_truncation = 21
    integer, parameter :: steady_steps = nint(86400.0_real64/time_step)
    real(real64), parameter :: wind_tolerance = 1.5_real64
    real(real64), parameter :: surface_pressure_tolerance = 100.0_real64
    type(harmonic_transform) :: transform
    type(dry_atmosphere_solver) :: solver
    real(real64), allocatable :: zeta(:, :, :), delta(:, :, :), temperature(:, :, :)
    real(real64), allocatable :: surface_pressure(:, :), u(:, :, :), v(:, :, :), initial_u(:, :, :)
    integer, allocatable :: nlon(:)
    real(real64) :: maximum_v, maximum_u_change, maximum_surface_pressure_change
    integer :: step, j

    call transform%init(steady_truncation)
    nlon = transform%get_nlon()
    call solver%init(steady_truncation, time_step)
    call solver%set_jablonowski_williamson_state(.false.)
    call solver%get_fields(zeta, delta, temperature, surface_pressure, initial_u, v)
    do step = 1, steady_steps
      call solver%advance()
    end do
    call solver%get_fields(zeta, delta, temperature, surface_pressure, u, v)

    maximum_v = 0.0_real64
    maximum_u_change = 0.0_real64
    maximum_surface_pressure_change = 0.0_real64
    do j = 1, size(nlon)
      maximum_v = max(maximum_v, maxval(abs(v(1:nlon(j), j, :))))
      maximum_u_change = max(maximum_u_change, &
        maxval(abs(u(1:nlon(j), j, :) - initial_u(1:nlon(j), j, :))))
      maximum_surface_pressure_change = max(maximum_surface_pressure_change, &
        maxval(abs(surface_pressure(1:nlon(j), j) - reference_surface_pressure)))
    end do
    if (maximum_v > wind_tolerance .or. maximum_u_change > wind_tolerance .or. &
        maximum_surface_pressure_change > surface_pressure_tolerance) then
      write (*, '(a,3es12.4)') 'max |v|, max |du|, max |dps| = ', maximum_v, maximum_u_change, &
        maximum_surface_pressure_change
      error stop 'Jablonowski-Williamson base state is not steady'
    end if
  end subroutine check_jablonowski_steady_state

  subroutine check_flat_terrain_balanced_state()
    ! Over flat terrain the Jablonowski-Williamson temperature is steady only with the
    ! rebalanced zonal wind used by the radiation case.  With the original wind instead,
    ! |v| reaches about 20 m/s and ps changes by about 9 hPa within a day.
    integer, parameter :: steady_truncation = 21
    integer, parameter :: steady_steps = nint(86400.0_real64/time_step)
    real(real64), parameter :: wind_tolerance = 1.5_real64
    real(real64), parameter :: surface_pressure_tolerance = 100.0_real64
    type(harmonic_transform) :: transform
    type(hybrid_sigma_coordinate) :: coordinate
    type(dry_atmosphere_solver) :: solver
    complex(real64), allocatable :: zeta_spectral(:, :, :), delta_spectral(:, :, :)
    complex(real64), allocatable :: temperature_spectral(:, :, :), log_ps_spectral(:, :)
    complex(real64), allocatable :: unused_surface_geopotential(:, :)
    real(real64), allocatable :: zeta(:, :, :), delta(:, :, :), temperature(:, :, :)
    real(real64), allocatable :: surface_pressure(:, :), u(:, :, :), v(:, :, :), initial_u(:, :, :)
    integer, allocatable :: nlon(:)
    real(real64) :: maximum_v, maximum_u_change, maximum_surface_pressure_change
    integer :: step, j

    call transform%init(steady_truncation)
    nlon = transform%get_nlon()
    call coordinate%init_default()
    call solver%init(steady_truncation, time_step)
    call jablonowski_williamson_initial_state(transform, steady_truncation, coordinate, &
                                              .false., zeta_spectral, delta_spectral, &
                                              temperature_spectral, log_ps_spectral, &
                                              unused_surface_geopotential, flat_terrain=.true.)
    if (maxval(abs(unused_surface_geopotential)) > 0.0_real64) then
      error stop 'flat-terrain Jablonowski-Williamson state has nonzero surface geopotential'
    end if
    call solver%set_initial_state(zeta_spectral, delta_spectral, temperature_spectral, log_ps_spectral)
    call solver%get_fields(zeta, delta, temperature, surface_pressure, initial_u, v)
    do j = 1, size(nlon)
      if (maxval(abs(initial_u(1:nlon(j), j, size(initial_u, 3)))) > 3.5_real64 .or. &
          maxval(initial_u(1:nlon(j), j, 1)) > 24.0_real64) then
        error stop 'flat-terrain balanced zonal wind has an unexpected magnitude'
      end if
    end do
    do step = 1, steady_steps
      call solver%advance()
    end do
    call solver%get_fields(zeta, delta, temperature, surface_pressure, u, v)

    maximum_v = 0.0_real64
    maximum_u_change = 0.0_real64
    maximum_surface_pressure_change = 0.0_real64
    do j = 1, size(nlon)
      maximum_v = max(maximum_v, maxval(abs(v(1:nlon(j), j, :))))
      maximum_u_change = max(maximum_u_change, &
        maxval(abs(u(1:nlon(j), j, :) - initial_u(1:nlon(j), j, :))))
      maximum_surface_pressure_change = max(maximum_surface_pressure_change, &
        maxval(abs(surface_pressure(1:nlon(j), j) - reference_surface_pressure)))
    end do
    if (maximum_v > wind_tolerance .or. maximum_u_change > wind_tolerance .or. &
        maximum_surface_pressure_change > surface_pressure_tolerance) then
      write (*, '(a,3es12.4)') 'max |v|, max |du|, max |dps| = ', maximum_v, maximum_u_change, &
        maximum_surface_pressure_change
      error stop 'flat-terrain balanced Jablonowski-Williamson state is not steady'
    end if
  end subroutine check_flat_terrain_balanced_state

  subroutine allocate_zero_state(levels, surface_pressure, delta, temperature)
    integer, intent(in) :: levels
    complex(real64), allocatable, intent(out) :: surface_pressure(:, :), delta(:, :, :), temperature(:, :, :)

    allocate (surface_pressure(0:truncation + 1, 0:truncation))
    allocate (delta(0:truncation + 1, 0:truncation, levels))
    allocate (temperature(0:truncation + 1, 0:truncation, levels))
    surface_pressure = 0.0_real64
    delta = 0.0_real64
    temperature = 0.0_real64
  end subroutine allocate_zero_state

end program check_dry_atmosphere
