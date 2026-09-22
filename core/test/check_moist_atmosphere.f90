!> Unit checks of the moist physics (docs/cases/moist.md, "checks after implementation")
!> and of the moist dynamics wiring.
program check_moist_atmosphere
  use iso_fortran_env, only: real64
  use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
  use harmonics, only: harmonic_transform
  use planet_parameters, only: earth_gravity, planet_config
  use dry_vertical_coordinate, only: hybrid_sigma_coordinate, reference_surface_pressure, dry_air_gas_constant, &
                                     dry_air_kappa
  use moist_thermodynamics, only: dry_air_specific_heat, latent_heat_of_condensation, water_vapor_gas_constant, &
                                  gas_constant_ratio, virtual_temperature_coefficient, saturation_vapor_pressure, &
                                  saturation_specific_humidity, saturation_specific_humidity_derivative, &
                                  equivalent_potential_temperature, lifting_condensation_level, &
                                  full_level_pressures, virtual_temperature
  use surface_exchange, only: surface_sensible_heat_flux, surface_evaporation_flux, lowest_full_level_pressure
  use large_scale_condensation, only: large_scale_condensation_tendency
  use moist_convection, only: moist_convective_adjustment_tendency, moist_convection_reference_profile
  use dry_convection, only: dry_convective_adjustment_tendency
  use dry_radiation, only: radiation_tendency, radiation_diagnostics
  use cloud_diagnostics, only: relative_humidity_cloud_cover, convective_cloud_cover, diagnose_cloud_cover
  use dry_physics_config, only: dry_model_physics_config, radiation_config, convection_config, &
                                moist_convection_config, condensation_config, cloud_config, &
                                radiation_surface_heat_capacity
  use dry_case_initial_conditions, only: moist_case_physics, slab_ocean_case_physics, radiation_case_planet, &
                                         set_radiation_case_state
  use dry_held_suarez, only: held_suarez_initial_state
  use dry_atmosphere, only: dry_atmosphere_solver
  use dry_state, only: dry_state_type, dry_tendency_type, allocate_dry_state, allocate_dry_surface_fields, &
                       allocate_dry_tendency
  use dry_tendency_workspace, only: dry_workspace_type
  use dry_tendency_evaluator, only: evaluate_dry_tendency
  use radiation_diagnostics_collector, only: radiation_case_diagnostics, radiation_monthly_means
  implicit none

  integer, parameter :: truncation = 5
  real(real64), parameter :: time_step = 1200.0_real64
  real(real64), parameter :: pressure_half(0:4) = [ &
    1000.0_real64, 20000.0_real64, 50000.0_real64, 80000.0_real64, 100000.0_real64]
  type(radiation_config) :: surface
  type(convection_config) :: dry_adjustment
  type(moist_convection_config) :: moist_adjustment
  type(condensation_config) :: condensation

  call check_saturation_relations()
  call check_lifting_condensation_level()
  call check_neutral_column_sensible_heat()
  call check_evaporation_flux_and_surface_budget()
  call check_large_scale_condensation()
  call check_moist_dry_adjustment()
  call check_moist_convection_stable_columns()
  call check_moist_convection_deep()
  call check_moist_convection_shallow()
  call check_cloud_diagnosis()
  call check_cloud_shortwave_reflection()
  call check_moist_case_configuration()
  call check_dry_atmosphere_unchanged_by_zero_humidity()
  call check_column_water_budget_closure()
  call check_moist_integration()
  write (*, '(a)') 'check_moist_atmosphere: all checks passed'

contains

  subroutine check_saturation_relations()
    real(real64) :: temperature, pressure, humidity, derivative, finite_difference, step
    real(real64), parameter :: temperatures(3) = [250.0_real64, 280.0_real64, 305.0_real64]
    integer :: i

    if (abs(gas_constant_ratio - dry_air_gas_constant/water_vapor_gas_constant) > 1.0e-15_real64 .or. &
        abs(virtual_temperature_coefficient - (water_vapor_gas_constant/dry_air_gas_constant - 1.0_real64)) > &
        1.0e-14_real64 .or. abs(dry_air_specific_heat - 1004.0_real64) > 0.0_real64 .or. &
        abs(latent_heat_of_condensation - 2.5e6_real64) > 0.0_real64) then
      error stop 'moist constants differ from the documented values'
    end if
    if (abs(saturation_vapor_pressure(273.16_real64) - 610.78_real64) > 1.0e-9_real64) then
      error stop 'saturation vapour pressure does not pass through (T_0, e_0)'
    end if
    pressure = 1.0e5_real64
    do i = 1, 3
      temperature = temperatures(i)
      humidity = saturation_specific_humidity(temperature, pressure)
      if (abs(humidity - gas_constant_ratio*saturation_vapor_pressure(temperature)/ &
              (pressure - (1.0_real64 - gas_constant_ratio)*saturation_vapor_pressure(temperature))) > 1.0e-15_real64) then
        error stop 'saturation specific humidity formula is incorrect'
      end if
      step = 1.0e-3_real64
      finite_difference = (saturation_specific_humidity(temperature + step, pressure) - &
                           saturation_specific_humidity(temperature - step, pressure))/(2.0_real64*step)
      derivative = saturation_specific_humidity_derivative(temperature, pressure)
      if (abs(derivative - finite_difference) > 1.0e-6_real64*derivative) then
        error stop 'saturation specific humidity derivative does not match a finite difference'
      end if
    end do
    ! Above the boiling point the layer counts as unsaturated.
    if (saturation_specific_humidity(400.0_real64, 1.0e5_real64) /= 1.0_real64 .or. &
        saturation_specific_humidity_derivative(400.0_real64, 1.0e5_real64) /= 0.0_real64) then
      error stop 'saturation specific humidity is not capped where p <= e_s'
    end if
    if (abs(virtual_temperature(300.0_real64, 0.01_real64) - (1.0_real64 + 0.01_real64*virtual_temperature_coefficient)* &
            300.0_real64) > 1.0e-12_real64 .or. virtual_temperature(300.0_real64, -0.01_real64) /= 300.0_real64) then
      error stop 'virtual temperature does not clip negative humidity'
    end if
  end subroutine check_saturation_relations

  subroutine check_lifting_condensation_level()
    real(real64) :: parcel_temperature, parcel_humidity, parcel_pressure, level_temperature, level_pressure

    parcel_temperature = 300.0_real64
    parcel_pressure = 1.0e5_real64
    parcel_humidity = 0.7_real64*saturation_specific_humidity(parcel_temperature, parcel_pressure)
    call lifting_condensation_level(parcel_temperature, parcel_humidity, parcel_pressure, level_temperature, level_pressure)
    if (level_temperature >= parcel_temperature .or. level_pressure >= parcel_pressure .or. &
        abs(level_pressure - parcel_pressure*(level_temperature/parcel_temperature)**(1.0_real64/dry_air_kappa)) > &
        1.0e-6_real64 .or. &
        abs(saturation_specific_humidity(level_temperature, level_pressure) - parcel_humidity) > 1.0e-4_real64*parcel_humidity) &
      then
      error stop 'lifting condensation level does not saturate the dry-adiabatically lifted parcel'
    end if
    ! A saturated parcel is already at its condensation level.
    parcel_humidity = saturation_specific_humidity(parcel_temperature, parcel_pressure)
    call lifting_condensation_level(parcel_temperature, parcel_humidity, parcel_pressure, level_temperature, level_pressure)
    if (level_temperature /= parcel_temperature .or. level_pressure /= parcel_pressure) then
      error stop 'saturated parcel was lifted to a condensation level'
    end if
  end subroutine check_lifting_condensation_level

  !> A dry-adiabatically neutral column exchanges no sensible heat with the surface.
  subroutine check_neutral_column_sensible_heat()
    real(real64) :: lowest_pressure, lowest_temperature, surface_temperature, flux, warm_flux

    lowest_pressure = lowest_full_level_pressure(pressure_half)
    surface_temperature = 290.0_real64
    lowest_temperature = surface_temperature*(lowest_pressure/pressure_half(4))**dry_air_kappa
    flux = surface_sensible_heat_flux(surface, pressure_half, lowest_temperature, surface_temperature, &
                                      3.0_real64, 4.0_real64)
    if (abs(flux) > 1.0e-9_real64) error stop 'neutral column has a non-zero sensible heat flux'
    warm_flux = surface_sensible_heat_flux(surface, pressure_half, lowest_temperature, surface_temperature + 1.0_real64, &
                                           3.0_real64, 4.0_real64)
    if (abs(warm_flux - lowest_pressure/(dry_air_gas_constant*lowest_temperature)*surface%dry_air_specific_heat* &
            surface%surface_exchange_coefficient*sqrt(9.0_real64 + 16.0_real64 + surface%gustiness_speed**2)) > &
        1.0e-9_real64) then
      error stop 'sensible heat flux does not use the lowest-level density and bulk coefficient'
    end if
  end subroutine check_neutral_column_sensible_heat

  !> E = beta rho_N C_E |U| [q_s(T_s, p_s) - q_N^+], and the surface budget loses exactly L E.
  subroutine check_evaporation_flux_and_surface_budget()
    real(real64), parameter :: temperature(4) = [220.0_real64, 250.0_real64, 275.0_real64, 288.0_real64]
    real(real64) :: lowest_pressure, expected, flux, dew_flux
    real(real64) :: temperature_tendency(4), surface_tendency, deep_tendency, reference_surface_tendency
    real(real64) :: incoming, reflected, outgoing
    type(radiation_config) :: ocean
    type(dry_model_physics_config) :: ocean_physics

    lowest_pressure = lowest_full_level_pressure(pressure_half)
    expected = lowest_pressure/(dry_air_gas_constant*temperature(4))*surface%surface_exchange_coefficient* &
      sqrt(9.0_real64 + 16.0_real64 + surface%gustiness_speed**2)* &
      (saturation_specific_humidity(290.0_real64, pressure_half(4)) - 0.005_real64)
    flux = surface_evaporation_flux(surface, 1.0_real64, pressure_half, temperature(4), 0.005_real64, 290.0_real64, &
                                    3.0_real64, 4.0_real64)
    if (abs(flux - expected) > 1.0e-12_real64*abs(expected) .or. flux <= 0.0_real64) then
      error stop 'evaporation flux formula is incorrect'
    end if
    if (abs(surface_evaporation_flux(surface, 0.25_real64, pressure_half, temperature(4), 0.005_real64, 290.0_real64, &
                                     3.0_real64, 4.0_real64) - 0.25_real64*flux) > 1.0e-12_real64*abs(flux)) then
      error stop 'evaporation flux does not scale with the surface wetness'
    end if
    ! Negative humidity is clipped; a saturated lowest level over a cold surface gives dew (E < 0), not zero.
    if (surface_evaporation_flux(surface, 1.0_real64, pressure_half, temperature(4), -0.1_real64, 290.0_real64, &
                                 3.0_real64, 4.0_real64) /= &
        surface_evaporation_flux(surface, 1.0_real64, pressure_half, temperature(4), 0.0_real64, 290.0_real64, &
                                 3.0_real64, 4.0_real64)) error stop 'evaporation does not clip negative humidity'
    dew_flux = surface_evaporation_flux(surface, 1.0_real64, pressure_half, temperature(4), 0.02_real64, 280.0_real64, &
                                        3.0_real64, 4.0_real64)
    if (dew_flux >= 0.0_real64) error stop 'dew (negative evaporation) is cut off'

    ocean_physics = slab_ocean_case_physics()
    ocean = ocean_physics%radiation
    call radiation_tendency(ocean, pressure_half, temperature, 290.0_real64, 0.0_real64, 3.0_real64, 4.0_real64, &
                            0.0_real64, 0.0_real64, 0.0_real64, temperature_tendency, reference_surface_tendency, &
                            deep_tendency, incoming, reflected, outgoing)
    call radiation_tendency(ocean, pressure_half, temperature, 290.0_real64, 0.0_real64, 3.0_real64, 4.0_real64, &
                            0.0_real64, 0.0_real64, 0.0_real64, temperature_tendency, surface_tendency, &
                            deep_tendency, incoming, reflected, outgoing, latent_heat_flux=100.0_real64)
    if (abs((reference_surface_tendency - surface_tendency)*radiation_surface_heat_capacity(ocean) - 100.0_real64) > &
        1.0e-6_real64) then
      error stop 'latent heat flux is not taken from the surface energy budget'
    end if
  end subroutine check_evaporation_flux_and_surface_budget

  !> After condensation every layer is saturated to 1e-4 q_s and c_p T + L q is unchanged per layer.
  subroutine check_large_scale_condensation()
    real(real64) :: temperature(4), humidity(4), full_level_pressure(4), delta_pressure(4)
    real(real64) :: temperature_tendency(4), humidity_tendency(4), precipitation, interval
    real(real64) :: new_temperature(4), new_humidity(4), saturation(4), enthalpy_change(4)
    integer :: k

    call full_level_pressures(pressure_half, full_level_pressure, delta_pressure)
    temperature = [220.0_real64, 255.0_real64, 285.0_real64, 300.0_real64]
    do k = 1, 4
      saturation(k) = saturation_specific_humidity(temperature(k), full_level_pressure(k))
    end do
    ! Unsaturated, mildly, strongly and very strongly supersaturated layers.
    humidity = saturation*[0.5_real64, 1.1_real64, 1.5_real64, 2.0_real64]
    interval = 2.0_real64*time_step
    call large_scale_condensation_tendency(condensation, pressure_half, temperature, humidity, interval, &
                                           temperature_tendency, humidity_tendency, precipitation)
    new_temperature = temperature + interval*temperature_tendency
    new_humidity = humidity + interval*humidity_tendency
    enthalpy_change = dry_air_specific_heat*interval*temperature_tendency + latent_heat_of_condensation*interval*humidity_tendency
    if (abs(temperature_tendency(1)) > 0.0_real64 .or. abs(humidity_tendency(1)) > 0.0_real64) then
      error stop 'large-scale condensation acted on an unsaturated layer'
    end if
    do k = 2, 4
      saturation(k) = saturation_specific_humidity(new_temperature(k), full_level_pressure(k))
      if (abs(new_humidity(k) - saturation(k)) > 1.0e-4_real64*saturation(k)) then
        error stop 'large-scale condensation does not end at saturation'
      end if
      if (new_humidity(k) >= humidity(k) .or. new_temperature(k) <= temperature(k)) then
        error stop 'large-scale condensation has the wrong sign'
      end if
    end do
    if (maxval(abs(enthalpy_change)) > 1.0e-9_real64*dry_air_specific_heat) then
      error stop 'large-scale condensation does not conserve layer moist enthalpy'
    end if
    if (abs(precipitation + sum(humidity_tendency*delta_pressure)/earth_gravity) > 1.0e-12_real64*precipitation .or. &
        precipitation <= 0.0_real64) then
      error stop 'large-scale precipitation does not equal the column water removed'
    end if
  end subroutine check_large_scale_condensation

  !> The moist dry adjustment mixes humidity within blocks, conserves column
  !> water, judges stability by virtual potential temperature, and reduces to the
  !> dry procedure at q = 0.
  subroutine check_moist_dry_adjustment()
    real(real64) :: full_level_pressure(4), delta_pressure(4), exner(4), temperature(4), humidity(4)
    real(real64) :: dry_tendency(4), temperature_tendency(4), humidity_tendency(4), zero_humidity(4)
    real(real64), parameter :: unstable_potential_temperature(4) = [340.0_real64, 280.0_real64, 300.0_real64, 250.0_real64]
    integer :: k

    call full_level_pressures(pressure_half, full_level_pressure, delta_pressure)
    exner = (full_level_pressure/reference_surface_pressure)**dry_air_kappa
    temperature = exner*unstable_potential_temperature
    humidity = [0.0001_real64, 0.002_real64, 0.008_real64, 0.012_real64]
    zero_humidity = 0.0_real64
    call dry_convective_adjustment_tendency(dry_adjustment, pressure_half, temperature, dry_tendency)
    call dry_convective_adjustment_tendency(dry_adjustment, pressure_half, temperature, temperature_tendency, &
                                            zero_humidity, humidity_tendency)
    if (any(temperature_tendency /= dry_tendency) .or. any(humidity_tendency /= 0.0_real64)) then
      error stop 'moist dry adjustment with q = 0 differs from the dry adjustment'
    end if
    call dry_convective_adjustment_tendency(dry_adjustment, pressure_half, temperature, temperature_tendency, &
                                            humidity, humidity_tendency)
    if (abs(sum(humidity_tendency*delta_pressure)) > 1.0e-12_real64*maxval(abs(humidity_tendency))*sum(delta_pressure) .or. &
        abs(sum(temperature_tendency*delta_pressure)) > 1.0e-10_real64) then
      error stop 'moist dry adjustment does not conserve column water and enthalpy'
    end if
    ! Layers 2 and 3 form the unstable block: humidity is mixed to one value there and untouched elsewhere.
    if (abs(humidity_tendency(1)) > 0.0_real64 .or. abs(humidity_tendency(4)) > 0.0_real64 .or. &
        abs((humidity(2) + dry_adjustment%adjustment_time*humidity_tendency(2)) - &
            (humidity(3) + dry_adjustment%adjustment_time*humidity_tendency(3))) > 1.0e-12_real64) then
      error stop 'moist dry adjustment does not mix humidity uniformly within the block'
    end if
    ! A column neutral in theta but with more vapour below is unstable in theta_v and is adjusted.
    temperature = exner*300.0_real64
    call dry_convective_adjustment_tendency(dry_adjustment, pressure_half, temperature, temperature_tendency, &
                                            humidity, humidity_tendency)
    if (maxval(abs(humidity_tendency)) <= 0.0_real64) then
      error stop 'moist dry adjustment ignores the virtual temperature effect of humidity'
    end if
    do k = 1, 4
      if (.not. ieee_is_finite(temperature_tendency(k))) error stop 'moist dry adjustment produced a non-finite value'
    end do
  end subroutine check_moist_dry_adjustment

  !> No buoyant layer, or a reference profile colder than the column (P_T <= 0), gives zero tendencies.
  subroutine check_moist_convection_stable_columns()
    real(real64) :: full_level_pressure(4), delta_pressure(4), temperature(4), humidity(4)
    real(real64) :: temperature_tendency(4), humidity_tendency(4), precipitation
    real(real64) :: reference_temperature(4), parcel_humidity(4), reference_humidity(4)
    integer :: k

    call full_level_pressures(pressure_half, full_level_pressure, delta_pressure)
    ! Very stable (nearly isothermal) column: the lifted parcel is never buoyant.
    temperature = [230.0_real64, 270.0_real64, 285.0_real64, 288.0_real64]
    do k = 1, 4
      humidity(k) = 0.5_real64*saturation_specific_humidity(temperature(k), full_level_pressure(k))
    end do
    call moist_convective_adjustment_tendency(moist_adjustment, pressure_half, temperature, humidity, &
                                              temperature_tendency, humidity_tendency, precipitation)
    if (any(temperature_tendency /= 0.0_real64) .or. any(humidity_tendency /= 0.0_real64) .or. precipitation /= 0.0_real64) &
      then
      error stop 'moist convection acted on a column without buoyancy'
    end if
    ! Dry parcel (q_N = 0): the reference is the dry adiabat, buoyant nowhere in a dry-stable column.
    humidity = 0.0_real64
    call moist_convection_reference_profile(moist_adjustment, full_level_pressure, temperature, humidity, &
                                            reference_temperature, parcel_humidity, reference_humidity)
    do k = 1, 4
      if (abs(reference_temperature(k) - temperature(4)*(full_level_pressure(k)/full_level_pressure(4))**dry_air_kappa) > &
          1.0e-9_real64 .or. parcel_humidity(k) /= 0.0_real64) then
        error stop 'dry parcel reference profile is not the dry adiabat'
      end if
    end do
    ! Buoyant only in one thin upper layer of an otherwise much warmer column: P_T <= 0, no adjustment.
    temperature = [190.0_real64, 262.0_real64, 289.0_real64, 300.0_real64]
    do k = 1, 4
      humidity(k) = 0.5_real64*saturation_specific_humidity(temperature(k), full_level_pressure(k))
    end do
    humidity(4) = 0.9_real64*saturation_specific_humidity(temperature(4), full_level_pressure(4))
    call moist_convection_reference_profile(moist_adjustment, full_level_pressure, temperature, humidity, &
                                            reference_temperature, parcel_humidity, reference_humidity)
    call moist_convective_adjustment_tendency(moist_adjustment, pressure_half, temperature, humidity, &
                                              temperature_tendency, humidity_tendency, precipitation)
    if (reference_temperature(1) > temperature(1) .and. reference_temperature(2) < temperature(2)) then
      ! Level 1 is buoyant but the free convection level search starts at the surface and stops at
      ! the first buoyant level; here levels 3 and 2 are not buoyant, so k_f = 1 and the layer is 1..4.
      if (sum((reference_temperature - temperature)*delta_pressure) <= 0.0_real64) then
        if (any(temperature_tendency /= 0.0_real64) .or. any(humidity_tendency /= 0.0_real64)) then
          error stop 'moist convection acted on a column with P_T <= 0'
        end if
      end if
    end if
  end subroutine check_moist_convection_stable_columns

  !> Deep convection: column moist enthalpy is conserved and the precipitation is the water removed.
  subroutine check_moist_convection_deep()
    real(real64) :: full_level_pressure(4), delta_pressure(4), temperature(4), humidity(4)
    real(real64) :: temperature_tendency(4), humidity_tendency(4), precipitation, enthalpy_change
    real(real64) :: reference_temperature(4), parcel_humidity(4), reference_humidity(4)
    integer :: k

    call full_level_pressures(pressure_half, full_level_pressure, delta_pressure)
    ! Conditionally unstable moist tropical column: a warm, humid lowest layer under a cool troposphere.
    temperature = [185.0_real64, 260.0_real64, 289.0_real64, 302.0_real64]
    do k = 1, 3
      humidity(k) = 0.8_real64*saturation_specific_humidity(temperature(k), full_level_pressure(k))
    end do
    humidity(4) = 0.95_real64*saturation_specific_humidity(temperature(4), full_level_pressure(4))
    call moist_convection_reference_profile(moist_adjustment, full_level_pressure, temperature, humidity, &
                                            reference_temperature, parcel_humidity, reference_humidity)
    if (reference_temperature(3) <= temperature(3)) error stop 'deep-convection test column is not buoyant'
    do k = 1, 3
      if (parcel_humidity(k) /= saturation_specific_humidity(reference_temperature(k), full_level_pressure(k)) .and. &
          parcel_humidity(k) /= humidity(4)) error stop 'parcel humidity is neither saturated nor conserved'
      if (abs(reference_humidity(k) - moist_adjustment%reference_relative_humidity* &
              saturation_specific_humidity(reference_temperature(k), full_level_pressure(k))) > 1.0e-15_real64) then
        error stop 'reference humidity is not RH_c q_s(T_ref)'
      end if
    end do
    call moist_convective_adjustment_tendency(moist_adjustment, pressure_half, temperature, humidity, &
                                              temperature_tendency, humidity_tendency, precipitation)
    if (precipitation <= 0.0_real64) error stop 'deep convection produced no precipitation'
    enthalpy_change = sum((dry_air_specific_heat*temperature_tendency + latent_heat_of_condensation*humidity_tendency)* &
                          delta_pressure)
    if (abs(enthalpy_change) > 1.0e-9_real64*dry_air_specific_heat*maxval(abs(temperature_tendency))*sum(delta_pressure)) then
      error stop 'deep convection does not conserve column moist enthalpy'
    end if
    if (abs(precipitation + sum(humidity_tendency*delta_pressure)/earth_gravity) > 1.0e-12_real64*precipitation) then
      error stop 'convective precipitation does not equal the column water removed'
    end if
    if (sum(temperature_tendency*delta_pressure) <= 0.0_real64) error stop 'deep convection cools the column'
  end subroutine check_moist_convection_deep

  !> Shallow convection: no precipitation, and both column water and column enthalpy are unchanged.
  subroutine check_moist_convection_shallow()
    real(real64) :: full_level_pressure(4), delta_pressure(4), temperature(4), humidity(4)
    real(real64) :: temperature_tendency(4), humidity_tendency(4), precipitation
    real(real64) :: reference_temperature(4), parcel_humidity(4), reference_humidity(4)
    integer :: k

    call full_level_pressures(pressure_half, full_level_pressure, delta_pressure)
    ! Warm humid surface parcel under a cool but very dry troposphere: the reference humidity exceeds
    ! the column humidity (P_q <= 0) while the parcel is buoyant (P_T > 0).
    temperature = [215.0_real64, 245.0_real64, 268.0_real64, 302.0_real64]
    do k = 1, 3
      humidity(k) = 0.05_real64*saturation_specific_humidity(temperature(k), full_level_pressure(k))
    end do
    humidity(4) = 0.95_real64*saturation_specific_humidity(temperature(4), full_level_pressure(4))
    call moist_convection_reference_profile(moist_adjustment, full_level_pressure, temperature, humidity, &
                                            reference_temperature, parcel_humidity, reference_humidity)
    if (sum((reference_humidity - humidity)*delta_pressure) <= 0.0_real64 .or. &
        sum((reference_temperature - temperature)*delta_pressure) <= 0.0_real64) then
      error stop 'shallow-convection test column is not of the shallow type'
    end if
    call moist_convective_adjustment_tendency(moist_adjustment, pressure_half, temperature, humidity, &
                                              temperature_tendency, humidity_tendency, precipitation)
    if (precipitation /= 0.0_real64) error stop 'shallow convection produced precipitation'
    if (maxval(abs(humidity_tendency)) <= 0.0_real64) error stop 'shallow convection did nothing'
    if (abs(sum(humidity_tendency*delta_pressure)) > 1.0e-12_real64*maxval(abs(humidity_tendency))*sum(delta_pressure)) then
      error stop 'shallow convection changes the column water'
    end if
    if (abs(sum(temperature_tendency*delta_pressure)) > 1.0e-9_real64*maxval(abs(temperature_tendency))*sum(delta_pressure)) &
      then
      error stop 'shallow convection changes the column enthalpy'
    end if
    ! The reference humidity below the top of the shallow layer never exceeds RH_c q_s(T_ref), so the
    ! adjusted column is not supersaturated with respect to the reference temperature.
    do k = 1, 4
      if (humidity(k) + moist_adjustment%adjustment_time*humidity_tendency(k) > &
          max(humidity(k), reference_humidity(k)) + 1.0e-15_real64) then
        error stop 'shallow convection overshoots the reference humidity'
      end if
    end do
  end subroutine check_moist_convection_shallow

  !> Cloud cover from the column-maximum relative humidity and the convective
  !> precipitation (docs/tendency/cloud.md).
  subroutine check_cloud_diagnosis()
    type(cloud_config) :: clouds
    real(real64), parameter :: mm_per_day = 1.0_real64/86400.0_real64
    real(real64), parameter :: temperature(4) = [220.0_real64, 250.0_real64, 275.0_real64, 288.0_real64]
    real(real64) :: full_pressure(4), delta_pressure(4), saturation(4), humidity(4), cover
    integer :: k

    if (abs(clouds%critical_relative_humidity - 0.8_real64) > 0.0_real64 .or. &
        abs(clouds%convective_intercept - 0.245_real64) > 0.0_real64 .or. &
        abs(clouds%convective_slope - 0.125_real64) > 0.0_real64 .or. &
        abs(clouds%convective_reference_precipitation - mm_per_day) > 1.0e-20_real64 .or. &
        abs(clouds%convective_maximum_cover - 0.8_real64) > 0.0_real64) then
      error stop 'cloud defaults differ from docs/tendency/cloud.md'
    end if
    if (relative_humidity_cloud_cover(clouds, 0.5_real64) /= 0.0_real64 .or. &
        relative_humidity_cloud_cover(clouds, 0.8_real64) /= 0.0_real64 .or. &
        abs(relative_humidity_cloud_cover(clouds, 0.9_real64) - 0.25_real64) > 1.0e-14_real64 .or. &
        abs(relative_humidity_cloud_cover(clouds, 1.0_real64) - 1.0_real64) > 1.0e-14_real64 .or. &
        abs(relative_humidity_cloud_cover(clouds, 1.2_real64) - 1.0_real64) > 1.0e-14_real64) then
      error stop 'relative-humidity cloud cover is incorrect'
    end if
    if (convective_cloud_cover(clouds, 0.0_real64) /= 0.0_real64 .or. &
        convective_cloud_cover(clouds, 0.1_real64*mm_per_day) /= 0.0_real64 .or. &
        abs(convective_cloud_cover(clouds, mm_per_day) - 0.245_real64) > 1.0e-14_real64 .or. &
        abs(convective_cloud_cover(clouds, 10.0_real64*mm_per_day) - (0.245_real64 + 0.125_real64*log(10.0_real64))) > &
        1.0e-14_real64 .or. &
        abs(convective_cloud_cover(clouds, 100.0_real64*mm_per_day) - 0.8_real64) > 0.0_real64 .or. &
        convective_cloud_cover(clouds, 5.0_real64*mm_per_day) <= convective_cloud_cover(clouds, 2.0_real64*mm_per_day)) then
      error stop 'convective cloud cover is incorrect'
    end if
    call full_level_pressures(pressure_half, full_pressure, delta_pressure)
    do k = 1, 4
      saturation(k) = saturation_specific_humidity(temperature(k), full_pressure(k))
    end do
    humidity = 0.5_real64*saturation
    if (diagnose_cloud_cover(clouds, full_pressure, temperature, humidity, 0.0_real64) /= 0.0_real64) then
      error stop 'an unsaturated column without convective rain has cloud'
    end if
    humidity(3) = saturation(3)
    if (abs(diagnose_cloud_cover(clouds, full_pressure, temperature, humidity, 0.0_real64) - 1.0_real64) > 1.0e-14_real64) then
      error stop 'a saturated layer does not give full cloud cover'
    end if
    humidity = 0.9_real64*saturation
    humidity(1) = -0.5_real64
    cover = diagnose_cloud_cover(clouds, full_pressure, temperature, humidity, 10.0_real64*mm_per_day)
    if (abs(cover - max(0.25_real64, 0.245_real64 + 0.125_real64*log(10.0_real64))) > 1.0e-14_real64) then
      error stop 'cloud cover is not the larger of the humidity and convective covers'
    end if
    humidity = -0.5_real64
    if (diagnose_cloud_cover(clouds, full_pressure, temperature, humidity, 0.0_real64) /= 0.0_real64) then
      error stop 'negative humidity produced cloud'
    end if
  end subroutine check_cloud_diagnosis

  !> The cloud reflects C alpha_c of the downward shortwave once below the ozone layer
  !> (docs/tendency/shortwave-radiation.md): C = 0 reproduces the cloud-free column bit for
  !> bit, the atmosphere is not heated by clouds, and the column energy budget closes.
  subroutine check_cloud_shortwave_reflection()
    real(real64), parameter :: temperature(4) = [220.0_real64, 250.0_real64, 275.0_real64, 288.0_real64]
    real(real64) :: clear_tendency(4), cloudy_tendency(4), clear_surface, cloudy_surface, deep
    real(real64) :: clear_incoming, clear_reflected, clear_outgoing, incoming, reflected, outgoing
    real(real64) :: surface_shortwave, expected_reflected, total_energy_tendency, cover
    type(dry_model_physics_config) :: physics
    type(radiation_config) :: ocean
    integer :: quarter

    physics = moist_case_physics()
    ocean = physics%radiation
    if (abs(ocean%cloud_shortwave_albedo - 0.43_real64) > 0.0_real64 .or. &
        abs(ocean%surface_shortwave_albedo - 0.06_real64) > 0.0_real64) then
      error stop 'moist case albedos differ from docs/tendency/shortwave-radiation.md'
    end if
    call radiation_tendency(ocean, pressure_half, temperature, 290.0_real64, 0.0_real64, 3.0_real64, 4.0_real64, &
                            0.0_real64, 0.0_real64, 0.0_real64, clear_tendency, clear_surface, deep, &
                            clear_incoming, clear_reflected, clear_outgoing, latent_heat_flux=50.0_real64)
    call radiation_tendency(ocean, pressure_half, temperature, 290.0_real64, 0.0_real64, 3.0_real64, 4.0_real64, &
                            0.0_real64, 0.0_real64, 0.0_real64, cloudy_tendency, cloudy_surface, deep, &
                            incoming, reflected, outgoing, latent_heat_flux=50.0_real64, cloud_cover=0.0_real64)
    if (any(cloudy_tendency /= clear_tendency) .or. cloudy_surface /= clear_surface .or. &
        reflected /= clear_reflected .or. incoming /= clear_incoming .or. outgoing /= clear_outgoing) then
      error stop 'zero cloud cover changes the radiation column'
    end if
    surface_shortwave = clear_reflected/ocean%surface_shortwave_albedo
    do quarter = 1, 4
      cover = 0.25_real64*real(quarter, real64)
      call radiation_tendency(ocean, pressure_half, temperature, 290.0_real64, 0.0_real64, 3.0_real64, 4.0_real64, &
                              0.0_real64, 0.0_real64, 0.0_real64, cloudy_tendency, cloudy_surface, deep, &
                              incoming, reflected, outgoing, latent_heat_flux=50.0_real64, cloud_cover=cover)
      expected_reflected = (cover*ocean%cloud_shortwave_albedo + &
        (1.0_real64 - cover*ocean%cloud_shortwave_albedo)*ocean%surface_shortwave_albedo)*surface_shortwave
      if (abs(reflected - expected_reflected) > 1.0e-12_real64*expected_reflected .or. &
          any(cloudy_tendency /= clear_tendency) .or. cloudy_surface >= clear_surface .or. &
          incoming /= clear_incoming .or. outgoing /= clear_outgoing) then
        error stop 'cloud reflection of the downward shortwave is incorrect'
      end if
      total_energy_tendency = sum(ocean%dry_air_specific_heat*(pressure_half(1:4) - pressure_half(0:3))/ &
        earth_gravity*cloudy_tendency) + radiation_surface_heat_capacity(ocean)*cloudy_surface + 50.0_real64
      if (abs(total_energy_tendency - (incoming - reflected - outgoing)) > 1.0e-9_real64) then
        error stop 'cloudy radiation column does not conserve energy'
      end if
    end do
  end subroutine check_cloud_shortwave_reflection

  subroutine check_moist_case_configuration()
    type(dry_model_physics_config) :: physics
    type(radiation_config) :: seasonal_radiation

    physics = moist_case_physics()
    if (.not. (physics%moisture%enabled .and. physics%evaporation%enabled .and. physics%moist_convection%enabled .and. &
               physics%condensation%enabled .and. physics%convection%enabled .and. physics%radiation%enabled .and. &
               physics%radiation%slab_ocean_enabled .and. physics%surface_friction%enabled .and. &
               physics%rayleigh_friction%enabled .and. physics%cloud%enabled) .or. physics%held_suarez%enabled .or. &
        physics%evaporation%surface_wetness /= 1.0_real64 .or. &
        physics%radiation%surface_shortwave_albedo /= seasonal_radiation%ocean_shortwave_albedo .or. &
        abs(physics%radiation%surface_shortwave_albedo - 0.06_real64) > 0.0_real64 .or. &
        physics%radiation%axial_tilt /= seasonal_radiation%axial_tilt .or. &
        abs(physics%radiation%axial_tilt - 23.4_real64*acos(-1.0_real64)/180.0_real64) > 1.0e-15_real64 .or. &
        abs(physics%moisture%initial_relative_humidity - 0.7_real64) > 0.0_real64 .or. &
        abs(physics%moisture%initial_humidity_top_pressure - 2.0e4_real64) > 0.0_real64 .or. &
        abs(physics%moist_convection%adjustment_time - 7200.0_real64) > 0.0_real64 .or. &
        abs(physics%moist_convection%reference_relative_humidity - 0.7_real64) > 0.0_real64) then
      error stop 'moist case configuration does not match docs/cases/moist.md'
    end if
  end subroutine check_moist_case_configuration

  !> Enabling the humidity variable with q = 0 leaves the dry tendencies unchanged.
  !> The longwave radiation sees the prognostic humidity once it is enabled, so the
  !> dry reference humidity it would otherwise use is set to zero for the comparison.
  subroutine check_dry_atmosphere_unchanged_by_zero_humidity()
    type(harmonic_transform) :: transform
    type(hybrid_sigma_coordinate) :: coordinate
    type(dry_workspace_type) :: workspace
    type(dry_state_type) :: state
    type(dry_tendency_type) :: dry, moist
    type(dry_model_physics_config) :: physics
    complex(real64), allocatable :: zeta(:, :, :), delta(:, :, :), temperature(:, :, :)
    complex(real64), allocatable :: log_ps(:, :), surface_geopotential(:, :)
    real(real64), allocatable :: lowest_grid(:, :)
    real(real64) :: speed_dry, speed_moist
    integer :: levels

    call transform%init(truncation)
    call coordinate%init_default()
    levels = coordinate%number_of_levels
    call held_suarez_initial_state(transform, truncation, coordinate, zeta, delta, temperature, log_ps, &
                                   surface_geopotential)
    call allocate_dry_state(state, truncation, levels)
    state%zeta = zeta
    state%delta = delta
    state%temperature = temperature
    state%log_surface_pressure = log_ps
    call transform%spectral_to_grid(temperature(:, :, levels), lowest_grid)
    call allocate_dry_surface_fields(state, size(lowest_grid, 1), size(lowest_grid, 2))
    state%surface_temperature = lowest_grid
    state%deep_temperature = lowest_grid
    call workspace%initialize(transform, truncation, levels)
    call allocate_dry_tendency(dry, truncation, levels, workspace%nx, workspace%ny)
    call allocate_dry_tendency(moist, truncation, levels, workspace%nx, workspace%ny)
    physics = slab_ocean_case_physics()
    physics%radiation%longwave_reference_surface_humidity = 0.0_real64
    call evaluate_dry_tendency(transform, truncation, coordinate, planet_config(), state, state, surface_geopotential, &
                               physics, workspace, 0.0_real64, 2.0_real64*time_step, dry, speed_dry)
    physics%moisture%enabled = .true.
    call evaluate_dry_tendency(transform, truncation, coordinate, planet_config(), state, state, surface_geopotential, &
                               physics, workspace, 0.0_real64, 2.0_real64*time_step, moist, speed_moist)
    if (any(moist%zeta /= dry%zeta) .or. any(moist%delta /= dry%delta) .or. any(moist%temperature /= dry%temperature) .or. &
        any(moist%log_surface_pressure /= dry%log_surface_pressure) .or. &
        any(moist%surface_temperature /= dry%surface_temperature) .or. speed_moist /= speed_dry) then
      error stop 'carrying zero humidity changes the dry tendencies'
    end if
    if (any(moist%specific_humidity /= 0.0_real64) .or. any(dry%specific_humidity /= 0.0_real64)) then
      error stop 'a dry atmosphere has a humidity tendency'
    end if
  end subroutine check_dry_atmosphere_unchanged_by_zero_humidity

  !> On a resting atmosphere the humidity tendency has no advective part, so its
  !> mass-weighted global integral must equal <E> - <P> exactly.
  subroutine check_column_water_budget_closure()
    type(harmonic_transform) :: transform
    type(hybrid_sigma_coordinate) :: coordinate
    type(dry_workspace_type) :: workspace
    type(dry_state_type) :: state
    type(dry_tendency_type) :: rhs
    type(dry_model_physics_config) :: physics
    type(radiation_diagnostics) :: sample
    complex(real64), allocatable :: zeta(:, :, :), delta(:, :, :), temperature(:, :, :)
    complex(real64), allocatable :: log_ps(:, :), surface_geopotential(:, :), humidity_spectral(:, :)
    real(real64), allocatable :: tendency_grid(:, :), weights(:), humidity_grid(:, :), full_level_pressure(:)
    real(real64), allocatable :: delta_pressure(:), pressure_half_column(:)
    integer, allocatable :: nlon(:)
    real(real64) :: speed, column_tendency, area_weight, budget
    integer :: levels, i, j, k

    call transform%init(truncation)
    nlon = transform%get_nlon()
    weights = transform%get_gaussian_weights()
    call coordinate%init_default()
    levels = coordinate%number_of_levels
    call held_suarez_initial_state(transform, truncation, coordinate, zeta, delta, temperature, log_ps, &
                                   surface_geopotential)
    call allocate_dry_state(state, truncation, levels)
    state%zeta = zeta
    state%delta = delta
    state%temperature = temperature
    state%log_surface_pressure = log_ps
    ! A warm ocean under a uniform 264 K atmosphere drives evaporation and convection;
    ! saturating two mid-levels drives large-scale condensation as well.
    call transform%allocate_field(humidity_grid)
    call allocate_dry_surface_fields(state, size(humidity_grid, 1), size(humidity_grid, 2))
    state%surface_temperature = 300.0_real64
    state%deep_temperature = 300.0_real64
    allocate (full_level_pressure(levels), delta_pressure(levels))
    pressure_half_column = coordinate%reference_p_half
    call full_level_pressures(pressure_half_column, full_level_pressure, delta_pressure)
    do k = 1, levels
      humidity_grid = 0.3_real64*saturation_specific_humidity(264.0_real64, full_level_pressure(k))
      if (k == 8 .or. k == 9) humidity_grid = 1.3_real64*saturation_specific_humidity(264.0_real64, full_level_pressure(k))
      call transform%grid_to_spectral(humidity_grid, humidity_spectral)
      state%specific_humidity(:, :, k) = humidity_spectral
    end do
    call workspace%initialize(transform, truncation, levels)
    call allocate_dry_tendency(rhs, truncation, levels, workspace%nx, workspace%ny)
    physics = moist_case_physics()
    call evaluate_dry_tendency(transform, truncation, coordinate, planet_config(), state, state, surface_geopotential, &
                               physics, workspace, 0.0_real64, 2.0_real64*time_step, rhs, speed, sample)
    if (speed /= 0.0_real64) error stop 'water budget test atmosphere is not at rest'
    if (sample%mean_evaporation <= 0.0_real64 .or. sample%mean_large_scale_precipitation <= 0.0_real64) then
      error stop 'water budget test does not exercise evaporation and condensation'
    end if
    budget = 0.0_real64
    do k = 1, levels
      call transform%spectral_to_grid(rhs%specific_humidity(:, :, k), tendency_grid)
      do j = 1, size(nlon)
        area_weight = 0.5_real64*weights(j)/real(nlon(j), real64)
        do i = 1, nlon(j)
          column_tendency = tendency_grid(i, j)*delta_pressure(k)/earth_gravity
          budget = budget + area_weight*column_tendency
        end do
      end do
    end do
    if (abs(budget - (sample%mean_evaporation - sample%mean_convective_precipitation - &
                      sample%mean_large_scale_precipitation)) > 1.0e-9_real64*sample%mean_evaporation) then
      write (*, '(a,3es16.8)') 'dW/dt, E, P = ', budget, sample%mean_evaporation, &
        sample%mean_convective_precipitation + sample%mean_large_scale_precipitation
      error stop 'column water tendency does not close as E - P'
    end if
    if (sample%mean_precipitable_water <= 0.0_real64 .or. &
        abs(sample%mean_precipitable_water - sample%mean_signed_column_water - sample%mean_negative_column_water) > &
        1.0e-12_real64*sample%mean_precipitable_water) then
      error stop 'column water diagnostics are inconsistent'
    end if
    if (sample%maximum_wind_speed /= 0.0_real64 .or. sample%maximum_wind_level < 1 .or. &
        sample%maximum_wind_level > levels .or. sample%maximum_wind_eta <= 0.0_real64) then
      error stop 'maximum wind diagnostic of a resting atmosphere is incorrect'
    end if
  end subroutine check_column_water_budget_closure

  !> The moist case starts at the configured relative humidity.  A short
  !> integration from a dry start stays finite, moistens from the ocean, and its
  !> daily-mean water budget is consistent with the change of column water.
  subroutine check_moist_integration()
    type(dry_atmosphere_solver) :: solver, dry_start_solver
    type(harmonic_transform) :: transform
    type(dry_model_physics_config) :: physics
    type(hybrid_sigma_coordinate) :: coordinate
    type(radiation_case_diagnostics) :: diagnostics
    type(radiation_diagnostics) :: sample, means
    type(radiation_monthly_means) :: monthly
    real(real64), allocatable :: zeta(:, :, :), delta(:, :, :), temperature(:, :, :), humidity(:, :, :)
    real(real64), allocatable :: surface_pressure(:, :), u(:, :, :), v(:, :, :), ocean(:, :), deep(:, :)
    real(real64), allocatable :: cloud_cover(:, :)
    integer, allocatable :: nlon(:)
    real(real64) :: first_water, last_water, accumulated_source, time_step_used
    real(real64), allocatable :: column_half(:), full_pressure(:), delta_pressure(:), weights(:)
    real(real64), allocatable :: expected(:, :)
    real(real64) :: ring_expected, ring_humidity, level_expected, level_humidity, level_maximum
    integer :: step, i, j, k, levels
    integer, parameter :: steps = 12

    call transform%init(truncation)
    nlon = transform%get_nlon()
    call solver%init(truncation, time_step)
    physics = moist_case_physics()
    call set_radiation_case_state(solver, transform, physics, radiation_case_planet(physics))
    call solver%get_fields(zeta, delta, temperature, surface_pressure, u, v, specific_humidity=humidity)
    ! q^0 = RH q_s(T_k, p_k) below the top pressure and zero above, up to the spectral truncation.  The
    ! truncation preserves the Gaussian-quadrature global mean of every level exactly.
    coordinate = solver%get_coordinate()
    levels = coordinate%number_of_levels
    weights = transform%get_gaussian_weights()
    allocate (column_half(0:levels), full_pressure(levels), delta_pressure(levels))
    allocate (expected(size(humidity, 1), size(humidity, 2)))
    do k = 1, levels
      expected = 0.0_real64
      level_expected = 0.0_real64
      level_humidity = 0.0_real64
      do j = 1, size(humidity, 2)
        ring_expected = 0.0_real64
        ring_humidity = 0.0_real64
        do i = 1, nlon(j)
          column_half = coordinate%a_half + coordinate%b_half*surface_pressure(i, j)
          call full_level_pressures(column_half, full_pressure, delta_pressure)
          if (full_pressure(k) >= physics%moisture%initial_humidity_top_pressure) then
            expected(i, j) = physics%moisture%initial_relative_humidity* &
              saturation_specific_humidity(temperature(i, j, k), full_pressure(k))
          end if
          ring_expected = ring_expected + expected(i, j)
          ring_humidity = ring_humidity + humidity(i, j, k)
        end do
        level_expected = level_expected + 0.5_real64*weights(j)*ring_expected/real(nlon(j), real64)
        level_humidity = level_humidity + 0.5_real64*weights(j)*ring_humidity/real(nlon(j), real64)
      end do
      if (abs(level_humidity - level_expected) > 1.0e-9_real64*level_expected + 1.0e-15_real64) then
        error stop 'global-mean initial humidity does not match the initial relative humidity'
      end if
      ! The truncation error is absolute and largest at the cold poles, so compare with the level maximum.
      level_maximum = maxval(expected)
      if (level_maximum == 0.0_real64) then
        if (maxval(abs(humidity(:, :, k))) > 1.0e-12_real64) error stop 'moist case does not start with a dry stratosphere'
      else
        do j = 1, size(humidity, 2)
          if (any(abs(humidity(1:nlon(j), j, k) - expected(1:nlon(j), j)) > 0.25_real64*level_maximum)) then
            error stop 'moist case does not start near the initial relative humidity'
          end if
        end do
      end if
    end do

    ! The budget closure is checked from a dry start, where the advection and filter error of the
    ! strongly varying initial humidity field does not mask the small accumulated E - P.
    physics%moisture%initial_relative_humidity = 0.0_real64
    call dry_start_solver%init(truncation, time_step)
    call set_radiation_case_state(dry_start_solver, transform, physics, radiation_case_planet(physics))
    call dry_start_solver%get_fields(zeta, delta, temperature, surface_pressure, u, v, specific_humidity=humidity)
    if (maxval(abs(humidity)) > 0.0_real64) error stop 'zero initial relative humidity does not start dry'
    call diagnostics%reset()
    accumulated_source = 0.0_real64
    first_water = 0.0_real64
    last_water = 0.0_real64
    time_step_used = time_step
    do step = 1, steps
      call dry_start_solver%advance()
      call dry_start_solver%take_latest_diagnostics(sample)
      if (step == 1) first_water = sample%mean_signed_column_water
      last_water = sample%mean_signed_column_water
      ! The first two advances are the dt/4 and dt/2 startup; afterwards each advance moves the state by dt.
      if (step > 1) accumulated_source = accumulated_source + time_step*(sample%mean_evaporation - &
        sample%mean_convective_precipitation - sample%mean_large_scale_precipitation)
      call diagnostics%add(sample)
    end do
    call diagnostics%take_daily(means)
    call diagnostics%take_monthly(monthly)
    call dry_start_solver%get_fields(zeta, delta, temperature, surface_pressure, u, v, surface_temperature=ocean, &
                           deep_temperature=deep, specific_humidity=humidity)
    if (.not. all(ieee_is_finite(humidity)) .or. .not. all(ieee_is_finite(temperature)) .or. &
        .not. all(ieee_is_finite(ocean)) .or. .not. all(ieee_is_finite(u))) then
      error stop 'moist integration produced a non-finite value'
    end if
    ! Spectral truncation is allowed to make q locally negative in dry regions;
    ! require that evaporation has nevertheless moistened the lowest level.
    if (maxval(humidity(:, :, size(humidity, 3))) <= 0.0_real64) then
      error stop 'the ocean did not moisten the lowest level'
    end if
    if (means%mean_evaporation <= 0.0_real64 .or. means%mean_latent_heat_flux <= 0.0_real64 .or. &
        means%mean_precipitable_water <= 0.0_real64 .or. means%maximum_wind_speed <= 0.0_real64 .or. &
        means%mean_cloud_cover < 0.0_real64 .or. means%mean_cloud_cover > 1.0_real64 .or. &
        means%maximum_wind_level < 1 .or. abs(means%maximum_wind_latitude_degrees) > 90.0_real64 .or. &
        means%maximum_wind_longitude_degrees < 0.0_real64 .or. means%maximum_wind_longitude_degrees >= 360.0_real64) then
      error stop 'moist daily diagnostics are incorrect'
    end if
    if (.not. allocated(monthly%precipitation) .or. .not. allocated(monthly%eddy_vq) .or. &
        .not. allocated(monthly%cloud_cover) .or. &
        .not. all(ieee_is_finite(monthly%precipitable_water)) .or. .not. all(ieee_is_finite(monthly%zonal_humidity))) then
      error stop 'moist monthly diagnostics are missing or non-finite'
    end if
    if (any(monthly%cloud_cover < 0.0_real64) .or. any(monthly%cloud_cover > 1.0_real64)) then
      error stop 'monthly cloud cover is outside [0,1]'
    end if
    call dry_start_solver%get_cloud_cover(cloud_cover)
    if (any(shape(cloud_cover) /= shape(ocean)) .or. any(cloud_cover < 0.0_real64) .or. any(cloud_cover > 1.0_real64)) then
      error stop 'solver cloud cover snapshot is missing or outside [0,1]'
    end if
    ! The column water sampled at the start of each step changes by the accumulated E - P up to the
    ! advection and filter error, which is small for the nearly balanced initial state over 4 hours.
    if (abs((last_water - first_water) - accumulated_source) > 0.1_real64*abs(accumulated_source)) then
      write (*, '(a,3es16.8)') 'dW, int(E-P), first W = ', last_water - first_water, accumulated_source, first_water
      error stop 'moist integration water budget does not close within tolerance'
    end if
  end subroutine check_moist_integration

end program check_moist_atmosphere
