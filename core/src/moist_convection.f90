!> Simplified Betts-Miller moist convective adjustment after Frierson (2007)
!> (docs/tendency/moist-convective-adjustment.md).
!>
!> The reference temperature is the moist adiabat of a parcel lifted from the
!> lowest level, the reference humidity a fixed relative humidity along it.  Only
!> columns in which the parcel is buoyant somewhere are adjusted.  The reference
!> profile is shifted so that column moist enthalpy c_p T + L q is conserved;
!> the water removed is convective precipitation.  Shallow (non-precipitating)
!> convection is limited to the depth over which the humidity change integrates
!> to zero and redistributes enthalpy and water vertically.
!>
!> The *_from_levels entry points take the full-level pressures, layer
!> thicknesses and Exner function (p/p_0)^kappa that the caller already has, so
!> that a column evaluated by several processes computes them once.
module moist_convection
  use iso_fortran_env, only: real64
  use planet_parameters, only: earth_gravity
  use dry_vertical_coordinate, only: dry_air_kappa, reference_surface_pressure
  use dry_physics_config, only: moist_convection_config
  use moist_thermodynamics, only: dry_air_specific_heat, latent_heat_of_condensation, &
                                  virtual_temperature_coefficient, saturation_specific_humidity, &
                                  saturation_specific_humidity_and_derivative, &
                                  equivalent_potential_temperature, &
                                  lifting_condensation_level, full_level_pressures
  implicit none
  private
  public :: moist_convective_adjustment_tendency, moist_convection_reference_profile
  public :: moist_convective_adjustment_from_levels, moist_convection_reference_profile_from_levels

  integer, parameter :: reference_maximum_iterations = 20
  real(real64), parameter :: reference_temperature_tolerance = 1.0e-3_real64
  !> L/c_p, the warming per unit condensate.
  real(real64), parameter :: heating_per_condensate = latent_heat_of_condensation/dry_air_specific_heat

contains

  !> Reference temperature and parcel humidity of the parcel lifted from level N,
  !> and the reference humidity RH_c q_s(T_ref, p).
  pure subroutine moist_convection_reference_profile(config, full_level_pressure, temperature, humidity, &
                                                     reference_temperature, parcel_humidity, reference_humidity)
    type(moist_convection_config), intent(in) :: config
    real(real64), intent(in) :: full_level_pressure(:), temperature(:), humidity(:)
    real(real64), intent(out) :: reference_temperature(:), parcel_humidity(:), reference_humidity(:)
    real(real64) :: exner(size(temperature))

    exner = (full_level_pressure/reference_surface_pressure)**dry_air_kappa
    call moist_convection_reference_profile_from_levels(config, full_level_pressure, exner, temperature, humidity, &
                                                        reference_temperature, parcel_humidity, reference_humidity)
  end subroutine moist_convection_reference_profile

  !> As above with the Exner function (p_k/p_0)^kappa of every level supplied.
  !>
  !> Along the dry part of the ascent the saturation humidity evaluated for the
  !> saturation test is also the reference humidity, and along the saturated part
  !> the parcel humidity is; so q_s is evaluated once per level.
  pure subroutine moist_convection_reference_profile_from_levels(config, full_level_pressure, exner, temperature, &
                                                                 humidity, reference_temperature, parcel_humidity, &
                                                                 reference_humidity)
    type(moist_convection_config), intent(in) :: config
    real(real64), intent(in) :: full_level_pressure(:), exner(:), temperature(:), humidity(:)
    real(real64), intent(out) :: reference_temperature(:), parcel_humidity(:), reference_humidity(:)
    real(real64) :: parcel_temperature, parcel_water, parcel_exner, dry_temperature, dry_saturation
    real(real64) :: level_temperature, level_pressure, log_theta_e
    logical :: saturated
    integer :: levels, k

    levels = size(temperature)
    parcel_temperature = temperature(levels)
    parcel_water = max(humidity(levels), 0.0_real64)
    parcel_exner = exner(levels)
    reference_temperature(levels) = parcel_temperature
    parcel_humidity(levels) = parcel_water
    reference_humidity(levels) = config%reference_relative_humidity* &
      saturation_specific_humidity(parcel_temperature, full_level_pressure(levels))
    saturated = .false.
    log_theta_e = 0.0_real64
    dry_saturation = 0.0_real64
    do k = levels - 1, 1, -1
      dry_temperature = parcel_temperature*exner(k)/parcel_exner
      if (.not. saturated) then
        dry_saturation = saturation_specific_humidity(dry_temperature, full_level_pressure(k))
        if (parcel_water > dry_saturation) then
          saturated = .true.
          call lifting_condensation_level(parcel_temperature, parcel_water, full_level_pressure(levels), &
                                          level_temperature, level_pressure)
          log_theta_e = log(equivalent_potential_temperature(level_temperature, parcel_water, level_pressure))
        end if
      end if
      if (saturated) then
        ! ln theta_e = ln T - ln Pi + L q_s/(c_p T): the pressure factor is folded into the target.
        reference_temperature(k) = saturated_reference_temperature(log_theta_e + log(exner(k)), &
                                                                   full_level_pressure(k), dry_temperature, &
                                                                   reference_temperature(k + 1))
        parcel_humidity(k) = saturation_specific_humidity(reference_temperature(k), full_level_pressure(k))
        reference_humidity(k) = config%reference_relative_humidity*parcel_humidity(k)
      else
        reference_temperature(k) = dry_temperature
        parcel_humidity(k) = parcel_water
        reference_humidity(k) = config%reference_relative_humidity*dry_saturation
      end if
    end do
  end subroutine moist_convection_reference_profile_from_levels

  !> Temperature at which the saturated parcel at `pressure` has the given reduced
  !> ln theta_e (ln theta_e + ln Pi), found by safeguarded Newton iteration in
  !> [lower_bound, upper_bound].  Each iterate costs one exp and one log.
  pure real(real64) function saturated_reference_temperature(reduced_target, pressure, lower_bound, &
                                                             upper_bound) result(temperature)
    real(real64), intent(in) :: reduced_target, pressure, lower_bound, upper_bound
    real(real64) :: lower, upper, residual, derivative, lower_residual, unused, change, trial
    integer :: iteration

    lower = min(lower_bound, upper_bound)
    upper = max(lower_bound, upper_bound)
    call saturated_residual(upper, pressure, reduced_target, residual, derivative)
    if (residual <= 0.0_real64) then
      temperature = upper
      return
    end if
    call saturated_residual(lower, pressure, reduced_target, lower_residual, unused)
    if (lower_residual >= 0.0_real64) then
      temperature = lower
      return
    end if
    temperature = upper
    do iteration = 1, reference_maximum_iterations
      if (residual > 0.0_real64) then
        upper = temperature
      else
        lower = temperature
      end if
      change = -residual/derivative
      trial = temperature + change
      if (trial <= lower .or. trial >= upper) then
        trial = 0.5_real64*(lower + upper)
        change = trial - temperature
      end if
      temperature = trial
      if (abs(change) < reference_temperature_tolerance .or. upper - lower < reference_temperature_tolerance) return
      call saturated_residual(temperature, pressure, reduced_target, residual, derivative)
    end do
    temperature = 0.5_real64*(lower + upper)
  end function saturated_reference_temperature

  !> ln T + L q_s(T, p)/(c_p T) - reduced_target and its temperature derivative
  !>   1/T + L/(c_p T) (dq_s/dT - q_s/T), which is always positive.
  pure subroutine saturated_residual(temperature, pressure, reduced_target, residual, derivative)
    real(real64), intent(in) :: temperature, pressure, reduced_target
    real(real64), intent(out) :: residual, derivative
    real(real64) :: humidity, humidity_derivative

    call saturation_specific_humidity_and_derivative(temperature, pressure, humidity, humidity_derivative)
    residual = log(temperature) + heating_per_condensate*humidity/temperature - reduced_target
    derivative = 1.0_real64/temperature + heating_per_condensate/temperature* &
      (humidity_derivative - humidity/temperature)
  end subroutine saturated_residual

  !> Convective tendencies of one column and its precipitation in kg m^-2 s^-1.
  !> `humidity` is the (non-negative) specific humidity the process sees.
  pure subroutine moist_convective_adjustment_tendency(config, pressure_half, temperature, humidity, &
                                                       temperature_tendency, humidity_tendency, precipitation)
    type(moist_convection_config), intent(in) :: config
    real(real64), intent(in) :: pressure_half(0:), temperature(:), humidity(:)
    real(real64), intent(out) :: temperature_tendency(:), humidity_tendency(:), precipitation
    real(real64) :: full_level_pressure(size(temperature)), delta_pressure(size(temperature))
    real(real64) :: exner(size(temperature))

    temperature_tendency = 0.0_real64
    humidity_tendency = 0.0_real64
    precipitation = 0.0_real64
    if (size(temperature) < 2) return
    call full_level_pressures(pressure_half, full_level_pressure, delta_pressure)
    exner = (full_level_pressure/reference_surface_pressure)**dry_air_kappa
    call moist_convective_adjustment_from_levels(config, full_level_pressure, delta_pressure, exner, temperature, &
                                                 humidity, temperature_tendency, humidity_tendency, precipitation)
  end subroutine moist_convective_adjustment_tendency

  !> As above with the full-level pressures, thicknesses and Exner function supplied.
  pure subroutine moist_convective_adjustment_from_levels(config, full_level_pressure, delta_pressure, exner, &
                                                          temperature, humidity, temperature_tendency, &
                                                          humidity_tendency, precipitation)
    type(moist_convection_config), intent(in) :: config
    real(real64), intent(in) :: full_level_pressure(:), delta_pressure(:), exner(:), temperature(:), humidity(:)
    real(real64), intent(out) :: temperature_tendency(:), humidity_tendency(:), precipitation
    integer :: levels, k, free_convection_level, convection_top, zero_crossing
    real(real64) :: reference_temperature(size(temperature)), parcel_humidity(size(temperature))
    real(real64) :: reference_humidity(size(temperature))
    real(real64) :: temperature_change(size(temperature)), humidity_change(size(temperature))
    real(real64) :: clipped_humidity(size(temperature))
    real(real64) :: mass_sum, temperature_sum, humidity_sum, temperature_shift, cumulative, previous_cumulative
    real(real64) :: scale, precipitation_from_humidity, precipitation_from_temperature
    logical :: buoyant

    levels = size(temperature)
    temperature_tendency = 0.0_real64
    humidity_tendency = 0.0_real64
    precipitation = 0.0_real64
    if (levels < 2) return
    clipped_humidity = max(humidity, 0.0_real64)
    call moist_convection_reference_profile_from_levels(config, full_level_pressure, exner, temperature, &
                                                        clipped_humidity, reference_temperature, parcel_humidity, &
                                                        reference_humidity)

    ! Level of free convection: the first level above the surface where the parcel is buoyant
    ! (judged by virtual temperature), and the top of the contiguous buoyant layer above it.
    free_convection_level = 0
    do k = levels - 1, 1, -1
      if (parcel_is_buoyant(k)) then
        free_convection_level = k
        exit
      end if
    end do
    if (free_convection_level == 0) return
    convection_top = free_convection_level
    do k = free_convection_level - 1, 1, -1
      buoyant = parcel_is_buoyant(k)
      if (.not. buoyant) exit
      convection_top = k
    end do

    temperature_change = 0.0_real64
    humidity_change = 0.0_real64
    mass_sum = 0.0_real64
    temperature_sum = 0.0_real64
    humidity_sum = 0.0_real64
    do k = convection_top, levels
      temperature_change(k) = reference_temperature(k) - temperature(k)
      humidity_change(k) = reference_humidity(k) - clipped_humidity(k)
      mass_sum = mass_sum + delta_pressure(k)
      temperature_sum = temperature_sum + temperature_change(k)*delta_pressure(k)
      humidity_sum = humidity_sum + humidity_change(k)*delta_pressure(k)
    end do
    precipitation_from_humidity = -humidity_sum/(earth_gravity*config%adjustment_time)
    precipitation_from_temperature = dry_air_specific_heat/latent_heat_of_condensation* &
      temperature_sum/(earth_gravity*config%adjustment_time)
    if (precipitation_from_temperature <= 0.0_real64) return

    if (precipitation_from_humidity > 0.0_real64) then
      ! Deep convection: shift the reference temperature so that column moist enthalpy is conserved.
      temperature_shift = (-latent_heat_of_condensation/dry_air_specific_heat*humidity_sum - temperature_sum)/mass_sum
      do k = convection_top, levels
        temperature_change(k) = temperature_change(k) + temperature_shift
      end do
      precipitation = precipitation_from_humidity
    else
      ! Shallow convection: restrict the depth so that the column humidity change vanishes.
      zero_crossing = levels
      cumulative = 0.0_real64
      previous_cumulative = 0.0_real64
      do k = levels, convection_top, -1
        previous_cumulative = cumulative
        cumulative = cumulative + humidity_change(k)*delta_pressure(k)
        if (cumulative >= 0.0_real64) then
          zero_crossing = k
          exit
        end if
      end do
      if (zero_crossing == levels) return
      scale = -previous_cumulative/(humidity_change(zero_crossing)*delta_pressure(zero_crossing))
      humidity_change(zero_crossing) = scale*humidity_change(zero_crossing)
      temperature_change(zero_crossing) = scale*temperature_change(zero_crossing)
      do k = convection_top, zero_crossing - 1
        temperature_change(k) = 0.0_real64
        humidity_change(k) = 0.0_real64
      end do
      mass_sum = 0.0_real64
      temperature_sum = 0.0_real64
      do k = zero_crossing, levels
        mass_sum = mass_sum + delta_pressure(k)
        temperature_sum = temperature_sum + temperature_change(k)*delta_pressure(k)
      end do
      temperature_shift = -temperature_sum/mass_sum
      do k = zero_crossing, levels
        temperature_change(k) = temperature_change(k) + temperature_shift
      end do
      convection_top = zero_crossing
      precipitation = 0.0_real64
    end if

    do k = convection_top, levels
      temperature_tendency(k) = temperature_change(k)/config%adjustment_time
      humidity_tendency(k) = humidity_change(k)/config%adjustment_time
    end do

  contains

    pure logical function parcel_is_buoyant(level)
      integer, intent(in) :: level
      parcel_is_buoyant = (1.0_real64 + virtual_temperature_coefficient*parcel_humidity(level))* &
        reference_temperature(level) > &
        (1.0_real64 + virtual_temperature_coefficient*clipped_humidity(level))*temperature(level)
    end function parcel_is_buoyant

  end subroutine moist_convective_adjustment_from_levels

end module moist_convection
