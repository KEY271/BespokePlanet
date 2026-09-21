!> Grid-scale condensation of supersaturated water vapour
!> (docs/tendency/large-scale-condensation.md).
!>
!> In every supersaturated layer the condensate C_k that brings the layer
!> exactly to saturation after its latent heating is found by Newton iteration,
!> and removed within one leapfrog step.  The condensate leaves the column as
!> precipitation; c_p T + L q of each layer is unchanged.
module large_scale_condensation
  use iso_fortran_env, only: real64
  use planet_parameters, only: earth_gravity
  use dry_physics_config, only: condensation_config
  use moist_thermodynamics, only: dry_air_specific_heat, latent_heat_of_condensation, &
                                  saturation_specific_humidity_and_derivative, full_level_pressures
  implicit none
  private
  public :: large_scale_condensation_tendency, large_scale_condensation_from_levels, saturation_condensate

contains

  !> Condensate C >= 0 with q - C = q_s(T + L C/c_p, p) for a layer at (T, q, p); zero when unsaturated.
  !> q_s and dq_s/dT of each iterate come from one evaluation of the saturation vapour pressure.
  pure real(real64) function saturation_condensate(config, temperature, humidity, pressure) result(condensate)
    type(condensation_config), intent(in) :: config
    real(real64), intent(in) :: temperature, humidity, pressure
    real(real64) :: saturation, warmed_saturation, derivative, heating_per_condensate, warmed_temperature, residual
    integer :: iteration

    condensate = 0.0_real64
    call saturation_specific_humidity_and_derivative(temperature, pressure, saturation, derivative)
    if (humidity <= saturation) return
    heating_per_condensate = latent_heat_of_condensation/dry_air_specific_heat
    residual = humidity - saturation
    do iteration = 1, config%maximum_iterations
      condensate = condensate + residual/(1.0_real64 + heating_per_condensate*derivative)
      warmed_temperature = temperature + heating_per_condensate*condensate
      call saturation_specific_humidity_and_derivative(warmed_temperature, pressure, warmed_saturation, derivative)
      residual = humidity - condensate - warmed_saturation
      if (abs(residual) <= config%saturation_tolerance*saturation) exit
    end do
  end function saturation_condensate

  !> Tendencies that remove the supersaturation of one column over `interval`
  !> (the leapfrog advance width 2 dt) and the resulting precipitation in kg m^-2 s^-1.
  pure subroutine large_scale_condensation_tendency(config, pressure_half, temperature, humidity, interval, &
                                                    temperature_tendency, humidity_tendency, precipitation)
    type(condensation_config), intent(in) :: config
    real(real64), intent(in) :: pressure_half(0:), temperature(:), humidity(:), interval
    real(real64), intent(out) :: temperature_tendency(:), humidity_tendency(:), precipitation
    real(real64) :: full_level_pressure(size(temperature)), delta_pressure(size(temperature))

    call full_level_pressures(pressure_half, full_level_pressure, delta_pressure)
    call large_scale_condensation_from_levels(config, full_level_pressure, delta_pressure, temperature, humidity, &
                                              interval, temperature_tendency, humidity_tendency, precipitation)
  end subroutine large_scale_condensation_tendency

  !> As above with the full-level pressures and layer thicknesses supplied.
  pure subroutine large_scale_condensation_from_levels(config, full_level_pressure, delta_pressure, temperature, &
                                                       humidity, interval, temperature_tendency, humidity_tendency, &
                                                       precipitation)
    type(condensation_config), intent(in) :: config
    real(real64), intent(in) :: full_level_pressure(:), delta_pressure(:), temperature(:), humidity(:), interval
    real(real64), intent(out) :: temperature_tendency(:), humidity_tendency(:), precipitation
    real(real64) :: condensate
    integer :: k

    precipitation = 0.0_real64
    do k = 1, size(temperature)
      condensate = saturation_condensate(config, temperature(k), max(humidity(k), 0.0_real64), full_level_pressure(k))
      humidity_tendency(k) = -condensate/interval
      temperature_tendency(k) = latent_heat_of_condensation/dry_air_specific_heat*condensate/interval
      precipitation = precipitation + condensate*delta_pressure(k)/(earth_gravity*interval)
    end do
  end subroutine large_scale_condensation_from_levels

end module large_scale_condensation
