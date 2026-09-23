!> Grid-scale condensation of supersaturated water vapour
!> (docs/tendency/large-scale-condensation.md), with the optional snowfall of
!> docs/tendency/snow.md.
!>
!> In every supersaturated layer the condensate C_k that brings the layer
!> exactly to saturation after its latent heating is found by Newton iteration,
!> and removed within one leapfrog step.  The condensate leaves the column as
!> precipitation; c_p T + L q of each layer is unchanged.
!>
!> With snow enabled, the condensate of a layer colder than the formation
!> temperature is snow and releases L_v + L_f.  The snow falls from the top
!> down; a layer warmer than the atmospheric melting temperature melts it with
!> its sensible heat above that temperature, at most down to that temperature.
module large_scale_condensation
  use iso_fortran_env, only: real64
  use planet_parameters, only: earth_gravity
  use dry_physics_config, only: condensation_config, snow_config
  use moist_thermodynamics, only: dry_air_specific_heat, latent_heat_of_condensation, &
                                  saturation_specific_humidity_and_derivative, full_level_pressures
  implicit none
  private
  public :: large_scale_condensation_tendency, large_scale_condensation_from_levels, saturation_condensate

contains

  !> Condensate C >= 0 with q - C = q_s(T + L C/c_p, p) for a layer at (T, q, p); zero when unsaturated.
  !> q_s and dq_s/dT of each iterate come from one evaluation of the saturation vapour pressure.
  !> L is the latent heat of condensation unless latent_heat supplies another (L_v + L_f for snow).
  pure real(real64) function saturation_condensate(config, temperature, humidity, pressure, latent_heat) &
    result(condensate)
    type(condensation_config), intent(in) :: config
    real(real64), intent(in) :: temperature, humidity, pressure
    real(real64), intent(in), optional :: latent_heat
    real(real64) :: saturation, warmed_saturation, derivative, heating_per_condensate, warmed_temperature, residual
    integer :: iteration

    condensate = 0.0_real64
    call saturation_specific_humidity_and_derivative(temperature, pressure, saturation, derivative)
    if (humidity <= saturation) return
    if (present(latent_heat)) then
      heating_per_condensate = latent_heat/dry_air_specific_heat
    else
      heating_per_condensate = latent_heat_of_condensation/dry_air_specific_heat
    end if
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
                                                    temperature_tendency, humidity_tendency, precipitation, &
                                                    snow, snowfall, snow_melt)
    type(condensation_config), intent(in) :: config
    real(real64), intent(in) :: pressure_half(0:), temperature(:), humidity(:), interval
    real(real64), intent(out) :: temperature_tendency(:), humidity_tendency(:), precipitation
    type(snow_config), intent(in), optional :: snow
    real(real64), intent(out), optional :: snowfall, snow_melt
    real(real64) :: full_level_pressure(size(temperature)), delta_pressure(size(temperature))

    call full_level_pressures(pressure_half, full_level_pressure, delta_pressure)
    call large_scale_condensation_from_levels(config, full_level_pressure, delta_pressure, temperature, humidity, &
                                              interval, temperature_tendency, humidity_tendency, precipitation, &
                                              snow, snowfall, snow_melt)
  end subroutine large_scale_condensation_tendency

  !> As above with the full-level pressures and layer thicknesses supplied.
  !> `precipitation` is the whole large-scale precipitation reaching the surface
  !> (rain and snow), `snowfall` its snow part and `snow_melt` the snow melted in
  !> the column, all in kg m^-2 s^-1.  Without an enabled snow configuration
  !> every condensate is rain and both snow fluxes are zero.
  pure subroutine large_scale_condensation_from_levels(config, full_level_pressure, delta_pressure, temperature, &
                                                       humidity, interval, temperature_tendency, humidity_tendency, &
                                                       precipitation, snow, snowfall, snow_melt)
    type(condensation_config), intent(in) :: config
    real(real64), intent(in) :: full_level_pressure(:), delta_pressure(:), temperature(:), humidity(:), interval
    real(real64), intent(out) :: temperature_tendency(:), humidity_tendency(:), precipitation
    type(snow_config), intent(in), optional :: snow
    real(real64), intent(out), optional :: snowfall, snow_melt
    real(real64) :: condensate, latent_heat, falling_snow, melted_snow, melting_capacity, column_melt
    real(real64) :: layer_mass_rate, warmed_temperature
    logical :: snowing, forms_snow
    integer :: k

    precipitation = 0.0_real64
    falling_snow = 0.0_real64
    column_melt = 0.0_real64
    snowing = .false.
    if (present(snow)) snowing = snow%enabled
    if (.not. snowing) then
      do k = 1, size(temperature)
        condensate = saturation_condensate(config, temperature(k), max(humidity(k), 0.0_real64), &
                                           full_level_pressure(k))
        humidity_tendency(k) = -condensate/interval
        temperature_tendency(k) = latent_heat_of_condensation/dry_air_specific_heat*condensate/interval
        precipitation = precipitation + condensate*delta_pressure(k)/(earth_gravity*interval)
      end do
    else
      ! Level 1 is the top: the snow flux leaving each layer falls into the next.
      do k = 1, size(temperature)
        forms_snow = temperature(k) < snow%formation_temperature
        if (forms_snow) then
          latent_heat = latent_heat_of_condensation + snow%latent_heat_of_fusion
        else
          latent_heat = latent_heat_of_condensation
        end if
        condensate = saturation_condensate(config, temperature(k), max(humidity(k), 0.0_real64), &
                                           full_level_pressure(k), latent_heat)
        layer_mass_rate = delta_pressure(k)/(earth_gravity*interval)
        humidity_tendency(k) = -condensate/interval
        temperature_tendency(k) = latent_heat/dry_air_specific_heat*condensate/interval
        precipitation = precipitation + condensate*layer_mass_rate
        if (forms_snow) falling_snow = falling_snow + condensate*layer_mass_rate
        ! The sensible heat above the melting temperature, after the layer's own condensation.
        warmed_temperature = temperature(k) + latent_heat/dry_air_specific_heat*condensate
        melting_capacity = dry_air_specific_heat*max(0.0_real64, warmed_temperature - &
          snow%atmospheric_melting_temperature)*layer_mass_rate/snow%latent_heat_of_fusion
        melted_snow = min(falling_snow, melting_capacity)
        if (melted_snow > 0.0_real64) then
          falling_snow = falling_snow - melted_snow
          column_melt = column_melt + melted_snow
          temperature_tendency(k) = temperature_tendency(k) - &
            snow%latent_heat_of_fusion/dry_air_specific_heat*melted_snow/(layer_mass_rate*interval)
        end if
      end do
    end if
    if (present(snowfall)) snowfall = falling_snow
    if (present(snow_melt)) snow_melt = column_melt
  end subroutine large_scale_condensation_from_levels

end module large_scale_condensation
