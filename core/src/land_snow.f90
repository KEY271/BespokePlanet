!> Land snowpack of docs/tendency/snow.md: snow cover, its surface effects, and
!> the melting with the heat of the land surface above the melting temperature.
!>
!> The snowpack S is a water-equivalent mass per unit land area (kg m^-2) and
!> has no heat capacity; every temperature is in K.
module land_snow
  use iso_fortran_env, only: real64
  use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
  use dry_physics_config, only: snow_config, dry_model_physics_config
  implicit none
  private
  public :: validate_snow_config, snow_cover_fraction, advance_land_snow

contains

  !> Checks the snow coefficients and that the processes snow relies on are active.
  subroutine validate_snow_config(physics)
    type(dry_model_physics_config), intent(in) :: physics
    type(snow_config) :: config
    config = physics%snow
    if (.not. all(ieee_is_finite([config%formation_temperature, config%atmospheric_melting_temperature, &
        config%surface_melting_temperature, config%masking_water_equivalent, config%albedo_increase, &
        config%latent_heat_of_fusion, config%initial_water_equivalent]))) error stop 'nonfinite snow configuration'
    if (min(config%formation_temperature, config%atmospheric_melting_temperature, &
            config%surface_melting_temperature) <= 0.0_real64 .or. &
        config%formation_temperature > config%atmospheric_melting_temperature .or. &
        config%masking_water_equivalent <= 0.0_real64 .or. config%latent_heat_of_fusion <= 0.0_real64 .or. &
        config%initial_water_equivalent < 0.0_real64 .or. config%albedo_increase < 0.0_real64 .or. &
        physics%radiation%land_shortwave_albedo + config%albedo_increase > 1.0_real64) &
      error stop 'invalid snow configuration'
    if (.not. physics%moisture%enabled .or. .not. physics%condensation%enabled .or. &
        .not. physics%radiation%enabled .or. &
        .not. (physics%radiation%land_sea_mixing_enabled .or. physics%sea_ice%enabled)) &
      error stop 'snow requires moisture, large-scale condensation and the tiled radiative surface'
    if (physics%radiation%land_sea_mixing_enabled .and. .not. physics%bucket%enabled) &
      error stop 'snow on land requires the land bucket for its meltwater'
  end subroutine validate_snow_config

  !> f = S/(S + S_0) in [0, 1); zero for an empty (or negative) snowpack.
  pure real(real64) function snow_cover_fraction(config, snow_water) result(fraction)
    type(snow_config), intent(in) :: config
    real(real64), intent(in) :: snow_water
    if (snow_water <= 0.0_real64) then
      fraction = 0.0_real64
    else
      fraction = snow_water/(snow_water + config%masking_water_equivalent)
    end if
  end function snow_cover_fraction

  !> One physical update over `interval` of the land temperature and the snowpack.
  !> The candidate T_L* = T_L + interval*land_rhs and S* = S + interval*snowfall
  !> are formed first; if T_L* exceeds the surface melting temperature, the heat
  !> C_L (T_L* - T_m) melts at most S*, and T_L* drops by the latent heat used.
  !> On return land_rhs includes the melt cooling, snow_rhs = (S_new - S)/interval
  !> and melt is the melted water per land area and time, kg m^-2 s^-1.
  pure subroutine advance_land_snow(config, heat_capacity, interval, land_temperature, snow_water, snowfall, &
                                    land_rhs, snow_rhs, melt)
    type(snow_config), intent(in) :: config
    real(real64), intent(in) :: heat_capacity, interval, land_temperature, snow_water, snowfall
    real(real64), intent(inout) :: land_rhs
    real(real64), intent(out) :: snow_rhs, melt
    real(real64) :: candidate_temperature, candidate_snow, melted

    candidate_snow = max(snow_water, 0.0_real64) + interval*max(snowfall, 0.0_real64)
    candidate_temperature = land_temperature + interval*land_rhs
    melted = 0.0_real64
    if (candidate_snow > 0.0_real64 .and. candidate_temperature > config%surface_melting_temperature) then
      melted = min(candidate_snow, heat_capacity*(candidate_temperature - config%surface_melting_temperature)/ &
                                   config%latent_heat_of_fusion)
      land_rhs = land_rhs - config%latent_heat_of_fusion*melted/(heat_capacity*interval)
    end if
    melt = melted/interval
    snow_rhs = (candidate_snow - melted - snow_water)/interval
  end subroutine advance_land_snow

end module land_snow
