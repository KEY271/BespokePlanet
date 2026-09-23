!> Land, open water and diagnostic sea ice share one atmospheric column.
!> With snow enabled the land snowpack brightens the land, scales its sensible
!> heat by 1 - f and melts with the land heat above the melting temperature;
!> snow falling on the ocean melts with the heat of the mixed layer
!> (docs/tendency/snow.md).
module surface_tiles
  use iso_fortran_env, only: real64
  use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
  use dry_physics_config, only: radiation_config, sea_ice_config, snow_config
  use dry_vertical_coordinate, only: dry_air_kappa
  use dry_radiation, only: radiation_downward_column, radiation_upward_column
  use surface_exchange, only: lowest_full_level_pressure, surface_transfer_mass_flux
  use sea_ice, only: solve_ice_surface, advance_sea_ice, sea_ice_budget
  use land_snow, only: snow_cover_fraction, advance_land_snow
  implicit none
  private
  public :: tiled_surface_tendency
contains
  subroutine tiled_surface_tendency(config, ice, pressure, air, u, v, mu, longitude, time, interval, &
                                    land, tl, td, to, area, volume, land_latent, ocean_latent, cloud, &
                                    air_rhs, land_rhs, deep_rhs, ocean_rhs, area_rhs, volume_rhs, &
                                    ice_temperature, incoming, reflected, outgoing, budget, humidity, surface_residual, &
                                    ocean_heat_convergence, snow, snow_water, snowfall, snow_rhs, snow_melt)
    type(radiation_config), intent(in) :: config
    type(sea_ice_config), intent(in) :: ice
    real(real64), intent(in) :: pressure(0:), air(:), u, v, mu, longitude, time, interval
    real(real64), intent(in) :: land, tl, td, to, area, volume, land_latent, ocean_latent, cloud
    real(real64), intent(in), optional :: humidity(:)
    real(real64), intent(out), optional :: surface_residual
    !> Prescribed Q flux per ocean area (docs/tendency/q-flux.md); zero when absent.
    real(real64), intent(in), optional :: ocean_heat_convergence
    !> Snow configuration, land snowpack S (kg m^-2 of land) and the large-scale
    !> snowfall reaching the surface (kg m^-2 s^-1).  snow_rhs is the tendency of S
    !> and snow_melt the melted water per land area (kg m^-2 s^-1).
    type(snow_config), intent(in), optional :: snow
    real(real64), intent(in), optional :: snow_water, snowfall
    real(real64), intent(out), optional :: snow_rhs, snow_melt
    real(real64), intent(out) :: air_rhs(:), land_rhs, deep_rhs, ocean_rhs, area_rhs, volume_rhs
    real(real64), intent(out) :: ice_temperature, incoming, reflected, outgoing
    type(sea_ice_budget), intent(out) :: budget
    real(real64) :: transmission(size(air)), emission(size(air)), downward(0:size(air)), sw_rhs(size(air))
    real(real64) :: sw, sw_surface, lower_pressure, exchange, ta, sigma, co, wl, wo, wi
    real(real64) :: hl, ho, hi, fl, fo, ground_heat, conduction, melting, surface_upward
    real(real64) :: new_temperature, new_area, new_volume, albedo, convergence
    real(real64) :: land_albedo, cover, snow_tendency, melt
    logical :: snowing
    integer :: levels

    levels = size(air)
    if (present(surface_residual)) surface_residual = 0.0_real64
    convergence = 0.0_real64
    if (present(ocean_heat_convergence)) convergence = ocean_heat_convergence
    if (.not. ieee_is_finite(convergence)) error stop 'nonfinite ocean heat convergence'
    if (land >= 1.0_real64 .and. convergence /= 0.0_real64) error stop 'Q flux on a pure land point'
    snowing = .false.
    if (present(snow)) snowing = snow%enabled
    if (snowing .and. .not. (present(snow_water) .and. present(snowfall))) &
      error stop 'snow requires the snowpack and the snowfall'
    if (snowing) then
      if (.not. all(ieee_is_finite([snow_water, snowfall])) .or. snow_water < 0.0_real64 .or. &
          snowfall < 0.0_real64) error stop 'invalid snowpack or snowfall'
    end if
    snow_tendency = 0.0_real64
    melt = 0.0_real64
    cover = 0.0_real64
    land_albedo = config%land_shortwave_albedo
    if (snowing .and. land > 0.0_real64) then
      cover = snow_cover_fraction(snow, snow_water)
      land_albedo = land_albedo + snow%albedo_increase*cover
    end if
    ! Snow falling on open water or ice melts at once with the heat of the whole mixed layer.
    if (snowing .and. land < 1.0_real64) convergence = convergence - snow%latent_heat_of_fusion*snowfall
    if (.not. all(ieee_is_finite([land, interval, tl, td, to, area, volume, cloud, land_latent, ocean_latent]))) &
      error stop 'nonfinite surface tile input'
    if (land > 0.0_real64) then
      if (.not. all(ieee_is_finite([config%surface_heat_capacity, config%deep_ground_heat_capacity, &
          config%ground_exchange_coefficient])) .or. config%surface_heat_capacity <= 0.0_real64 .or. &
          config%deep_ground_heat_capacity <= 0.0_real64 .or. config%ground_exchange_coefficient < 0.0_real64) &
        error stop 'invalid land heat capacity or exchange'
    end if
    if (land < 1.0_real64) then
      if (.not. all(ieee_is_finite([config%seawater_density, config%seawater_specific_heat, config%slab_ocean_depth])) .or. &
          min(config%seawater_density, config%seawater_specific_heat, config%slab_ocean_depth) <= 0.0_real64) &
        error stop 'invalid ocean heat capacity'
    end if
    if (land < 0.0_real64 .or. land > 1.0_real64 .or. interval <= 0.0_real64) &
      error stop 'invalid surface tile geometry or interval'
    if (cloud < 0.0_real64 .or. cloud > 1.0_real64) error stop 'invalid cloud fraction'
    if ((land > 0.0_real64 .and. (tl <= 0.0_real64 .or. td <= 0.0_real64)) .or. &
        (land < 1.0_real64 .and. to <= 0.0_real64)) error stop 'invalid surface tile temperature'
    call radiation_downward_column(config, pressure, air, mu, longitude, time, transmission, emission, &
                                   downward, sw, incoming, sw_rhs, humidity)
    sw_surface = (1.0_real64 - cloud*config%cloud_shortwave_albedo)*sw
    lower_pressure = lowest_full_level_pressure(pressure)
    exchange = surface_transfer_mass_flux(config, lower_pressure, air(levels), u, v)*config%dry_air_specific_heat
    ta = air(levels)*(pressure(levels)/lower_pressure)**dry_air_kappa
    sigma = config%stefan_boltzmann_constant
    co = config%seawater_density*config%seawater_specific_heat*config%slab_ocean_depth
    wl = land
    wo = 1.0_real64 - land
    wi = 0.0_real64
    if (ice%enabled) then
      wi = wo*area
      wo = wo*(1.0_real64 - area)
    end if
    hl = 0.0_real64
    ho = 0.0_real64
    hi = 0.0_real64
    fo = 0.0_real64
    land_rhs = 0.0_real64
    deep_rhs = 0.0_real64
    ocean_rhs = 0.0_real64
    area_rhs = 0.0_real64
    volume_rhs = 0.0_real64
    ice_temperature = 0.0_real64
    conduction = 0.0_real64
    melting = 0.0_real64
    surface_upward = 0.0_real64
    budget = sea_ice_budget()
    if (wl > 0.0_real64) then
      hl = exchange*(tl - ta)
      if (snowing) hl = (1.0_real64 - cover)*hl
      fl = (1.0_real64 - land_albedo)*sw_surface + downward(levels) - sigma*tl**4 - hl - land_latent
      ground_heat = config%ground_exchange_coefficient*(tl - td)
      land_rhs = (fl - ground_heat)/config%surface_heat_capacity
      deep_rhs = ground_heat/config%deep_ground_heat_capacity
      if (snowing) call advance_land_snow(snow, config%surface_heat_capacity, interval, tl, snow_water, snowfall, &
                                          land_rhs, snow_tendency, melt)
      surface_upward = wl*sigma*tl**4
    end if
    if (land < 1.0_real64) then
      ho = exchange*(to - ta)
      fo = (1.0_real64 - config%ocean_shortwave_albedo)*sw_surface + downward(levels) - sigma*to**4 - ho - ocean_latent
      surface_upward = surface_upward + wo*sigma*to**4
      if (wi > 0.0_real64) then
        call solve_ice_surface(ice, volume/area, (1.0_real64 - ice%albedo)*sw_surface + downward(levels), &
                              sigma, exchange, ta, ice_temperature, conduction, melting, surface_residual)
        hi = exchange*(ice_temperature - ta)
        surface_upward = surface_upward + wi*sigma*ice_temperature**4
      end if
      if (ice%enabled) then
        call advance_sea_ice(ice, co, interval, to, area, volume, fo, conduction, melting, &
                             new_temperature, new_area, new_volume, budget, convergence)
        ocean_rhs = (new_temperature - to)/interval
        area_rhs = (new_area - area)/interval
        volume_rhs = (new_volume - volume)/interval
      else
        ocean_rhs = (fo + convergence)/co
      end if
    end if
    call radiation_upward_column(config, pressure, transmission, emission, downward, surface_upward, air_rhs, outgoing)
    air_rhs = air_rhs + sw_rhs
    air_rhs(levels) = air_rhs(levels) + config%gravity_acceleration/(config%dry_air_specific_heat* &
      (pressure(levels) - pressure(levels - 1)))*(wl*hl + wo*ho + wi*hi)
    albedo = wl*land_albedo + wo*config%ocean_shortwave_albedo + wi*ice%albedo
    reflected = cloud*config%cloud_shortwave_albedo*sw + albedo*sw_surface
    if (present(snow_rhs)) snow_rhs = snow_tendency
    if (present(snow_melt)) snow_melt = melt
  end subroutine tiled_surface_tendency
end module surface_tiles
