!> Radiative heating of the atmosphere together with the selected ground or
!> slab-ocean surface energy budget, including the sensible heat exchange and
!> the latent heat flux recorded by the evaporation tendency.
!>
!> The column fluxes that the case diagnostics aggregate are recorded in the
!> workspace; dry_tendency_diagnostics assembles the sample once every tendency
!> has been added.  This module accumulates nothing over time and writes no files.
module dry_radiation_tendency
  use iso_fortran_env, only: real64
  use dry_physics_config, only: radiation_config, sea_ice_config, snow_config, radiation_scheme_band
  use cloud_diagnostics, only: cloud_layers
  use band_radiation, only: band_column_fluxes, band_longwave_subbands
  use surface_tiles, only: tiled_surface_tendency
  use sea_ice, only: sea_ice_budget, solve_ice_surface
  use dry_vertical_coordinate, only: dry_air_kappa
  use surface_exchange, only: lowest_full_level_pressure, surface_transfer_mass_flux
  use moist_thermodynamics, only: latent_heat_of_condensation
  use dry_radiation, only: radiation_tendency, reference_layer_humidity, radiation_downward_column, &
                           radiation_band_downward_column
  use dry_tendency_workspace, only: dry_workspace_type
  implicit none
  private
  public :: add_dry_radiation_tendency, diagnose_surface_tiles

contains

  !> Radiation and surface exchange use one consistent RAW-filtered previous-time column.
  !> With prognostic water vapour the longwave optical depth sees the previous-time
  !> specific humidity and the shortwave reflection the cloud cover diagnosed by the
  !> convection tendency of this evaluation (zero unless the cloud diagnosis is
  !> enabled); otherwise the fixed reference humidity stands in for the water vapour.
  !> With snow enabled the tiled surface also advances the land snowpack
  !> (docs/tendency/snow.md); its tendency and melt are recorded in the workspace.
  subroutine add_dry_radiation_tendency(config, ice, snow, interval, moisture_enabled, transform_mu, workspace)
    type(radiation_config), intent(in) :: config
    type(sea_ice_config), intent(in) :: ice
    type(snow_config), intent(in) :: snow
    real(real64), intent(in) :: interval
    logical, intent(in) :: moisture_enabled
    real(real64), intent(in) :: transform_mu(:)
    type(dry_workspace_type), intent(inout) :: workspace
    integer :: i, j, levels, k
    real(real64) :: longitude, incoming_shortwave, reflected_shortwave, outgoing_longwave
    real(real64) :: temperature_contribution(workspace%number_of_levels)
    real(real64) :: surface_contribution, deep_contribution
    type(sea_ice_budget) :: budget
    real(real64) :: humidity(workspace%number_of_levels)
    type(cloud_layers) :: clouds
    type(band_column_fluxes) :: fluxes

    levels = workspace%number_of_levels
    !$omp parallel do default(shared) schedule(dynamic, 2) &
    !$omp   private(i, j, k, longitude, incoming_shortwave, reflected_shortwave, outgoing_longwave) &
    !$omp   private(temperature_contribution, surface_contribution, deep_contribution, budget, humidity) &
    !$omp   private(clouds, fluxes)
    do j = 1, workspace%ny
      do i = 1, workspace%ring_nlon(j)
        longitude = 2.0_real64*acos(-1.0_real64)*real(i - 1, real64)/real(workspace%ring_nlon(j), real64)
        clouds = column_cloud_layers(workspace, i, j)
        fluxes = band_column_fluxes()
        if (workspace%surface_tiles_enabled) then
          ! Dry columns use the same prescribed humidity profile as the legacy radiation path.
          if (moisture_enabled) then
            humidity = workspace%previous_humidity_grid(i, j, :)
          else
            do k = 1, levels
              humidity(k) = reference_layer_humidity(config, workspace%previous_pressure_half(i, j, k - 1), &
                workspace%previous_pressure_half(i, j, k), workspace%previous_pressure_half(i, j, levels))
            end do
          end if
          call tiled_surface_tendency(config, ice, workspace%previous_pressure_half(i, j, :), &
            workspace%previous_temperature_grid(i, j, :), workspace%previous_u(i, j, levels), &
            workspace%previous_v(i, j, levels), transform_mu(j), longitude, workspace%evaluation_time, interval, &
            workspace%land_fraction(i, j), workspace%previous_land_temperature(i, j), &
            workspace%previous_deep_temperature_grid(i, j), workspace%previous_ocean_temperature(i, j), &
            workspace%previous_sea_ice_fraction(i, j), workspace%previous_sea_ice_volume(i, j), &
            latent_heat_of_condensation*workspace%land_evaporation(i, j), &
            latent_heat_of_condensation*workspace%ocean_evaporation(i, j), workspace%cloud_cover(i, j), &
            temperature_contribution, workspace%forcing_land_temperature(i, j), &
            workspace%forcing_deep_temperature(i, j), workspace%forcing_ocean_temperature(i, j), &
            workspace%forcing_sea_ice_fraction(i, j), workspace%forcing_sea_ice_volume(i, j), &
            workspace%sea_ice_temperature(i, j), workspace%incoming_shortwave(i, j), &
            workspace%reflected_shortwave(i, j), workspace%outgoing_longwave(i, j), budget, humidity, &
            workspace%ice_surface_residual(i, j), workspace%ocean_q_flux(i, j), snow, &
            workspace%previous_snow_water(i, j), workspace%snowfall(i, j), &
            workspace%forcing_snow_water(i, j), workspace%snow_melt(i, j), clouds, fluxes)
          call record_band_fluxes(config, fluxes, workspace, i, j)
          workspace%ice_energy_residual(i, j) = budget%energy_residual
          ! Pure ocean keeps no snowpack, so only land points have a snow budget.
          if (snow%enabled .and. workspace%land_fraction(i, j) > 0.0_real64) then
            workspace%snow_budget_residual(i, j) = workspace%forcing_snow_water(i, j) - &
              (workspace%snowfall(i, j) - workspace%snow_melt(i, j))
          end if
          workspace%forcing_temperature(i, j, :) = workspace%forcing_temperature(i, j, :) + temperature_contribution
          cycle
        end if
        if (moisture_enabled) then
          if (config%land_sea_mixing_enabled) then
            call radiation_tendency(config, workspace%previous_pressure_half(i, j, :), &
              workspace%previous_temperature_grid(i, j, :), &
              workspace%previous_surface_temperature_grid(i, j), &
              workspace%previous_deep_temperature_grid(i, j), &
              workspace%previous_u(i, j, levels), workspace%previous_v(i, j, levels), transform_mu(j), &
              longitude, workspace%evaluation_time, temperature_contribution, surface_contribution, &
              deep_contribution, incoming_shortwave, reflected_shortwave, outgoing_longwave, &
              latent_heat_flux=workspace%latent_heat_flux(i, j), &
              specific_humidity=workspace%previous_humidity_grid(i, j, :), &
              cloud_cover=workspace%cloud_cover(i, j), land_fraction=workspace%land_fraction(i, j), &
              clouds=clouds, band_fluxes=fluxes)
          else
            call radiation_tendency(config, workspace%previous_pressure_half(i, j, :), &
            workspace%previous_temperature_grid(i, j, :), &
            workspace%previous_surface_temperature_grid(i, j), &
            workspace%previous_deep_temperature_grid(i, j), &
            workspace%previous_u(i, j, levels), workspace%previous_v(i, j, levels), transform_mu(j), &
            longitude, workspace%evaluation_time, temperature_contribution, surface_contribution, &
            deep_contribution, incoming_shortwave, reflected_shortwave, outgoing_longwave, &
            latent_heat_flux=workspace%latent_heat_flux(i, j), &
            specific_humidity=workspace%previous_humidity_grid(i, j, :), cloud_cover=workspace%cloud_cover(i, j), &
            clouds=clouds, band_fluxes=fluxes)
          end if
        else
          if (config%land_sea_mixing_enabled) then
            call radiation_tendency(config, workspace%previous_pressure_half(i, j, :), &
              workspace%previous_temperature_grid(i, j, :), &
              workspace%previous_surface_temperature_grid(i, j), &
              workspace%previous_deep_temperature_grid(i, j), &
              workspace%previous_u(i, j, levels), workspace%previous_v(i, j, levels), transform_mu(j), &
              longitude, workspace%evaluation_time, temperature_contribution, surface_contribution, &
              deep_contribution, incoming_shortwave, reflected_shortwave, outgoing_longwave, &
              latent_heat_flux=workspace%latent_heat_flux(i, j), &
              land_fraction=workspace%land_fraction(i, j), band_fluxes=fluxes)
          else
            call radiation_tendency(config, workspace%previous_pressure_half(i, j, :), &
            workspace%previous_temperature_grid(i, j, :), &
            workspace%previous_surface_temperature_grid(i, j), &
            workspace%previous_deep_temperature_grid(i, j), &
            workspace%previous_u(i, j, levels), workspace%previous_v(i, j, levels), transform_mu(j), &
            longitude, workspace%evaluation_time, temperature_contribution, surface_contribution, &
            deep_contribution, incoming_shortwave, reflected_shortwave, outgoing_longwave, &
            latent_heat_flux=workspace%latent_heat_flux(i, j), band_fluxes=fluxes)
          end if
        end if
        workspace%forcing_temperature(i, j, :) = workspace%forcing_temperature(i, j, :) + temperature_contribution
        workspace%forcing_surface_temperature(i, j) = workspace%forcing_surface_temperature(i, j) + &
                                                      surface_contribution
        workspace%forcing_deep_temperature(i, j) = workspace%forcing_deep_temperature(i, j) + deep_contribution
        workspace%incoming_shortwave(i, j) = incoming_shortwave
        workspace%reflected_shortwave(i, j) = reflected_shortwave
        workspace%outgoing_longwave(i, j) = outgoing_longwave
        call record_band_fluxes(config, fluxes, workspace, i, j)
      end do
    end do
    !$omp end parallel do
  end subroutine add_dry_radiation_tendency

  !> Cloud sub-columns of one grid point from the latest cloud diagnosis.
  pure function column_cloud_layers(workspace, i, j) result(clouds)
    type(dry_workspace_type), intent(in) :: workspace
    integer, intent(in) :: i, j
    type(cloud_layers) :: clouds

    clouds%large_scale_fraction = workspace%large_scale_cloud_fraction(i, j)
    clouds%convective_fraction = workspace%convective_cloud_fraction(i, j)
    clouds%large_scale_level = workspace%large_scale_cloud_level(i, j)
    clouds%convective_level = workspace%convective_cloud_level(i, j)
  end function column_cloud_layers

  !> Keeps the band-radiation column fluxes of one grid point for the diagnostics.
  subroutine record_band_fluxes(config, fluxes, workspace, i, j)
    type(radiation_config), intent(in) :: config
    type(band_column_fluxes), intent(in) :: fluxes
    type(dry_workspace_type), intent(inout) :: workspace
    integer, intent(in) :: i, j

    if (config%scheme /= radiation_scheme_band) return
    workspace%clear_reflected_shortwave(i, j) = fluxes%clear_reflected_shortwave
    workspace%clear_outgoing_longwave(i, j) = fluxes%clear_outgoing_longwave
    workspace%window_outgoing_longwave(i, j) = fluxes%window_outgoing_longwave
    workspace%surface_incident_shortwave(i, j) = fluxes%surface_incident_shortwave
    workspace%atmospheric_shortwave_absorption(i, j) = fluxes%atmospheric_shortwave_absorption
    workspace%surface_downward_longwave(i, j) = fluxes%surface_downward_longwave
    workspace%surface_upward_longwave(i, j) = fluxes%surface_upward_longwave
  end subroutine record_band_fluxes

  !> Current-state temperatures for state output; the actual coupling fluxes stay untouched.
  subroutine diagnose_surface_tiles(config, ice, moist, mu, workspace)
    type(radiation_config), intent(in) :: config
    type(sea_ice_config), intent(in) :: ice
    logical, intent(in) :: moist
    real(real64), intent(in) :: mu(:)
    type(dry_workspace_type), intent(inout) :: workspace
    real(real64) :: transmission(workspace%number_of_levels), emission(workspace%number_of_levels)
    real(real64) :: downward(0:workspace%number_of_levels), sw_rhs(workspace%number_of_levels)
    real(real64) :: humidity(workspace%number_of_levels), sw, incoming, longitude, p, exchange, ta, qc, melt
    real(real64) :: land, area, land_albedo, surface_albedo
    real(real64) :: band_transmission(band_longwave_subbands, workspace%number_of_levels)
    real(real64) :: band_source(band_longwave_subbands, workspace%number_of_levels)
    type(band_column_fluxes) :: fluxes
    integer :: i, j, k, levels
    if (.not. workspace%surface_tiles_enabled) return
    levels = workspace%number_of_levels
    workspace%sea_ice_temperature = 0.0_real64
    workspace%sea_ice_thickness = 0.0_real64
    workspace%surface_temperature_grid = 0.0_real64
    !$omp parallel do default(shared) schedule(dynamic, 2) &
    !$omp private(i, j, k, transmission, emission, downward, sw_rhs, humidity, sw, incoming, longitude, &
    !$omp         p, exchange, ta, qc, melt, land, area, land_albedo, surface_albedo, band_transmission, &
    !$omp         band_source, fluxes)
    do j = 1, workspace%ny
      do i = 1, workspace%ring_nlon(j)
        land = workspace%land_fraction(i, j)
        area = workspace%sea_ice_fraction(i, j)
        if (ice%enabled .and. area > 0.0_real64 .and. land < 1.0_real64) then
          longitude = 2.0_real64*acos(-1.0_real64)*real(i - 1, real64)/real(workspace%ring_nlon(j), real64)
          if (moist) then
            humidity = workspace%humidity_grid(i, j, :)
          else
            do k = 1, levels
              humidity(k) = reference_layer_humidity(config, workspace%pressure_half(i, j, k - 1), &
                workspace%pressure_half(i, j, k), workspace%pressure_half(i, j, levels))
            end do
          end if
          if (config%scheme == radiation_scheme_band) then
            land_albedo = config%land_shortwave_albedo
            if (workspace%snow_enabled .and. workspace%snow_water(i, j) > 0.0_real64) land_albedo = land_albedo + &
              workspace%snow_albedo_increase*workspace%snow_water(i, j)/ &
              (workspace%snow_water(i, j) + workspace%snow_masking_water_equivalent)
            surface_albedo = land*land_albedo + (1.0_real64 - land)*((1.0_real64 - area)*config%ocean_shortwave_albedo + &
              area*ice%albedo)
            call radiation_band_downward_column(config, workspace%pressure_half(i, j, :), &
              workspace%temperature_grid(i, j, :), mu(j), longitude, workspace%evaluation_time, surface_albedo, &
              column_cloud_layers(workspace, i, j), band_transmission, band_source, downward, sw_rhs, fluxes, humidity)
            sw = fluxes%surface_incident_shortwave
          else
            call radiation_downward_column(config, workspace%pressure_half(i, j, :), &
              workspace%temperature_grid(i, j, :), mu(j), longitude, workspace%evaluation_time, &
              transmission, emission, downward, sw, incoming, sw_rhs, humidity)
            sw = sw*(1.0_real64 - workspace%cloud_cover(i, j)*config%cloud_shortwave_albedo)
          end if
          p = lowest_full_level_pressure(workspace%pressure_half(i, j, :))
          exchange = surface_transfer_mass_flux(config, p, workspace%temperature_grid(i, j, levels), &
            workspace%u(i, j, levels), workspace%v(i, j, levels))*config%dry_air_specific_heat
          ta = workspace%temperature_grid(i, j, levels)*(workspace%pressure_half(i, j, levels)/p)**dry_air_kappa
          workspace%sea_ice_thickness(i, j) = workspace%sea_ice_volume(i, j)/area
          call solve_ice_surface(ice, workspace%sea_ice_thickness(i, j), (1.0_real64 - ice%albedo)*sw + &
            downward(levels), config%stefan_boltzmann_constant, exchange, ta, workspace%sea_ice_temperature(i, j), qc, melt)
        end if
        workspace%surface_temperature_grid(i, j) = land*workspace%land_temperature(i, j) + &
          (1.0_real64 - land)*((1.0_real64 - area)*workspace%ocean_temperature(i, j) + &
                               area*workspace%sea_ice_temperature(i, j))
      end do
    end do
    !$omp end parallel do
  end subroutine diagnose_surface_tiles
end module dry_radiation_tendency
