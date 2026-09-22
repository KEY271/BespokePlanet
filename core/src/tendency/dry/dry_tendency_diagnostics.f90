!> Assembles the instantaneous diagnostic sample of one tendency evaluation
!> from the workspace, after every tendency has been added.
!>
!> State means are taken from the current state on the grid; the fluxes are
!> those the physical tendencies just evaluated on the previous time level.
!> The moist quantities are filled only when the atmosphere carries moisture.
module dry_tendency_diagnostics
  use iso_fortran_env, only: real64
  use planet_parameters, only: earth_gravity
  use dry_radiation, only: radiation_diagnostics
  use dry_tendency_workspace, only: dry_workspace_type
  implicit none
  private
  public :: collect_dry_diagnostics

contains

  subroutine collect_dry_diagnostics(transform_mu, workspace, diagnostics)
    real(real64), intent(in) :: transform_mu(:)
    type(dry_workspace_type), intent(in) :: workspace
    type(radiation_diagnostics), intent(out) :: diagnostics
    integer :: i, j, k, levels
    real(real64) :: area_weight, ring_weight, atmospheric_mass, temperature_mass_sum, kinetic_energy_mass_sum
    real(real64) :: column_water, signed_column_water, negative_column_water, speed_squared, maximum_speed_squared
    real(real64) :: land_area, ocean_area, land_weight, ocean_weight, total_precipitation
    real(real64), parameter :: degrees = 180.0_real64/acos(-1.0_real64)

    levels = workspace%number_of_levels
    diagnostics%time_seconds = workspace%evaluation_time
    diagnostics%surface_temperature = workspace%surface_temperature_grid
    diagnostics%deep_temperature = workspace%deep_temperature_grid
    diagnostics%surface_pressure = workspace%ps
    diagnostics%surface_water = workspace%surface_water
    diagnostics%surface_wetness = workspace%surface_wetness
    diagnostics%runoff = workspace%runoff
    allocate (diagnostics%zonal_temperature(workspace%ny, levels))
    allocate (diagnostics%zonal_u(workspace%ny, levels), diagnostics%zonal_v(workspace%ny, levels))
    allocate (diagnostics%zonal_uv(workspace%ny, levels), diagnostics%zonal_vt(workspace%ny, levels))
    diagnostics%zonal_temperature = 0.0_real64
    diagnostics%zonal_u = 0.0_real64
    diagnostics%zonal_v = 0.0_real64
    diagnostics%zonal_uv = 0.0_real64
    diagnostics%zonal_vt = 0.0_real64
    if (workspace%moisture_grids_ready) then
      allocate (diagnostics%precipitation(workspace%nx, workspace%ny))
      allocate (diagnostics%precipitable_water(workspace%nx, workspace%ny))
      allocate (diagnostics%zonal_humidity(workspace%ny, levels), diagnostics%zonal_vq(workspace%ny, levels))
      diagnostics%evaporation = workspace%evaporation
      diagnostics%precipitation = workspace%convective_precipitation + workspace%large_scale_precipitation
      diagnostics%cloud_cover = workspace%cloud_cover
      diagnostics%precipitable_water = 0.0_real64
      diagnostics%zonal_humidity = 0.0_real64
      diagnostics%zonal_vq = 0.0_real64
    end if
    atmospheric_mass = 0.0_real64
    temperature_mass_sum = 0.0_real64
    kinetic_energy_mass_sum = 0.0_real64
    maximum_speed_squared = -1.0_real64
    land_area = 0.0_real64
    ocean_area = 0.0_real64
    do j = 1, workspace%ny
      ring_weight = 1.0_real64/real(workspace%ring_nlon(j), real64)
      area_weight = 0.5_real64*workspace%gaussian_weights(j)*ring_weight
      do i = 1, workspace%ring_nlon(j)
        diagnostics%mean_surface_temperature = diagnostics%mean_surface_temperature + &
          area_weight*workspace%surface_temperature_grid(i, j)
        diagnostics%mean_deep_temperature = diagnostics%mean_deep_temperature + &
          area_weight*workspace%deep_temperature_grid(i, j)
        land_weight = area_weight*workspace%land_fraction(i, j)
        ocean_weight = area_weight*(1.0_real64 - workspace%land_fraction(i, j))
        land_area = land_area + land_weight
        ocean_area = ocean_area + ocean_weight
        diagnostics%mean_land_surface_temperature = diagnostics%mean_land_surface_temperature + &
          land_weight*workspace%surface_temperature_grid(i, j)
        diagnostics%mean_ocean_surface_temperature = diagnostics%mean_ocean_surface_temperature + &
          ocean_weight*workspace%surface_temperature_grid(i, j)
        diagnostics%mean_surface_pressure = diagnostics%mean_surface_pressure + area_weight*workspace%ps(i, j)
        diagnostics%mean_incoming_shortwave = diagnostics%mean_incoming_shortwave + &
          area_weight*workspace%incoming_shortwave(i, j)
        diagnostics%mean_reflected_shortwave = diagnostics%mean_reflected_shortwave + &
          area_weight*workspace%reflected_shortwave(i, j)
        diagnostics%mean_outgoing_longwave = diagnostics%mean_outgoing_longwave + &
          area_weight*workspace%outgoing_longwave(i, j)
        column_water = 0.0_real64
        signed_column_water = 0.0_real64
        negative_column_water = 0.0_real64
        do k = 1, levels
          diagnostics%zonal_temperature(j, k) = diagnostics%zonal_temperature(j, k) + &
            ring_weight*workspace%temperature_grid(i, j, k)
          diagnostics%zonal_u(j, k) = diagnostics%zonal_u(j, k) + ring_weight*workspace%u(i, j, k)
          diagnostics%zonal_v(j, k) = diagnostics%zonal_v(j, k) + ring_weight*workspace%v(i, j, k)
          diagnostics%zonal_uv(j, k) = diagnostics%zonal_uv(j, k) + &
            ring_weight*workspace%u(i, j, k)*workspace%v(i, j, k)
          diagnostics%zonal_vt(j, k) = diagnostics%zonal_vt(j, k) + &
            ring_weight*workspace%v(i, j, k)*workspace%temperature_grid(i, j, k)
          atmospheric_mass = atmospheric_mass + area_weight*workspace%delta_p(i, j, k)
          temperature_mass_sum = temperature_mass_sum + &
            area_weight*workspace%delta_p(i, j, k)*workspace%temperature_grid(i, j, k)
          speed_squared = workspace%u(i, j, k)**2 + workspace%v(i, j, k)**2
          kinetic_energy_mass_sum = kinetic_energy_mass_sum + &
            area_weight*workspace%delta_p(i, j, k)*0.5_real64*speed_squared
          ! The first grid point reaching the maximum in scan order is the one recorded.
          if (speed_squared > maximum_speed_squared) then
            maximum_speed_squared = speed_squared
            diagnostics%maximum_wind_longitude_degrees = 360.0_real64*real(i - 1, real64)*ring_weight
            diagnostics%maximum_wind_latitude_degrees = asin(transform_mu(j))*degrees
            diagnostics%maximum_wind_level = k
            diagnostics%maximum_wind_eta = workspace%full_level_eta(k)
          end if
          if (workspace%moisture_grids_ready) then
            diagnostics%zonal_humidity(j, k) = diagnostics%zonal_humidity(j, k) + &
              ring_weight*workspace%humidity_grid(i, j, k)
            diagnostics%zonal_vq(j, k) = diagnostics%zonal_vq(j, k) + &
              ring_weight*workspace%v(i, j, k)*workspace%humidity_grid(i, j, k)
            column_water = column_water + max(workspace%humidity_grid(i, j, k), 0.0_real64)*workspace%delta_p(i, j, k)
            signed_column_water = signed_column_water + workspace%humidity_grid(i, j, k)*workspace%delta_p(i, j, k)
            negative_column_water = negative_column_water - &
              min(workspace%humidity_grid(i, j, k), 0.0_real64)*workspace%delta_p(i, j, k)
          end if
        end do
        if (workspace%moisture_grids_ready) then
          diagnostics%precipitable_water(i, j) = column_water/earth_gravity
          diagnostics%mean_precipitable_water = diagnostics%mean_precipitable_water + &
            area_weight*column_water/earth_gravity
          diagnostics%mean_signed_column_water = diagnostics%mean_signed_column_water + &
            area_weight*signed_column_water/earth_gravity
          diagnostics%mean_negative_column_water = diagnostics%mean_negative_column_water + &
            area_weight*negative_column_water/earth_gravity
          diagnostics%mean_convective_precipitation = diagnostics%mean_convective_precipitation + &
            area_weight*workspace%convective_precipitation(i, j)
          diagnostics%mean_large_scale_precipitation = diagnostics%mean_large_scale_precipitation + &
            area_weight*workspace%large_scale_precipitation(i, j)
          diagnostics%mean_evaporation = diagnostics%mean_evaporation + area_weight*workspace%evaporation(i, j)
          diagnostics%mean_latent_heat_flux = diagnostics%mean_latent_heat_flux + &
            area_weight*workspace%latent_heat_flux(i, j)
          diagnostics%mean_cloud_cover = diagnostics%mean_cloud_cover + area_weight*workspace%cloud_cover(i, j)
          total_precipitation = workspace%convective_precipitation(i, j) + &
            workspace%large_scale_precipitation(i, j)
          diagnostics%mean_land_precipitation = diagnostics%mean_land_precipitation + &
            land_weight*total_precipitation
          diagnostics%mean_ocean_precipitation = diagnostics%mean_ocean_precipitation + &
            ocean_weight*total_precipitation
          diagnostics%mean_land_evaporation = diagnostics%mean_land_evaporation + &
            land_weight*workspace%land_evaporation(i, j)
          diagnostics%mean_ocean_evaporation = diagnostics%mean_ocean_evaporation + &
            ocean_weight*workspace%ocean_evaporation(i, j)
          diagnostics%mean_surface_water = diagnostics%mean_surface_water + &
            land_weight*workspace%surface_water(i, j)
          diagnostics%mean_surface_wetness = diagnostics%mean_surface_wetness + &
            land_weight*workspace%surface_wetness(i, j)
          diagnostics%mean_runoff = diagnostics%mean_runoff + land_weight*workspace%runoff(i, j)
          if (workspace%surface_wetness(i, j) < workspace%bucket_dry_threshold_fraction) then
            diagnostics%dry_land_fraction = diagnostics%dry_land_fraction + land_weight
          end if
          diagnostics%mean_water_budget_residual = diagnostics%mean_water_budget_residual + &
            land_weight*workspace%water_budget_residual(i, j)
          diagnostics%maximum_water_budget_residual = max(diagnostics%maximum_water_budget_residual, &
            abs(workspace%water_budget_residual(i, j)))
        end if
      end do
    end do
    diagnostics%mean_atmospheric_temperature = temperature_mass_sum/atmospheric_mass
    diagnostics%mean_kinetic_energy = kinetic_energy_mass_sum/atmospheric_mass
    if (land_area > 0.0_real64) then
      diagnostics%mean_land_surface_temperature = diagnostics%mean_land_surface_temperature/land_area
      diagnostics%mean_land_precipitation = diagnostics%mean_land_precipitation/land_area
      diagnostics%mean_land_evaporation = diagnostics%mean_land_evaporation/land_area
      diagnostics%mean_surface_water = diagnostics%mean_surface_water/land_area
      diagnostics%mean_surface_wetness = diagnostics%mean_surface_wetness/land_area
      diagnostics%mean_runoff = diagnostics%mean_runoff/land_area
      diagnostics%dry_land_fraction = diagnostics%dry_land_fraction/land_area
      diagnostics%mean_water_budget_residual = diagnostics%mean_water_budget_residual/land_area
    end if
    if (ocean_area > 0.0_real64) then
      diagnostics%mean_ocean_surface_temperature = diagnostics%mean_ocean_surface_temperature/ocean_area
      diagnostics%mean_ocean_precipitation = diagnostics%mean_ocean_precipitation/ocean_area
      diagnostics%mean_ocean_evaporation = diagnostics%mean_ocean_evaporation/ocean_area
    end if
    diagnostics%maximum_wind_speed = sqrt(max(maximum_speed_squared, 0.0_real64))
  end subroutine collect_dry_diagnostics

end module dry_tendency_diagnostics
