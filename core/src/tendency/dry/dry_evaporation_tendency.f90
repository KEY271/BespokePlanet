!> Bulk evaporation from the surface into the lowest model level
!> (docs/tendency/evaporation.md).
!>
!> The water vapour source g E/dp_N is added to the lowest-level humidity, and
!> the latent heat flux L E is recorded in the workspace so that the radiation
!> tendency, evaluated afterwards, takes exactly the same value from the surface
!> energy budget.  Evaluated on the RAW-filtered previous time level with the
!> wind and density the sensible heat flux uses.
module dry_evaporation_tendency
  use iso_fortran_env, only: real64
  use planet_parameters, only: earth_gravity
  use dry_physics_config, only: radiation_config, evaporation_config, bucket_config
  use moist_thermodynamics, only: latent_heat_of_condensation
  use surface_exchange, only: surface_evaporation_flux
  use land_bucket, only: bucket_wetness, limit_bucket_evaporation
  use dry_tendency_workspace, only: dry_workspace_type
  implicit none
  private
  public :: add_dry_evaporation_tendency

contains

  subroutine add_dry_evaporation_tendency(config, surface, bucket, interval, workspace)
    type(evaporation_config), intent(in) :: config
    type(radiation_config), intent(in) :: surface
    type(bucket_config), intent(in) :: bucket
    real(real64), intent(in) :: interval
    type(dry_workspace_type), intent(inout) :: workspace
    integer :: i, j, levels
    real(real64) :: evaporation, land_evaporation, ocean_evaporation, potential_evaporation
    real(real64) :: lowest_thickness, surface_wetness, land_fraction

    levels = workspace%number_of_levels
    !$omp parallel do default(shared) &
    !$omp private(i, j, evaporation, land_evaporation, ocean_evaporation, potential_evaporation, &
    !$omp         lowest_thickness, surface_wetness, land_fraction) schedule(dynamic, 2)
    do j = 1, workspace%ny
      do i = 1, workspace%ring_nlon(j)
        land_fraction = workspace%land_fraction(i, j)
        if (workspace%surface_tiles_enabled) then
          land_evaporation = 0.0_real64
          ocean_evaporation = 0.0_real64
          surface_wetness = 0.0_real64
          if (land_fraction > 0.0_real64) then
            potential_evaporation = surface_evaporation_flux(surface, 1.0_real64, &
              workspace%previous_pressure_half(i, j, :), workspace%previous_temperature_grid(i, j, levels), &
              workspace%previous_humidity_grid(i, j, levels), workspace%previous_land_temperature(i, j), &
              workspace%previous_u(i, j, levels), workspace%previous_v(i, j, levels))
            if (bucket%enabled) then
              land_evaporation = limit_bucket_evaporation(potential_evaporation, &
                workspace%previous_surface_water(i, j), bucket%capacity, interval)
              surface_wetness = bucket_wetness(workspace%previous_surface_water(i, j), bucket%capacity)
            else
              surface_wetness = config%land_surface_wetness
              land_evaporation = surface_wetness*potential_evaporation
            end if
          end if
          if (land_fraction < 1.0_real64) then
            ocean_evaporation = surface_evaporation_flux(surface, config%ocean_surface_wetness, &
              workspace%previous_pressure_half(i, j, :), workspace%previous_temperature_grid(i, j, levels), &
              workspace%previous_humidity_grid(i, j, levels), workspace%previous_ocean_temperature(i, j), &
              workspace%previous_u(i, j, levels), workspace%previous_v(i, j, levels))
          end if
          evaporation = land_fraction*land_evaporation + (1.0_real64 - land_fraction)* &
            (1.0_real64 - workspace%previous_sea_ice_fraction(i, j))*ocean_evaporation
        else if (surface%land_sea_mixing_enabled .and. bucket%enabled) then
          potential_evaporation = surface_evaporation_flux(surface, 1.0_real64, &
            workspace%previous_pressure_half(i, j, :), workspace%previous_temperature_grid(i, j, levels), &
            workspace%previous_humidity_grid(i, j, levels), workspace%previous_surface_temperature_grid(i, j), &
            workspace%previous_u(i, j, levels), workspace%previous_v(i, j, levels))
          land_evaporation = limit_bucket_evaporation(potential_evaporation, &
            workspace%previous_surface_water(i, j), bucket%capacity, interval)
          ocean_evaporation = config%ocean_surface_wetness*potential_evaporation
          surface_wetness = bucket_wetness(workspace%previous_surface_water(i, j), bucket%capacity)
          evaporation = land_fraction*land_evaporation + (1.0_real64 - land_fraction)*ocean_evaporation
        else if (surface%land_sea_mixing_enabled) then
          surface_wetness = workspace%land_fraction(i, j)*config%land_surface_wetness + &
            (1.0_real64 - workspace%land_fraction(i, j))*config%ocean_surface_wetness
          evaporation = surface_evaporation_flux(surface, surface_wetness, &
            workspace%previous_pressure_half(i, j, :), workspace%previous_temperature_grid(i, j, levels), &
            workspace%previous_humidity_grid(i, j, levels), workspace%previous_surface_temperature_grid(i, j), &
            workspace%previous_u(i, j, levels), workspace%previous_v(i, j, levels))
          land_evaporation = evaporation
          ocean_evaporation = evaporation
        else
          surface_wetness = config%surface_wetness
          evaporation = surface_evaporation_flux(surface, surface_wetness, &
            workspace%previous_pressure_half(i, j, :), workspace%previous_temperature_grid(i, j, levels), &
            workspace%previous_humidity_grid(i, j, levels), workspace%previous_surface_temperature_grid(i, j), &
            workspace%previous_u(i, j, levels), workspace%previous_v(i, j, levels))
          land_evaporation = evaporation
          ocean_evaporation = evaporation
        end if
        lowest_thickness = workspace%previous_pressure_half(i, j, levels) - &
                           workspace%previous_pressure_half(i, j, levels - 1)
        workspace%forcing_humidity(i, j, levels) = workspace%forcing_humidity(i, j, levels) + &
          earth_gravity/lowest_thickness*evaporation
        workspace%evaporation(i, j) = evaporation
        workspace%land_evaporation(i, j) = land_evaporation
        workspace%ocean_evaporation(i, j) = ocean_evaporation
        workspace%latent_heat_flux(i, j) = latent_heat_of_condensation*evaporation
        workspace%surface_wetness(i, j) = surface_wetness
      end do
    end do
    !$omp end parallel do
  end subroutine add_dry_evaporation_tendency

end module dry_evaporation_tendency
