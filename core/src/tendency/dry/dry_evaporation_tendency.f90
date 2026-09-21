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
  use dry_physics_config, only: radiation_config, evaporation_config
  use moist_thermodynamics, only: latent_heat_of_condensation
  use surface_exchange, only: surface_evaporation_flux
  use dry_tendency_workspace, only: dry_workspace_type
  implicit none
  private
  public :: add_dry_evaporation_tendency

contains

  subroutine add_dry_evaporation_tendency(config, surface, workspace)
    type(evaporation_config), intent(in) :: config
    type(radiation_config), intent(in) :: surface
    type(dry_workspace_type), intent(inout) :: workspace
    integer :: i, j, levels
    real(real64) :: evaporation, lowest_thickness, surface_wetness

    levels = workspace%number_of_levels
    !$omp parallel do default(shared) private(i, j, evaporation, lowest_thickness, surface_wetness) schedule(dynamic, 2)
    do j = 1, workspace%ny
      do i = 1, workspace%ring_nlon(j)
        if (surface%land_sea_mixing_enabled) then
          surface_wetness = workspace%land_fraction(i, j)*config%land_surface_wetness + &
            (1.0_real64 - workspace%land_fraction(i, j))*config%ocean_surface_wetness
        else
          surface_wetness = config%surface_wetness
        end if
        evaporation = surface_evaporation_flux(surface, surface_wetness, &
          workspace%previous_pressure_half(i, j, :), workspace%previous_temperature_grid(i, j, levels), &
          workspace%previous_humidity_grid(i, j, levels), workspace%previous_surface_temperature_grid(i, j), &
          workspace%previous_u(i, j, levels), workspace%previous_v(i, j, levels))
        lowest_thickness = workspace%previous_pressure_half(i, j, levels) - &
                           workspace%previous_pressure_half(i, j, levels - 1)
        workspace%forcing_humidity(i, j, levels) = workspace%forcing_humidity(i, j, levels) + &
          earth_gravity/lowest_thickness*evaporation
        workspace%evaporation(i, j) = evaporation
        workspace%latent_heat_flux(i, j) = latent_heat_of_condensation*evaporation
      end do
    end do
    !$omp end parallel do
  end subroutine add_dry_evaporation_tendency

end module dry_evaporation_tendency
