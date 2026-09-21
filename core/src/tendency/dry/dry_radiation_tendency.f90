!> Radiative heating of the atmosphere together with the selected ground or
!> slab-ocean surface energy budget, including the sensible heat exchange and
!> the latent heat flux recorded by the evaporation tendency.
!>
!> The column fluxes that the case diagnostics aggregate are recorded in the
!> workspace; dry_tendency_diagnostics assembles the sample once every tendency
!> has been added.  This module accumulates nothing over time and writes no files.
module dry_radiation_tendency
  use iso_fortran_env, only: real64
  use dry_physics_config, only: radiation_config
  use dry_radiation, only: radiation_tendency
  use dry_tendency_workspace, only: dry_workspace_type
  implicit none
  private
  public :: add_dry_radiation_tendency

contains

  !> Radiation and surface exchange use one consistent RAW-filtered previous-time column.
  subroutine add_dry_radiation_tendency(config, transform_mu, workspace)
    type(radiation_config), intent(in) :: config
    real(real64), intent(in) :: transform_mu(:)
    type(dry_workspace_type), intent(inout) :: workspace
    integer :: i, j, levels
    real(real64) :: longitude, incoming_shortwave, reflected_shortwave, outgoing_longwave
    real(real64) :: temperature_contribution(workspace%number_of_levels)
    real(real64) :: surface_contribution, deep_contribution

    levels = workspace%number_of_levels
    do j = 1, workspace%ny
      do i = 1, workspace%ring_nlon(j)
        longitude = 2.0_real64*acos(-1.0_real64)*real(i - 1, real64)/real(workspace%ring_nlon(j), real64)
        call radiation_tendency(config, workspace%previous_pressure_half(i, j, :), &
          workspace%previous_temperature_grid(i, j, :), &
          workspace%previous_surface_temperature_grid(i, j), &
          workspace%previous_deep_temperature_grid(i, j), &
          workspace%previous_u(i, j, levels), workspace%previous_v(i, j, levels), transform_mu(j), &
          longitude, workspace%evaluation_time, temperature_contribution, surface_contribution, &
          deep_contribution, incoming_shortwave, reflected_shortwave, outgoing_longwave, &
          latent_heat_flux=workspace%latent_heat_flux(i, j))
        workspace%forcing_temperature(i, j, :) = workspace%forcing_temperature(i, j, :) + temperature_contribution
        workspace%forcing_surface_temperature(i, j) = workspace%forcing_surface_temperature(i, j) + &
                                                      surface_contribution
        workspace%forcing_deep_temperature(i, j) = workspace%forcing_deep_temperature(i, j) + deep_contribution
        workspace%incoming_shortwave(i, j) = incoming_shortwave
        workspace%reflected_shortwave(i, j) = reflected_shortwave
        workspace%outgoing_longwave(i, j) = outgoing_longwave
      end do
    end do
  end subroutine add_dry_radiation_tendency

end module dry_radiation_tendency
