!> Held-Suarez Newtonian relaxation of temperature.
!>
!> The boundary-layer drag that accompanies it in the original Held-Suarez
!> formulation is a separate tendency (dry_surface_friction_tendency), because
!> the radiation case applies the drag without this thermal relaxation.
module dry_held_suarez_tendency
  use iso_fortran_env, only: real64
  use harmonics, only: harmonic_transform
  use dry_physics_config, only: held_suarez_config
  use dry_held_suarez, only: held_suarez_thermal_relaxation
  use dry_tendency_workspace, only: dry_workspace_type
  implicit none
  private
  public :: add_dry_held_suarez_tendency

contains

  !> Evaluated on the RAW-filtered previous time level.
  subroutine add_dry_held_suarez_tendency(config, transform, workspace)
    type(held_suarez_config), intent(in) :: config
    type(harmonic_transform), intent(in) :: transform
    type(dry_workspace_type), intent(inout) :: workspace

    call apply_thermal_relaxation(config, workspace%nx, workspace%ny, workspace%number_of_levels, &
      workspace%ring_nlon, transform%mu, workspace%previous_full_level_pressure, &
      workspace%previous_ps, workspace%previous_temperature_grid, workspace%forcing_temperature)
  end subroutine add_dry_held_suarez_tendency

  subroutine apply_thermal_relaxation(config, nx, ny, levels, ring_nlon, mu, full_level_pressure, ps, &
                                      temperature, forcing_temperature)
    type(held_suarez_config), intent(in) :: config
    integer, intent(in) :: nx, ny, levels, ring_nlon(ny)
    real(real64), intent(in) :: mu(ny)
    real(real64), intent(in) :: full_level_pressure(nx, ny, levels), ps(nx, ny)
    real(real64), intent(in) :: temperature(nx, ny, levels)
    real(real64), intent(inout) :: forcing_temperature(nx, ny, levels)
    integer :: i, j, k
    real(real64) :: relaxation

    do k = 1, levels
      do j = 1, ny
        do i = 1, ring_nlon(j)
          call held_suarez_thermal_relaxation(config, mu(j), full_level_pressure(i, j, k), ps(i, j), &
                                              temperature(i, j, k), relaxation)
          forcing_temperature(i, j, k) = forcing_temperature(i, j, k) + relaxation
        end do
      end do
    end do
  end subroutine apply_thermal_relaxation

end module dry_held_suarez_tendency
