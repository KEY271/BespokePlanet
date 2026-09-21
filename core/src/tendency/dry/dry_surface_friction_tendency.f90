!> Held-Suarez boundary-layer Rayleigh drag on the horizontal wind.
!>
!> This is an independent momentum tendency: the radiation case applies it
!> without the Held-Suarez thermal relaxation, which is why it is a module of
!> its own rather than part of dry_held_suarez_tendency.
module dry_surface_friction_tendency
  use iso_fortran_env, only: real64
  use dry_physics_config, only: surface_friction_config
  use dry_held_suarez, only: held_suarez_friction
  use dry_tendency_workspace, only: dry_workspace_type
  implicit none
  private
  public :: add_dry_surface_friction_tendency

contains

  !> Evaluated on the RAW-filtered previous time level, like every prescribed
  !> physical tendency in this model.
  subroutine add_dry_surface_friction_tendency(config, workspace)
    type(surface_friction_config), intent(in) :: config
    type(dry_workspace_type), intent(inout) :: workspace

    call apply_surface_friction(config, workspace%nx, workspace%ny, workspace%number_of_levels, &
      workspace%ring_nlon, workspace%previous_full_level_pressure, workspace%previous_ps, &
      workspace%previous_u, workspace%previous_v, workspace%forcing_u, workspace%forcing_v)
  end subroutine add_dry_surface_friction_tendency

  subroutine apply_surface_friction(config, nx, ny, levels, ring_nlon, full_level_pressure, ps, u, v, &
                                    forcing_u, forcing_v)
    type(surface_friction_config), intent(in) :: config
    integer, intent(in) :: nx, ny, levels, ring_nlon(ny)
    real(real64), intent(in) :: full_level_pressure(nx, ny, levels), ps(nx, ny)
    real(real64), intent(in) :: u(nx, ny, levels), v(nx, ny, levels)
    real(real64), intent(inout) :: forcing_u(nx, ny, levels), forcing_v(nx, ny, levels)
    integer :: i, j, k
    real(real64) :: drag_u, drag_v

    !$omp parallel do default(shared) private(i, j, k, drag_u, drag_v) schedule(static)
    do k = 1, levels
      do j = 1, ny
        do i = 1, ring_nlon(j)
          call held_suarez_friction(config, full_level_pressure(i, j, k), ps(i, j), &
                                    u(i, j, k), v(i, j, k), drag_u, drag_v)
          forcing_u(i, j, k) = forcing_u(i, j, k) + drag_u
          forcing_v(i, j, k) = forcing_v(i, j, k) + drag_v
        end do
      end do
    end do
    !$omp end parallel do
  end subroutine apply_surface_friction

end module dry_surface_friction_tendency
