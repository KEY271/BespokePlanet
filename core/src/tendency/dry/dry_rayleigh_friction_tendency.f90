!> Rayleigh damping of the horizontal wind in the topmost model levels, used as
!> a sponge layer by the radiation case.
module dry_rayleigh_friction_tendency
  use iso_fortran_env, only: real64
  use dry_physics_config, only: rayleigh_friction_config
  use dry_tendency_workspace, only: dry_workspace_type
  implicit none
  private
  public :: dry_rayleigh_friction_rate, add_dry_rayleigh_friction, add_dry_rayleigh_friction_tendency

contains

  !> Damping rate at one model level: the configured rate within the sponge, zero below it.
  pure real(real64) function dry_rayleigh_friction_rate(config, level) result(rate)
    type(rayleigh_friction_config), intent(in) :: config
    integer, intent(in) :: level

    rate = 0.0_real64
    if (level >= 1 .and. level <= config%top_levels) rate = config%rate
  end function dry_rayleigh_friction_rate

  pure subroutine add_dry_rayleigh_friction(rate, u, v, momentum_u_rhs, momentum_v_rhs)
    real(real64), intent(in) :: rate, u, v
    real(real64), intent(inout) :: momentum_u_rhs, momentum_v_rhs
    momentum_u_rhs = momentum_u_rhs - rate*u
    momentum_v_rhs = momentum_v_rhs - rate*v
  end subroutine add_dry_rayleigh_friction

  !> Evaluated on the RAW-filtered previous time level.  The damping rate is
  !> zero below the sponge layer, so the lower levels are untouched.
  subroutine add_dry_rayleigh_friction_tendency(config, workspace)
    type(rayleigh_friction_config), intent(in) :: config
    type(dry_workspace_type), intent(inout) :: workspace

    call apply_rayleigh_friction(config, workspace%nx, workspace%ny, workspace%number_of_levels, &
      workspace%ring_nlon, workspace%previous_u, workspace%previous_v, &
      workspace%forcing_u, workspace%forcing_v)
  end subroutine add_dry_rayleigh_friction_tendency

  subroutine apply_rayleigh_friction(config, nx, ny, levels, ring_nlon, u, v, forcing_u, forcing_v)
    type(rayleigh_friction_config), intent(in) :: config
    integer, intent(in) :: nx, ny, levels, ring_nlon(ny)
    real(real64), intent(in) :: u(nx, ny, levels), v(nx, ny, levels)
    real(real64), intent(inout) :: forcing_u(nx, ny, levels), forcing_v(nx, ny, levels)
    integer :: i, j, k
    real(real64) :: rate

    do k = 1, levels
      rate = dry_rayleigh_friction_rate(config, k)
      do j = 1, ny
        do i = 1, ring_nlon(j)
          call add_dry_rayleigh_friction(rate, u(i, j, k), v(i, j, k), &
                                         forcing_u(i, j, k), forcing_v(i, j, k))
        end do
      end do
    end do
  end subroutine apply_rayleigh_friction

end module dry_rayleigh_friction_tendency
