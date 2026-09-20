!> Dry convective adjustment: relaxes statically unstable columns towards a
!> neutral profile while conserving column enthalpy.
module dry_convection_tendency
  use iso_fortran_env, only: real64
  use dry_physics_config, only: convection_config
  use dry_convection, only: dry_convective_adjustment_tendency
  use dry_tendency_workspace, only: dry_workspace_type
  implicit none
  private
  public :: add_dry_convection_column_tendency, add_dry_convection_tendency

contains

  subroutine add_dry_convection_column_tendency(config, pressure_half, temperature, temperature_rhs)
    type(convection_config), intent(in) :: config
    real(real64), intent(in) :: pressure_half(0:), temperature(:)
    real(real64), intent(inout) :: temperature_rhs(:)
    real(real64) :: contribution(size(temperature))
    call dry_convective_adjustment_tendency(config, pressure_half, temperature, contribution)
    temperature_rhs = temperature_rhs + contribution
  end subroutine add_dry_convection_column_tendency

  !> Evaluated column by column on the RAW-filtered previous time level.
  subroutine add_dry_convection_tendency(config, workspace)
    type(convection_config), intent(in) :: config
    type(dry_workspace_type), intent(inout) :: workspace
    integer :: i, j

    do j = 1, workspace%ny
      do i = 1, workspace%ring_nlon(j)
        call add_dry_convection_column_tendency(config, workspace%previous_pressure_half(i, j, :), &
          workspace%previous_temperature_grid(i, j, :), workspace%forcing_temperature(i, j, :))
      end do
    end do
  end subroutine add_dry_convection_tendency

end module dry_convection_tendency
