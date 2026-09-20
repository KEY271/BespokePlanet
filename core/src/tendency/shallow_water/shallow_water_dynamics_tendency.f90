module shallow_water_dynamics_tendency
  use harmonics, only: harmonic_transform
  use shallow_water_state, only: shallow_water_state_type, shallow_water_tendency_type
  use shallow_water_nonlinear, only: compute_shallow_water_nonlinear_tendency
  implicit none
  private
  public :: add_shallow_water_dynamics_tendency

contains

  subroutine add_shallow_water_dynamics_tendency(transform, truncation, state, rhs)
    type(harmonic_transform), intent(inout) :: transform
    integer, intent(in) :: truncation
    type(shallow_water_state_type), intent(in) :: state
    type(shallow_water_tendency_type), intent(inout) :: rhs
    complex(kind(state%zeta)), allocatable :: zeta(:, :), delta(:, :), eta(:, :)

    call compute_shallow_water_nonlinear_tendency(transform, truncation, &
      state%zeta, state%delta, state%eta, zeta, delta, eta)
    rhs%zeta = rhs%zeta + zeta
    rhs%delta = rhs%delta + delta
    rhs%eta = rhs%eta + eta
  end subroutine add_shallow_water_dynamics_tendency

end module shallow_water_dynamics_tendency
