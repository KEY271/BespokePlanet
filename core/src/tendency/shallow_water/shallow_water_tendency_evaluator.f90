module shallow_water_tendency_evaluator
  use harmonics, only: harmonic_transform
  use shallow_water_state, only: shallow_water_state_type, shallow_water_tendency_type, &
                                 zero_shallow_water_tendency
  use shallow_water_dynamics_tendency, only: add_shallow_water_dynamics_tendency
  implicit none
  private
  public :: evaluate_shallow_water_tendency

contains

  subroutine evaluate_shallow_water_tendency(transform, truncation, state, rhs)
    type(harmonic_transform), intent(inout) :: transform
    integer, intent(in) :: truncation
    type(shallow_water_state_type), intent(in) :: state
    type(shallow_water_tendency_type), intent(inout) :: rhs
    call zero_shallow_water_tendency(rhs)
    call add_shallow_water_dynamics_tendency(transform, truncation, state, rhs)
  end subroutine evaluate_shallow_water_tendency

end module shallow_water_tendency_evaluator
