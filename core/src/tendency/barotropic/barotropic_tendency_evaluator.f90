module barotropic_tendency_evaluator
  use harmonics, only: harmonic_transform
  use barotropic_state, only: barotropic_state_type, barotropic_tendency_type, &
                              zero_barotropic_tendency
  use barotropic_dynamics_tendency, only: add_barotropic_dynamics_tendency
  implicit none
  private
  public :: evaluate_barotropic_tendency

contains

  subroutine evaluate_barotropic_tendency(transform, truncation, state, rhs)
    type(harmonic_transform), intent(inout) :: transform
    integer, intent(in) :: truncation
    type(barotropic_state_type), intent(in) :: state
    type(barotropic_tendency_type), intent(inout) :: rhs
    call zero_barotropic_tendency(rhs)
    call add_barotropic_dynamics_tendency(transform, truncation, state, rhs)
  end subroutine evaluate_barotropic_tendency

end module barotropic_tendency_evaluator
