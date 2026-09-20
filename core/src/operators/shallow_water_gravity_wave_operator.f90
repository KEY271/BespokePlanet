module shallow_water_gravity_wave_operator
  use iso_fortran_env, only: real64
  use shallow_water_gravity_wave, only: solve_implicit_gravity_wave
  use shallow_water_state, only: shallow_water_state_type, shallow_water_tendency_type
  use shallow_water_config, only: shallow_water_equation_config
  implicit none
  private
  public :: solve_shallow_water_gravity_wave

contains

  subroutine solve_shallow_water_gravity_wave(truncation, centered_interval, equation, previous, rhs, candidate)
    integer, intent(in) :: truncation
    real(real64), intent(in) :: centered_interval
    type(shallow_water_equation_config), intent(in) :: equation
    type(shallow_water_state_type), intent(in) :: previous
    type(shallow_water_tendency_type), intent(in) :: rhs
    type(shallow_water_state_type), intent(inout) :: candidate
    call solve_implicit_gravity_wave(truncation, centered_interval, previous%delta, previous%eta, &
      rhs%delta, rhs%eta, candidate%delta, candidate%eta, equation)
  end subroutine solve_shallow_water_gravity_wave

end module shallow_water_gravity_wave_operator
