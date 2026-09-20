module dry_gravity_wave_operator
  use iso_fortran_env, only: real64
  use dry_gravity_wave, only: dry_gravity_wave_solver
  use dry_state, only: dry_state_type, dry_tendency_type
  implicit none
  private
  public :: solve_dry_gravity_wave

contains

  subroutine solve_dry_gravity_wave(operator, centered_interval, previous, current, rhs, candidate)
    type(dry_gravity_wave_solver), intent(in) :: operator
    real(real64), intent(in) :: centered_interval
    type(dry_state_type), intent(in) :: previous, current
    type(dry_tendency_type), intent(in) :: rhs
    type(dry_state_type), intent(inout) :: candidate
    call operator%solve(centered_interval, previous%log_surface_pressure, previous%delta, &
      previous%temperature, current%log_surface_pressure, current%delta, current%temperature, &
      rhs%log_surface_pressure, rhs%delta, rhs%temperature, candidate%log_surface_pressure, &
      candidate%delta, candidate%temperature)
  end subroutine solve_dry_gravity_wave

end module dry_gravity_wave_operator
