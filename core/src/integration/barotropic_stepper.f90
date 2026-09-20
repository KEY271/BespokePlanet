module barotropic_stepper
  use iso_fortran_env, only: real64
  use harmonics, only: harmonic_transform
  use barotropic_state, only: barotropic_state_type, barotropic_tendency_type, &
                              allocate_barotropic_state, allocate_barotropic_tendency, &
                              copy_barotropic_state, swap_barotropic_states, &
                              enforce_barotropic_constraints
  use barotropic_tendency_evaluator, only: evaluate_barotropic_tendency
  use raw_filter, only: apply_raw_filter
  use spectral_hyperdiffusion, only: apply_spectral_hyperdiffusion
  use numerics_config, only: model_numerics_config
  implicit none
  private

  !> Midpoint start-up followed by RAW-filtered leapfrog steps.  The work
  !> storage is allocated once and exchanged with the model state by swapping.
  type, public :: barotropic_stepper_type
    private
    integer :: truncation = -1
    type(barotropic_tendency_type) :: rhs
    type(barotropic_state_type) :: midpoint, candidate, next, filtered
  contains
    procedure, public :: initialize => initialize_barotropic_stepper
    procedure, public :: advance => advance_barotropic
  end type barotropic_stepper_type

contains

  subroutine initialize_barotropic_stepper(this, truncation)
    class(barotropic_stepper_type), intent(inout) :: this
    integer, intent(in) :: truncation

    this%truncation = truncation
    call allocate_barotropic_tendency(this%rhs, truncation)
    call allocate_barotropic_state(this%midpoint, truncation)
    call allocate_barotropic_state(this%candidate, truncation)
    call allocate_barotropic_state(this%next, truncation)
    call allocate_barotropic_state(this%filtered, truncation)
  end subroutine initialize_barotropic_stepper

  subroutine advance_barotropic(this, transform, numerics, step_number, previous_filtered, current)
    class(barotropic_stepper_type), intent(inout) :: this
    type(harmonic_transform), intent(inout) :: transform
    type(model_numerics_config), intent(in) :: numerics
    integer, intent(inout) :: step_number
    type(barotropic_state_type), intent(inout) :: previous_filtered, current
    integer :: truncation
    real(real64) :: time_step, order, timescale

    truncation = numerics%truncation
    if (this%truncation /= truncation) error stop 'barotropic stepper is not initialized for this truncation'
    time_step = numerics%time_step
    order = real(numerics%hyperdiffusion%order, real64)
    timescale = numerics%hyperdiffusion%timescale_seconds

    if (step_number == 0) then
      call evaluate_barotropic_tendency(transform, truncation, current, this%rhs)
      this%midpoint%zeta(:, :) = current%zeta + 0.5_real64*time_step*this%rhs%zeta
      call apply_spectral_hyperdiffusion(truncation, 0.5_real64*time_step, this%midpoint%zeta, timescale, order)
      call enforce_barotropic_constraints(this%midpoint%zeta, truncation)

      call evaluate_barotropic_tendency(transform, truncation, this%midpoint, this%rhs)
      this%next%zeta(:, :) = current%zeta + time_step*this%rhs%zeta
      call apply_spectral_hyperdiffusion(truncation, time_step, this%next%zeta, timescale, order)
      call copy_barotropic_state(current, previous_filtered)
    else
      call evaluate_barotropic_tendency(transform, truncation, current, this%rhs)
      this%candidate%zeta(:, :) = previous_filtered%zeta + 2.0_real64*time_step*this%rhs%zeta
      call apply_spectral_hyperdiffusion(truncation, 2.0_real64*time_step, this%candidate%zeta, timescale, order)
      call enforce_barotropic_constraints(this%candidate%zeta, truncation)
      call apply_raw_filter(previous_filtered%zeta, current%zeta, this%candidate%zeta, &
                            this%filtered%zeta, this%next%zeta, numerics%raw_filter)
      call swap_barotropic_states(previous_filtered, this%filtered)
      call enforce_barotropic_constraints(previous_filtered%zeta, truncation)
    end if

    call enforce_barotropic_constraints(this%next%zeta, truncation)
    call swap_barotropic_states(current, this%next)
    step_number = step_number + 1
  end subroutine advance_barotropic

end module barotropic_stepper
