module shallow_water_stepper
  use iso_fortran_env, only: real64
  use harmonics, only: harmonic_transform
  use shallow_water_state, only: shallow_water_state_type, shallow_water_tendency_type, &
                                 allocate_shallow_water_state, allocate_shallow_water_tendency, &
                                 copy_shallow_water_state, swap_shallow_water_states, &
                                 enforce_shallow_water_constraints
  use shallow_water_tendency_evaluator, only: evaluate_shallow_water_tendency
  use shallow_water_gravity_wave_operator, only: solve_shallow_water_gravity_wave
  use spectral_hyperdiffusion, only: apply_spectral_hyperdiffusion
  use raw_filter, only: apply_raw_filter
  use numerics_config, only: model_numerics_config
  use shallow_water_config, only: shallow_water_equation_config
  implicit none
  private

  !> Half-step/midpoint start-up followed by semi-implicit RAW-filtered leapfrog
  !> steps.  The work storage is allocated once and swapped with the model state.
  type, public :: shallow_water_stepper_type
    private
    integer :: truncation = -1
    type(shallow_water_tendency_type) :: rhs
    type(shallow_water_state_type) :: half, candidate, next, filtered
  contains
    procedure, public :: initialize => initialize_shallow_water_stepper
    procedure, public :: advance => advance_shallow_water
  end type shallow_water_stepper_type

contains

  subroutine initialize_shallow_water_stepper(this, truncation)
    class(shallow_water_stepper_type), intent(inout) :: this
    integer, intent(in) :: truncation

    this%truncation = truncation
    call allocate_shallow_water_tendency(this%rhs, truncation)
    call allocate_shallow_water_state(this%half, truncation)
    call allocate_shallow_water_state(this%candidate, truncation)
    call allocate_shallow_water_state(this%next, truncation)
    call allocate_shallow_water_state(this%filtered, truncation)
  end subroutine initialize_shallow_water_stepper

  subroutine advance_shallow_water(this, transform, numerics, equation, step_number, previous, current)
    class(shallow_water_stepper_type), intent(inout) :: this
    type(harmonic_transform), intent(inout) :: transform
    type(model_numerics_config), intent(in) :: numerics
    type(shallow_water_equation_config), intent(in) :: equation
    integer, intent(inout) :: step_number
    type(shallow_water_state_type), intent(inout) :: previous, current
    real(real64) :: time_step

    if (this%truncation /= numerics%truncation) then
      error stop 'shallow-water stepper is not initialized for this truncation'
    end if
    time_step = numerics%time_step

    if (step_number == 0) then
      ! The two start-up calls reproduce the half-step/midpoint initialization
      ! in docs/dynamics/shallow-water-equation.md.  The half-step is not RAW-filtered.
      call integration_step(transform, numerics, equation, 0.25_real64*time_step, current, current, .false., &
                            this%rhs, this%candidate, this%half, this%filtered)
      call integration_step(transform, numerics, equation, 0.5_real64*time_step, current, this%half, .true., &
                            this%rhs, this%candidate, this%next, this%filtered)
      call copy_shallow_water_state(current, previous)
    else
      call integration_step(transform, numerics, equation, time_step, previous, current, .true., &
                            this%rhs, this%candidate, this%next, this%filtered)
      call swap_shallow_water_states(previous, this%filtered)
    end if
    call swap_shallow_water_states(current, this%next)
    step_number = step_number + 1
  end subroutine advance_shallow_water

  subroutine integration_step(transform, numerics, equation, interval, previous, current, apply_raw, &
                              rhs, candidate, next, filtered)
    type(harmonic_transform), intent(inout) :: transform
    type(model_numerics_config), intent(in) :: numerics
    type(shallow_water_equation_config), intent(in) :: equation
    real(real64), intent(in) :: interval
    type(shallow_water_state_type), intent(in) :: previous, current
    logical, intent(in) :: apply_raw
    type(shallow_water_tendency_type), intent(inout) :: rhs
    type(shallow_water_state_type), intent(inout) :: candidate, next, filtered
    real(real64) :: centered_interval, order, timescale
    integer :: truncation

    truncation = numerics%truncation
    centered_interval = 2.0_real64*interval
    order = real(numerics%hyperdiffusion%order, real64)
    timescale = numerics%hyperdiffusion%timescale_seconds
    call evaluate_shallow_water_tendency(transform, truncation, current, rhs)
    candidate%zeta(:, :) = previous%zeta + centered_interval*rhs%zeta
    call solve_shallow_water_gravity_wave(truncation, centered_interval, equation, previous, rhs, candidate)
    call apply_spectral_hyperdiffusion(truncation, centered_interval, candidate%zeta, timescale, order)
    call apply_spectral_hyperdiffusion(truncation, centered_interval, candidate%delta, timescale, order)
    call enforce_shallow_water_constraints(candidate, truncation)

    if (apply_raw) then
      call apply_raw_filter(previous%zeta, current%zeta, candidate%zeta, filtered%zeta, next%zeta, &
                            numerics%raw_filter)
      call apply_raw_filter(previous%delta, current%delta, candidate%delta, filtered%delta, next%delta, &
                            numerics%raw_filter)
      call apply_raw_filter(previous%eta, current%eta, candidate%eta, filtered%eta, next%eta, &
                            numerics%raw_filter)
    else
      call copy_shallow_water_state(current, filtered)
      call copy_shallow_water_state(candidate, next)
    end if
    call enforce_shallow_water_constraints(filtered, truncation)
    call enforce_shallow_water_constraints(next, truncation)
  end subroutine integration_step

end module shallow_water_stepper
