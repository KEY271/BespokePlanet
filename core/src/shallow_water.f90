module shallow_water
  use iso_fortran_env, only: real64
  use harmonics, only: harmonic_transform
  use barotropic_vorticity, only: earth_radius
  use shallow_water_nonlinear, only: compute_shallow_water_nonlinear_tendency, &
                                     diagnose_shallow_water_velocity
  use shallow_water_gravity_wave, only: gravity_acceleration, mean_depth, &
                                        gravity_wave_implicitness, &
                                        solve_implicit_gravity_wave
  use spectral_hyperdiffusion, only: apply_spectral_hyperdiffusion
  use raw_filter, only: apply_raw_filter
  implicit none
  private

  public :: gravity_acceleration
  public :: mean_depth
  public :: gravity_wave_implicitness

  type, public :: shallow_water_solver
    private
    type(harmonic_transform) :: transform
    integer :: truncation = -1
    integer :: step_number = -1
    real(real64) :: dt = 0.0_real64
    complex(real64), allocatable :: previous_zeta(:, :), previous_delta(:, :), previous_eta(:, :)
    complex(real64), allocatable :: current_zeta(:, :), current_delta(:, :), current_eta(:, :)
  contains
    procedure, public :: init
    procedure, public :: set_initial_state
    procedure, public :: advance
    procedure, public :: get_fields
    procedure, public :: get_spectral_state
    procedure, public :: get_step
  end type shallow_water_solver

contains

  subroutine init(this, truncation, dt)
    class(shallow_water_solver), intent(inout) :: this
    integer, intent(in) :: truncation
    real(real64), intent(in) :: dt

    if (truncation < 1) error stop 'shallow_water_solver: truncation must be positive'
    if (dt <= 0.0_real64) error stop 'shallow_water_solver: dt must be positive'

    this%truncation = truncation
    this%dt = dt
    this%step_number = -1
    call this%transform%init(truncation)
    call allocate_state(this, this%previous_zeta, this%previous_delta, this%previous_eta)
    call allocate_state(this, this%current_zeta, this%current_delta, this%current_eta)
    this%previous_zeta = cmplx(0.0_real64, 0.0_real64, kind=real64)
    this%previous_delta = cmplx(0.0_real64, 0.0_real64, kind=real64)
    this%previous_eta = cmplx(0.0_real64, 0.0_real64, kind=real64)
    this%current_zeta = cmplx(0.0_real64, 0.0_real64, kind=real64)
    this%current_delta = cmplx(0.0_real64, 0.0_real64, kind=real64)
    this%current_eta = cmplx(0.0_real64, 0.0_real64, kind=real64)
  end subroutine init

  subroutine set_initial_state(this, zeta, delta, eta)
    class(shallow_water_solver), intent(inout) :: this
    complex(real64), intent(in) :: zeta(0:, 0:), delta(0:, 0:), eta(0:, 0:)

    call check_initialized(this)
    call check_state_shape(this, zeta, 'zeta')
    call check_state_shape(this, delta, 'delta')
    call check_state_shape(this, eta, 'eta')
    this%current_zeta = zeta(0:this%truncation + 1, 0:this%truncation)
    this%current_delta = delta(0:this%truncation + 1, 0:this%truncation)
    this%current_eta = eta(0:this%truncation + 1, 0:this%truncation)
    call enforce_spectral_constraints(this, this%current_zeta)
    call enforce_spectral_constraints(this, this%current_delta)
    call enforce_spectral_constraints(this, this%current_eta)
    this%previous_zeta = this%current_zeta
    this%previous_delta = this%current_delta
    this%previous_eta = this%current_eta
    this%step_number = 0
  end subroutine set_initial_state

  subroutine advance(this)
    class(shallow_water_solver), intent(inout) :: this
    complex(real64), allocatable :: next_zeta(:, :), next_delta(:, :), next_eta(:, :)
    complex(real64), allocatable :: filtered_zeta(:, :), filtered_delta(:, :), filtered_eta(:, :)
    complex(real64), allocatable :: half_zeta(:, :), half_delta(:, :), half_eta(:, :)

    call check_ready(this)
    if (this%step_number == 0) then
      ! The two start-up calls reproduce the half-step/midpoint initialization
      ! in docs/shallow-water-equation.md.  The half-step is not RAW-filtered.
      call integration_step(this, 0.25_real64*this%dt, &
                            this%current_zeta, this%current_delta, this%current_eta, &
                            this%current_zeta, this%current_delta, this%current_eta, .false., &
                            half_zeta, half_delta, half_eta, &
                            filtered_zeta, filtered_delta, filtered_eta)
      call integration_step(this, 0.5_real64*this%dt, &
                            this%current_zeta, this%current_delta, this%current_eta, &
                            half_zeta, half_delta, half_eta, .true., &
                            next_zeta, next_delta, next_eta, &
                            filtered_zeta, filtered_delta, filtered_eta)
      this%previous_zeta = this%current_zeta
      this%previous_delta = this%current_delta
      this%previous_eta = this%current_eta
    else
      call integration_step(this, this%dt, &
                            this%previous_zeta, this%previous_delta, this%previous_eta, &
                            this%current_zeta, this%current_delta, this%current_eta, .true., &
                            next_zeta, next_delta, next_eta, &
                            filtered_zeta, filtered_delta, filtered_eta)
      this%previous_zeta = filtered_zeta
      this%previous_delta = filtered_delta
      this%previous_eta = filtered_eta
    end if
    this%current_zeta = next_zeta
    this%current_delta = next_delta
    this%current_eta = next_eta
    this%step_number = this%step_number + 1
  end subroutine advance

  subroutine integration_step(this, interval, previous_zeta, previous_delta, previous_eta, &
                              current_zeta, current_delta, current_eta, apply_raw, &
                              next_zeta, next_delta, next_eta, &
                              filtered_zeta, filtered_delta, filtered_eta)
    class(shallow_water_solver), intent(inout) :: this
    real(real64), intent(in) :: interval
    complex(real64), intent(in) :: previous_zeta(0:, 0:), previous_delta(0:, 0:), previous_eta(0:, 0:)
    complex(real64), intent(in) :: current_zeta(0:, 0:), current_delta(0:, 0:), current_eta(0:, 0:)
    logical, intent(in) :: apply_raw
    complex(real64), allocatable, intent(out) :: next_zeta(:, :), next_delta(:, :), next_eta(:, :)
    complex(real64), allocatable, intent(out) :: filtered_zeta(:, :), filtered_delta(:, :), filtered_eta(:, :)
    complex(real64), allocatable :: rhs_zeta(:, :), rhs_delta(:, :), rhs_eta(:, :)
    complex(real64), allocatable :: candidate_zeta(:, :), candidate_delta(:, :), candidate_eta(:, :)
    integer :: T
    real(real64) :: centered_interval

    T = this%truncation
    centered_interval = 2.0_real64*interval
    call compute_shallow_water_nonlinear_tendency(this%transform, T, &
                                                  current_zeta, current_delta, current_eta, &
                                                  rhs_zeta, rhs_delta, rhs_eta)
    call allocate_state(this, candidate_zeta, candidate_delta, candidate_eta)
    call allocate_state(this, next_zeta, next_delta, next_eta)
    call allocate_state(this, filtered_zeta, filtered_delta, filtered_eta)
    candidate_zeta = previous_zeta + centered_interval*rhs_zeta
    call solve_implicit_gravity_wave(T, centered_interval, previous_delta, previous_eta, &
                                     rhs_delta, rhs_eta, candidate_delta, candidate_eta)
    call apply_spectral_hyperdiffusion(T, centered_interval, candidate_zeta)
    call apply_spectral_hyperdiffusion(T, centered_interval, candidate_delta)
    call enforce_spectral_constraints(this, candidate_zeta)
    call enforce_spectral_constraints(this, candidate_delta)
    call enforce_spectral_constraints(this, candidate_eta)

    if (apply_raw) then
      call apply_raw_filter(previous_zeta, current_zeta, candidate_zeta, filtered_zeta, next_zeta)
      call apply_raw_filter(previous_delta, current_delta, candidate_delta, filtered_delta, next_delta)
      call apply_raw_filter(previous_eta, current_eta, candidate_eta, filtered_eta, next_eta)
    else
      filtered_zeta = current_zeta
      filtered_delta = current_delta
      filtered_eta = current_eta
      next_zeta = candidate_zeta
      next_delta = candidate_delta
      next_eta = candidate_eta
    end if
    call enforce_spectral_constraints(this, filtered_zeta)
    call enforce_spectral_constraints(this, filtered_delta)
    call enforce_spectral_constraints(this, filtered_eta)
    call enforce_spectral_constraints(this, next_zeta)
    call enforce_spectral_constraints(this, next_delta)
    call enforce_spectral_constraints(this, next_eta)
  end subroutine integration_step

  subroutine get_fields(this, zeta, delta, eta, u, v, cfl)
    class(shallow_water_solver), intent(inout) :: this
    real(real64), allocatable, intent(out) :: zeta(:, :), delta(:, :), eta(:, :), u(:, :), v(:, :)
    real(real64), intent(out), optional :: cfl
    integer, allocatable :: nlon(:)
    integer :: j
    real(real64) :: maximum_speed

    call check_ready(this)
    call this%transform%spectral_to_grid(this%current_zeta, zeta)
    call this%transform%spectral_to_grid(this%current_delta, delta)
    call this%transform%spectral_to_grid(this%current_eta, eta)
    call diagnose_shallow_water_velocity(this%transform, this%truncation, &
                                         this%current_zeta, this%current_delta, u, v)
    if (present(cfl)) then
      maximum_speed = 0.0_real64
      nlon = this%transform%get_nlon()
      do j = 1, size(nlon)
        maximum_speed = max(maximum_speed, &
          maxval(sqrt(u(1:nlon(j), j)**2 + v(1:nlon(j), j)**2)))
      end do
      cfl = maximum_speed*this%dt/earth_radius* &
            sqrt(real(this%truncation*(this%truncation + 1), real64))
    end if
  end subroutine get_fields

  subroutine get_spectral_state(this, zeta, delta, eta)
    class(shallow_water_solver), intent(in) :: this
    complex(real64), allocatable, intent(out) :: zeta(:, :), delta(:, :), eta(:, :)

    call check_ready(this)
    zeta = this%current_zeta
    delta = this%current_delta
    eta = this%current_eta
  end subroutine get_spectral_state

  integer function get_step(this) result(step)
    class(shallow_water_solver), intent(in) :: this
    step = this%step_number
  end function get_step

  subroutine allocate_state(this, zeta, delta, eta)
    class(shallow_water_solver), intent(in) :: this
    complex(real64), allocatable, intent(out) :: zeta(:, :), delta(:, :), eta(:, :)
    integer :: T

    T = this%truncation
    allocate (zeta(0:T + 1, 0:T), delta(0:T + 1, 0:T), eta(0:T + 1, 0:T))
  end subroutine allocate_state

  subroutine enforce_spectral_constraints(this, field)
    class(shallow_water_solver), intent(in) :: this
    complex(real64), intent(inout) :: field(0:, 0:)
    integer :: T, n, m

    T = this%truncation
    field(0, 0) = cmplx(0.0_real64, 0.0_real64, kind=real64)
    field(T + 1, :) = cmplx(0.0_real64, 0.0_real64, kind=real64)
    do m = 1, T
      do n = 0, m - 1
        field(n, m) = cmplx(0.0_real64, 0.0_real64, kind=real64)
      end do
    end do
  end subroutine enforce_spectral_constraints

  subroutine check_initialized(this)
    class(shallow_water_solver), intent(in) :: this
    if (this%truncation < 1 .or. .not. allocated(this%current_zeta)) then
      error stop 'shallow_water_solver: call init before use'
    end if
  end subroutine check_initialized

  subroutine check_ready(this)
    class(shallow_water_solver), intent(in) :: this
    call check_initialized(this)
    if (this%step_number < 0) error stop 'shallow_water_solver: set initial state before integration'
  end subroutine check_ready

  subroutine check_state_shape(this, field, name)
    class(shallow_water_solver), intent(in) :: this
    complex(real64), intent(in) :: field(0:, 0:)
    character(*), intent(in) :: name
    if (ubound(field, 1) < this%truncation + 1 .or. ubound(field, 2) < this%truncation) then
      error stop 'shallow_water_solver: initial '//name//' has an inconsistent shape'
    end if
  end subroutine check_state_shape

end module shallow_water
