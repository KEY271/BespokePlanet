module shallow_water
  use iso_fortran_env, only: real64
  use harmonics, only: harmonic_transform
  use planet_parameters, only: earth_radius
  use spectral_vector_operators, only: diagnose_horizontal_velocity
  use shallow_water_state, only: shallow_water_state_type, allocate_shallow_water_state, &
                                 copy_shallow_water_state, enforce_shallow_water_constraints
  use shallow_water_stepper, only: shallow_water_stepper_type
  use numerics_config, only: model_numerics_config
  use shallow_water_config, only: shallow_water_equation_config
  implicit none
  private

  type, public :: shallow_water_solver
    private
    type(harmonic_transform) :: transform
    integer :: truncation = -1
    integer :: step_number = -1
    type(model_numerics_config) :: numerics
    type(shallow_water_equation_config) :: equation
    type(shallow_water_stepper_type) :: stepper
    type(shallow_water_state_type) :: previous
    type(shallow_water_state_type) :: current
  contains
    procedure, public :: init
    procedure, public :: init_with_config
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
    this%numerics = model_numerics_config()
    this%numerics%truncation = truncation
    this%numerics%time_step = dt
    this%step_number = -1
    call this%transform%init(truncation)
    call allocate_shallow_water_state(this%previous, truncation)
    call allocate_shallow_water_state(this%current, truncation)
    this%equation = shallow_water_equation_config()
    call this%stepper%initialize(truncation)
  end subroutine init

  subroutine init_with_config(this, config, equation)
    class(shallow_water_solver), intent(inout) :: this
    type(model_numerics_config), intent(in) :: config
    type(shallow_water_equation_config), intent(in), optional :: equation
    call this%init(config%truncation, config%time_step)
    this%numerics = config
    if (present(equation)) this%equation = equation
  end subroutine init_with_config

  subroutine set_initial_state(this, zeta, delta, eta)
    class(shallow_water_solver), intent(inout) :: this
    complex(real64), intent(in) :: zeta(0:, 0:), delta(0:, 0:), eta(0:, 0:)

    call check_initialized(this)
    call check_state_shape(this, zeta, 'zeta')
    call check_state_shape(this, delta, 'delta')
    call check_state_shape(this, eta, 'eta')
    this%current%zeta = zeta(0:this%truncation + 1, 0:this%truncation)
    this%current%delta = delta(0:this%truncation + 1, 0:this%truncation)
    this%current%eta = eta(0:this%truncation + 1, 0:this%truncation)
    call enforce_shallow_water_constraints(this%current, this%truncation)
    call copy_shallow_water_state(this%current, this%previous)
    this%step_number = 0
  end subroutine set_initial_state

  subroutine advance(this)
    class(shallow_water_solver), intent(inout) :: this
    call check_ready(this)
    call this%stepper%advance(this%transform, this%numerics, this%equation, this%step_number, &
                              this%previous, this%current)
  end subroutine advance

  subroutine get_fields(this, zeta, delta, eta, u, v, cfl)
    class(shallow_water_solver), intent(inout) :: this
    real(real64), allocatable, intent(out) :: zeta(:, :), delta(:, :), eta(:, :), u(:, :), v(:, :)
    real(real64), intent(out), optional :: cfl
    integer, allocatable :: nlon(:)
    integer :: j
    real(real64) :: maximum_speed

    call check_ready(this)
    call this%transform%spectral_to_grid(this%current%zeta, zeta)
    call this%transform%spectral_to_grid(this%current%delta, delta)
    call this%transform%spectral_to_grid(this%current%eta, eta)
    call diagnose_horizontal_velocity(this%transform, this%truncation, &
                                         this%current%zeta, this%current%delta, u, v)
    if (present(cfl)) then
      maximum_speed = 0.0_real64
      nlon = this%transform%get_nlon()
      do j = 1, size(nlon)
        maximum_speed = max(maximum_speed, &
          maxval(sqrt(u(1:nlon(j), j)**2 + v(1:nlon(j), j)**2)))
      end do
      cfl = maximum_speed*this%numerics%time_step/earth_radius* &
            sqrt(real(this%truncation*(this%truncation + 1), real64))
    end if
  end subroutine get_fields

  subroutine get_spectral_state(this, zeta, delta, eta)
    class(shallow_water_solver), intent(in) :: this
    complex(real64), allocatable, intent(out) :: zeta(:, :), delta(:, :), eta(:, :)

    call check_ready(this)
    zeta = this%current%zeta
    delta = this%current%delta
    eta = this%current%eta
  end subroutine get_spectral_state

  integer function get_step(this) result(step)
    class(shallow_water_solver), intent(in) :: this
    step = this%step_number
  end function get_step

  subroutine check_initialized(this)
    class(shallow_water_solver), intent(in) :: this
    if (this%truncation < 1 .or. .not. allocated(this%current%zeta)) then
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
