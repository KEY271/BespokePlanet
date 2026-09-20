module barotropic_vorticity
  use iso_fortran_env, only: real64
  use harmonics, only: harmonic_transform
  use planet_parameters, only: earth_radius
  use barotropic_state, only: barotropic_state_type, allocate_barotropic_state, &
                              copy_barotropic_state, enforce_barotropic_constraints
  use barotropic_stepper, only: barotropic_stepper_type
  use numerics_config, only: model_numerics_config
  implicit none
  private

  type, public :: barotropic_solver
    private
    type(harmonic_transform) :: transform
    integer :: truncation = -1
    integer :: step_number = -1
    type(model_numerics_config) :: numerics
    type(barotropic_stepper_type) :: stepper
    type(barotropic_state_type) :: previous_filtered
    type(barotropic_state_type) :: current
  contains
    procedure, public :: init
    procedure, public :: init_with_config
    procedure, public :: set_initial_vorticity
    procedure, public :: advance
    procedure, public :: get_fields
    procedure, public :: get_step
  end type barotropic_solver

contains

  subroutine init(this, truncation, dt)
    class(barotropic_solver), intent(inout) :: this
    integer, intent(in) :: truncation
    real(real64), intent(in) :: dt

    if (truncation < 1) error stop 'barotropic_solver: truncation must be positive'
    if (dt <= 0.0_real64) error stop 'barotropic_solver: dt must be positive'

    this%truncation = truncation
    this%numerics = model_numerics_config()
    this%numerics%truncation = truncation
    this%numerics%time_step = dt
    this%step_number = -1
    call this%transform%init(truncation)

    call allocate_barotropic_state(this%previous_filtered, truncation)
    call allocate_barotropic_state(this%current, truncation)
    call this%stepper%initialize(truncation)
  end subroutine init

  subroutine init_with_config(this, config)
    class(barotropic_solver), intent(inout) :: this
    type(model_numerics_config), intent(in) :: config
    call this%init(config%truncation, config%time_step)
    this%numerics = config
  end subroutine init_with_config

  subroutine set_initial_vorticity(this, zeta)
    class(barotropic_solver), intent(inout) :: this
    complex(real64), intent(in) :: zeta(0:, 0:)

    call check_initialized(this)
    if (ubound(zeta, 1) < this%truncation + 1 .or. &
        ubound(zeta, 2) < this%truncation) then
      error stop 'barotropic_solver: initial vorticity has an inconsistent shape'
    end if

    this%current%zeta = zeta(0:this%truncation + 1, 0:this%truncation)
    call enforce_barotropic_constraints(this%current%zeta, this%truncation)
    call copy_barotropic_state(this%current, this%previous_filtered)
    this%step_number = 0
  end subroutine set_initial_vorticity

  subroutine advance(this)
    class(barotropic_solver), intent(inout) :: this
    call check_ready(this)
    call this%stepper%advance(this%transform, this%numerics, this%step_number, &
                              this%previous_filtered, this%current)
  end subroutine advance

  subroutine get_fields(this, zeta, u, v, cfl)
    class(barotropic_solver), intent(inout) :: this
    real(real64), allocatable, intent(out) :: zeta(:, :), u(:, :), v(:, :)
    real(real64), intent(out), optional :: cfl
    complex(real64), allocatable :: psi(:, :)
    real(real64), allocatable :: dpsi_dlambda(:, :), dpsi_dphi(:, :)
    integer, allocatable :: nlon(:)
    integer :: j, n, m, T
    real(real64) :: cosphi, maximum_speed

    call check_ready(this)
    T = this%truncation
    allocate (psi(0:T + 1, 0:T))
    psi = cmplx(0.0_real64, 0.0_real64, kind=real64)
    do m = 0, T
      do n = max(1, m), T
        psi(n, m) = -earth_radius**2*this%current%zeta(n, m)/real(n*(n + 1), real64)
      end do
    end do

    call this%transform%spectral_to_grid(this%current%zeta, zeta)
    call this%transform%gradient_to_grid(psi, dpsi_dlambda, dpsi_dphi)
    call this%transform%allocate_field(u)
    call this%transform%allocate_field(v)
    u = 0.0_real64
    v = 0.0_real64
    nlon = this%transform%get_nlon()
    maximum_speed = 0.0_real64
    do j = 1, size(nlon)
      cosphi = sqrt(max(0.0_real64, 1.0_real64 - this%transform%mu(j)**2))
      u(1:nlon(j), j) = -dpsi_dphi(1:nlon(j), j)/earth_radius
      v(1:nlon(j), j) = dpsi_dlambda(1:nlon(j), j)/(earth_radius*cosphi)
      maximum_speed = max(maximum_speed, &
                          maxval(sqrt(u(1:nlon(j), j)**2 + v(1:nlon(j), j)**2)))
    end do

    if (present(cfl)) then
      cfl = maximum_speed*this%numerics%time_step/earth_radius*sqrt(real(T*(T + 1), real64))
    end if
  end subroutine get_fields

  integer function get_step(this) result(step)
    class(barotropic_solver), intent(in) :: this

    step = this%step_number
  end function get_step

  subroutine check_initialized(this)
    class(barotropic_solver), intent(in) :: this

    if (this%truncation < 1 .or. .not. allocated(this%current%zeta)) then
      error stop 'barotropic_solver: call init before use'
    end if
  end subroutine check_initialized

  subroutine check_ready(this)
    class(barotropic_solver), intent(in) :: this

    call check_initialized(this)
    if (this%step_number < 0) then
      error stop 'barotropic_solver: set initial vorticity before integration'
    end if
  end subroutine check_ready

end module barotropic_vorticity
