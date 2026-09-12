module barotropic_vorticity
  use iso_fortran_env, only: real64
  use harmonics, only: harmonic_transform
  implicit none
  private

  real(real64), parameter, public :: earth_radius = 6.371e6_real64
  real(real64), parameter, public :: rotation_rate = 7.2921159e-5_real64
  real(real64), parameter, public :: raw_filter_epsilon = 0.1_real64
  real(real64), parameter, public :: raw_filter_alpha = 0.53_real64
  integer, parameter, public :: hyperdiffusion_order = 4
  real(real64), parameter, public :: hyperdiffusion_timescale_seconds = &
    6.0_real64*3600.0_real64

  type, public :: barotropic_solver
    private
    type(harmonic_transform) :: transform
    integer :: truncation = -1
    integer :: step_number = -1
    real(real64) :: dt = 0.0_real64
    complex(real64), allocatable :: previous_filtered(:, :)
    complex(real64), allocatable :: current(:, :)
  contains
    procedure, public :: init
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
    this%dt = dt
    this%step_number = -1
    call this%transform%init(truncation)

    if (allocated(this%previous_filtered)) deallocate (this%previous_filtered)
    if (allocated(this%current)) deallocate (this%current)
    allocate (this%previous_filtered(0:truncation + 1, 0:truncation))
    allocate (this%current(0:truncation + 1, 0:truncation))
    this%previous_filtered = cmplx(0.0_real64, 0.0_real64, kind=real64)
    this%current = cmplx(0.0_real64, 0.0_real64, kind=real64)
  end subroutine init

  subroutine set_initial_vorticity(this, zeta)
    class(barotropic_solver), intent(inout) :: this
    complex(real64), intent(in) :: zeta(0:, 0:)

    call check_initialized(this)
    if (ubound(zeta, 1) < this%truncation + 1 .or. &
        ubound(zeta, 2) < this%truncation) then
      error stop 'barotropic_solver: initial vorticity has an inconsistent shape'
    end if

    this%current = zeta(0:this%truncation + 1, 0:this%truncation)
    call enforce_spectral_constraints(this, this%current)
    this%previous_filtered = this%current
    this%step_number = 0
  end subroutine set_initial_vorticity

  subroutine advance(this)
    class(barotropic_solver), intent(inout) :: this
    complex(real64), allocatable :: rhs(:, :), midpoint(:, :), next(:, :), delta(:, :)
    integer :: n, m, T
    real(real64) :: damping

    call check_ready(this)
    T = this%truncation
    allocate (midpoint(0:T + 1, 0:T), next(0:T + 1, 0:T))
    midpoint = cmplx(0.0_real64, 0.0_real64, kind=real64)
    next = cmplx(0.0_real64, 0.0_real64, kind=real64)

    if (this%step_number == 0) then
      call tendency(this, this%current, rhs)
      do m = 0, T
        do n = m, T
          damping = damping_rate(this, n)
          midpoint(n, m) = (this%current(n, m) + 0.5_real64*this%dt*rhs(n, m))/ &
                           (1.0_real64 + 0.5_real64*this%dt*damping)
        end do
      end do
      call enforce_spectral_constraints(this, midpoint)

      call tendency(this, midpoint, rhs)
      do m = 0, T
        do n = m, T
          damping = damping_rate(this, n)
          next(n, m) = (this%current(n, m) + this%dt*rhs(n, m))/ &
                       (1.0_real64 + this%dt*damping)
        end do
      end do
      this%previous_filtered = this%current
    else
      call tendency(this, this%current, rhs)
      do m = 0, T
        do n = m, T
          damping = damping_rate(this, n)
          next(n, m) = (this%previous_filtered(n, m) + 2.0_real64*this%dt*rhs(n, m))/ &
                       (1.0_real64 + 2.0_real64*this%dt*damping)
        end do
      end do
      call enforce_spectral_constraints(this, next)

      allocate (delta(0:T + 1, 0:T))
      delta = 0.5_real64*raw_filter_epsilon* &
              (this%previous_filtered - 2.0_real64*this%current + next)
      this%previous_filtered = this%current + raw_filter_alpha*delta
      next = next - (1.0_real64 - raw_filter_alpha)*delta
      call enforce_spectral_constraints(this, this%previous_filtered)
    end if

    call enforce_spectral_constraints(this, next)
    this%current = next
    this%step_number = this%step_number + 1
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
        psi(n, m) = -earth_radius**2*this%current(n, m)/real(n*(n + 1), real64)
      end do
    end do

    call this%transform%spectral_to_grid(this%current, zeta)
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
      cfl = maximum_speed*this%dt/earth_radius*sqrt(real(T*(T + 1), real64))
    end if
  end subroutine get_fields

  integer function get_step(this) result(step)
    class(barotropic_solver), intent(in) :: this

    step = this%step_number
  end function get_step

  subroutine tendency(this, zeta, rhs)
    class(barotropic_solver), intent(inout) :: this
    complex(real64), intent(in) :: zeta(0:, 0:)
    complex(real64), allocatable, intent(out) :: rhs(:, :)
    complex(real64), allocatable :: psi(:, :), flux_u_spec(:, :), flux_v_spec(:, :)
    real(real64), allocatable :: zeta_grid(:, :), u(:, :), v(:, :), q(:, :)
    real(real64), allocatable :: dpsi_dlambda(:, :), dpsi_dphi(:, :)
    real(real64), allocatable :: flux_u(:, :), flux_v(:, :), tendency_grid(:, :)
    real(real64), allocatable :: dflux_u_dlambda(:, :), unused_u(:, :)
    real(real64), allocatable :: unused_v(:, :), dflux_v_dphi(:, :)
    integer, allocatable :: nlon(:)
    integer :: j, n, m, T
    real(real64) :: cosphi

    T = this%truncation
    allocate (psi(0:T + 1, 0:T))
    psi = cmplx(0.0_real64, 0.0_real64, kind=real64)
    do m = 0, T
      do n = max(1, m), T
        psi(n, m) = -earth_radius**2*zeta(n, m)/real(n*(n + 1), real64)
      end do
    end do

    call this%transform%spectral_to_grid(zeta, zeta_grid)
    call this%transform%gradient_to_grid(psi, dpsi_dlambda, dpsi_dphi)
    call this%transform%allocate_field(u)
    call this%transform%allocate_field(v)
    call this%transform%allocate_field(q)
    call this%transform%allocate_field(flux_u)
    call this%transform%allocate_field(flux_v)
    u = 0.0_real64
    v = 0.0_real64
    q = 0.0_real64
    flux_u = 0.0_real64
    flux_v = 0.0_real64
    nlon = this%transform%get_nlon()

    do j = 1, size(nlon)
      cosphi = sqrt(max(0.0_real64, 1.0_real64 - this%transform%mu(j)**2))
      u(1:nlon(j), j) = -dpsi_dphi(1:nlon(j), j)/earth_radius
      v(1:nlon(j), j) = dpsi_dlambda(1:nlon(j), j)/(earth_radius*cosphi)
      q(1:nlon(j), j) = zeta_grid(1:nlon(j), j) + &
                         2.0_real64*rotation_rate*this%transform%mu(j)
      flux_u(1:nlon(j), j) = u(1:nlon(j), j)*q(1:nlon(j), j)
      flux_v(1:nlon(j), j) = v(1:nlon(j), j)*q(1:nlon(j), j)*cosphi
    end do

    call this%transform%grid_to_spectral(flux_u, flux_u_spec)
    call this%transform%grid_to_spectral(flux_v, flux_v_spec)
    call this%transform%gradient_to_grid(flux_u_spec, dflux_u_dlambda, unused_u)
    call this%transform%gradient_to_grid(flux_v_spec, unused_v, dflux_v_dphi)
    call this%transform%allocate_field(tendency_grid)
    tendency_grid = 0.0_real64
    do j = 1, size(nlon)
      cosphi = sqrt(max(0.0_real64, 1.0_real64 - this%transform%mu(j)**2))
      tendency_grid(1:nlon(j), j) = &
        -(dflux_u_dlambda(1:nlon(j), j) + dflux_v_dphi(1:nlon(j), j))/ &
        (earth_radius*cosphi)
    end do

    call this%transform%grid_to_spectral(tendency_grid, rhs)
    call enforce_spectral_constraints(this, rhs)
  end subroutine tendency

  pure real(real64) function damping_rate(this, n) result(rate)
    class(barotropic_solver), intent(in) :: this
    integer, intent(in) :: n
    real(real64) :: ratio

    ratio = real(n*(n + 1), real64)/ &
            real(this%truncation*(this%truncation + 1), real64)
    rate = ratio**hyperdiffusion_order/hyperdiffusion_timescale_seconds
  end function damping_rate

  subroutine enforce_spectral_constraints(this, field)
    class(barotropic_solver), intent(in) :: this
    complex(real64), intent(inout) :: field(0:, 0:)
    integer :: m, n, T

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
    class(barotropic_solver), intent(in) :: this

    if (this%truncation < 1 .or. .not. allocated(this%current)) then
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
