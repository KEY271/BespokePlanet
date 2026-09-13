module barotropic_vorticity
  use iso_fortran_env, only: real64
  use harmonics, only: harmonic_transform
  use raw_filter, only: raw_filter_epsilon, raw_filter_alpha, apply_raw_filter
  use spectral_hyperdiffusion, only: hyperdiffusion_order, &
                                       hyperdiffusion_timescale_seconds, &
                                       apply_spectral_hyperdiffusion
  implicit none
  private

  real(real64), parameter, public :: earth_radius = 6.371e6_real64
  real(real64), parameter, public :: rotation_rate = 7.2921159e-5_real64
  public :: raw_filter_epsilon, raw_filter_alpha
  public :: hyperdiffusion_order, hyperdiffusion_timescale_seconds

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
    complex(real64), allocatable :: rhs(:, :), midpoint(:, :), candidate(:, :), next(:, :)
    complex(real64), allocatable :: filtered_current(:, :)
    integer :: T

    call check_ready(this)
    T = this%truncation
    allocate (midpoint(0:T + 1, 0:T), candidate(0:T + 1, 0:T))
    midpoint = cmplx(0.0_real64, 0.0_real64, kind=real64)
    candidate = cmplx(0.0_real64, 0.0_real64, kind=real64)

    if (this%step_number == 0) then
      call tendency(this, this%current, rhs)
      midpoint = this%current + 0.5_real64*this%dt*rhs
      call apply_spectral_hyperdiffusion(T, 0.5_real64*this%dt, midpoint)
      call enforce_spectral_constraints(this, midpoint)

      call tendency(this, midpoint, rhs)
      allocate (next(0:T + 1, 0:T))
      next = this%current + this%dt*rhs
      call apply_spectral_hyperdiffusion(T, this%dt, next)
      this%previous_filtered = this%current
    else
      call tendency(this, this%current, rhs)
      candidate = this%previous_filtered + 2.0_real64*this%dt*rhs
      call apply_spectral_hyperdiffusion(T, 2.0_real64*this%dt, candidate)
      call enforce_spectral_constraints(this, candidate)
      call apply_raw_filter(this%previous_filtered, this%current, candidate, &
                            filtered_current, next)
      this%previous_filtered = filtered_current
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
