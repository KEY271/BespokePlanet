module harmonics
  use iso_c_binding, only: c_associated, c_int, c_null_ptr, c_ptr
  use iso_fortran_env, only: real64
  use fftw_interfaces, only: FFTW_ESTIMATE, FFTW_UNALIGNED, &
                             fftw_destroy_plan, fftw_execute_dft_c2r, fftw_execute_dft_r2c, &
                             fftw_plan_dft_c2r_1d, fftw_plan_dft_r2c_1d
  implicit none
  private

  type :: fft_plan
    integer :: length = 0
    type(c_ptr) :: forward = c_null_ptr
    type(c_ptr) :: backward = c_null_ptr
  contains
    final :: destroy_fft_plan
    final :: destroy_fft_plans
  end type fft_plan

  type, public :: harmonic_transform
    private
    integer, allocatable :: nlon(:)
    real(real64), allocatable :: w(:), pnm(:, :, :)
    type(fft_plan), allocatable :: fft_plans(:)
    integer :: current_T = -1
    real(real64), allocatable, public :: mu(:)
  contains
    procedure, public :: init
    procedure, public :: allocate_field
    procedure, public :: grid_to_spectral
    procedure, public :: spectral_to_grid
    procedure, public :: gradient
    procedure, public :: gradient_to_grid
    procedure, public :: get_nlon
  end type harmonic_transform

contains

  subroutine initialize_fft_plan(plan, n)
    type(fft_plan), intent(out) :: plan
    integer, intent(in) :: n
    real(real64), allocatable :: real_buffer(:)
    complex(real64), allocatable :: complex_buffer(:)
    integer(c_int) :: flags

    plan%length = n
    allocate (real_buffer(n), complex_buffer(n/2 + 1))
    flags = ior(FFTW_ESTIMATE, FFTW_UNALIGNED)
    plan%forward = fftw_plan_dft_r2c_1d(int(n, c_int), real_buffer, complex_buffer, flags)
    plan%backward = fftw_plan_dft_c2r_1d(int(n, c_int), complex_buffer, real_buffer, flags)

    if (.not. c_associated(plan%forward) .or. .not. c_associated(plan%backward)) then
      call destroy_fft_plan(plan)
      error stop "initialize_fft_plan: FFTW plan creation failed"
    end if
  end subroutine initialize_fft_plan

  subroutine execute_fft_forward(plan, values, spectrum)
    type(fft_plan), intent(in) :: plan
    real(real64), intent(inout) :: values(:)
    complex(real64), intent(out) :: spectrum(:)

    if (size(values) /= plan%length) error stop "execute_fft_forward: inconsistent transform length"
    if (size(spectrum) /= plan%length/2 + 1) error stop "execute_fft_forward: inconsistent spectrum length"

    call fftw_execute_dft_r2c(plan%forward, values, spectrum)
  end subroutine execute_fft_forward

  subroutine execute_fft_backward(plan, spectrum, values)
    type(fft_plan), intent(in) :: plan
    complex(real64), intent(inout) :: spectrum(:)
    real(real64), intent(out) :: values(:)

    if (size(spectrum) /= plan%length/2 + 1) error stop "execute_fft_backward: inconsistent spectrum length"
    if (size(values) /= plan%length) error stop "execute_fft_backward: inconsistent transform length"

    ! FFTW's c2r transform is unnormalized, matching the synthesis formula.
    call fftw_execute_dft_c2r(plan%backward, spectrum, values)
  end subroutine execute_fft_backward

  subroutine destroy_fft_plan(plan)
    type(fft_plan), intent(inout) :: plan

    if (c_associated(plan%forward)) call fftw_destroy_plan(plan%forward)
    if (c_associated(plan%backward)) call fftw_destroy_plan(plan%backward)
    plan%forward = c_null_ptr
    plan%backward = c_null_ptr
    plan%length = 0
  end subroutine destroy_fft_plan

  subroutine destroy_fft_plans(plans)
    type(fft_plan), intent(inout) :: plans(:)
    integer :: j

    do j = 1, size(plans)
      call destroy_fft_plan(plans(j))
    end do
  end subroutine destroy_fft_plans

  subroutine gauss_legendre(this)
    use lapack_interfaces, only: dstev

    class(harmonic_transform), intent(inout) :: this
    real(real64), allocatable :: e(:), z(:, :), work(:)
    integer :: T, n, k, info

    T = this%current_T
    n = 2*(T + 1)
    allocate (this%mu(n))
    allocate (this%w(n))
    allocate (e(max(1, n - 1)))
    allocate (z(n, n))
    allocate (work(max(1, 2*n - 2)))

    this%mu(:) = 0.0_real64
    do k = 1, n - 1
      e(k) = real(k, real64)/sqrt(4.0_real64*real(k, real64)**2 - 1.0_real64)
    end do

    call dstev('V', n, this%mu, e, z, n, work, info)

    if (info /= 0) error stop "DSTEV failed"

    this%w(:) = 2.0_real64*z(1, :)**2
  end subroutine gauss_legendre

  subroutine associated_legendre(this)
    class(harmonic_transform), intent(inout) :: this
    real(real64) :: x, s, anm, bnm
    integer :: T, j, n, m

    T = this%current_T
    allocate (this%pnm(2*(T + 1), 0:T + 1, 0:T))

    this%pnm = 0.0_real64
    do j = 1, 2*(T + 1)
      this%pnm(j, 0, 0) = 1.0_real64
      x = this%mu(j)
      s = sqrt(max(0.0_real64, 1.0_real64 - x*x))
      do m = 1, T
        this%pnm(j, m, m) = sqrt(real(2*m + 1, real64)/real(2*m, real64))*s*this%pnm(j, m - 1, m - 1)
      end do
      do m = 0, T
        this%pnm(j, m + 1, m) = sqrt(real(2*m + 3, real64))*x*this%pnm(j, m, m)
      end do
      do m = 0, T
        do n = m + 2, T + 1
          anm = sqrt(real(4*n*n - 1, real64)/real(n*n - m*m, real64))
          bnm = sqrt(real((2*n + 1)*((n - 1)*(n - 1) - m*m), real64)/real((2*n - 3)*(n*n - m*m), real64))
          this%pnm(j, n, m) = anm*x*this%pnm(j, n - 1, m) - bnm*this%pnm(j, n - 2, m)
        end do
      end do
    end do
  end subroutine associated_legendre

  subroutine init(this, T)
    class(harmonic_transform), intent(inout) :: this
    integer, intent(in) :: T
    integer :: G, j

    if (T < 0) error stop "init: T must be non-negative"

    if (allocated(this%mu)) deallocate (this%mu)
    if (allocated(this%w)) deallocate (this%w)
    if (allocated(this%pnm)) deallocate (this%pnm)
    if (allocated(this%nlon)) deallocate (this%nlon)
    if (allocated(this%fft_plans)) deallocate (this%fft_plans)

    this%current_T = T
    call gauss_legendre(this)
    call associated_legendre(this)

    G = T + 1
    allocate (this%nlon(2*G))
    allocate (this%fft_plans(G))
    do j = 1, G
      this%nlon(j) = 20 + 4*(j - 1)
      this%nlon(2*G + 1 - j) = this%nlon(j)
      call initialize_fft_plan(this%fft_plans(j), this%nlon(j))
    end do
  end subroutine init

  subroutine allocate_field(this, field)
    class(harmonic_transform), intent(in) :: this
    real(real64), allocatable, intent(out) :: field(:, :)
    integer :: T

    call check_transform_state(this)

    T = this%current_T
    allocate (field(4*(T + 1) + 16, 2*(T + 1)))
  end subroutine allocate_field

  subroutine grid_to_spectral(this, field, a)
    class(harmonic_transform), intent(in) :: this
    real(real64), intent(in) :: field(:, :)
    complex(real64), allocatable, intent(out) :: a(:, :)

    integer :: T, nlat, nlon_j, mmax_j, plan_index
    integer :: j, n, m
    real(real64), allocatable :: samples(:)
    complex(real64), allocatable :: fourier(:, :), spectrum(:)

    call check_transform_state(this)
    call check_field_shape(this, field)

    T = this%current_T
    nlat = 2*(T + 1)

    allocate (a(0:T + 1, 0:T))
    allocate (fourier(0:T, nlat))
    allocate (samples(maxval(this%nlon)))
    allocate (spectrum(maxval(this%nlon)/2 + 1))
    a = cmplx(0.0_real64, 0.0_real64, kind=real64)
    fourier = cmplx(0.0_real64, 0.0_real64, kind=real64)

    do j = 1, nlat
      nlon_j = this%nlon(j)
      plan_index = min(j, nlat + 1 - j)
      samples(1:nlon_j) = field(1:nlon_j, j)
      call execute_fft_forward(this%fft_plans(plan_index), samples(1:nlon_j), spectrum(1:nlon_j/2 + 1))
      mmax_j = min(T, nlon_j/2 - 1)
      fourier(0:mmax_j, j) = spectrum(1:mmax_j + 1)/real(nlon_j, real64)
    end do

    do m = 0, T
      do n = m, T
        do j = 1, nlat
          a(n, m) = a(n, m) + 0.5_real64*this%w(j)*this%pnm(j, n, m)*fourier(m, j)
        end do
      end do
    end do
  end subroutine grid_to_spectral

  subroutine spectral_to_grid(this, a, field)
    class(harmonic_transform), intent(in) :: this
    complex(real64), intent(in) :: a(0:, 0:)
    real(real64), allocatable, intent(out) :: field(:, :)

    integer :: T, nlat, nlon_j, mmax_j, plan_index
    integer :: j, n, m
    complex(real64) :: coefficient
    real(real64), allocatable :: samples(:)
    complex(real64), allocatable :: spectrum(:)

    call check_transform_state(this)
    call check_spectral_shape(this, a, "spectral_to_grid")

    T = this%current_T
    nlat = 2*(T + 1)

    call this%allocate_field(field)
    field = 0.0_real64
    allocate (samples(maxval(this%nlon)))
    allocate (spectrum(maxval(this%nlon)/2 + 1))

    do j = 1, nlat
      nlon_j = this%nlon(j)
      mmax_j = min(T, nlon_j/2 - 1)
      spectrum(1:nlon_j/2 + 1) = cmplx(0.0_real64, 0.0_real64, kind=real64)
      do m = 0, mmax_j
        coefficient = cmplx(0.0_real64, 0.0_real64, kind=real64)
        do n = m, T + 1
          coefficient = coefficient + a(n, m)*this%pnm(j, n, m)
        end do
        spectrum(m + 1) = coefficient
      end do

      plan_index = min(j, nlat + 1 - j)
      call execute_fft_backward(this%fft_plans(plan_index), spectrum(1:nlon_j/2 + 1), samples(1:nlon_j))
      field(1:nlon_j, j) = samples(1:nlon_j)
    end do
  end subroutine spectral_to_grid

  subroutine gradient(this, a, zonal, meridional)
    class(harmonic_transform), intent(in) :: this
    complex(real64), intent(in) :: a(0:, 0:)
    complex(real64), allocatable, intent(out) :: zonal(:, :), meridional(:, :)

    ! Spectral coefficients of the cos(latitude)-scaled spherical gradient:
    ! zonal = partial f / partial lambda,
    ! meridional = cos(latitude) * partial f / partial phi.
    integer :: T, n, m

    call check_transform_state(this)
    call check_spectral_shape(this, a, "gradient")

    T = this%current_T
    allocate (zonal(0:T + 1, 0:T), meridional(0:T + 1, 0:T))
    zonal = cmplx(0.0_real64, 0.0_real64, kind=real64)
    meridional = cmplx(0.0_real64, 0.0_real64, kind=real64)

    do m = 0, T
      do n = m, T
        zonal(n, m) = cmplx(0.0_real64, real(m, real64), kind=real64)*a(n, m)
        if (n > m) then
          meridional(n - 1, m) = meridional(n - 1, m) + &
                                 real(n + 1, real64)*legendre_epsilon(n, m)*a(n, m)
        end if
        meridional(n + 1, m) = meridional(n + 1, m) - &
                               real(n, real64)*legendre_epsilon(n + 1, m)*a(n, m)
      end do
    end do
  end subroutine gradient

  subroutine gradient_to_grid(this, a, dfdlambda, dfdphi)
    class(harmonic_transform) :: this
    complex(real64), intent(in) :: a(0:, 0:)
    real(real64), allocatable, intent(out) :: dfdlambda(:, :), dfdphi(:, :)
    complex(real64), allocatable :: zonal_spec(:, :), meridional_spec(:, :)
    integer :: j

    call this%gradient(a, zonal_spec, meridional_spec)

    call this%spectral_to_grid(zonal_spec, dfdlambda)
    call this%spectral_to_grid(meridional_spec, dfdphi)

    do j = 1, size(this%mu)
      dfdphi(1:this%nlon(j), j) = dfdphi(1:this%nlon(j), j)/sqrt(1.0_real64 - this%mu(j)**2)
    end do
  end subroutine

  subroutine check_transform_state(this)
    class(harmonic_transform), intent(in) :: this

    if (this%current_T < 0) then
      error stop "harmonics: call init(T) before using the transform"
    end if
    if (.not. allocated(this%w) .or. .not. allocated(this%pnm) .or. &
        .not. allocated(this%nlon) .or. .not. allocated(this%fft_plans)) then
      error stop "harmonics: transform tables are not initialized"
    end if
  end subroutine check_transform_state

  subroutine check_spectral_shape(this, a, procedure_name)
    class(harmonic_transform), intent(in) :: this
    complex(real64), intent(in) :: a(0:, 0:)
    character(*), intent(in) :: procedure_name
    integer :: T

    T = this%current_T
    if (ubound(a, 1) < T + 1 .or. ubound(a, 2) < T) then
      error stop procedure_name//": a must contain indices 0:T+1,0:T"
    end if
  end subroutine check_spectral_shape

  pure function legendre_epsilon(n, m) result(value)
    integer, intent(in) :: n, m
    real(real64) :: value
    real(real64) :: rn, rm

    if (n == 0) then
      value = 0.0_real64
      return
    end if

    rn = real(n, real64)
    rm = real(m, real64)
    value = sqrt((rn*rn - rm*rm)/(4.0_real64*rn*rn - 1.0_real64))
  end function legendre_epsilon

  function get_nlon(this) result(nlon)
    class(harmonic_transform), intent(in) :: this
    integer, allocatable :: nlon(:)

    if (.not. allocated(this%nlon)) then
      error stop "harmonics: call init(T) before accessing nlon"
    end if

    nlon = this%nlon
  end function get_nlon

  subroutine check_field_shape(this, field)
    class(harmonic_transform), intent(in) :: this
    real(real64), intent(in) :: field(:, :)
    integer :: T, nlat, max_nlon

    T = this%current_T
    nlat = 2*(T + 1)
    max_nlon = 4*(T + 1) + 16

    if (size(field, 1) /= max_nlon .or. size(field, 2) /= nlat) then
      error stop "harmonics: field has an inconsistent shape for T"
    end if
  end subroutine check_field_shape
end module harmonics
