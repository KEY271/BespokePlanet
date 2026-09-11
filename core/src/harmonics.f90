module harmonics
  use iso_fortran_env, only: real64
  implicit none
  private

  type, public :: harmonic_transform
    private
    integer, allocatable :: nlon(:)
    real(real64), allocatable :: mu(:), w(:), pnm(:, :, :)
    integer :: current_T = -1
  contains
    procedure, public :: init
    procedure, public :: allocate_field
    procedure, public :: field_to_a
    procedure, public :: a_to_field
    procedure, public :: get_nlon
  end type harmonic_transform

contains

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
    allocate (this%pnm(2*(T + 1), 0:T, 0:T))

    this%pnm = 0.0_real64
    do j = 1, 2*(T + 1)
      this%pnm(j, 0, 0) = 1.0_real64
      x = this%mu(j)
      s = sqrt(max(0.0_real64, 1.0_real64 - x*x))
      do m = 1, T
        this%pnm(j, m, m) = sqrt(real(2*m + 1, real64)/real(2*m, real64))*s*this%pnm(j, m - 1, m - 1)
      end do
      do m = 0, T - 1
        this%pnm(j, m + 1, m) = sqrt(real(2*m + 3, real64))*x*this%pnm(j, m, m)
      end do
      do m = 0, T
        do n = m + 2, T
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

    this%current_T = T
    call gauss_legendre(this)
    call associated_legendre(this)

    G = T + 1
    allocate (this%nlon(2*G))
    do j = 1, G
      this%nlon(j) = 20 + 4*(j - 1)
      this%nlon(2*G + 1 - j) = this%nlon(j)
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

  subroutine field_to_a(this, field, a)
    class(harmonic_transform), intent(in) :: this
    real(real64), intent(in) :: field(:, :)
    complex(real64), allocatable, intent(out) :: a(:, :)

    integer :: T, nlat, nlon_j, mmax_j
    integer :: j, k, n, m
    real(real64) :: pi, angle
    complex(real64) :: phase
    complex(real64), allocatable :: fm(:)

    call check_transform_state(this)
    call check_field_shape(this, field)

    T = this%current_T
    nlat = 2*(T + 1)
    pi = acos(-1.0_real64)

    allocate (a(0:T, 0:T))
    allocate (fm(nlat))
    a = cmplx(0.0_real64, 0.0_real64, kind=real64)

    do m = 0, T
      fm = cmplx(0.0_real64, 0.0_real64, kind=real64)

      do j = 1, nlat
        nlon_j = this%nlon(j)
        mmax_j = min(T, nlon_j/2 - 1)
        if (m > mmax_j) cycle

        do k = 1, nlon_j
          angle = 2.0_real64*pi*real(m*(k - 1), real64)/real(nlon_j, real64)
          phase = cmplx(cos(angle), -sin(angle), kind=real64)
          fm(j) = fm(j) + field(k, j)*phase
        end do
        fm(j) = fm(j)/real(nlon_j, real64)
      end do

      do n = m, T
        do j = 1, nlat
          a(n, m) = a(n, m) + 0.5_real64*this%w(j)*this%pnm(j, n, m)*fm(j)
        end do
      end do
    end do
  end subroutine field_to_a

  subroutine a_to_field(this, a, field)
    class(harmonic_transform), intent(in) :: this
    complex(real64), intent(in) :: a(0:, 0:)
    real(real64), allocatable, intent(out) :: field(:, :)

    integer :: T, nlat, nlon_j, mmax_j
    integer :: j, k, n, m
    real(real64) :: pi, angle
    complex(real64) :: phase
    complex(real64), allocatable :: fm(:)

    call check_transform_state(this)

    T = this%current_T
    if (ubound(a, 1) < T .or. ubound(a, 2) < T) then
      error stop "a_to_field: a must contain indices 0:T,0:T"
    end if

    nlat = 2*(T + 1)
    pi = acos(-1.0_real64)

    call this%allocate_field(field)
    field = 0.0_real64
    allocate (fm(nlat))

    do m = 0, T
      fm = cmplx(0.0_real64, 0.0_real64, kind=real64)

      do j = 1, nlat
        mmax_j = min(T, this%nlon(j)/2 - 1)
        if (m > mmax_j) cycle

        do n = m, T
          fm(j) = fm(j) + a(n, m)*this%pnm(j, n, m)
        end do
      end do

      do j = 1, nlat
        nlon_j = this%nlon(j)
        mmax_j = min(T, nlon_j/2 - 1)
        if (m > mmax_j) cycle

        do k = 1, nlon_j
          angle = 2.0_real64*pi*real(m*(k - 1), real64)/real(nlon_j, real64)

          if (m == 0) then
            field(k, j) = field(k, j) + real(fm(j), kind=real64)
          else
            phase = cmplx(cos(angle), sin(angle), kind=real64)
            field(k, j) = field(k, j) + 2.0_real64*real(fm(j)*phase, kind=real64)
          end if
        end do
      end do
    end do
  end subroutine a_to_field

  subroutine check_transform_state(this)
    class(harmonic_transform), intent(in) :: this

    if (this%current_T < 0) then
      error stop "harmonics: call init(T) before using the transform"
    end if
    if (.not. allocated(this%w) .or. .not. allocated(this%pnm) .or. .not. allocated(this%nlon)) then
      error stop "harmonics: transform tables are not initialized"
    end if
  end subroutine check_transform_state

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
