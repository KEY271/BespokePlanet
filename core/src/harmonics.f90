module harmonics
  use iso_fortran_env, only: real64
  implicit none
  private

  integer, allocatable, public :: nlon(:)
  real(real64), allocatable :: mu(:), w(:), pnm(:, :, :)
  integer :: current_T = -1

  public :: init_harmonics, allocate_field
  public :: field_to_a, a_to_field

contains

  subroutine gauss_legendre(T)
    use lapack_interfaces, only: dstev

    integer, intent(in) :: T
    real(real64), allocatable :: e(:), z(:, :), work(:)
    integer :: n, k, info

    n = 2*(T + 1)
    allocate (mu(n))
    allocate (w(n))
    allocate (e(max(1, n - 1)))
    allocate (z(n, n))
    allocate (work(max(1, 2*n - 2)))

    mu(:) = 0.0_real64
    do k = 1, n - 1
      e(k) = real(k, real64)/sqrt(4.0_real64*real(k, real64)**2 - 1.0_real64)
    end do

    call dstev('V', n, mu, e, z, n, work, info)

    if (info /= 0) error stop "DSTEV failed"

    w(:) = 2.0_real64*z(1, :)**2
  end subroutine gauss_legendre

  subroutine associated_legendre(T)
    integer, intent(in) :: T
    real(real64) :: x, s, anm, bnm
    integer :: j, n, m

    allocate (pnm(2*(T + 1), 0:T, 0:T))

    pnm = 0.0_real64
    do j = 1, 2*(T + 1)
      pnm(j, 0, 0) = 1.0_real64
      x = mu(j)
      s = sqrt(max(0.0_real64, 1.0_real64 - x*x))
      do m = 1, T
        pnm(j, m, m) = sqrt(real(2*m + 1, real64)/real(2*m, real64))*s*pnm(j, m - 1, m - 1)
      end do
      do m = 0, T - 1
        pnm(j, m + 1, m) = sqrt(real(2*m + 3, real64))*x*pnm(j, m, m)
      end do
      do m = 0, T
        do n = m + 2, T
          anm = sqrt(real(4*n*n - 1, real64)/real(n*n - m*m, real64))
          bnm = sqrt(real((2*n + 1)*((n - 1)*(n - 1) - m*m), real64)/real((2*n - 3)*(n*n - m*m), real64))
          pnm(j, n, m) = anm*x*pnm(j, n - 1, m) - bnm*pnm(j, n - 2, m)
        end do
      end do
    end do
  end subroutine associated_legendre

  subroutine init_harmonics(T)
    integer, intent(in) :: T
    integer :: G, j

    if (T < 0) error stop "init_harmonics: T must be non-negative"

    if (allocated(mu)) deallocate (mu)
    if (allocated(w)) deallocate (w)
    if (allocated(pnm)) deallocate (pnm)
    if (allocated(nlon)) deallocate (nlon)

    call gauss_legendre(T)
    call associated_legendre(T)

    G = T + 1
    allocate (nlon(2*G))
    do j = 1, G
      nlon(j) = 20 + 4*(j - 1)
      nlon(2*G + 1 - j) = nlon(j)
    end do

    current_T = T
  end subroutine init_harmonics

  subroutine allocate_field(T, field)
    integer, intent(in) :: T
    real(real64), allocatable, intent(out) :: field(:, :)

    if (T < 0) error stop "allocate_field: T must be non-negative"
    allocate (field(4*(T + 1) + 16, 2*(T + 1)))
  end subroutine allocate_field

  subroutine field_to_a(T, field, a)
    integer, intent(in) :: T
    real(real64), intent(in) :: field(:, :)
    complex(real64), allocatable, intent(out) :: a(:, :)

    integer :: nlat, nlon_j, mmax_j
    integer :: j, k, n, m
    real(real64) :: pi, angle
    complex(real64) :: phase
    complex(real64), allocatable :: fm(:)

    call check_transform_state(T)
    call check_field_shape(T, field)

    nlat = 2*(T + 1)
    pi = acos(-1.0_real64)

    allocate (a(0:T, 0:T))
    allocate (fm(nlat))
    a = cmplx(0.0_real64, 0.0_real64, kind=real64)

    do m = 0, T
      fm = cmplx(0.0_real64, 0.0_real64, kind=real64)

      do j = 1, nlat
        nlon_j = nlon(j)
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
          a(n, m) = a(n, m) + 0.5_real64*w(j)*pnm(j, n, m)*fm(j)
        end do
      end do
    end do
  end subroutine field_to_a

  subroutine a_to_field(T, a, field)
    integer, intent(in) :: T
    complex(real64), intent(in) :: a(0:, 0:)
    real(real64), allocatable, intent(out) :: field(:, :)

    integer :: nlat, nlon_j, mmax_j
    integer :: j, k, n, m
    real(real64) :: pi, angle
    complex(real64) :: phase
    complex(real64), allocatable :: fm(:)

    call check_transform_state(T)

    if (ubound(a, 1) < T .or. ubound(a, 2) < T) then
      error stop "a_to_field: a must contain indices 0:T,0:T"
    end if

    nlat = 2*(T + 1)
    pi = acos(-1.0_real64)

    call allocate_field(T, field)
    field = 0.0_real64
    allocate (fm(nlat))

    do m = 0, T
      fm = cmplx(0.0_real64, 0.0_real64, kind=real64)

      do j = 1, nlat
        mmax_j = min(T, nlon(j)/2 - 1)
        if (m > mmax_j) cycle

        do n = m, T
          fm(j) = fm(j) + a(n, m)*pnm(j, n, m)
        end do
      end do

      do j = 1, nlat
        nlon_j = nlon(j)
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

  subroutine check_transform_state(T)
    integer, intent(in) :: T

    if (current_T /= T) then
      error stop "harmonics: call init_harmonics(T) before transforming"
    end if
    if (.not. allocated(w) .or. .not. allocated(pnm) .or. .not. allocated(nlon)) then
      error stop "harmonics: transform tables are not initialized"
    end if
  end subroutine check_transform_state

  subroutine check_field_shape(T, field)
    integer, intent(in) :: T
    real(real64), intent(in) :: field(:, :)
    integer :: nlat, max_nlon

    nlat = 2*(T + 1)
    max_nlon = 4*(T + 1) + 16

    if (size(field, 1) /= max_nlon .or. size(field, 2) /= nlat) then
      error stop "harmonics: field has an inconsistent shape for T"
    end if
  end subroutine check_field_shape
end module harmonics
