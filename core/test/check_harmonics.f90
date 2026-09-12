program check_harmonics
  use iso_fortran_env, only: real64
  use harmonics, only: harmonic_transform
  implicit none

  integer, parameter :: T = 63
  real(real64), parameter :: tolerance = 1.0e-10_real64
  type(harmonic_transform) :: transform
  integer, allocatable :: nlon(:)
  integer :: mmax_all

  call transform%init(T)
  nlon = transform%get_nlon()
  mmax_all = minval(nlon)/2 - 1

  call check_constant()
  call check_mode(1, 0, cmplx(1.0_real64, 0.0_real64, kind=real64))
  call check_mode(2, 1, cmplx(0.75_real64, 0.0_real64, kind=real64))
  call check_mode(3, 2, cmplx(0.0_real64, -0.5_real64, kind=real64))
  call check_mode(8, 3, cmplx(0.25_real64, 0.5_real64, kind=real64))
  call check_mode(T, 0, cmplx(1.0_real64, 0.0_real64, kind=real64))
  call check_mode(T, mmax_all, cmplx(-0.25_real64, 0.5_real64, kind=real64))

contains

  subroutine check_constant()
    real(real64), allocatable :: field(:, :), field2(:, :)
    complex(real64), allocatable :: a(:, :)
    real(real64) :: max_err, rel_rms_err

    call transform%allocate_field(field)
    field = 1.0_real64

    call transform%field_to_a(field, a)
    call transform%a_to_field(a, field2)
    call calc_field_error(field, field2, max_err, rel_rms_err)

    write (*, '(a,2(a,es12.4))') 'constant:', ' max error = ', max_err, &
      ' relative RMS error = ', rel_rms_err
    if (max_err > tolerance .or. rel_rms_err > tolerance) then
      error stop 'constant round-trip failed'
    end if
  end subroutine check_constant

  subroutine check_mode(n_mode, m_mode, amplitude)
    integer, intent(in) :: n_mode, m_mode
    complex(real64), intent(in) :: amplitude
    real(real64), allocatable :: field(:, :), field2(:, :)
    complex(real64), allocatable :: a(:, :), a2(:, :)
    real(real64) :: coefficient_err, field_err, rel_field_err

    allocate (a(0:T, 0:T))
    a = cmplx(0.0_real64, 0.0_real64, kind=real64)
    a(n_mode, m_mode) = amplitude

    call transform%a_to_field(a, field)
    call transform%field_to_a(field, a2)
    call transform%a_to_field(a2, field2)

    coefficient_err = maxval(abs(a2 - a))
    call calc_field_error(field, field2, field_err, rel_field_err)

    write (*, '(a,i0,a,i0,a,3(a,es12.4))') 'mode (n=', n_mode, ', m=', m_mode, '):', &
      ' coefficient error = ', coefficient_err, ' max field error = ', field_err, &
      ' relative RMS error = ', rel_field_err
    if (coefficient_err > tolerance .or. field_err > tolerance .or. rel_field_err > tolerance) then
      error stop 'spherical harmonic mode round-trip failed'
    end if
  end subroutine check_mode

  subroutine calc_field_error(field, field2, max_err, rel_rms_err)
    real(real64), intent(in) :: field(:, :), field2(:, :)
    real(real64), intent(out) :: max_err, rel_rms_err
    integer :: j, k, npts
    real(real64) :: err, sum_err2, sum_ref2

    max_err = 0.0_real64
    sum_err2 = 0.0_real64
    sum_ref2 = 0.0_real64
    npts = 0

    do j = 1, 2*(T + 1)
      do k = 1, nlon(j)
        err = field2(k, j) - field(k, j)
        max_err = max(max_err, abs(err))
        sum_err2 = sum_err2 + err**2
        sum_ref2 = sum_ref2 + field(k, j)**2
        npts = npts + 1
      end do
    end do

    if (sum_ref2 > 0.0_real64) then
      rel_rms_err = sqrt(sum_err2/sum_ref2)
    else
      rel_rms_err = sqrt(sum_err2/real(npts, real64))
    end if
  end subroutine calc_field_error
end program check_harmonics
