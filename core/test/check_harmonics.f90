program check_harmonics
  use iso_fortran_env, only: real64
  use harmonics, only: init_harmonics, allocate_field, field_to_a, a_to_field, nlon
  implicit none

  integer, parameter :: T = 63
  real(real64), allocatable :: field(:, :), field2(:, :)
  complex(real64), allocatable :: a(:, :)
  integer :: j, k, npts
  real(real64) :: err, max_err, rms_err, rel_rms_err, sum_err2, sum_ref2

  call allocate_field(T, field)
  call init_harmonics(T)

  field = 1.0_real64
  call field_to_a(T, field, a)
  call a_to_field(T, a, field2)
  call calc_error(field, field2)

contains
  subroutine calc_error(f, f2)
    real(real64), intent(in) :: f(:, :), f2(:, :)

    max_err = 0.0_real64
    sum_err2 = 0.0_real64
    sum_ref2 = 0.0_real64
    npts = 0

    do j = 1, 2*(T + 1)
      do k = 1, nlon(j)
        err = f2(k, j) - f(k, j)
        max_err = max(max_err, abs(err))
        sum_err2 = sum_err2 + err**2
        sum_ref2 = sum_ref2 + f(k, j)**2
        npts = npts + 1
      end do
    end do

    rms_err = sqrt(sum_err2/real(npts, real64))

    if (sum_ref2 > 0.0_real64) then
      rel_rms_err = sqrt(sum_err2/sum_ref2)
    else
      rel_rms_err = 0.0_real64
    end if

    write (*, '(a,es24.16)') 'max error          = ', max_err
    write (*, '(a,es24.16)') 'RMS error          = ', rms_err
    write (*, '(a,es24.16)') 'relative RMS error = ', rel_rms_err
  end subroutine calc_error
end program check_harmonics
