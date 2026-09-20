module field_binary_writer
  use iso_fortran_env, only: real64
  use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
  implicit none
  private
  public :: write_field, write_rectangular_field, write_spectral_field
  public :: check_finite, check_spectral_finite, check_rectangular_finite

contains

  subroutine write_field(path, nlon, field)
    character(*), intent(in) :: path
    integer, intent(in) :: nlon(:)
    real(real64), intent(in) :: field(:, :)
    integer :: unit, j
    open (newunit=unit, file=path, status='replace', access='stream', &
          form='unformatted', action='write', convert='little_endian')
    do j = 1, size(nlon)
      write (unit) field(1:nlon(j), j)
    end do
    close (unit)
  end subroutine write_field

  subroutine write_rectangular_field(path, field)
    character(*), intent(in) :: path
    real(real64), intent(in) :: field(:, :)
    integer :: unit, k
    open (newunit=unit, file=path, status='replace', access='stream', &
          form='unformatted', action='write', convert='little_endian')
    do k = 1, size(field, 2)
      write (unit) field(:, k)
    end do
    close (unit)
  end subroutine write_rectangular_field

  subroutine write_spectral_field(path, field)
    character(*), intent(in) :: path
    complex(real64), intent(in) :: field(0:, 0:)
    integer :: unit, n, m
    open (newunit=unit, file=path, status='replace', access='stream', &
          form='unformatted', action='write', convert='little_endian')
    do m = 0, ubound(field, 2)
      do n = 0, ubound(field, 1) - 1
        write (unit) real(field(n, m), real64), aimag(field(n, m))
      end do
    end do
    close (unit)
  end subroutine write_spectral_field

  subroutine check_finite(name, nlon, field)
    character(*), intent(in) :: name
    integer, intent(in) :: nlon(:)
    real(real64), intent(in) :: field(:, :)
    integer :: j
    do j = 1, size(nlon)
      if (.not. all(ieee_is_finite(field(1:nlon(j), j)))) then
        error stop 'non-finite value in output field '//name
      end if
    end do
  end subroutine check_finite

  subroutine check_spectral_finite(name, field)
    character(*), intent(in) :: name
    complex(real64), intent(in) :: field(0:, 0:)
    if (.not. all(ieee_is_finite(real(field, real64))) .or. &
        .not. all(ieee_is_finite(aimag(field)))) then
      error stop 'non-finite value in output field '//name
    end if
  end subroutine check_spectral_finite

  subroutine check_rectangular_finite(name, field)
    character(*), intent(in) :: name
    real(real64), intent(in) :: field(:, :)
    if (.not. all(ieee_is_finite(field))) then
      error stop 'non-finite value in output field '//name
    end if
  end subroutine check_rectangular_finite

end module field_binary_writer
