!> Low-level JSON fragments shared by case metadata writers.  This module knows
!> nothing about equations, variables, or case schemas.
module json_writer
  use iso_fortran_env, only: real64
  implicit none
  private
  public :: write_inline_real_values, write_real_array, write_integer_array, write_ring_offsets

contains

  !> Writes comma-separated values followed by "],"; the caller writes the key and "[".
  subroutine write_inline_real_values(unit, values)
    integer, intent(in) :: unit
    real(real64), intent(in) :: values(:)
    integer :: i

    do i = 1, size(values)
      write (unit, '(es24.16e3)', advance='no') values(i)
      if (i < size(values)) write (unit, '(a)', advance='no') ','
    end do
    write (unit, '(a)') '],'
  end subroutine write_inline_real_values

  subroutine write_real_array(unit, name, values, trailing_comma)
    integer, intent(in) :: unit
    character(*), intent(in) :: name
    real(real64), intent(in) :: values(:)
    logical, intent(in) :: trailing_comma
    integer :: j

    write (unit, '(3a)', advance='no') '    "', trim(name), '": ['
    do j = 1, size(values)
      write (unit, '(es24.16e3)', advance='no') values(j)
      if (j < size(values)) write (unit, '(a)', advance='no') ','
    end do
    if (trailing_comma) then
      write (unit, '(a)') '],'
    else
      write (unit, '(a)') ']'
    end if
  end subroutine write_real_array

  subroutine write_integer_array(unit, name, values, trailing_comma)
    integer, intent(in) :: unit
    character(*), intent(in) :: name
    integer, intent(in) :: values(:)
    logical, intent(in) :: trailing_comma
    integer :: j

    write (unit, '(3a)', advance='no') '    "', trim(name), '": ['
    do j = 1, size(values)
      write (unit, '(i0)', advance='no') values(j)
      if (j < size(values)) write (unit, '(a)', advance='no') ','
    end do
    if (trailing_comma) then
      write (unit, '(a)') '],'
    else
      write (unit, '(a)') ']'
    end if
  end subroutine write_integer_array

  !> Cumulative ring offsets of a ring-major array, starting at 0.
  subroutine write_ring_offsets(unit, nlon)
    integer, intent(in) :: unit
    integer, intent(in) :: nlon(:)
    integer :: j, offset

    write (unit, '(a)', advance='no') '    "ring_offsets": [0'
    offset = 0
    do j = 1, size(nlon)
      offset = offset + nlon(j)
      write (unit, '(a,i0)', advance='no') ',', offset
    end do
    write (unit, '(a)') ']'
  end subroutine write_ring_offsets

end module json_writer
