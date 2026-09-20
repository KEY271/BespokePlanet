!> Low-level CSV output: a header line and appended rows.  Field formatting is
!> fixed here so that every case writes numbers the same way; the module knows
!> nothing about what the columns mean.
module csv_writer
  use iso_fortran_env, only: real64
  implicit none
  private
  public :: write_csv_header, append_csv_row, csv_real, csv_integer

contains

  !> Creates (or truncates) the file and writes the header line.
  subroutine write_csv_header(path, header)
    character(*), intent(in) :: path, header
    integer :: unit

    open (newunit=unit, file=path, status='replace', action='write')
    write (unit, '(a)') header
    close (unit)
  end subroutine write_csv_header

  !> Appends one already comma-joined row to an existing file.
  subroutine append_csv_row(path, row)
    character(*), intent(in) :: path, row
    integer :: unit

    open (newunit=unit, file=path, status='old', position='append', action='write')
    write (unit, '(a)') row
    close (unit)
  end subroutine append_csv_row

  !> Fixed-width scientific field (es24.16e3), the format used by every numeric column.
  function csv_real(value) result(text)
    real(real64), intent(in) :: value
    character(len=24) :: text

    write (text, '(es24.16e3)') value
  end function csv_real

  function csv_integer(value) result(text)
    integer, intent(in) :: value
    character(len=:), allocatable :: text
    character(len=16) :: buffer

    write (buffer, '(i0)') value
    text = trim(buffer)
  end function csv_integer

end module csv_writer
