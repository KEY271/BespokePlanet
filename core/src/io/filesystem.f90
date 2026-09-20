module filesystem
  implicit none
  private
  public :: make_directory

contains

  subroutine make_directory(path)
    character(*), intent(in) :: path
    integer :: exit_status
    call execute_command_line('mkdir -p "'//trim(path)//'"', exitstat=exit_status)
    if (exit_status /= 0) error stop 'failed to create output directory'
  end subroutine make_directory

end module filesystem
