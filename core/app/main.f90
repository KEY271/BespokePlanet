program main
  use case_runners, only: run_requested_cases
  implicit none

  character(len=64) :: requested_case
  integer :: argument_count

  argument_count = command_argument_count()
  if (argument_count == 0) then
    requested_case = 'dry'
  else if (argument_count == 1) then
    call get_command_argument(1, requested_case)
  else
    call run_requested_cases('help')
    error stop 'too many command-line arguments'
  end if

  call run_requested_cases(trim(requested_case))
end program main
