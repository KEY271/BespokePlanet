!> Resources and helpers shared by every case runner: the spectral transform and
!> output root of one command invocation, the common resolution and time step,
!> and the progress log.  No case-specific configuration lives here.
module case_runtime
  use iso_fortran_env, only: real64, int64
  use harmonics, only: harmonic_transform
  use numerics_config, only: model_numerics_config
  use filesystem, only: make_directory
  implicit none
  private

  integer, parameter, public :: T = 63
  real(real64), parameter, public :: dt = 1200.0_real64
  real(real64), parameter, public :: duration = 10.0_real64*24.0_real64*3600.0_real64
  integer, parameter, public :: output_interval_steps = 16
  !> Progress is logged once per simulated day (and at the final step).
  integer, parameter :: log_interval_steps = max(1, nint(86400.0_real64/dt))

  !> Shared, read-mostly resources of one command invocation.  Each case derives its
  !> own configuration from `numerics` instead of modifying it.
  type, public :: case_context
    logical :: initialized = .false.
    type(harmonic_transform) :: transform
    integer, allocatable :: nlon(:)
    character(len=:), allocatable :: output_root
    type(model_numerics_config) :: numerics
  end type case_context

  public :: ensure_context, write_case_header, is_log_step, write_progress, elapsed_seconds

contains

  !> Builds the transform and output directory on first use, so that usage
  !> messages do not initialize the model or touch the file system.
  subroutine ensure_context(context)
    type(case_context), intent(inout) :: context

    if (context%initialized) return
    call context%transform%init(T)
    context%nlon = context%transform%get_nlon()
    context%output_root = find_output_root()
    call make_directory(context%output_root)
    context%numerics = model_numerics_config()
    context%numerics%truncation = T
    context%numerics%time_step = dt
    context%initialized = .true.
  end subroutine ensure_context

  !> Printed once when a case starts.
  subroutine write_case_header(case_name, case_duration, time_step)
    character(*), intent(in) :: case_name
    real(real64), intent(in) :: case_duration, time_step

    write (*, '(3a,f0.1,a,i0,a,f0.1,a)') '== ', trim(case_name), ': ', case_duration/86400.0_real64, &
      ' days, ', nint(case_duration/time_step), ' steps, dt = ', time_step, ' s'
  end subroutine write_case_header

  logical function is_log_step(step, number_of_steps)
    integer, intent(in) :: step, number_of_steps

    is_log_step = mod(step, log_interval_steps) == 0 .or. step == number_of_steps
  end function is_log_step

  !> One progress line: simulated day, CFL, wall time so far and a linear estimate of the
  !> remaining wall time (elapsed per completed step times the steps still to go).
  subroutine write_progress(step, number_of_steps, completed_steps, cfl_label, cfl, start_count, time_step)
    integer, intent(in) :: step, number_of_steps, completed_steps
    character(*), intent(in) :: cfl_label
    real(real64), intent(in) :: cfl
    real(real64), intent(in) :: time_step
    integer(int64), intent(in) :: start_count
    character(len=:), allocatable :: remaining_text
    real(real64) :: elapsed

    elapsed = elapsed_seconds(start_count)
    if (completed_steps >= number_of_steps) then
      remaining_text = 'done'
    else if (step > 0) then
      ! The first step also carries the setup cost, so no estimate is made before day 1.
      remaining_text = 'remaining ~'//format_duration(elapsed/real(completed_steps, real64)* &
                                                      real(number_of_steps - completed_steps, real64))
    else
      remaining_text = 'remaining --'
    end if
    write (*, '(a,f7.1,3a,f5.3,4a)') '  day ', real(step, real64)*time_step/86400.0_real64, &
      '  ', cfl_label, ' = ', cfl, '  elapsed ', format_duration(elapsed), '  ', remaining_text
  end subroutine write_progress

  function format_duration(seconds) result(text)
    real(real64), intent(in) :: seconds
    character(len=:), allocatable :: text
    character(len=32) :: buffer
    integer :: total

    if (seconds < 59.95_real64) then
      write (buffer, '(f4.1,a)') seconds, ' s'
    else
      total = nint(seconds)
      if (total >= 3600) then
        write (buffer, '(i0,a,i2.2,a,i2.2,a)') total/3600, 'h', mod(total, 3600)/60, 'm', mod(total, 60), 's'
      else
        write (buffer, '(i0,a,i2.2,a)') total/60, 'm', mod(total, 60), 's'
      end if
    end if
    text = trim(adjustl(buffer))
  end function format_duration

  real(real64) function elapsed_seconds(start_count) result(seconds)
    integer(int64), intent(in) :: start_count
    integer(int64) :: now, rate, maximum

    call system_clock(now, rate, maximum)
    if (now >= start_count) then
      seconds = real(now - start_count, real64)/real(rate, real64)
    else
      seconds = real(maximum - start_count + now + 1_int64, real64)/real(rate, real64)
    end if
  end function elapsed_seconds

  function find_output_root() result(path)
    character(len=:), allocatable :: path
    logical :: exists

    inquire (file='docs/shallow-water-equation.md', exist=exists)
    if (exists) then
      path = 'output'
      return
    end if
    inquire (file='../docs/shallow-water-equation.md', exist=exists)
    if (exists) then
      path = '../output'
    else
      path = 'output'
    end if
  end function find_output_root

end module case_runtime
