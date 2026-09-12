program main
  use iso_fortran_env, only: real64, int64
  use harmonics, only: harmonic_transform
  use barotropic_vorticity, only: barotropic_solver
  use barotropic_initial_conditions, only: single_harmonic_vorticity, &
                                           rossby_haurwitz_vorticity, &
                                           random_low_wavenumber_vorticity, &
                                           single_mode_degree, single_mode_order, &
                                           prescribed_maximum_vorticity, &
                                           rossby_haurwitz_wave_number, &
                                           rossby_haurwitz_omega, &
                                           rossby_haurwitz_wave_amplitude, &
                                           random_minimum_degree, random_maximum_degree
  use field_output, only: make_directory, write_snapshot, write_run_metadata
  implicit none

  integer, parameter :: T = 63
  real(real64), parameter :: dt = 450.0_real64
  real(real64), parameter :: duration = 24.0_real64*3600.0_real64
  integer(int64), parameter :: random_seed_value = 20260913_int64
  type(harmonic_transform) :: transform
  integer, allocatable :: nlon(:)
  character(len=:), allocatable :: output_root

  call transform%init(T)
  nlon = transform%get_nlon()
  output_root = find_output_root()
  call make_directory(output_root)

  call run_case('single_harmonic', 1)
  call run_case('rossby_haurwitz_r4', 2)
  call run_case('random_n8_n12_seed_20260913', 3)

contains

  subroutine run_case(case_name, initial_condition)
    character(*), intent(in) :: case_name
    integer, intent(in) :: initial_condition
    type(barotropic_solver) :: solver
    complex(real64), allocatable :: initial_zeta(:, :)
    real(real64), allocatable :: zeta(:, :), u(:, :), v(:, :)
    character(len=:), allocatable :: case_directory, initial_condition_json
    integer :: step, number_of_steps
    integer(int64) :: start_count, end_count, clock_rate, clock_max
    real(real64) :: cfl, maximum_cfl, elapsed_wall_seconds

    call system_clock(start_count, clock_rate, clock_max)
    select case (initial_condition)
    case (1)
      call single_harmonic_vorticity(transform, T, initial_zeta)
    case (2)
      call rossby_haurwitz_vorticity(transform, T, initial_zeta)
    case (3)
      call random_low_wavenumber_vorticity(transform, T, random_seed_value, initial_zeta)
    case default
      error stop 'unknown initial condition'
    end select

    case_directory = output_root//'/'//case_name
    call make_directory(case_directory)
    call solver%init(T, dt)
    call solver%set_initial_vorticity(initial_zeta)
    number_of_steps = nint(duration/dt)
    maximum_cfl = 0.0_real64

    do step = 0, number_of_steps
      call solver%get_fields(zeta, u, v, cfl)
      maximum_cfl = max(maximum_cfl, cfl)
      call write_snapshot(case_directory, step, nlon, zeta, u, v)
      if (mod(step, 12) == 0 .or. step == number_of_steps) then
        write (*, '(2a,i0,a,f6.3)') trim(case_name), ': step ', step, ', CFL = ', cfl
      end if
      if (step < number_of_steps) call solver%advance()
    end do
    call system_clock(end_count)
    if (end_count >= start_count) then
      elapsed_wall_seconds = real(end_count - start_count, real64)/real(clock_rate, real64)
    else
      elapsed_wall_seconds = real(clock_max - start_count + end_count + 1_int64, real64)/ &
                             real(clock_rate, real64)
    end if
    initial_condition_json = make_initial_condition_json(initial_condition)
    call write_run_metadata(case_directory, case_name, initial_condition_json, T, dt, &
                            duration, number_of_steps, maximum_cfl, elapsed_wall_seconds, &
                            nlon, transform%mu)
  end subroutine run_case

  function make_initial_condition_json(initial_condition) result(json)
    integer, intent(in) :: initial_condition
    character(len=:), allocatable :: json
    character(len=512) :: buffer

    select case (initial_condition)
    case (1)
      write (buffer, '(a,i0,a,i0,a,es24.16e3,a)') &
        '{"type":"single_spherical_harmonic","degree_n":', single_mode_degree, &
        ',"zonal_order_m":', single_mode_order, &
        ',"maximum_absolute_vorticity_s-1":', prescribed_maximum_vorticity, '}'
    case (2)
      write (buffer, '(a,i0,a,es24.16e3,a,es24.16e3,a)') &
        '{"type":"rossby_haurwitz","R":', rossby_haurwitz_wave_number, &
        ',"omega_s-1":', rossby_haurwitz_omega, &
        ',"wave_amplitude_s-1":', rossby_haurwitz_wave_amplitude, '}'
    case (3)
      write (buffer, '(a,i0,a,i0,a,i0,a,es24.16e3,a)') &
        '{"type":"random_low_wavenumber_vorticity","minimum_degree_n":', &
        random_minimum_degree, ',"maximum_degree_n":', random_maximum_degree, &
        ',"seed":', random_seed_value, ',"maximum_absolute_vorticity_s-1":', &
        prescribed_maximum_vorticity, '}'
    case default
      error stop 'unknown initial condition metadata'
    end select
    json = trim(buffer)
  end function make_initial_condition_json

  function find_output_root() result(path)
    character(len=:), allocatable :: path
    logical :: exists

    inquire (file='docs/barotropic-vorticity-equation.md', exist=exists)
    if (exists) then
      path = 'output'
      return
    end if
    inquire (file='../docs/barotropic-vorticity-equation.md', exist=exists)
    if (exists) then
      path = '../output'
    else
      path = 'output'
    end if
  end function find_output_root
end program main
