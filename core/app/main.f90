program main
  use iso_fortran_env, only: real64, int64
  use harmonics, only: harmonic_transform
  use barotropic_vorticity, only: barotropic_solver
  use shallow_water, only: shallow_water_solver, gravity_acceleration, mean_depth, &
                           gravity_wave_implicitness
  use shallow_water_initial_conditions, only: isolated_height_mountain, &
                                                single_harmonic_height, &
                                                mountain_height_metres, &
                                                mountain_longitude_degrees, &
                                                mountain_latitude_degrees, &
                                                mountain_angular_radius_degrees, &
                                                height_mode_degree, height_mode_order, &
                                                height_mode_amplitude_metres
  use dry_atmosphere, only: dry_atmosphere_solver
  use barotropic_initial_conditions, only: single_harmonic_vorticity, &
                                           rossby_haurwitz_vorticity, &
                                           random_low_wavenumber_vorticity, &
                                           single_mode_degree, single_mode_order, &
                                           prescribed_maximum_vorticity, &
                                           rossby_haurwitz_wave_number, &
                                           rossby_haurwitz_omega, &
                                           rossby_haurwitz_wave_amplitude, &
                                           random_minimum_degree, random_maximum_degree
  use field_output, only: make_directory, write_snapshot, write_run_metadata, &
                          write_shallow_water_snapshot, write_shallow_water_metadata, &
                          write_dry_snapshot, write_dry_metadata
  implicit none

  integer, parameter :: T = 63
  real(real64), parameter :: dt = 900.0_real64
  real(real64), parameter :: duration = 10.0_real64*24.0_real64*3600.0_real64
  integer, parameter :: output_interval_steps = 4
  integer(int64), parameter :: random_seed_value = 20260913_int64
  type(harmonic_transform) :: transform
  integer, allocatable :: nlon(:)
  character(len=:), allocatable :: output_root
  character(len=64) :: equation_argument
  integer :: argument_count

  call transform%init(T)
  nlon = transform%get_nlon()
  output_root = find_output_root()
  call make_directory(output_root)

  argument_count = command_argument_count()
  if (argument_count == 0) then
    equation_argument = 'dry'
  else if (argument_count == 1) then
    call get_command_argument(1, equation_argument)
  else
    call print_usage()
    error stop 'too many command-line arguments'
  end if

  select case (trim(equation_argument))
  case ('shallow-water', '--shallow-water', 'shallow_water', 'swe')
    call run_shallow_water_case('shallow_water_mountain', 1)
    call run_shallow_water_case('shallow_water_single_harmonic', 2)
  case ('barotropic', '--barotropic', 'barotropic-vorticity', 'bve')
    call run_barotropic_case('single_harmonic', 1)
    call run_barotropic_case('rossby_haurwitz_r4', 2)
    call run_barotropic_case('random_n8_n12_seed_20260913', 3)
  case ('dry', '--dry', 'dry-atmosphere')
    call run_dry_case()
  case ('all')
    call run_dry_case()
    call run_shallow_water_case('shallow_water_mountain', 1)
    call run_shallow_water_case('shallow_water_single_harmonic', 2)
    call run_barotropic_case('single_harmonic', 1)
    call run_barotropic_case('rossby_haurwitz_r4', 2)
    call run_barotropic_case('random_n8_n12_seed_20260913', 3)
  case ('--help', '-h', 'help')
    call print_usage()
  case default
    call print_usage()
    error stop 'unknown equation argument'
  end select

contains

  subroutine run_barotropic_case(case_name, initial_condition)
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
      if (mod(step, output_interval_steps) == 0 .or. step == number_of_steps) then
        call write_snapshot(case_directory, step, nlon, zeta, u, v)
      end if
      if (mod(step, 24) == 0 .or. step == number_of_steps) then
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
                            duration, number_of_steps, output_interval_steps, maximum_cfl, &
                            elapsed_wall_seconds, nlon, transform%mu)
  end subroutine run_barotropic_case

  subroutine run_shallow_water_case(case_name, initial_condition)
    character(*), intent(in) :: case_name
    integer, intent(in) :: initial_condition
    type(shallow_water_solver) :: solver
    complex(real64), allocatable :: initial_zeta(:, :), initial_delta(:, :), initial_eta(:, :)
    complex(real64), allocatable :: zeta_spectral(:, :), delta_spectral(:, :), eta_spectral(:, :)
    real(real64), allocatable :: zeta(:, :), delta(:, :), eta(:, :), u(:, :), v(:, :)
    character(len=:), allocatable :: case_directory, initial_condition_json
    integer :: step, number_of_steps
    integer(int64) :: start_count, end_count, clock_rate, clock_max
    real(real64) :: cfl, maximum_cfl, elapsed_wall_seconds

    call system_clock(start_count, clock_rate, clock_max)
    select case (initial_condition)
    case (1)
      call isolated_height_mountain(transform, T, initial_zeta, initial_delta, initial_eta)
    case (2)
      call single_harmonic_height(transform, T, initial_zeta, initial_delta, initial_eta)
    case default
      error stop 'unknown shallow-water initial condition'
    end select

    case_directory = output_root//'/'//case_name
    call make_directory(case_directory)
    call solver%init(T, dt)
    call solver%set_initial_state(initial_zeta, initial_delta, initial_eta)
    number_of_steps = nint(duration/dt)
    maximum_cfl = 0.0_real64

    do step = 0, number_of_steps
      call solver%get_fields(zeta, delta, eta, u, v, cfl)
      maximum_cfl = max(maximum_cfl, cfl)
      if (mod(step, output_interval_steps) == 0 .or. step == number_of_steps) then
        call solver%get_spectral_state(zeta_spectral, delta_spectral, eta_spectral)
        call write_shallow_water_snapshot(case_directory, step, nlon, &
                                          zeta_spectral, delta_spectral, eta_spectral, &
                                          zeta, delta, eta, u, v)
      end if
      if (mod(step, 24) == 0 .or. step == number_of_steps) then
        write (*, '(2a,i0,a,f6.3)') trim(case_name), ': step ', step, ', advective CFL = ', cfl
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
    initial_condition_json = make_shallow_water_initial_condition_json(initial_condition)
    call write_shallow_water_metadata(case_directory, case_name, initial_condition_json, &
                                      T, dt, duration, number_of_steps, output_interval_steps, &
                                      maximum_cfl, elapsed_wall_seconds, nlon, transform%mu, &
                                      gravity_acceleration, mean_depth, gravity_wave_implicitness)
  end subroutine run_shallow_water_case

  subroutine run_dry_case()
    type(dry_atmosphere_solver) :: solver
    real(real64), allocatable :: zeta(:, :, :), delta(:, :, :), temperature(:, :, :)
    real(real64), allocatable :: surface_pressure(:, :), u(:, :, :), v(:, :, :)
    real(real64), allocatable :: pressure_half(:), delta_pressure(:), layer_l(:), alpha(:), reference_temperature(:)
    real(real64), allocatable :: a_half(:), b_half(:)
    character(len=:), allocatable :: case_directory
    integer :: step, number_of_steps
    integer(int64) :: start_count, end_count, clock_rate, clock_max
    real(real64) :: cfl, maximum_cfl, elapsed_wall_seconds

    call system_clock(start_count, clock_rate, clock_max)
    case_directory = output_root//'/dry_jablonowski_williamson'
    call make_directory(case_directory)
    call solver%init(T, dt)
    call solver%set_jablonowski_williamson_state(.true.)
    number_of_steps = nint(duration/dt)
    maximum_cfl = 0.0_real64
    do step = 0, number_of_steps
      call solver%get_fields(zeta, delta, temperature, surface_pressure, u, v, cfl)
      maximum_cfl = max(maximum_cfl, cfl)
      if (mod(step, output_interval_steps) == 0 .or. step == number_of_steps) then
        call write_dry_snapshot(case_directory, step, nlon, zeta, delta, temperature, surface_pressure, u, v)
      end if
      if (mod(step, 24) == 0 .or. step == number_of_steps) then
        write (*, '(a,i0,a,f6.3)') 'dry_jablonowski_williamson: step ', step, ', advective CFL = ', cfl
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
    call solver%get_reference_atmosphere(pressure_half, delta_pressure, layer_l, alpha, reference_temperature, &
                                         a_half, b_half)
    call write_dry_metadata(case_directory, T, dt, duration, number_of_steps, output_interval_steps, &
                            maximum_cfl, elapsed_wall_seconds, nlon, transform%mu, &
                            pressure_half, delta_pressure, layer_l, alpha, reference_temperature, a_half, b_half)
  end subroutine run_dry_case

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

  function make_shallow_water_initial_condition_json(initial_condition) result(json)
    integer, intent(in) :: initial_condition
    character(len=:), allocatable :: json
    character(len=768) :: buffer

    select case (initial_condition)
    case (1)
      write (buffer, '(a,es24.16e3,a,es24.16e3,a,es24.16e3,a,es24.16e3,a)') &
        '{"type":"isolated_gaussian_height_mountain","height_m":', mountain_height_metres, &
        ',"longitude_degrees":', mountain_longitude_degrees, &
        ',"latitude_degrees":', mountain_latitude_degrees, &
        ',"angular_radius_degrees":', mountain_angular_radius_degrees, '}'
    case (2)
      write (buffer, '(a,i0,a,i0,a,es24.16e3,a)') &
        '{"type":"single_spherical_harmonic_height","degree_n":', height_mode_degree, &
        ',"zonal_order_m":', height_mode_order, &
        ',"maximum_absolute_height_m":', height_mode_amplitude_metres, '}'
    case default
      error stop 'unknown shallow-water initial condition metadata'
    end select
    json = trim(buffer)
  end function make_shallow_water_initial_condition_json

  subroutine print_usage()
    write (*, '(a)') 'Usage: core [shallow-water|barotropic|dry|all]'
    write (*, '(a)') '  shallow-water:           run the mountain and single-harmonic height cases'
    write (*, '(a)') '  barotropic:             run the three barotropic-vorticity cases'
    write (*, '(a)') '  dry (default):          run the 10-day Jablonowski-Williamson dry-atmosphere case'
    write (*, '(a)') '  all:                    run every case'
  end subroutine print_usage

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
end program main
