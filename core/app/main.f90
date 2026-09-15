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
                          write_dry_snapshot, write_dry_metadata, &
                          initialize_radiation_daily_output, append_radiation_daily_output, &
                          write_radiation_monthly_output, write_radiation_yearly_snapshot, &
                          write_radiation_metadata
  use dry_radiation, only: solar_day, days_per_month, months_per_year, days_per_year, orbital_period
  implicit none

  integer, parameter :: T = 63
  real(real64), parameter :: dt = 900.0_real64
  real(real64), parameter :: duration = 10.0_real64*24.0_real64*3600.0_real64
  integer, parameter :: output_interval_steps = 16
  real(real64), parameter :: held_suarez_duration = 200.0_real64*24.0_real64*3600.0_real64
  integer, parameter :: held_suarez_output_interval_steps = nint(5.0_real64*24.0_real64*3600.0_real64/dt)
  real(real64), parameter :: radiation_time_step = 1200.0_real64
  integer, parameter :: radiation_number_of_years = 5
  real(real64), parameter :: radiation_duration = real(radiation_number_of_years, real64)*orbital_period
  !> Progress is logged once per simulated day (and at the final step).
  integer, parameter :: log_interval_steps = max(1, nint(86400.0_real64/dt))
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
    call run_dry_case('dry_jablonowski_williamson_steady', .false.)
    call run_dry_case('dry_jablonowski_williamson_perturbed', .true.)
  case ('held-suarez', '--held-suarez', 'held_suarez')
    call run_dry_case('dry_held_suarez', .false., held_suarez=.true.)
  case ('radiation', '--radiation', 'uniform-radiation', '--uniform-radiation', 'uniform_radiation')
    call run_radiation_case()
  case ('all')
    call run_dry_case('dry_jablonowski_williamson_steady', .false.)
    call run_dry_case('dry_jablonowski_williamson_perturbed', .true.)
    call run_shallow_water_case('shallow_water_mountain', 1)
    call run_shallow_water_case('shallow_water_single_harmonic', 2)
    call run_barotropic_case('single_harmonic', 1)
    call run_barotropic_case('rossby_haurwitz_r4', 2)
    call run_barotropic_case('random_n8_n12_seed_20260913', 3)
    call run_dry_case('dry_held_suarez', .false., held_suarez=.true.)
    call run_radiation_case()
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
    integer(int64) :: start_count
    real(real64) :: cfl, maximum_cfl, elapsed_wall_seconds

    call system_clock(start_count)
    call write_case_header(case_name, duration, dt)
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
      if (is_log_step(step, number_of_steps)) then
        call write_progress(step, number_of_steps, step, 'CFL', cfl, start_count, dt)
      end if
      if (step < number_of_steps) call solver%advance()
    end do
    elapsed_wall_seconds = elapsed_seconds(start_count)
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
    integer(int64) :: start_count
    real(real64) :: cfl, maximum_cfl, elapsed_wall_seconds

    call system_clock(start_count)
    call write_case_header(case_name, duration, dt)
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
      if (is_log_step(step, number_of_steps)) then
        call write_progress(step, number_of_steps, step, 'advective CFL', cfl, start_count, dt)
      end if
      if (step < number_of_steps) call solver%advance()
    end do
    elapsed_wall_seconds = elapsed_seconds(start_count)
    initial_condition_json = make_shallow_water_initial_condition_json(initial_condition)
    call write_shallow_water_metadata(case_directory, case_name, initial_condition_json, &
                                      T, dt, duration, number_of_steps, output_interval_steps, &
                                      maximum_cfl, elapsed_wall_seconds, nlon, transform%mu, &
                                      gravity_acceleration, mean_depth, gravity_wave_implicitness)
  end subroutine run_shallow_water_case

  subroutine run_dry_case(case_name, include_perturbation, held_suarez)
    character(*), intent(in) :: case_name
    !> .false. keeps the balanced zonal base state (should stay steady);
    !> .true. adds the localized wind perturbation that triggers the baroclinic wave.
    logical, intent(in) :: include_perturbation
    logical, intent(in), optional :: held_suarez
    type(dry_atmosphere_solver) :: solver
    real(real64), allocatable :: zeta(:, :, :), delta(:, :, :), temperature(:, :, :)
    real(real64), allocatable :: surface_pressure(:, :), u(:, :, :), v(:, :, :)
    real(real64), allocatable :: pressure_half(:), delta_pressure(:), layer_l(:), alpha(:), reference_temperature(:)
    real(real64), allocatable :: a_half(:), b_half(:)
    character(len=:), allocatable :: case_directory, initial_condition
    integer :: step, number_of_steps, snapshot_interval_steps
    integer(int64) :: start_count
    real(real64) :: cfl, maximum_cfl, elapsed_wall_seconds, case_duration
    logical :: use_held_suarez

    call system_clock(start_count)
    use_held_suarez = .false.
    if (present(held_suarez)) use_held_suarez = held_suarez
    if (use_held_suarez) then
      initial_condition = 'Held-Suarez resting isothermal atmosphere with thermal and Rayleigh forcing'
      case_duration = held_suarez_duration
      snapshot_interval_steps = held_suarez_output_interval_steps
    else if (include_perturbation) then
      initial_condition = 'Jablonowski-Williamson with localized wind perturbation'
      case_duration = duration
      snapshot_interval_steps = output_interval_steps
    else
      initial_condition = 'Jablonowski-Williamson balanced base state without perturbation'
      case_duration = duration
      snapshot_interval_steps = output_interval_steps
    end if
    call write_case_header(case_name, case_duration, dt)
    case_directory = output_root//'/'//case_name
    call make_directory(case_directory)
    call solver%init(T, dt)
    if (use_held_suarez) then
      call solver%set_held_suarez_state()
    else
      call solver%set_jablonowski_williamson_state(include_perturbation)
    end if
    number_of_steps = nint(case_duration/dt)
    maximum_cfl = 0.0_real64
    do step = 0, number_of_steps
      ! Grid fields are synthesized only for snapshots.  On every other step the CFL of
      ! the current state comes from the grid winds that advance already builds.
      if (mod(step, snapshot_interval_steps) == 0 .or. step == number_of_steps) then
        call solver%get_fields(zeta, delta, temperature, surface_pressure, u, v, cfl)
        maximum_cfl = max(maximum_cfl, cfl)
        call write_dry_snapshot(case_directory, step, nlon, zeta, delta, temperature, surface_pressure, u, v)
      end if
      if (step < number_of_steps) then
        call solver%advance()
        cfl = solver%get_last_advance_cfl()
        maximum_cfl = max(maximum_cfl, cfl)
      end if
      ! Logged after advance, so the elapsed time already covers step + 1 completed steps.
      if (is_log_step(step, number_of_steps)) then
        call write_progress(step, number_of_steps, min(step + 1, number_of_steps), &
                            'advective CFL', cfl, start_count, dt)
      end if
    end do
    elapsed_wall_seconds = elapsed_seconds(start_count)
    call solver%get_reference_atmosphere(pressure_half, delta_pressure, layer_l, alpha, reference_temperature, &
                                         a_half, b_half)
    call write_dry_metadata(case_directory, initial_condition, T, dt, case_duration, number_of_steps, &
                            snapshot_interval_steps, &
                            maximum_cfl, elapsed_wall_seconds, nlon, transform%mu, &
                            pressure_half, delta_pressure, layer_l, alpha, reference_temperature, a_half, b_half)
  end subroutine run_dry_case

  subroutine run_radiation_case()
    type(dry_atmosphere_solver) :: solver
    real(real64), allocatable :: monthly_surface_temperature(:, :), monthly_surface_pressure(:, :)
    real(real64), allocatable :: zonal_temperature(:, :), zonal_u(:, :), zonal_v(:, :)
    real(real64), allocatable :: eddy_uv(:, :), eddy_vt(:, :)
    real(real64), allocatable :: pressure_half(:), delta_pressure(:), layer_l(:), alpha(:)
    real(real64), allocatable :: reference_temperature(:), a_half(:), b_half(:)
    character(len=:), allocatable :: case_directory
    integer :: completed_step, number_of_steps, daily_interval_steps
    integer :: month, month_boundary_step, year, year_boundary_step
    integer(int64) :: start_count
    real(real64) :: cfl, maximum_cfl, elapsed_wall_seconds, diagnostic_time
    real(real64) :: mean_atmospheric_temperature, mean_surface_temperature, mean_deep_temperature
    real(real64) :: mean_kinetic_energy, mean_surface_pressure
    real(real64) :: mean_incoming_shortwave, mean_reflected_shortwave, mean_outgoing_longwave

    call system_clock(start_count)
    call write_case_header('dry_radiation', radiation_duration, radiation_time_step)
    case_directory = output_root//'/dry_radiation'
    call make_directory(case_directory)
    call initialize_radiation_daily_output(case_directory)
    call solver%init(T, radiation_time_step)
    call solver%set_radiation_state()
    call write_current_radiation_snapshot(solver, case_directory, 1, nlon)
    number_of_steps = nint(radiation_duration/radiation_time_step)
    daily_interval_steps = nint(solar_day/radiation_time_step)
    if (abs(real(daily_interval_steps, real64)*radiation_time_step - solar_day) > 1.0e-12_real64) then
      error stop 'radiation time step must divide the solar day exactly'
    end if
    maximum_cfl = 0.0_real64
    month = 1
    month_boundary_step = month*days_per_month*daily_interval_steps
    year = 2
    year_boundary_step = (year - 1)*days_per_year*daily_interval_steps

    do completed_step = 1, number_of_steps
      call solver%advance()
      cfl = solver%get_last_advance_cfl()
      maximum_cfl = max(maximum_cfl, cfl)

      ! Each row is the mean over the steps of one solar day, stamped with the day's start time.
      if (mod(completed_step, daily_interval_steps) == 0) then
        call solver%take_radiation_daily_means(diagnostic_time, mean_atmospheric_temperature, &
          mean_surface_temperature, mean_deep_temperature, mean_kinetic_energy, mean_surface_pressure, &
          mean_incoming_shortwave, mean_reflected_shortwave, mean_outgoing_longwave)
        call append_radiation_daily_output(case_directory, diagnostic_time, mean_atmospheric_temperature, &
          mean_surface_temperature, mean_deep_temperature, mean_kinetic_energy, mean_surface_pressure, &
          mean_incoming_shortwave, mean_reflected_shortwave, mean_outgoing_longwave)
      end if

      if (completed_step == month_boundary_step) then
        call solver%take_radiation_monthly_means(monthly_surface_temperature, monthly_surface_pressure, &
          zonal_temperature, zonal_u, zonal_v, eddy_uv, eddy_vt)
        call write_radiation_monthly_output(case_directory, month, nlon, monthly_surface_temperature, &
          monthly_surface_pressure, zonal_temperature, zonal_u, zonal_v, eddy_uv, eddy_vt)
        month = month + 1
        if (month <= radiation_number_of_years*months_per_year) then
          month_boundary_step = month*days_per_month*daily_interval_steps
        end if
      end if

      if (completed_step == year_boundary_step) then
        call write_current_radiation_snapshot(solver, case_directory, year, nlon)
        year = year + 1
        if (year <= radiation_number_of_years + 1) then
          year_boundary_step = (year - 1)*days_per_year*daily_interval_steps
        end if
      end if

      if (mod(completed_step, daily_interval_steps) == 0 .or. completed_step == number_of_steps) then
        call write_progress(completed_step, number_of_steps, completed_step, &
                            'advective CFL', cfl, start_count, radiation_time_step)
      end if
    end do

    elapsed_wall_seconds = elapsed_seconds(start_count)
    call solver%get_reference_atmosphere(pressure_half, delta_pressure, layer_l, alpha, &
                                         reference_temperature, a_half, b_half)
    call write_radiation_metadata(case_directory, T, radiation_time_step, radiation_duration, number_of_steps, &
                                  maximum_cfl, elapsed_wall_seconds, nlon, transform%mu, pressure_half, &
                                  delta_pressure, layer_l, alpha, reference_temperature, a_half, b_half)
  end subroutine run_radiation_case

  subroutine write_current_radiation_snapshot(solver, case_directory, year, ring_nlon)
    type(dry_atmosphere_solver), intent(inout) :: solver
    character(*), intent(in) :: case_directory
    integer, intent(in) :: year, ring_nlon(:)
    complex(real64), allocatable :: zeta_spectral(:, :, :), delta_spectral(:, :, :)
    complex(real64), allocatable :: temperature_spectral(:, :, :), log_ps_spectral(:, :)
    real(real64), allocatable :: zeta(:, :, :), delta(:, :, :), temperature(:, :, :)
    real(real64), allocatable :: surface_pressure(:, :), log_surface_pressure(:, :)
    real(real64), allocatable :: u(:, :, :), v(:, :, :)
    real(real64), allocatable :: surface_temperature(:, :), deep_temperature(:, :)

    call solver%get_spectral_state(zeta_spectral, delta_spectral, temperature_spectral, log_ps_spectral)
    call solver%get_fields(zeta, delta, temperature, surface_pressure, u, v, &
                           surface_temperature=surface_temperature, deep_temperature=deep_temperature)
    log_surface_pressure = log(surface_pressure)
    call write_radiation_yearly_snapshot(case_directory, year, ring_nlon, &
      zeta_spectral, delta_spectral, temperature_spectral, log_ps_spectral, &
      zeta, delta, temperature, u, v, log_surface_pressure, surface_temperature, deep_temperature)
  end subroutine write_current_radiation_snapshot

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

  subroutine print_usage()
    write (*, '(a)') 'Usage: core [shallow-water|barotropic|dry|held-suarez|radiation|all]'
    write (*, '(a)') '  shallow-water:           run the mountain and single-harmonic height cases'
    write (*, '(a)') '  barotropic:             run the three barotropic-vorticity cases'
    write (*, '(a)') '  dry (default):          run the 10-day Jablonowski-Williamson dry-atmosphere cases'
    write (*, '(a)') '                          (steady base state and localized wind perturbation)'
    write (*, '(a)') '  held-suarez:            run the 200-day forced dry-atmosphere case (output every 5 days)'
    write (*, '(a)') '  radiation:              run the 5-year diurnal/seasonal radiation case'
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
