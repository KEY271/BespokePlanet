!> The barotropic-vorticity cases: initial condition, run length and snapshot schedule.
module barotropic_case
  use iso_fortran_env, only: real64, int64
  use planet_parameters, only: earth_radius, earth_rotation_rate
  use numerics_config, only: model_numerics_config
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
  use filesystem, only: make_directory
  use barotropic_case_output, only: write_barotropic_snapshot, write_barotropic_metadata
  use case_runtime, only: case_context, ensure_context, T, duration, output_interval_steps, &
                          write_case_header, is_log_step, write_progress, elapsed_seconds
  implicit none
  private
  public :: run_barotropic_case

  integer(int64), parameter :: random_seed_value = 20260913_int64

contains

  subroutine run_barotropic_case(context, case_name, initial_condition)
    type(case_context), intent(inout) :: context
    character(*), intent(in) :: case_name
    integer, intent(in) :: initial_condition
    type(barotropic_solver) :: solver
    type(model_numerics_config) :: numerics
    complex(real64), allocatable :: initial_zeta(:, :)
    real(real64), allocatable :: zeta(:, :), u(:, :), v(:, :)
    character(len=:), allocatable :: case_directory, initial_condition_json
    integer :: step, number_of_steps
    integer(int64) :: start_count
    real(real64) :: cfl, maximum_cfl, elapsed_wall_seconds

    call ensure_context(context)
    numerics = context%numerics
    call system_clock(start_count)
    call write_case_header(case_name, duration, numerics%time_step)
    select case (initial_condition)
    case (1)
      call single_harmonic_vorticity(context%transform, T, initial_zeta)
    case (2)
      call rossby_haurwitz_vorticity(context%transform, T, initial_zeta)
    case (3)
      call random_low_wavenumber_vorticity(context%transform, T, random_seed_value, initial_zeta)
    case default
      error stop 'unknown initial condition'
    end select

    case_directory = context%output_root//'/'//case_name
    call make_directory(case_directory)
    call solver%init_with_config(numerics)
    call solver%set_initial_vorticity(initial_zeta)
    number_of_steps = nint(duration/numerics%time_step)
    maximum_cfl = 0.0_real64

    do step = 0, number_of_steps
      call solver%get_fields(zeta, u, v, cfl)
      maximum_cfl = max(maximum_cfl, cfl)
      if (mod(step, output_interval_steps) == 0 .or. step == number_of_steps) then
        call write_barotropic_snapshot(case_directory, step, context%nlon, zeta, u, v)
      end if
      if (is_log_step(step, number_of_steps)) then
        call write_progress(step, number_of_steps, step, 'CFL', cfl, start_count, numerics%time_step)
      end if
      if (step < number_of_steps) call solver%advance()
    end do
    elapsed_wall_seconds = elapsed_seconds(start_count)
    initial_condition_json = make_initial_condition_json(initial_condition)
    call write_barotropic_metadata(case_directory, case_name, initial_condition_json, numerics, &
                                   duration, number_of_steps, output_interval_steps, maximum_cfl, &
                                   elapsed_wall_seconds, context%nlon, context%transform%mu, &
                                   earth_radius, earth_rotation_rate)
  end subroutine run_barotropic_case

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

end module barotropic_case
