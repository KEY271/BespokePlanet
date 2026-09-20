!> The shallow-water cases: initial condition, run length and snapshot schedule.
module shallow_water_case
  use iso_fortran_env, only: real64, int64
  use planet_parameters, only: earth_radius, earth_rotation_rate
  use numerics_config, only: model_numerics_config
  use shallow_water_config, only: shallow_water_equation_config
  use shallow_water, only: shallow_water_solver
  use shallow_water_initial_conditions, only: isolated_height_mountain, &
                                                single_harmonic_height, &
                                                mountain_height_metres, &
                                                mountain_longitude_degrees, &
                                                mountain_latitude_degrees, &
                                                mountain_angular_radius_degrees, &
                                                height_mode_degree, height_mode_order, &
                                                height_mode_amplitude_metres
  use filesystem, only: make_directory
  use shallow_water_case_output, only: write_shallow_water_snapshot, write_shallow_water_metadata
  use case_runtime, only: case_context, ensure_context, T, duration, output_interval_steps, &
                          write_case_header, is_log_step, write_progress, elapsed_seconds
  implicit none
  private
  public :: run_shallow_water_case

contains

  subroutine run_shallow_water_case(context, case_name, initial_condition)
    type(case_context), intent(inout) :: context
    character(*), intent(in) :: case_name
    integer, intent(in) :: initial_condition
    type(shallow_water_solver) :: solver
    type(model_numerics_config) :: numerics
    type(shallow_water_equation_config) :: equation
    complex(real64), allocatable :: initial_zeta(:, :), initial_delta(:, :), initial_eta(:, :)
    complex(real64), allocatable :: zeta_spectral(:, :), delta_spectral(:, :), eta_spectral(:, :)
    real(real64), allocatable :: zeta(:, :), delta(:, :), eta(:, :), u(:, :), v(:, :)
    character(len=:), allocatable :: case_directory, initial_condition_json
    integer :: step, number_of_steps
    integer(int64) :: start_count
    real(real64) :: cfl, maximum_cfl, elapsed_wall_seconds

    call ensure_context(context)
    numerics = context%numerics
    equation = shallow_water_equation_config()
    call system_clock(start_count)
    call write_case_header(case_name, duration, numerics%time_step)
    select case (initial_condition)
    case (1)
      call isolated_height_mountain(context%transform, T, initial_zeta, initial_delta, initial_eta)
    case (2)
      call single_harmonic_height(context%transform, T, initial_zeta, initial_delta, initial_eta)
    case default
      error stop 'unknown shallow-water initial condition'
    end select

    case_directory = context%output_root//'/'//case_name
    call make_directory(case_directory)
    call solver%init_with_config(numerics, equation)
    call solver%set_initial_state(initial_zeta, initial_delta, initial_eta)
    number_of_steps = nint(duration/numerics%time_step)
    maximum_cfl = 0.0_real64

    do step = 0, number_of_steps
      call solver%get_fields(zeta, delta, eta, u, v, cfl)
      maximum_cfl = max(maximum_cfl, cfl)
      if (mod(step, output_interval_steps) == 0 .or. step == number_of_steps) then
        call solver%get_spectral_state(zeta_spectral, delta_spectral, eta_spectral)
        call write_shallow_water_snapshot(case_directory, step, context%nlon, &
                                          zeta_spectral, delta_spectral, eta_spectral, &
                                          zeta, delta, eta, u, v)
      end if
      if (is_log_step(step, number_of_steps)) then
        call write_progress(step, number_of_steps, step, 'advective CFL', cfl, start_count, numerics%time_step)
      end if
      if (step < number_of_steps) call solver%advance()
    end do
    elapsed_wall_seconds = elapsed_seconds(start_count)
    initial_condition_json = make_shallow_water_initial_condition_json(initial_condition)
    call write_shallow_water_metadata(case_directory, case_name, initial_condition_json, &
                                      numerics, equation, duration, number_of_steps, output_interval_steps, &
                                      maximum_cfl, elapsed_wall_seconds, context%nlon, context%transform%mu, &
                                      earth_radius, earth_rotation_rate)
  end subroutine run_shallow_water_case

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

end module shallow_water_case
