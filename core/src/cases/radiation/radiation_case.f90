!> The dry radiation case: five calendar years with daily global means, monthly
!> zonal means and yearly snapshots.  The case chooses the physics, the planet,
!> and the calendar boundaries at which each output is written.
module radiation_case
  use iso_fortran_env, only: real64, int64
  use planet_parameters, only: planet_config
  use numerics_config, only: model_numerics_config
  use dry_physics_config, only: dry_model_physics_config, radiation_days_per_year, radiation_orbital_period
  use dry_atmosphere, only: dry_atmosphere_solver
  use dry_case_initial_conditions, only: radiation_case_physics, radiation_case_planet, &
                                         set_radiation_case_state
  use dry_radiation, only: radiation_diagnostics, radiation_calendar_date
  use radiation_diagnostics_collector, only: radiation_case_diagnostics
  use filesystem, only: make_directory
  use radiation_case_output, only: initialize_radiation_daily_output, append_radiation_daily_output, &
                                    write_radiation_monthly_output, write_radiation_yearly_snapshot, &
                                    write_radiation_metadata
  use case_runtime, only: case_context, ensure_context, dt, write_case_header, write_progress, &
                          elapsed_seconds
  implicit none
  private
  public :: run_radiation_case

  real(real64), parameter :: radiation_time_step = dt
  integer, parameter :: radiation_number_of_years = 5

contains

  subroutine run_radiation_case(context)
    type(case_context), intent(inout) :: context
    type(dry_atmosphere_solver) :: solver
    type(model_numerics_config) :: numerics
    type(dry_model_physics_config) :: physics
    type(planet_config) :: planet
    type(radiation_case_diagnostics) :: diagnostics
    type(radiation_diagnostics) :: sample, daily_means
    real(real64), allocatable :: monthly_surface_temperature(:, :), monthly_surface_pressure(:, :)
    real(real64), allocatable :: zonal_temperature(:, :), zonal_u(:, :), zonal_v(:, :)
    real(real64), allocatable :: eddy_uv(:, :), eddy_vt(:, :)
    real(real64), allocatable :: pressure_half(:), delta_pressure(:), layer_l(:), alpha(:)
    real(real64), allocatable :: reference_temperature(:), a_half(:), b_half(:)
    character(len=:), allocatable :: case_directory
    integer :: completed_step, number_of_steps, daily_interval_steps
    integer :: month, month_boundary_step, year, year_boundary_step
    integer :: calendar_year, calendar_month, calendar_day, days_per_month, days_per_year
    integer(int64) :: start_count
    real(real64) :: cfl, maximum_cfl, elapsed_wall_seconds, seconds_of_day, solar_day, radiation_duration

    call ensure_context(context)
    numerics = context%numerics
    numerics%time_step = radiation_time_step
    physics = radiation_case_physics()
    planet = radiation_case_planet(physics)
    solar_day = physics%radiation%solar_day
    days_per_month = physics%radiation%days_per_month
    days_per_year = radiation_days_per_year(physics%radiation)
    radiation_duration = real(radiation_number_of_years, real64)*radiation_orbital_period(physics%radiation)
    call system_clock(start_count)
    call write_case_header('dry_radiation', radiation_duration, numerics%time_step)
    case_directory = context%output_root//'/dry_radiation'
    call make_directory(case_directory)
    call initialize_radiation_daily_output(case_directory)
    call solver%init_with_config(numerics)
    call set_radiation_case_state(solver, context%transform, physics, planet)
    call diagnostics%reset()
    call write_current_radiation_snapshot(solver, case_directory, 1, context%nlon)
    number_of_steps = nint(radiation_duration/numerics%time_step)
    daily_interval_steps = nint(solar_day/numerics%time_step)
    if (abs(real(daily_interval_steps, real64)*numerics%time_step - solar_day) > 1.0e-12_real64) then
      error stop 'radiation time step must divide the solar day exactly'
    end if
    maximum_cfl = 0.0_real64
    month = 1
    month_boundary_step = month*days_per_month*daily_interval_steps
    year = 2
    year_boundary_step = (year - 1)*days_per_year*daily_interval_steps

    do completed_step = 1, number_of_steps
      call solver%advance()
      call solver%take_latest_diagnostics(sample)
      call diagnostics%add(sample)
      cfl = solver%get_last_advance_cfl()
      maximum_cfl = max(maximum_cfl, cfl)

      ! Each row is the mean over the steps of one solar day, stamped with the day's start time.
      if (mod(completed_step, daily_interval_steps) == 0) then
        call diagnostics%take_daily(daily_means)
        call radiation_calendar_date(physics%radiation, daily_means%time_seconds, calendar_year, &
                                     calendar_month, calendar_day, seconds_of_day)
        call append_radiation_daily_output(case_directory, daily_means, daily_means%time_seconds/solar_day, &
                                           calendar_year, calendar_month, calendar_day)
      end if

      if (completed_step == month_boundary_step) then
        call diagnostics%take_monthly(monthly_surface_temperature, monthly_surface_pressure, &
          zonal_temperature, zonal_u, zonal_v, eddy_uv, eddy_vt)
        call write_radiation_monthly_output(case_directory, month, context%nlon, monthly_surface_temperature, &
          monthly_surface_pressure, zonal_temperature, zonal_u, zonal_v, eddy_uv, eddy_vt)
        month = month + 1
        if (month <= radiation_number_of_years*physics%radiation%months_per_year) then
          month_boundary_step = month*days_per_month*daily_interval_steps
        end if
      end if

      if (completed_step == year_boundary_step) then
        call write_current_radiation_snapshot(solver, case_directory, year, context%nlon)
        year = year + 1
        if (year <= radiation_number_of_years + 1) then
          year_boundary_step = (year - 1)*days_per_year*daily_interval_steps
        end if
      end if

      if (mod(completed_step, daily_interval_steps) == 0 .or. completed_step == number_of_steps) then
        call write_progress(completed_step, number_of_steps, completed_step, &
                            'advective CFL', cfl, start_count, numerics%time_step)
      end if
    end do

    elapsed_wall_seconds = elapsed_seconds(start_count)
    call solver%get_reference_atmosphere(pressure_half, delta_pressure, layer_l, alpha, &
                                         reference_temperature, a_half, b_half)
    call write_radiation_metadata(case_directory, physics%radiation, physics%convection, planet, &
                                  numerics%truncation, numerics%time_step, radiation_duration, number_of_steps, &
                                  maximum_cfl, elapsed_wall_seconds, context%nlon, context%transform%mu, &
                                  pressure_half, delta_pressure, layer_l, alpha, reference_temperature, &
                                  a_half, b_half)
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

end module radiation_case
