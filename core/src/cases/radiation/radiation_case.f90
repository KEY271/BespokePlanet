!> The dry radiation, slab-ocean and moist slab-ocean cases: five calendar
!> years with daily global means, monthly zonal means and yearly snapshots.
!> Each case chooses the surface physics, planet, resolution and calendar
!> boundaries for its output.
module radiation_case
  use iso_fortran_env, only: real64, int64
  use harmonics, only: harmonic_transform
  use planet_parameters, only: planet_config
  use numerics_config, only: model_numerics_config
  use dry_physics_config, only: dry_model_physics_config, radiation_days_per_year, radiation_orbital_period
  use dry_atmosphere, only: dry_atmosphere_solver
  use dry_case_initial_conditions, only: radiation_case_physics, slab_ocean_case_physics, moist_case_physics, &
                                         land_sea_case_physics, radiation_case_planet, set_radiation_case_state, &
                                         set_land_sea_case_state
  use topography, only: topography_config, topography_diagnostics, generate_topography
  use earth_topography, only: earth_topography_config, earth_topography_diagnostics, generate_earth_topography
  use field_binary_writer, only: write_field
  use ocean_q_flux, only: q_flux_diagnostics
  use dry_radiation, only: radiation_diagnostics, radiation_calendar_date
  use radiation_diagnostics_collector, only: radiation_case_diagnostics, radiation_monthly_means
  use filesystem, only: make_directory
  use radiation_case_output, only: radiation_output_options, initialize_radiation_daily_output, &
                                    append_radiation_daily_output, write_radiation_monthly_output, &
                                    write_radiation_yearly_snapshot, write_radiation_metadata
  use case_runtime, only: case_context, ensure_context, dt, write_case_header, write_progress, &
                          elapsed_seconds
  implicit none
  private
  public :: run_radiation_case, run_slab_ocean_case, run_moist_case, run_land_case, run_land_t63_case
  public :: run_land_earth_case, run_land_earth_t63_case

  real(real64), parameter :: radiation_time_step = dt
  integer, parameter :: radiation_number_of_years = 5
  !> The moist case is first run at a lower resolution than the other cases (docs/cases/moist.md).
  integer, parameter, public :: moist_case_truncation = 31

  integer, parameter :: ground_variant = 1, slab_ocean_variant = 2, moist_variant = 3, land_variant = 4
  integer, parameter :: land_earth_variant = 5

contains

  subroutine run_radiation_case(context)
    type(case_context), intent(inout) :: context
    call ensure_context(context)
    call run_radiative_surface_case(context, context%transform, context%nlon, ground_variant)
  end subroutine run_radiation_case

  subroutine run_slab_ocean_case(context)
    type(case_context), intent(inout) :: context
    call ensure_context(context)
    call run_radiative_surface_case(context, context%transform, context%nlon, slab_ocean_variant)
  end subroutine run_slab_ocean_case

  !> The moist case builds its initial state and output grid on its own, coarser transform.
  subroutine run_moist_case(context)
    type(case_context), intent(inout) :: context
    type(harmonic_transform) :: transform
    integer, allocatable :: nlon(:)

    call ensure_context(context)
    call transform%init(moist_case_truncation)
    nlon = transform%get_nlon()
    call run_radiative_surface_case(context, transform, nlon, moist_variant)
  end subroutine run_moist_case

  subroutine run_land_case(context)
    type(case_context), intent(inout) :: context
    call run_land_case_at_truncation(context, 31)
  end subroutine run_land_case

  subroutine run_land_t63_case(context)
    type(case_context), intent(inout) :: context
    call run_land_case_at_truncation(context, 63)
  end subroutine run_land_t63_case

  subroutine run_land_earth_case(context)
    type(case_context), intent(inout) :: context
    call run_land_case_at_truncation(context, 31, land_earth_variant)
  end subroutine run_land_earth_case

  subroutine run_land_earth_t63_case(context)
    type(case_context), intent(inout) :: context
    call run_land_case_at_truncation(context, 63, land_earth_variant)
  end subroutine run_land_earth_t63_case

  subroutine run_land_case_at_truncation(context, truncation, variant)
    type(case_context), intent(inout) :: context
    integer, intent(in) :: truncation
    integer, intent(in), optional :: variant
    type(harmonic_transform) :: transform
    integer, allocatable :: nlon(:)
    integer :: chosen_variant

    chosen_variant = land_variant
    if (present(variant)) chosen_variant = variant
    call ensure_context(context)
    call transform%init(truncation)
    nlon = transform%get_nlon()
    call run_radiative_surface_case(context, transform, nlon, chosen_variant, truncation)
  end subroutine run_land_case_at_truncation

  subroutine run_radiative_surface_case(context, transform, nlon, variant, requested_truncation)
    type(case_context), intent(inout) :: context
    type(harmonic_transform), intent(inout) :: transform
    integer, intent(in) :: nlon(:), variant
    integer, intent(in), optional :: requested_truncation
    type(dry_atmosphere_solver) :: solver
    type(model_numerics_config) :: numerics
    type(dry_model_physics_config) :: physics
    type(planet_config) :: planet
    type(radiation_case_diagnostics) :: diagnostics
    type(radiation_diagnostics) :: sample, daily_means
    type(radiation_monthly_means) :: monthly_means
    type(radiation_output_options) :: options
    type(topography_config) :: terrain
    type(topography_diagnostics) :: terrain_diagnostics
    type(earth_topography_config) :: earth_terrain
    type(earth_topography_diagnostics) :: earth_terrain_diagnostics
    real(real64), allocatable :: land_fraction(:, :), analytic_height(:, :), truncated_height(:, :)
    real(real64), allocatable :: q_flux(:, :)
    type(q_flux_diagnostics) :: q_flux_summary
    complex(real64), allocatable :: surface_geopotential(:, :)
    real(real64), allocatable :: pressure_half(:), delta_pressure(:), layer_l(:), alpha(:)
    real(real64), allocatable :: reference_temperature(:), a_half(:), b_half(:)
    character(len=:), allocatable :: case_directory, case_name
    integer :: completed_step, number_of_steps, daily_interval_steps
    integer :: month, month_boundary_step, year, year_boundary_step
    integer :: calendar_year, calendar_month, calendar_day, days_per_month, days_per_year
    integer(int64) :: start_count
    real(real64) :: cfl, maximum_cfl, elapsed_wall_seconds, seconds_of_day, solar_day, radiation_duration

    numerics = context%numerics
    numerics%time_step = radiation_time_step
    select case (variant)
    case (ground_variant)
      physics = radiation_case_physics()
      case_name = 'dry_radiation'
      options%include_deep_temperature = .true.
      options%include_moisture = .false.
    case (slab_ocean_variant)
      physics = slab_ocean_case_physics()
      case_name = 'dry_slab_ocean'
      options%include_deep_temperature = .false.
      options%include_moisture = .false.
    case (moist_variant)
      physics = moist_case_physics()
      case_name = 'moist_slab_ocean'
      numerics%truncation = moist_case_truncation
      options%include_deep_temperature = .false.
      options%include_moisture = .true.
    case (land_variant, land_earth_variant)
      if (.not. present(requested_truncation)) error stop 'land case truncation is required'
      physics = land_sea_case_physics()
      numerics%truncation = requested_truncation
      if (requested_truncation == 31) then
        case_name = 'moist_land_sea_t31'
      else if (requested_truncation == 63) then
        case_name = 'moist_land_sea_t63'
      else
        error stop 'land case supports only T31 and T63'
      end if
      if (variant == land_earth_variant) case_name = case_name(1:14)//'_earth'//case_name(15:)
      options%include_deep_temperature = .true.
      options%include_moisture = .true.
      options%include_land_sea = .true.
    case default
      error stop 'unknown radiative surface case variant'
    end select
    options%include_surface_tiles = physics%radiation%land_sea_mixing_enabled .or. physics%sea_ice%enabled
    options%include_sea_ice = physics%sea_ice%enabled
    planet = radiation_case_planet(physics)
    solar_day = physics%radiation%solar_day
    days_per_month = physics%radiation%days_per_month
    days_per_year = radiation_days_per_year(physics%radiation)
    radiation_duration = real(radiation_number_of_years, real64)*radiation_orbital_period(physics%radiation)
    call system_clock(start_count)
    call write_case_header(case_name, radiation_duration, numerics%time_step)
    case_directory = context%output_root//'/'//case_name
    call make_directory(case_directory)
    call initialize_radiation_daily_output(case_directory, options)
    call solver%init_with_config(numerics)
    if (variant == land_variant .or. variant == land_earth_variant) then
      if (variant == land_variant) then
        call generate_topography(transform, terrain, land_fraction, analytic_height, surface_geopotential, &
                                 truncated_height, terrain_diagnostics)
      else
        call generate_earth_topography(transform, earth_terrain, land_fraction, analytic_height, &
                                       surface_geopotential, truncated_height, earth_terrain_diagnostics)
      end if
      call set_land_sea_case_state(solver, transform, physics, planet, surface_geopotential, land_fraction)
      call write_field(trim(case_directory)//'/land_fraction.bin', nlon, land_fraction)
      call write_field(trim(case_directory)//'/surface_height.bin', nlon, truncated_height)
    else
      call set_radiation_case_state(solver, transform, physics, planet)
    end if
    call solver%get_ocean_q_flux(q_flux, q_flux_summary)
    if (physics%q_flux%enabled) call write_field(trim(case_directory)//'/ocean_q_flux.bin', nlon, q_flux)
    call diagnostics%reset()
    call write_current_radiation_snapshot(solver, case_directory, 1, nlon, options)
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
                                           calendar_year, calendar_month, calendar_day, options)
      end if

      if (completed_step == month_boundary_step) then
        call diagnostics%take_monthly(monthly_means)
        call write_radiation_monthly_output(case_directory, month, nlon, monthly_means, options)
        month = month + 1
        if (month <= radiation_number_of_years*physics%radiation%months_per_year) then
          month_boundary_step = month*days_per_month*daily_interval_steps
        end if
      end if

      if (completed_step == year_boundary_step) then
        call write_current_radiation_snapshot(solver, case_directory, year, nlon, options)
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
    if (variant == land_variant) then
      call write_radiation_metadata(case_directory, case_name, physics, planet, &
                                    numerics%truncation, numerics%time_step, radiation_duration, number_of_steps, &
                                    maximum_cfl, elapsed_wall_seconds, nlon, transform%mu, &
                                    pressure_half, delta_pressure, layer_l, alpha, reference_temperature, &
                                    a_half, b_half, options, terrain, terrain_diagnostics, q_flux=q_flux_summary)
    else if (variant == land_earth_variant) then
      call write_radiation_metadata(case_directory, case_name, physics, planet, &
                                    numerics%truncation, numerics%time_step, radiation_duration, number_of_steps, &
                                    maximum_cfl, elapsed_wall_seconds, nlon, transform%mu, &
                                    pressure_half, delta_pressure, layer_l, alpha, reference_temperature, &
                                    a_half, b_half, options, q_flux=q_flux_summary, earth_terrain=earth_terrain, &
                                    earth_terrain_diagnostics=earth_terrain_diagnostics)
    else
      call write_radiation_metadata(case_directory, case_name, physics, planet, &
                                    numerics%truncation, numerics%time_step, radiation_duration, number_of_steps, &
                                    maximum_cfl, elapsed_wall_seconds, nlon, transform%mu, &
                                    pressure_half, delta_pressure, layer_l, alpha, reference_temperature, &
                                    a_half, b_half, options, q_flux=q_flux_summary)
    end if
  end subroutine run_radiative_surface_case

  subroutine write_current_radiation_snapshot(solver, case_directory, year, ring_nlon, options)
    type(dry_atmosphere_solver), intent(inout) :: solver
    character(*), intent(in) :: case_directory
    integer, intent(in) :: year, ring_nlon(:)
    type(radiation_output_options), intent(in) :: options
    complex(real64), allocatable :: zeta_spectral(:, :, :), delta_spectral(:, :, :)
    complex(real64), allocatable :: temperature_spectral(:, :, :), log_ps_spectral(:, :)
    complex(real64), allocatable :: humidity_spectral(:, :, :)
    real(real64), allocatable :: zeta(:, :, :), delta(:, :, :), temperature(:, :, :)
    real(real64), allocatable :: surface_pressure(:, :), log_surface_pressure(:, :)
    real(real64), allocatable :: u(:, :, :), v(:, :, :), humidity(:, :, :)
    real(real64), allocatable :: surface_temperature(:, :), deep_temperature(:, :), cloud_cover(:, :)
    real(real64), allocatable :: surface_water(:, :)
    type(radiation_monthly_means) :: tiles

    if (options%include_moisture) then
      call solver%get_spectral_state(zeta_spectral, delta_spectral, temperature_spectral, log_ps_spectral, &
                                     specific_humidity=humidity_spectral)
      call solver%get_fields(zeta, delta, temperature, surface_pressure, u, v, &
                             surface_temperature=surface_temperature, deep_temperature=deep_temperature, &
                             specific_humidity=humidity, surface_water=surface_water, &
                             land_temperature=tiles%land_temperature, ocean_temperature=tiles%ocean_temperature, &
                             sea_ice_fraction=tiles%sea_ice_fraction, sea_ice_volume=tiles%sea_ice_volume, &
                             sea_ice_temperature=tiles%sea_ice_temperature, sea_ice_thickness=tiles%sea_ice_thickness)
      call solver%get_cloud_cover(cloud_cover)
      log_surface_pressure = log(surface_pressure)
      call write_radiation_yearly_snapshot(case_directory, year, ring_nlon, &
        zeta_spectral, delta_spectral, temperature_spectral, log_ps_spectral, &
        zeta, delta, temperature, u, v, log_surface_pressure, surface_temperature, deep_temperature, &
        options, humidity_spectral=humidity_spectral, humidity=humidity, time_seconds=solver%get_time(), &
        step=solver%get_step(), cloud_cover=cloud_cover, surface_water=surface_water, tiles=tiles)
    else
      call solver%get_spectral_state(zeta_spectral, delta_spectral, temperature_spectral, log_ps_spectral)
      call solver%get_fields(zeta, delta, temperature, surface_pressure, u, v, &
                             surface_temperature=surface_temperature, deep_temperature=deep_temperature, &
                             surface_water=surface_water, &
                             land_temperature=tiles%land_temperature, ocean_temperature=tiles%ocean_temperature, &
                             sea_ice_fraction=tiles%sea_ice_fraction, sea_ice_volume=tiles%sea_ice_volume, &
                             sea_ice_temperature=tiles%sea_ice_temperature, sea_ice_thickness=tiles%sea_ice_thickness)
      log_surface_pressure = log(surface_pressure)
      call write_radiation_yearly_snapshot(case_directory, year, ring_nlon, &
        zeta_spectral, delta_spectral, temperature_spectral, log_ps_spectral, &
        zeta, delta, temperature, u, v, log_surface_pressure, surface_temperature, deep_temperature, options, &
        surface_water=surface_water, tiles=tiles)
    end if
  end subroutine write_current_radiation_snapshot

end module radiation_case
