!> Output files and metadata schema shared by the dry radiation, slab-ocean and
!> moist slab-ocean cases.
module radiation_case_output
  use iso_fortran_env, only: real64
  use field_binary_writer, only: write_field, write_rectangular_field, write_spectral_field, &
                                 check_finite, check_spectral_finite, check_rectangular_finite
  use json_writer, only: write_real_array, write_integer_array, write_ring_offsets
  use csv_writer, only: write_csv_header, append_csv_row, csv_real, csv_integer
  use dry_case_output, only: write_dry_reference_atmosphere
  use planet_parameters, only: planet_config
  use dry_physics_config, only: dry_model_physics_config, radiation_config, radiation_days_per_year, &
                                radiation_orbital_period, radiation_surface_heat_capacity
  use dry_radiation, only: radiation_diagnostics
  use radiation_diagnostics_collector, only: radiation_monthly_means
  use moist_thermodynamics, only: water_vapor_gas_constant, latent_heat_of_condensation, &
                                  reference_saturation_vapor_pressure, saturation_reference_temperature, &
                                  dry_air_specific_heat
  use topography, only: topography_config, topography_diagnostics
  implicit none
  private

  public :: radiation_output_options
  public :: initialize_radiation_daily_output, append_radiation_daily_output
  public :: write_radiation_monthly_output, write_radiation_yearly_snapshot
  public :: write_radiation_metadata

  !> Which optional fields a case writes: the deep ground temperature exists in
  !> the ground and mixed land--sea cases, and water fields in moist cases.
  type :: radiation_output_options
    logical :: include_deep_temperature = .true.
    logical :: include_moisture = .false.
    logical :: include_land_sea = .false.
  end type radiation_output_options

  !> Precipitation and evaporation are aggregated in kg m^-2 s^-1 and written in mm day^-1.
  real(real64), parameter :: mm_per_day = 86400.0_real64

contains

  subroutine initialize_radiation_daily_output(case_directory, options)
    character(*), intent(in) :: case_directory
    type(radiation_output_options), intent(in) :: options
    character(len=:), allocatable :: header

    header = 'time_seconds,simulation_day,calendar_year,calendar_month,calendar_day,mean_atmospheric_temperature_k,'
    if (options%include_deep_temperature) then
      header = header//'mean_surface_temperature_k,mean_deep_temperature_k,'
    else
      header = header//'mean_ocean_temperature_k,'
    end if
    header = header//'mean_kinetic_energy_j_kg-1,mean_surface_pressure_pa,incoming_shortwave_w_m-2,'// &
             'reflected_shortwave_w_m-2,outgoing_longwave_w_m-2'
    if (options%include_moisture) then
      header = header//',precipitation_mm_day-1,convective_precipitation_mm_day-1,'// &
               'large_scale_precipitation_mm_day-1,evaporation_mm_day-1,latent_heat_flux_w_m-2,'// &
               'precipitable_water_kg_m-2,signed_column_water_kg_m-2,negative_column_water_kg_m-2,'// &
               'maximum_wind_speed_m_s-1,maximum_wind_longitude_deg,maximum_wind_latitude_deg,'// &
               'maximum_wind_level,maximum_wind_eta'
    end if
    if (options%include_land_sea) then
      header = header//',mean_land_surface_temperature_k,mean_ocean_surface_temperature_k,'// &
        'land_precipitation_mm_day-1,ocean_precipitation_mm_day-1,'// &
        'land_evaporation_mm_day-1,ocean_evaporation_mm_day-1'
    end if
    call write_csv_header(trim(case_directory)//'/daily_global.csv', header)
  end subroutine initialize_radiation_daily_output

  !> Appends one row of daily means; the calendar date of the day start is supplied by the case runner.
  subroutine append_radiation_daily_output(case_directory, means, simulation_day, &
                                           calendar_year, calendar_month, calendar_day, options)
    character(*), intent(in) :: case_directory
    type(radiation_diagnostics), intent(in) :: means
    real(real64), intent(in) :: simulation_day
    integer, intent(in) :: calendar_year, calendar_month, calendar_day
    type(radiation_output_options), intent(in) :: options
    character(len=:), allocatable :: row

    row = csv_real(means%time_seconds)//','//csv_real(simulation_day)//','// &
          csv_integer(calendar_year)//','//csv_integer(calendar_month)//','//csv_integer(calendar_day)//','// &
          csv_real(means%mean_atmospheric_temperature)//','//csv_real(means%mean_surface_temperature)//','
    if (options%include_deep_temperature) row = row//csv_real(means%mean_deep_temperature)//','
    row = row//csv_real(means%mean_kinetic_energy)//','//csv_real(means%mean_surface_pressure)//','// &
          csv_real(means%mean_incoming_shortwave)//','//csv_real(means%mean_reflected_shortwave)//','// &
          csv_real(means%mean_outgoing_longwave)
    if (options%include_moisture) then
      row = row//','// &
            csv_real(mm_per_day*(means%mean_convective_precipitation + means%mean_large_scale_precipitation))//','// &
            csv_real(mm_per_day*means%mean_convective_precipitation)//','// &
            csv_real(mm_per_day*means%mean_large_scale_precipitation)//','// &
            csv_real(mm_per_day*means%mean_evaporation)//','//csv_real(means%mean_latent_heat_flux)//','// &
            csv_real(means%mean_precipitable_water)//','//csv_real(means%mean_signed_column_water)//','// &
            csv_real(means%mean_negative_column_water)//','//csv_real(means%maximum_wind_speed)//','// &
            csv_real(means%maximum_wind_longitude_degrees)//','//csv_real(means%maximum_wind_latitude_degrees)//','// &
            csv_integer(means%maximum_wind_level)//','//csv_real(means%maximum_wind_eta)
    end if
    if (options%include_land_sea) then
      row = row//','//csv_real(means%mean_land_surface_temperature)//','// &
        csv_real(means%mean_ocean_surface_temperature)//','// &
        csv_real(mm_per_day*means%mean_land_precipitation)//','// &
        csv_real(mm_per_day*means%mean_ocean_precipitation)//','// &
        csv_real(mm_per_day*means%mean_land_evaporation)//','// &
        csv_real(mm_per_day*means%mean_ocean_evaporation)
    end if
    call append_csv_row(trim(case_directory)//'/daily_global.csv', row)
  end subroutine append_radiation_daily_output

  subroutine write_radiation_monthly_output(case_directory, month, nlon, means, options)
    character(*), intent(in) :: case_directory
    integer, intent(in) :: month
    integer, intent(in) :: nlon(:)
    type(radiation_monthly_means), intent(in) :: means
    type(radiation_output_options), intent(in) :: options
    character(len=4) :: month_text

    if (month < 1 .or. month > 9999) error stop 'radiation output month is outside the supported range'
    if (size(means%zonal_temperature, 1) /= size(nlon) .or. &
        any(shape(means%zonal_temperature) /= shape(means%zonal_u)) .or. &
        any(shape(means%zonal_temperature) /= shape(means%zonal_v)) .or. &
        any(shape(means%zonal_temperature) /= shape(means%eddy_uv)) .or. &
        any(shape(means%zonal_temperature) /= shape(means%eddy_vt))) then
      error stop 'radiation monthly zonal fields have inconsistent shapes'
    end if
    call check_finite('monthly_surface_temperature', nlon, means%surface_temperature)
    call check_finite('monthly_surface_pressure', nlon, means%surface_pressure)
    call check_rectangular_finite('monthly_zonal_temperature', means%zonal_temperature)
    call check_rectangular_finite('monthly_zonal_u', means%zonal_u)
    call check_rectangular_finite('monthly_zonal_v', means%zonal_v)
    call check_rectangular_finite('monthly_eddy_uv', means%eddy_uv)
    call check_rectangular_finite('monthly_eddy_vt', means%eddy_vt)
    write (month_text, '(i4.4)') month
    call write_field(trim(case_directory)//'/monthly_surface_temperature_m'//month_text//'.bin', &
                     nlon, means%surface_temperature)
    call write_field(trim(case_directory)//'/monthly_surface_pressure_m'//month_text//'.bin', &
                     nlon, means%surface_pressure)
    if (options%include_deep_temperature) then
      call check_finite('monthly_deep_temperature', nlon, means%deep_temperature)
      call write_field(trim(case_directory)//'/monthly_deep_temperature_m'//month_text//'.bin', &
                       nlon, means%deep_temperature)
    end if
    call write_rectangular_field(trim(case_directory)//'/monthly_zonal_temperature_m'//month_text//'.bin', &
                                 means%zonal_temperature)
    call write_rectangular_field(trim(case_directory)//'/monthly_zonal_u_m'//month_text//'.bin', means%zonal_u)
    call write_rectangular_field(trim(case_directory)//'/monthly_zonal_v_m'//month_text//'.bin', means%zonal_v)
    call write_rectangular_field(trim(case_directory)//'/monthly_eddy_uv_m'//month_text//'.bin', means%eddy_uv)
    call write_rectangular_field(trim(case_directory)//'/monthly_eddy_vt_m'//month_text//'.bin', means%eddy_vt)
    if (options%include_moisture) then
      if (.not. allocated(means%precipitation) .or. .not. allocated(means%evaporation) .or. &
          .not. allocated(means%precipitable_water) .or. .not. allocated(means%zonal_humidity) .or. &
          .not. allocated(means%eddy_vq)) error stop 'moist monthly means are missing'
      call check_finite('monthly_precipitation', nlon, means%precipitation)
      call check_finite('monthly_evaporation', nlon, means%evaporation)
      call check_finite('monthly_precipitable_water', nlon, means%precipitable_water)
      call check_rectangular_finite('monthly_zonal_humidity', means%zonal_humidity)
      call check_rectangular_finite('monthly_eddy_vq', means%eddy_vq)
      call write_field(trim(case_directory)//'/monthly_precipitation_m'//month_text//'.bin', &
                       nlon, mm_per_day*means%precipitation)
      call write_field(trim(case_directory)//'/monthly_evaporation_m'//month_text//'.bin', &
                       nlon, mm_per_day*means%evaporation)
      call write_field(trim(case_directory)//'/monthly_precipitable_water_m'//month_text//'.bin', &
                       nlon, means%precipitable_water)
      call write_rectangular_field(trim(case_directory)//'/monthly_zonal_humidity_m'//month_text//'.bin', &
                                   means%zonal_humidity)
      call write_rectangular_field(trim(case_directory)//'/monthly_eddy_vq_m'//month_text//'.bin', means%eddy_vq)
    end if
  end subroutine write_radiation_monthly_output

  subroutine write_radiation_yearly_snapshot(case_directory, year, nlon, &
                                              zeta_spectral, delta_spectral, temperature_spectral, &
                                              log_surface_pressure_spectral, zeta, delta, temperature, &
                                              u, v, log_surface_pressure, surface_temperature, deep_temperature, &
                                              options, humidity_spectral, humidity, surface_temperature_spectral, &
                                              time_seconds, step, deep_temperature_spectral)
    character(*), intent(in) :: case_directory
    integer, intent(in) :: year
    integer, intent(in) :: nlon(:)
    complex(real64), intent(in) :: zeta_spectral(0:, 0:, :), delta_spectral(0:, 0:, :)
    complex(real64), intent(in) :: temperature_spectral(0:, 0:, :)
    complex(real64), intent(in) :: log_surface_pressure_spectral(0:, 0:)
    real(real64), intent(in) :: zeta(:, :, :), delta(:, :, :), temperature(:, :, :)
    real(real64), intent(in) :: u(:, :, :), v(:, :, :), log_surface_pressure(:, :)
    real(real64), intent(in) :: surface_temperature(:, :), deep_temperature(:, :)
    type(radiation_output_options), intent(in) :: options
    !> Moist-case additions: spectral and grid specific humidity, spectral surface
    !> (ocean) temperature and the model time of the snapshot.
    complex(real64), intent(in), optional :: humidity_spectral(0:, 0:, :), surface_temperature_spectral(0:, 0:)
    complex(real64), intent(in), optional :: deep_temperature_spectral(0:, 0:)
    real(real64), intent(in), optional :: humidity(:, :, :), time_seconds
    integer, intent(in), optional :: step
    character(len=4) :: year_text
    character(len=2) :: level_text
    integer :: k, unit

    if (year < 1 .or. year > 9999) error stop 'radiation snapshot year is outside the supported range'
    if (size(zeta, 3) /= size(delta, 3) .or. size(zeta, 3) /= size(temperature, 3) .or. &
        size(zeta, 3) /= size(u, 3) .or. size(zeta, 3) /= size(v, 3) .or. &
        size(zeta, 3) /= size(zeta_spectral, 3) .or. &
        size(zeta, 3) /= size(delta_spectral, 3) .or. &
        size(zeta, 3) /= size(temperature_spectral, 3)) then
      error stop 'radiation yearly snapshot fields have different level counts'
    end if
    if (options%include_moisture) then
      if (.not. (present(humidity_spectral) .and. present(humidity) .and. &
                 present(surface_temperature_spectral) .and. present(time_seconds) .and. present(step))) then
        error stop 'moist yearly snapshot requires the humidity, ocean temperature and time'
      end if
      if (size(humidity, 3) /= size(zeta, 3) .or. size(humidity_spectral, 3) /= size(zeta, 3)) then
        error stop 'moist yearly snapshot humidity has a different level count'
      end if
    end if
    write (year_text, '(i4.4)') year
    call check_finite('yearly_log_surface_pressure', nlon, log_surface_pressure)
    call check_finite('yearly_surface_temperature', nlon, surface_temperature)
    if (options%include_deep_temperature) call check_finite('yearly_deep_temperature', nlon, deep_temperature)
    call check_spectral_finite('yearly_log_surface_pressure_spectral', log_surface_pressure_spectral)
    call write_field(trim(case_directory)//'/yearly_log_surface_pressure_y'//year_text//'.bin', &
                     nlon, log_surface_pressure)
    call write_field(trim(case_directory)//'/yearly_surface_temperature_y'//year_text//'.bin', &
                     nlon, surface_temperature)
    if (options%include_deep_temperature) then
      call write_field(trim(case_directory)//'/yearly_deep_temperature_y'//year_text//'.bin', &
                       nlon, deep_temperature)
      if (present(deep_temperature_spectral)) then
        call check_spectral_finite('yearly_deep_temperature_spectral', deep_temperature_spectral)
        call write_spectral_field(trim(case_directory)//'/yearly_deep_temperature_spectral_y'// &
                                  year_text//'.bin', deep_temperature_spectral)
      end if
    end if
    call write_spectral_field(trim(case_directory)//'/yearly_log_surface_pressure_spectral_y'// &
                              year_text//'.bin', log_surface_pressure_spectral)
    if (options%include_moisture) then
      call check_spectral_finite('yearly_surface_temperature_spectral', surface_temperature_spectral)
      call write_spectral_field(trim(case_directory)//'/yearly_surface_temperature_spectral_y'// &
                                year_text//'.bin', surface_temperature_spectral)
      open (newunit=unit, file=trim(case_directory)//'/yearly_time_y'//year_text//'.json', status='replace', &
            action='write')
      write (unit, '(a)') '{'
      write (unit, '(a,es24.16e3,a)') '  "time_seconds": ', time_seconds, ','
      write (unit, '(a,i0)') '  "step": ', step
      write (unit, '(a)') '}'
      close (unit)
    end if
    do k = 1, size(zeta, 3)
      write (level_text, '(i2.2)') k
      call check_finite('yearly_zeta', nlon, zeta(:, :, k))
      call check_finite('yearly_delta', nlon, delta(:, :, k))
      call check_finite('yearly_temperature', nlon, temperature(:, :, k))
      call check_finite('yearly_u', nlon, u(:, :, k))
      call check_finite('yearly_v', nlon, v(:, :, k))
      call check_spectral_finite('yearly_zeta_spectral', zeta_spectral(:, :, k))
      call check_spectral_finite('yearly_delta_spectral', delta_spectral(:, :, k))
      call check_spectral_finite('yearly_temperature_spectral', temperature_spectral(:, :, k))
      call write_field(trim(case_directory)//'/yearly_zeta_y'//year_text//'_l'//level_text//'.bin', &
                       nlon, zeta(:, :, k))
      call write_field(trim(case_directory)//'/yearly_delta_y'//year_text//'_l'//level_text//'.bin', &
                       nlon, delta(:, :, k))
      call write_field(trim(case_directory)//'/yearly_temperature_y'//year_text//'_l'//level_text//'.bin', &
                       nlon, temperature(:, :, k))
      call write_field(trim(case_directory)//'/yearly_u_y'//year_text//'_l'//level_text//'.bin', nlon, u(:, :, k))
      call write_field(trim(case_directory)//'/yearly_v_y'//year_text//'_l'//level_text//'.bin', nlon, v(:, :, k))
      call write_spectral_field(trim(case_directory)//'/yearly_zeta_spectral_y'//year_text//'_l'// &
                                level_text//'.bin', zeta_spectral(:, :, k))
      call write_spectral_field(trim(case_directory)//'/yearly_delta_spectral_y'//year_text//'_l'// &
                                level_text//'.bin', delta_spectral(:, :, k))
      call write_spectral_field(trim(case_directory)//'/yearly_temperature_spectral_y'//year_text//'_l'// &
                                level_text//'.bin', temperature_spectral(:, :, k))
      if (options%include_moisture) then
        call check_finite('yearly_specific_humidity', nlon, humidity(:, :, k))
        call check_spectral_finite('yearly_specific_humidity_spectral', humidity_spectral(:, :, k))
        call write_field(trim(case_directory)//'/yearly_specific_humidity_y'//year_text//'_l'//level_text//'.bin', &
                         nlon, humidity(:, :, k))
        call write_spectral_field(trim(case_directory)//'/yearly_specific_humidity_spectral_y'//year_text//'_l'// &
                                  level_text//'.bin', humidity_spectral(:, :, k))
      end if
    end do
  end subroutine write_radiation_yearly_snapshot

  !> Metadata of one run.  The physics and planet configuration are the values
  !> the case runner integrated with; this module imports no case constants.
  subroutine write_radiation_metadata(case_directory, case_name, physics, planet, truncation, time_step, &
                                      duration, number_of_steps, maximum_cfl, elapsed_wall_seconds, nlon, mu, &
                                      pressure_half, delta_pressure, layer_l, alpha, reference_temperature, &
                                      a_half, b_half, options, terrain, terrain_diagnostics)
    character(*), intent(in) :: case_directory, case_name
    type(dry_model_physics_config), intent(in) :: physics
    type(planet_config), intent(in) :: planet
    integer, intent(in) :: truncation, number_of_steps
    real(real64), intent(in) :: time_step, duration, maximum_cfl, elapsed_wall_seconds
    integer, intent(in) :: nlon(:)
    real(real64), intent(in) :: mu(:), pressure_half(0:), delta_pressure(:), layer_l(:), alpha(:)
    real(real64), intent(in) :: reference_temperature(:), a_half(0:), b_half(0:)
    type(radiation_output_options), intent(in) :: options
    type(topography_config), intent(in), optional :: terrain
    type(topography_diagnostics), intent(in), optional :: terrain_diagnostics
    type(radiation_config) :: radiation
    integer :: unit
    real(real64) :: orbital_period

    radiation = physics%radiation
    orbital_period = radiation_orbital_period(radiation)
    open (newunit=unit, file=trim(case_directory)//'/metadata.json', status='replace', action='write')
    write (unit, '(a)') '{'
    write (unit, '(a)') '  "schema_version": 3,'
    write (unit, '(a)') '  "case_name": "'//trim(case_name)//'",'
    if (physics%moisture%enabled) then
      write (unit, '(a)') '  "equation": "moist_hydrostatic_atmosphere",'
    else
      write (unit, '(a)') '  "equation": "dry_hydrostatic_atmosphere",'
    end if
    if (radiation%land_sea_mixing_enabled) then
      write (unit, '(a)') '  "initial_condition": '// &
        '"pressure-coordinate Jablonowski-Williamson basic state over analytic truncated terrain, mixed land and ocean",'
    else if (physics%moisture%enabled) then
      write (unit, '(a)') '  "initial_condition": '// &
        '"unperturbed Jablonowski-Williamson temperature, balanced wind, flat terrain, slab ocean, '// &
        'q = initial_relative_humidity * q_s(T, p) at and below initial_humidity_top_pressure, dry above",'
    else if (radiation%slab_ocean_enabled) then
      write (unit, '(a)') '  "initial_condition": '// &
        '"unperturbed Jablonowski-Williamson temperature, balanced wind, flat terrain, slab ocean",'
    else
      write (unit, '(a)') '  "initial_condition": '// &
        '"unperturbed Jablonowski-Williamson temperature, balanced wind, flat terrain, two-layer ground",'
    end if
    write (unit, '(a)') '  "simulation": {'
    write (unit, '(a,es24.16e3,a)') '    "duration_seconds": ', duration, ','
    write (unit, '(a,es24.16e3,a)') '    "time_step_seconds": ', time_step, ','
    write (unit, '(a,i0,a)') '    "number_of_steps": ', number_of_steps, ','
    write (unit, '(a,es24.16e3,a)') &
      '    "actual_end_time_seconds": ', real(number_of_steps, real64)*time_step, ','
    write (unit, '(a)') '    "start_season": "northern_spring_equinox",'
    write (unit, '(a)') '    "start_calendar_time": "0001-04-01 00:00:00",'
    write (unit, '(a)') '    "longitude_zero_local_solar_time_at_start": "noon"'
    write (unit, '(a)') '  },'
    write (unit, '(a)') '  "calendar": {'
    write (unit, '(a)') '    "type": "360_day",'
    write (unit, '(a,es24.16e3,a)') '    "solar_day_seconds": ', radiation%solar_day, ','
    write (unit, '(a,i0,a)') '    "days_per_month": ', radiation%days_per_month, ','
    write (unit, '(a,i0,a)') '    "months_per_year": ', radiation%months_per_year, ','
    write (unit, '(a,i0)') '    "days_per_year": ', radiation_days_per_year(radiation)
    write (unit, '(a)') '  },'
    write (unit, '(a)') '  "radiation": {'
    write (unit, '(a,es24.16e3,a)') '    "solar_constant_w_m-2": ', radiation%solar_constant, ','
    write (unit, '(a,es24.16e3,a)') '    "surface_shortwave_albedo": ', radiation%surface_shortwave_albedo, ','
    write (unit, '(a,es24.16e3,a)') '    "axial_tilt_radians": ', radiation%axial_tilt, ','
    write (unit, '(a,es24.16e3,a)') '    "orbital_period_seconds": ', orbital_period, ','
    write (unit, '(a,es24.16e3,a)') '    "rotation_rate_rad_s": ', planet%rotation_rate, ','
    write (unit, '(a)') '    "longwave_optical_depth": "grey; d tau/dp = (a mu + b q)/p_0 plus ozone '// &
      '(Byrne & O''Gorman 2013 as in Isca)",'
    write (unit, '(a,es24.16e3,a)') &
      '    "longwave_well_mixed_optical_depth_a": ', radiation%longwave_well_mixed_optical_depth, ','
    write (unit, '(a,es24.16e3,a)') &
      '    "longwave_water_vapor_optical_depth_b": ', radiation%longwave_water_vapor_optical_depth, ','
    write (unit, '(a,es24.16e3,a)') &
      '    "longwave_well_mixed_scaling_mu": ', radiation%longwave_well_mixed_scaling, ','
    write (unit, '(a,es24.16e3,a)') &
      '    "longwave_reference_pressure_pa": ', radiation%longwave_reference_pressure, ','
    if (.not. physics%moisture%enabled) then
      write (unit, '(a,es24.16e3,a)') &
        '    "longwave_reference_surface_humidity_q0": ', radiation%longwave_reference_surface_humidity, ','
    end if
    write (unit, '(a,es24.16e3,a)') &
      '    "stefan_boltzmann_constant_w_m-2_k-4": ', radiation%stefan_boltzmann_constant, ','
    write (unit, '(a,es24.16e3,a)') '    "dry_air_specific_heat_j_kg-1_k-1": ', radiation%dry_air_specific_heat, ','
    write (unit, '(a,es24.16e3,a)') '    "gravity_acceleration_m_s-2": ', radiation%gravity_acceleration, ','
    write (unit, '(a)') '    "shortwave_atmosphere": "prescribed ozone absorption; no scattering",'
    write (unit, '(a,es24.16e3,a)') &
      '    "ultraviolet_shortwave_fraction": ', radiation%ultraviolet_shortwave_fraction, ','
    write (unit, '(a,es24.16e3,a)') &
      '    "ozone_shortwave_optical_depth": ', radiation%ozone_shortwave_optical_depth, ','
    write (unit, '(a,es24.16e3,a)') &
      '    "ozone_longwave_optical_depth": ', radiation%ozone_longwave_optical_depth, ','
    write (unit, '(a,es24.16e3,a)') &
      '    "ozone_pressure_lower_bound_pa": ', radiation%ozone_pressure_lower_bound, ','
    write (unit, '(a,es24.16e3,a)') &
      '    "ozone_pressure_upper_bound_pa": ', radiation%ozone_pressure_upper_bound, ','
    write (unit, '(a,es24.16e3,a)') '    "ozone_peak_pressure_pa": ', radiation%ozone_peak_pressure, ','
    write (unit, '(a,es24.16e3,a)') '    "ozone_log_pressure_width": ', radiation%ozone_log_pressure_width, ','
    write (unit, '(a)') '    "ozone_optical_path": "vertical; no solar-zenith slant correction",'
    write (unit, '(a)') '    "shortwave_forcing": "instantaneous local solar zenith angle",'
    if (physics%moisture%enabled) then
      write (unit, '(a)') '    "water_vapor_feedback": "grey longwave optical depth follows the prognostic '// &
        'specific humidity (previous time level, clipped at zero)",'
    else
      write (unit, '(a)') '    "water_vapor_feedback": "none; fixed reference humidity q_ref = q_0 (p/p_s)^3",'
    end if
    write (unit, '(a)') '    "orbit_eccentricity": 0'
    write (unit, '(a)') '  },'
    write (unit, '(a)') '  "surface_exchange": {'
    write (unit, '(a)') '    "sensible_heat": '// &
      '"bulk; lowest level extrapolated dry-adiabatically to the surface pressure; rho_N = p_N/(R T_N)",'
    write (unit, '(a,es24.16e3,a)') &
      '    "exchange_coefficient": ', radiation%surface_exchange_coefficient, ','
    write (unit, '(a,es24.16e3)') '    "gustiness_speed_m_s-1": ', radiation%gustiness_speed
    write (unit, '(a)') '  },'
    write (unit, '(a)') '  "dry_convective_adjustment": {'
    if (physics%moisture%enabled) then
      write (unit, '(a)') '    "method": "enthalpy-conserving pool adjacent violators on virtual potential temperature; '// &
        'humidity mixed within each block",'
    else
      write (unit, '(a)') '    "method": "enthalpy-conserving pool adjacent violators",'
    end if
    write (unit, '(a,es24.16e3)') &
      '    "relaxation_time_seconds": ', physics%convection%adjustment_time
    write (unit, '(a)') '  },'
    if (physics%moisture%enabled) then
      write (unit, '(a)') '  "moisture": {'
      write (unit, '(a)') '    "prognostic_variable": "specific_humidity",'
      write (unit, '(a)') '    "dynamics": "virtual temperature in hydrostatic, pressure-gradient and adiabatic terms; '// &
        'surface pressure unaffected by water vapour",'
      write (unit, '(a)') '    "negative_humidity": "allowed on the grid; clipped to zero for physics and T_v only",'
      write (unit, '(a,es24.16e3,a)') '    "water_vapor_gas_constant_j_kg-1_k-1": ', water_vapor_gas_constant, ','
      write (unit, '(a,es24.16e3,a)') '    "latent_heat_j_kg-1": ', latent_heat_of_condensation, ','
      write (unit, '(a,es24.16e3,a)') '    "dry_air_specific_heat_j_kg-1_k-1": ', dry_air_specific_heat, ','
      write (unit, '(a,es24.16e3,a)') &
        '    "reference_saturation_vapor_pressure_pa": ', reference_saturation_vapor_pressure, ','
      write (unit, '(a,es24.16e3,a)') '    "reference_temperature_k": ', saturation_reference_temperature, ','
      write (unit, '(a,es24.16e3,a)') '    "evaporation_surface_wetness": ', physics%evaporation%surface_wetness, ','
      write (unit, '(a,es24.16e3,a)') &
        '    "initial_relative_humidity": ', physics%moisture%initial_relative_humidity, ','
      write (unit, '(a,es24.16e3,a)') &
        '    "initial_humidity_top_pressure_pa": ', physics%moisture%initial_humidity_top_pressure, ','
      write (unit, '(a)') '    "moist_convective_adjustment": {'
      write (unit, '(a)') '      "method": "simplified Betts-Miller (Frierson 2007), enthalpy-conserving shift, '// &
        'shallower non-precipitating scheme",'
      write (unit, '(a,es24.16e3,a)') &
        '      "relaxation_time_seconds": ', physics%moist_convection%adjustment_time, ','
      write (unit, '(a,es24.16e3)') &
        '      "reference_relative_humidity": ', physics%moist_convection%reference_relative_humidity
      write (unit, '(a)') '    },'
      write (unit, '(a)') '    "large_scale_condensation": {'
      write (unit, '(a)') '      "method": "Newton solve to saturation within one leapfrog step; condensate falls out",'
      write (unit, '(a,es24.16e3)') &
        '      "saturation_tolerance": ', physics%condensation%saturation_tolerance
      write (unit, '(a)') '    },'
      write (unit, '(a)') '    "physics_order": "dry adjustment, moist adjustment, condensation on provisional fields"'
      write (unit, '(a)') '  },'
    end if
    if (radiation%land_sea_mixing_enabled) then
      if (.not. present(terrain) .or. .not. present(terrain_diagnostics)) then
        error stop 'land-sea metadata requires terrain configuration and diagnostics'
      end if
      write (unit, '(a)') '  "land_sea_surface": {'
      write (unit, '(a)') '    "mixing": "linear in land_fraction",'
      write (unit, '(a,es24.16e3,a)') &
        '    "land_surface_heat_capacity_j_m-2_k-1": ', radiation%surface_heat_capacity, ','
      write (unit, '(a,es24.16e3,a)') '    "ocean_depth_m": ', radiation%slab_ocean_depth, ','
      write (unit, '(a,es24.16e3,a)') '    "deep_heat_capacity_j_m-2_k-1": ', &
        radiation%deep_ground_heat_capacity, ','
      write (unit, '(a,es24.16e3,a)') '    "land_ground_exchange_w_m-2_k-1": ', &
        radiation%ground_exchange_coefficient, ','
      write (unit, '(a,es24.16e3,a)') '    "land_albedo": ', radiation%land_shortwave_albedo, ','
      write (unit, '(a,es24.16e3,a)') '    "ocean_albedo": ', radiation%ocean_shortwave_albedo, ','
      write (unit, '(a,es24.16e3,a)') '    "land_wetness": ', physics%evaporation%land_surface_wetness, ','
      write (unit, '(a,es24.16e3)') '    "ocean_wetness": ', physics%evaporation%ocean_surface_wetness
      write (unit, '(a)') '  },'
      write (unit, '(a)') '  "topography": {'
      write (unit, '(a)') '    "shapes": "three ellipse continents, one south polar cap, two mountain ranges",'
      write (unit, '(a,es24.16e3,a)') '    "coast_width_degrees": ', terrain%coast_width_degrees, ','
      write (unit, '(a,es24.16e3,a)') '    "base_height_m": ', terrain%base_height_metres, ','
      write (unit, '(a,es24.16e3,4(",",es24.16e3),a)') '    "continent_a_lon_lat_a_b_theta_deg": [', &
        terrain%continent_a%longitude_degrees, terrain%continent_a%latitude_degrees, &
        terrain%continent_a%semi_axis_east_degrees, terrain%continent_a%semi_axis_north_degrees, &
        terrain%continent_a%orientation_degrees, '],'
      write (unit, '(a,es24.16e3,4(",",es24.16e3),a)') '    "continent_b_lon_lat_a_b_theta_deg": [', &
        terrain%continent_b%longitude_degrees, terrain%continent_b%latitude_degrees, &
        terrain%continent_b%semi_axis_east_degrees, terrain%continent_b%semi_axis_north_degrees, &
        terrain%continent_b%orientation_degrees, '],'
      write (unit, '(a,es24.16e3,4(",",es24.16e3),a)') '    "continent_c_lon_lat_a_b_theta_deg": [', &
        terrain%continent_c%longitude_degrees, terrain%continent_c%latitude_degrees, &
        terrain%continent_c%semi_axis_east_degrees, terrain%continent_c%semi_axis_north_degrees, &
        terrain%continent_c%orientation_degrees, '],'
      write (unit, '(a,es24.16e3,a,i0,a)') '    "south_polar_cap_edge_deg_sign": [', &
        terrain%south_polar_cap%edge_latitude_degrees, ',', terrain%south_polar_cap%hemisphere_sign, '],'
      write (unit, '(a,es24.16e3,5(",",es24.16e3),a)') '    "mountain_a_lon_lat_length_width_theta_height": [', &
        terrain%mountain_a%longitude_degrees, terrain%mountain_a%latitude_degrees, &
        terrain%mountain_a%length_degrees, terrain%mountain_a%half_width_degrees, &
        terrain%mountain_a%orientation_degrees, terrain%mountain_a%height_metres, '],'
      write (unit, '(a,es24.16e3,5(",",es24.16e3),a)') '    "mountain_b_lon_lat_length_width_theta_height": [', &
        terrain%mountain_b%longitude_degrees, terrain%mountain_b%latitude_degrees, &
        terrain%mountain_b%length_degrees, terrain%mountain_b%half_width_degrees, &
        terrain%mountain_b%orientation_degrees, terrain%mountain_b%height_metres, '],'
      write (unit, '(a,es24.16e3,a)') '    "global_land_fraction": ', &
        terrain_diagnostics%global_land_fraction, ','
      write (unit, '(a,es24.16e3,a)') '    "minimum_truncated_height_m": ', &
        terrain_diagnostics%minimum_truncated_height_metres, ','
      write (unit, '(a,es24.16e3,a)') '    "maximum_truncated_height_m": ', &
        terrain_diagnostics%maximum_truncated_height_metres, ','
      write (unit, '(a,es24.16e3)') '    "height_rms_error_m": ', &
        terrain_diagnostics%height_rms_error_metres
      write (unit, '(a)') '  },'
    else if (radiation%slab_ocean_enabled) then
      write (unit, '(a)') '  "slab_ocean": {'
      write (unit, '(a,i0,a)') '    "number_of_layers": ', 1, ','
      write (unit, '(a,es24.16e3,a)') '    "depth_m": ', radiation%slab_ocean_depth, ','
      write (unit, '(a,es24.16e3,a)') '    "density_kg_m-3": ', radiation%seawater_density, ','
      write (unit, '(a,es24.16e3,a)') &
        '    "specific_heat_j_kg-1_k-1": ', radiation%seawater_specific_heat, ','
      write (unit, '(a,es24.16e3,a)') &
        '    "heat_capacity_j_m-2_k-1": ', radiation_surface_heat_capacity(radiation), ','
      write (unit, '(a)') '    "ocean_heat_transport": "none",'
      write (unit, '(a)') '    "bottom_heat_flux": "none",'
      if (physics%evaporation%enabled) then
        write (unit, '(a)') '    "latent_heat_flux": "subtracted from the ocean budget",'
      else
        write (unit, '(a)') '    "latent_heat_flux": "none",'
      end if
      write (unit, '(a,es24.16e3,a)') &
        '    "surface_exchange_coefficient": ', radiation%surface_exchange_coefficient, ','
      write (unit, '(a,es24.16e3)') '    "gustiness_speed_m_s-1": ', radiation%gustiness_speed
      write (unit, '(a)') '  },'
    else
      write (unit, '(a)') '  "ground": {'
      write (unit, '(a,es24.16e3,a)') '    "surface_heat_capacity_j_m-2_k-1": ', radiation%surface_heat_capacity, ','
      write (unit, '(a,es24.16e3,a)') '    "deep_heat_capacity_j_m-2_k-1": ', radiation%deep_ground_heat_capacity, ','
      write (unit, '(a,es24.16e3,a)') &
        '    "exchange_coefficient_w_m-2_k-1": ', radiation%ground_exchange_coefficient, ','
      write (unit, '(a,es24.16e3,a)') &
        '    "surface_exchange_coefficient": ', radiation%surface_exchange_coefficient, ','
      write (unit, '(a,es24.16e3)') '    "gustiness_speed_m_s-1": ', radiation%gustiness_speed
      write (unit, '(a)') '  },'
    end if
    write (unit, '(a)') '  "numerics": {'
    write (unit, '(a,i0,a)') '    "spectral_truncation": ', truncation, ','
    write (unit, '(a)') '    "time_integrator": "semi-implicit RAW-filtered leapfrog",'
    write (unit, '(a)') '    "physics_tendency_state": "RAW-filtered previous time level",'
    if (physics%rayleigh_friction%enabled) then
      write (unit, '(a,i0,a)') '    "upper_rayleigh_friction_levels": ', physics%rayleigh_friction%top_levels, ','
    else
      write (unit, '(a)') '    "upper_rayleigh_friction_levels": 0,'
    end if
    write (unit, '(a,es24.16e3,a)') '    "maximum_advective_cfl": ', maximum_cfl, ','
    write (unit, '(a,es24.16e3)') '    "elapsed_wall_seconds": ', elapsed_wall_seconds
    write (unit, '(a)') '  },'
    call write_dry_reference_atmosphere(unit, pressure_half, delta_pressure, layer_l, alpha, &
                                        reference_temperature, a_half, b_half)
    write (unit, '(a)') '  "grid": {'
    write (unit, '(a)') '    "type": "octahedral_gaussian",'
    write (unit, '(a)') '    "ring_order": "south_to_north",'
    call write_real_array(unit, 'mu', mu, .true.)
    call write_integer_array(unit, 'nlon', nlon, .true.)
    call write_ring_offsets(unit, nlon)
    write (unit, '(a)') '  },'
    write (unit, '(a)') '  "output": {'
    write (unit, '(a)') '    "daily_global": "daily_global.csv",'
    write (unit, '(a)') '    "daily_sampling": "online mean of every model step in each solar day; time columns give the day start",'
    write (unit, '(a)') '    "atmospheric_global_means": "area and pressure-mass weighted",'
    if (options%include_moisture) then
      write (unit, '(a)') '    "daily_maximum_wind": "largest grid-point wind speed over all levels and steps of the day, '// &
        'with its longitude, latitude (degrees), level and eta",'
      write (unit, '(a)') '    "water_flux_units": "mm day^-1 in files (kg m^-2 s^-1 times 86400)",'
      write (unit, '(a)') '    "column_water": "precipitable_water uses max(q,0); signed_column_water uses q; '// &
        'negative_column_water = precipitable_water - signed_column_water",'
    end if
    write (unit, '(a,es24.16e3,a)') &
      '    "monthly_period_seconds": ', real(radiation%days_per_month, real64)*radiation%solar_day, ','
    write (unit, '(a)') '    "monthly_sampling": "online mean of every model step in each 30-solar-day interval",'
    write (unit, '(a)') '    "monthly_index_origin": "m0001 is 0001-04-01 through 0001-04-30",'
    write (unit, '(a,i0,a)') '    "monthly_count": ', &
      nint(duration/orbital_period)*radiation%months_per_year, ','
    write (unit, '(a,i0,a)') '    "yearly_snapshot_count": ', 1 + nint(duration/orbital_period), ','
    write (unit, '(a)') '    "yearly_sampling": "instantaneous at every April 1 00:00, including initialization",'
    write (unit, '(a)') '    "monthly_grid_layout": "flat ring-major; longitude varies fastest",'
    write (unit, '(a)') '    "monthly_zonal_layout": "[level][latitude], latitude varies fastest",'
    write (unit, '(a)') '    "monthly_eddy_definition": '// &
      '"monthly zonal mean of product minus product of monthly zonal means",'
    write (unit, '(a)') '    "monthly_surface_temperature": "monthly_surface_temperature_m{month:04d}.bin",'
    write (unit, '(a)') '    "monthly_surface_pressure": "monthly_surface_pressure_m{month:04d}.bin",'
    if (options%include_deep_temperature) then
      write (unit, '(a)') '    "monthly_deep_temperature": "monthly_deep_temperature_m{month:04d}.bin",'
    end if
    if (options%include_moisture) then
      write (unit, '(a)') '    "monthly_precipitation": "monthly_precipitation_m{month:04d}.bin",'
      write (unit, '(a)') '    "monthly_evaporation": "monthly_evaporation_m{month:04d}.bin",'
      write (unit, '(a)') '    "monthly_precipitable_water": "monthly_precipitable_water_m{month:04d}.bin",'
    end if
    if (options%include_land_sea) then
      write (unit, '(a)') '    "static_land_fraction": "land_fraction.bin",'
      write (unit, '(a)') '    "static_surface_height": "surface_height.bin",'
    end if
    write (unit, '(a)') '    "monthly_zonal_fields": "monthly_{name}_m{month:04d}.bin",'
    write (unit, '(a)') '    "yearly_grid_level_fields": '// &
      '"yearly_{name}_y{year:04d}_l{level:02d}.bin",'
    write (unit, '(a)') '    "yearly_grid_surface_fields": "yearly_{name}_y{year:04d}.bin",'
    write (unit, '(a)') '    "yearly_spectral_level_fields": '// &
      '"yearly_{name}_spectral_y{year:04d}_l{level:02d}.bin",'
    write (unit, '(a)') '    "yearly_spectral_surface_field": '// &
      '"yearly_log_surface_pressure_spectral_y{year:04d}.bin",'
    if (options%include_moisture) then
      write (unit, '(a)') '    "yearly_spectral_surface_temperature": '// &
        '"yearly_surface_temperature_spectral_y{year:04d}.bin",'
      write (unit, '(a)') '    "yearly_time": "yearly_time_y{year:04d}.json",'
    end if
    if (options%include_deep_temperature) then
      write (unit, '(a)') '    "yearly_spectral_deep_temperature": '// &
        '"yearly_deep_temperature_spectral_y{year:04d}.bin",'
    end if
    write (unit, '(a)') '    "binary_dtype": "float64 little-endian",'
    write (unit, '(a)') '    "spectral_dtype": "complex128 as interleaved float64 real,imag",'
    write (unit, '(a)') '    "spectral_layout": "rectangular (n,m); n=0:T varies fastest, then m=0:T",'
    write (unit, '(a)') '    "units": {'
    write (unit, '(a)') '      "temperature": "K", "surface_pressure": "Pa",'
    write (unit, '(a)') '      "zeta": "s^-1", "delta": "s^-1", "u": "m s^-1", "v": "m s^-1",'
    if (options%include_moisture) then
      write (unit, '(a)') '      "eddy_uv": "m2 s^-2", "eddy_vt": "K m s^-1",'
      write (unit, '(a)') '      "specific_humidity": "kg kg^-1", "precipitation": "mm day^-1",'
      write (unit, '(a)') '      "evaporation": "mm day^-1", "precipitable_water": "kg m^-2", "eddy_vq": "m s^-1"'
    else
      write (unit, '(a)') '      "eddy_uv": "m2 s^-2", "eddy_vt": "K m s^-1"'
    end if
    write (unit, '(a)') '    }'
    write (unit, '(a)') '  }'
    write (unit, '(a)') '}'
    close (unit)
  end subroutine write_radiation_metadata

end module radiation_case_output
