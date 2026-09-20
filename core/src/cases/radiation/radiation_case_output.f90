!> Output files and metadata schema shared by the dry radiation and slab-ocean cases.
module radiation_case_output
  use iso_fortran_env, only: real64
  use field_binary_writer, only: write_field, write_rectangular_field, write_spectral_field, &
                                 check_finite, check_spectral_finite, check_rectangular_finite
  use json_writer, only: write_real_array, write_integer_array, write_ring_offsets
  use csv_writer, only: write_csv_header, append_csv_row, csv_real, csv_integer
  use dry_case_output, only: write_dry_reference_atmosphere
  use planet_parameters, only: planet_config
  use dry_physics_config, only: radiation_config, convection_config, radiation_days_per_year, &
                                radiation_orbital_period, radiation_surface_heat_capacity
  use dry_radiation, only: radiation_diagnostics
  implicit none
  private

  public :: initialize_radiation_daily_output, append_radiation_daily_output
  public :: write_radiation_monthly_output, write_radiation_yearly_snapshot
  public :: write_radiation_metadata

contains

  subroutine initialize_radiation_daily_output(case_directory, include_deep_temperature)
    character(*), intent(in) :: case_directory
    logical, intent(in) :: include_deep_temperature

    if (include_deep_temperature) then
      call write_csv_header(trim(case_directory)//'/daily_global.csv', &
        'time_seconds,simulation_day,calendar_year,calendar_month,calendar_day,'// &
        'mean_atmospheric_temperature_k,'// &
        'mean_surface_temperature_k,mean_deep_temperature_k,mean_kinetic_energy_j_kg-1,'// &
        'mean_surface_pressure_pa,incoming_shortwave_w_m-2,reflected_shortwave_w_m-2,'// &
        'outgoing_longwave_w_m-2')
    else
      call write_csv_header(trim(case_directory)//'/daily_global.csv', &
        'time_seconds,simulation_day,calendar_year,calendar_month,calendar_day,'// &
        'mean_atmospheric_temperature_k,mean_ocean_temperature_k,mean_kinetic_energy_j_kg-1,'// &
        'mean_surface_pressure_pa,incoming_shortwave_w_m-2,reflected_shortwave_w_m-2,'// &
        'outgoing_longwave_w_m-2')
    end if
  end subroutine initialize_radiation_daily_output

  !> Appends one row of daily means; the calendar date of the day start is supplied by the case runner.
  subroutine append_radiation_daily_output(case_directory, means, simulation_day, &
                                           calendar_year, calendar_month, calendar_day, include_deep_temperature)
    character(*), intent(in) :: case_directory
    type(radiation_diagnostics), intent(in) :: means
    real(real64), intent(in) :: simulation_day
    integer, intent(in) :: calendar_year, calendar_month, calendar_day
    logical, intent(in) :: include_deep_temperature

    if (include_deep_temperature) then
      call append_csv_row(trim(case_directory)//'/daily_global.csv', &
        csv_real(means%time_seconds)//','//csv_real(simulation_day)//','// &
        csv_integer(calendar_year)//','//csv_integer(calendar_month)//','//csv_integer(calendar_day)//','// &
        csv_real(means%mean_atmospheric_temperature)//','//csv_real(means%mean_surface_temperature)//','// &
        csv_real(means%mean_deep_temperature)//','//csv_real(means%mean_kinetic_energy)//','// &
        csv_real(means%mean_surface_pressure)//','//csv_real(means%mean_incoming_shortwave)//','// &
        csv_real(means%mean_reflected_shortwave)//','//csv_real(means%mean_outgoing_longwave))
    else
      call append_csv_row(trim(case_directory)//'/daily_global.csv', &
        csv_real(means%time_seconds)//','//csv_real(simulation_day)//','// &
        csv_integer(calendar_year)//','//csv_integer(calendar_month)//','//csv_integer(calendar_day)//','// &
        csv_real(means%mean_atmospheric_temperature)//','//csv_real(means%mean_surface_temperature)//','// &
        csv_real(means%mean_kinetic_energy)//','//csv_real(means%mean_surface_pressure)//','// &
        csv_real(means%mean_incoming_shortwave)//','//csv_real(means%mean_reflected_shortwave)//','// &
        csv_real(means%mean_outgoing_longwave))
    end if
  end subroutine append_radiation_daily_output

  subroutine write_radiation_monthly_output(case_directory, month, nlon, surface_temperature, &
                                            surface_pressure, zonal_temperature, zonal_u, zonal_v, &
                                            eddy_uv, eddy_vt)
    character(*), intent(in) :: case_directory
    integer, intent(in) :: month
    integer, intent(in) :: nlon(:)
    real(real64), intent(in) :: surface_temperature(:, :), surface_pressure(:, :)
    real(real64), intent(in) :: zonal_temperature(:, :), zonal_u(:, :), zonal_v(:, :)
    real(real64), intent(in) :: eddy_uv(:, :), eddy_vt(:, :)
    character(len=4) :: month_text

    if (month < 1 .or. month > 9999) error stop 'radiation output month is outside the supported range'
    if (size(zonal_temperature, 1) /= size(nlon) .or. &
        any(shape(zonal_temperature) /= shape(zonal_u)) .or. &
        any(shape(zonal_temperature) /= shape(zonal_v)) .or. &
        any(shape(zonal_temperature) /= shape(eddy_uv)) .or. &
        any(shape(zonal_temperature) /= shape(eddy_vt))) then
      error stop 'radiation monthly zonal fields have inconsistent shapes'
    end if
    call check_finite('monthly_surface_temperature', nlon, surface_temperature)
    call check_finite('monthly_surface_pressure', nlon, surface_pressure)
    call check_rectangular_finite('monthly_zonal_temperature', zonal_temperature)
    call check_rectangular_finite('monthly_zonal_u', zonal_u)
    call check_rectangular_finite('monthly_zonal_v', zonal_v)
    call check_rectangular_finite('monthly_eddy_uv', eddy_uv)
    call check_rectangular_finite('monthly_eddy_vt', eddy_vt)
    write (month_text, '(i4.4)') month
    call write_field(trim(case_directory)//'/monthly_surface_temperature_m'//month_text//'.bin', &
                     nlon, surface_temperature)
    call write_field(trim(case_directory)//'/monthly_surface_pressure_m'//month_text//'.bin', &
                     nlon, surface_pressure)
    call write_rectangular_field(trim(case_directory)//'/monthly_zonal_temperature_m'//month_text//'.bin', &
                                 zonal_temperature)
    call write_rectangular_field(trim(case_directory)//'/monthly_zonal_u_m'//month_text//'.bin', zonal_u)
    call write_rectangular_field(trim(case_directory)//'/monthly_zonal_v_m'//month_text//'.bin', zonal_v)
    call write_rectangular_field(trim(case_directory)//'/monthly_eddy_uv_m'//month_text//'.bin', eddy_uv)
    call write_rectangular_field(trim(case_directory)//'/monthly_eddy_vt_m'//month_text//'.bin', eddy_vt)
  end subroutine write_radiation_monthly_output

  subroutine write_radiation_yearly_snapshot(case_directory, year, nlon, &
                                              zeta_spectral, delta_spectral, temperature_spectral, &
                                              log_surface_pressure_spectral, zeta, delta, temperature, &
                                              u, v, log_surface_pressure, surface_temperature, deep_temperature, &
                                              write_deep_temperature)
    character(*), intent(in) :: case_directory
    integer, intent(in) :: year
    integer, intent(in) :: nlon(:)
    complex(real64), intent(in) :: zeta_spectral(0:, 0:, :), delta_spectral(0:, 0:, :)
    complex(real64), intent(in) :: temperature_spectral(0:, 0:, :)
    complex(real64), intent(in) :: log_surface_pressure_spectral(0:, 0:)
    real(real64), intent(in) :: zeta(:, :, :), delta(:, :, :), temperature(:, :, :)
    real(real64), intent(in) :: u(:, :, :), v(:, :, :), log_surface_pressure(:, :)
    real(real64), intent(in) :: surface_temperature(:, :), deep_temperature(:, :)
    logical, intent(in) :: write_deep_temperature
    character(len=4) :: year_text
    character(len=2) :: level_text
    integer :: k

    if (year < 1 .or. year > 9999) error stop 'radiation snapshot year is outside the supported range'
    if (size(zeta, 3) /= size(delta, 3) .or. size(zeta, 3) /= size(temperature, 3) .or. &
        size(zeta, 3) /= size(u, 3) .or. size(zeta, 3) /= size(v, 3) .or. &
        size(zeta, 3) /= size(zeta_spectral, 3) .or. &
        size(zeta, 3) /= size(delta_spectral, 3) .or. &
        size(zeta, 3) /= size(temperature_spectral, 3)) then
      error stop 'radiation yearly snapshot fields have different level counts'
    end if
    write (year_text, '(i4.4)') year
    call check_finite('yearly_log_surface_pressure', nlon, log_surface_pressure)
    call check_finite('yearly_surface_temperature', nlon, surface_temperature)
    if (write_deep_temperature) call check_finite('yearly_deep_temperature', nlon, deep_temperature)
    call check_spectral_finite('yearly_log_surface_pressure_spectral', log_surface_pressure_spectral)
    call write_field(trim(case_directory)//'/yearly_log_surface_pressure_y'//year_text//'.bin', &
                     nlon, log_surface_pressure)
    call write_field(trim(case_directory)//'/yearly_surface_temperature_y'//year_text//'.bin', &
                     nlon, surface_temperature)
    if (write_deep_temperature) then
      call write_field(trim(case_directory)//'/yearly_deep_temperature_y'//year_text//'.bin', &
                       nlon, deep_temperature)
    end if
    call write_spectral_field(trim(case_directory)//'/yearly_log_surface_pressure_spectral_y'// &
                              year_text//'.bin', log_surface_pressure_spectral)
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
    end do
  end subroutine write_radiation_yearly_snapshot

  !> Metadata of one run.  The radiation, convection and planet configuration are
  !> the values the case runner integrated with; this module imports no constants.
  subroutine write_radiation_metadata(case_directory, case_name, radiation, convection, planet, truncation, time_step, &
                                      duration, number_of_steps, maximum_cfl, elapsed_wall_seconds, nlon, mu, &
                                      pressure_half, delta_pressure, layer_l, alpha, reference_temperature, &
                                      a_half, b_half)
    character(*), intent(in) :: case_directory, case_name
    type(radiation_config), intent(in) :: radiation
    type(convection_config), intent(in) :: convection
    type(planet_config), intent(in) :: planet
    integer, intent(in) :: truncation, number_of_steps
    real(real64), intent(in) :: time_step, duration, maximum_cfl, elapsed_wall_seconds
    integer, intent(in) :: nlon(:)
    real(real64), intent(in) :: mu(:), pressure_half(0:), delta_pressure(:), layer_l(:), alpha(:)
    real(real64), intent(in) :: reference_temperature(:), a_half(0:), b_half(0:)
    integer :: unit
    real(real64) :: orbital_period

    orbital_period = radiation_orbital_period(radiation)
    open (newunit=unit, file=trim(case_directory)//'/metadata.json', status='replace', action='write')
    write (unit, '(a)') '{'
    write (unit, '(a)') '  "schema_version": 2,'
    write (unit, '(a)') '  "case_name": "'//trim(case_name)//'",'
    write (unit, '(a)') '  "equation": "dry_hydrostatic_atmosphere",'
    if (radiation%slab_ocean_enabled) then
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
    write (unit, '(a,es24.16e3,a)') &
      '    "longwave_surface_optical_depth": ', radiation%longwave_surface_optical_depth, ','
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
    write (unit, '(a)') '    "orbit_eccentricity": 0'
    write (unit, '(a)') '  },'
    write (unit, '(a)') '  "dry_convective_adjustment": {'
    write (unit, '(a)') '    "method": "enthalpy-conserving pool adjacent violators",'
    write (unit, '(a,es24.16e3)') &
      '    "relaxation_time_seconds": ', convection%adjustment_time
    write (unit, '(a)') '  },'
    if (radiation%slab_ocean_enabled) then
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
    write (unit, '(a)') '    "monthly_zonal_fields": "monthly_{name}_m{month:04d}.bin",'
    write (unit, '(a)') '    "yearly_grid_level_fields": '// &
      '"yearly_{name}_y{year:04d}_l{level:02d}.bin",'
    write (unit, '(a)') '    "yearly_grid_surface_fields": "yearly_{name}_y{year:04d}.bin",'
    write (unit, '(a)') '    "yearly_spectral_level_fields": '// &
      '"yearly_{name}_spectral_y{year:04d}_l{level:02d}.bin",'
    write (unit, '(a)') '    "yearly_spectral_surface_field": '// &
      '"yearly_log_surface_pressure_spectral_y{year:04d}.bin",'
    write (unit, '(a)') '    "binary_dtype": "float64 little-endian",'
    write (unit, '(a)') '    "spectral_dtype": "complex128 as interleaved float64 real,imag",'
    write (unit, '(a)') '    "spectral_layout": "rectangular (n,m); n=0:T varies fastest, then m=0:T",'
    write (unit, '(a)') '    "units": {'
    write (unit, '(a)') '      "temperature": "K", "surface_pressure": "Pa",'
    write (unit, '(a)') '      "zeta": "s^-1", "delta": "s^-1", "u": "m s^-1", "v": "m s^-1",'
    write (unit, '(a)') '      "eddy_uv": "m2 s^-2", "eddy_vt": "K m s^-1"'
    write (unit, '(a)') '    }'
    write (unit, '(a)') '  }'
    write (unit, '(a)') '}'
    close (unit)
  end subroutine write_radiation_metadata

end module radiation_case_output
