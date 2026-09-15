module field_output
  use iso_fortran_env, only: real64
  use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
  use barotropic_vorticity, only: earth_radius, rotation_rate
  use raw_filter, only: raw_filter_epsilon, raw_filter_alpha
  use spectral_hyperdiffusion, only: hyperdiffusion_order, &
                                       hyperdiffusion_timescale_seconds
  use dry_radiation, only: solar_constant, surface_shortwave_albedo, axial_tilt, orbital_period, &
                           planetary_rotation_rate, longwave_surface_optical_depth, &
                           stefan_boltzmann_constant, dry_air_specific_heat, dry_gravity_acceleration, &
                           surface_heat_capacity, deep_ground_heat_capacity, ground_exchange_coefficient, &
                           surface_exchange_coefficient, gustiness_speed, solar_day, days_per_month, &
                           months_per_year, days_per_year, radiation_calendar_date
  implicit none
  private

  public :: make_directory
  public :: write_snapshot
  public :: write_run_metadata
  public :: write_shallow_water_snapshot
  public :: write_shallow_water_metadata
  public :: write_dry_snapshot
  public :: write_dry_metadata
  public :: initialize_radiation_daily_output
  public :: append_radiation_daily_output
  public :: write_radiation_monthly_output
  public :: write_radiation_yearly_snapshot
  public :: write_radiation_metadata

contains

  subroutine make_directory(path)
    character(*), intent(in) :: path
    integer :: exit_status

    call execute_command_line('mkdir -p "'//trim(path)//'"', exitstat=exit_status)
    if (exit_status /= 0) error stop 'failed to create output directory'
  end subroutine make_directory

  subroutine write_snapshot(case_directory, step, nlon, zeta, u, v)
    character(*), intent(in) :: case_directory
    integer, intent(in) :: step
    integer, intent(in) :: nlon(:)
    real(real64), intent(in) :: zeta(:, :), u(:, :), v(:, :)
    character(len=5) :: step_text

    if (step < 0 .or. step > 99999) error stop 'output step is outside the supported range'
    call check_finite('zeta', nlon, zeta)
    call check_finite('u', nlon, u)
    call check_finite('v', nlon, v)
    write (step_text, '(i5.5)') step
    call write_field(trim(case_directory)//'/zeta_'//step_text//'.bin', nlon, zeta)
    call write_field(trim(case_directory)//'/u_'//step_text//'.bin', nlon, u)
    call write_field(trim(case_directory)//'/v_'//step_text//'.bin', nlon, v)
  end subroutine write_snapshot

  subroutine write_shallow_water_snapshot(case_directory, step, nlon, &
                                          zeta_spectral, delta_spectral, eta_spectral, &
                                          zeta, delta, eta, u, v)
    character(*), intent(in) :: case_directory
    integer, intent(in) :: step
    integer, intent(in) :: nlon(:)
    complex(real64), intent(in) :: zeta_spectral(0:, 0:), delta_spectral(0:, 0:), eta_spectral(0:, 0:)
    real(real64), intent(in) :: zeta(:, :), delta(:, :), eta(:, :), u(:, :), v(:, :)
    character(len=5) :: step_text

    if (step < 0 .or. step > 99999) error stop 'output step is outside the supported range'
    call check_finite('zeta', nlon, zeta)
    call check_finite('delta', nlon, delta)
    call check_finite('eta', nlon, eta)
    call check_finite('u', nlon, u)
    call check_finite('v', nlon, v)
    call check_spectral_finite('zeta_spectral', zeta_spectral)
    call check_spectral_finite('delta_spectral', delta_spectral)
    call check_spectral_finite('eta_spectral', eta_spectral)
    write (step_text, '(i5.5)') step
    call write_field(trim(case_directory)//'/zeta_'//step_text//'.bin', nlon, zeta)
    call write_field(trim(case_directory)//'/delta_'//step_text//'.bin', nlon, delta)
    call write_field(trim(case_directory)//'/eta_'//step_text//'.bin', nlon, eta)
    call write_field(trim(case_directory)//'/u_'//step_text//'.bin', nlon, u)
    call write_field(trim(case_directory)//'/v_'//step_text//'.bin', nlon, v)
    call write_spectral_field(trim(case_directory)//'/zeta_spectral_'//step_text//'.bin', zeta_spectral)
    call write_spectral_field(trim(case_directory)//'/delta_spectral_'//step_text//'.bin', delta_spectral)
    call write_spectral_field(trim(case_directory)//'/eta_spectral_'//step_text//'.bin', eta_spectral)
  end subroutine write_shallow_water_snapshot

  subroutine write_dry_snapshot(case_directory, step, nlon, zeta, delta, temperature, &
                                surface_pressure, u, v, surface_temperature, deep_temperature)
    character(*), intent(in) :: case_directory
    integer, intent(in) :: step
    integer, intent(in) :: nlon(:)
    real(real64), intent(in) :: zeta(:, :, :), delta(:, :, :), temperature(:, :, :)
    real(real64), intent(in) :: surface_pressure(:, :), u(:, :, :), v(:, :, :)
    real(real64), intent(in), optional :: surface_temperature(:, :), deep_temperature(:, :)
    character(len=5) :: step_text
    character(len=2) :: level_text
    integer :: k

    if (step < 0 .or. step > 99999) error stop 'output step is outside the supported range'
    if (size(zeta, 3) /= size(delta, 3) .or. size(zeta, 3) /= size(temperature, 3) .or. &
        size(zeta, 3) /= size(u, 3) .or. size(zeta, 3) /= size(v, 3)) then
      error stop 'dry output fields have different level counts'
    end if
    call check_finite('surface_pressure', nlon, surface_pressure)
    write (step_text, '(i5.5)') step
    call write_field(trim(case_directory)//'/surface_pressure_'//step_text//'.bin', nlon, surface_pressure)
    if (present(surface_temperature) .neqv. present(deep_temperature)) then
      error stop 'dry output ground temperatures must be supplied together'
    end if
    if (present(surface_temperature)) then
      call check_finite('surface_temperature', nlon, surface_temperature)
      call check_finite('deep_temperature', nlon, deep_temperature)
      call write_field(trim(case_directory)//'/surface_temperature_'//step_text//'.bin', nlon, surface_temperature)
      call write_field(trim(case_directory)//'/deep_temperature_'//step_text//'.bin', nlon, deep_temperature)
    end if
    do k = 1, size(zeta, 3)
      write (level_text, '(i2.2)') k
      call check_finite('zeta', nlon, zeta(:, :, k))
      call check_finite('delta', nlon, delta(:, :, k))
      call check_finite('temperature', nlon, temperature(:, :, k))
      call check_finite('u', nlon, u(:, :, k))
      call check_finite('v', nlon, v(:, :, k))
      call write_field(trim(case_directory)//'/zeta_l'//level_text//'_'//step_text//'.bin', nlon, zeta(:, :, k))
      call write_field(trim(case_directory)//'/delta_l'//level_text//'_'//step_text//'.bin', nlon, delta(:, :, k))
      call write_field(trim(case_directory)//'/temperature_l'//level_text//'_'//step_text//'.bin', &
                       nlon, temperature(:, :, k))
      call write_field(trim(case_directory)//'/u_l'//level_text//'_'//step_text//'.bin', nlon, u(:, :, k))
      call write_field(trim(case_directory)//'/v_l'//level_text//'_'//step_text//'.bin', nlon, v(:, :, k))
    end do
  end subroutine write_dry_snapshot

  subroutine initialize_radiation_daily_output(case_directory)
    character(*), intent(in) :: case_directory
    integer :: unit

    open (newunit=unit, file=trim(case_directory)//'/daily_global.csv', status='replace', action='write')
    write (unit, '(a)') 'time_seconds,simulation_day,calendar_year,calendar_month,calendar_day,'// &
      'mean_atmospheric_temperature_k,'// &
      'mean_surface_temperature_k,mean_deep_temperature_k,mean_kinetic_energy_j_kg-1,'// &
      'mean_surface_pressure_pa,incoming_shortwave_w_m-2,reflected_shortwave_w_m-2,'// &
      'outgoing_longwave_w_m-2'
    close (unit)
  end subroutine initialize_radiation_daily_output

  subroutine append_radiation_daily_output(case_directory, time_seconds, mean_atmospheric_temperature, &
                                           mean_surface_temperature, mean_deep_temperature, &
                                           mean_kinetic_energy, mean_surface_pressure, &
                                           mean_incoming_shortwave, mean_reflected_shortwave, &
                                           mean_outgoing_longwave)
    character(*), intent(in) :: case_directory
    real(real64), intent(in) :: time_seconds, mean_atmospheric_temperature
    real(real64), intent(in) :: mean_surface_temperature, mean_deep_temperature
    real(real64), intent(in) :: mean_kinetic_energy, mean_surface_pressure
    real(real64), intent(in) :: mean_incoming_shortwave, mean_reflected_shortwave, mean_outgoing_longwave
    integer :: unit, calendar_year, calendar_month, calendar_day
    real(real64) :: seconds_of_day

    call radiation_calendar_date(time_seconds, calendar_year, calendar_month, calendar_day, seconds_of_day)
    open (newunit=unit, file=trim(case_directory)//'/daily_global.csv', status='old', &
          position='append', action='write')
    write (unit, '(es24.16e3,",",es24.16e3,3(",",i0),8(",",es24.16e3))') &
      time_seconds, time_seconds/solar_day, calendar_year, calendar_month, calendar_day, &
      mean_atmospheric_temperature, mean_surface_temperature, mean_deep_temperature, &
      mean_kinetic_energy, mean_surface_pressure, mean_incoming_shortwave, &
      mean_reflected_shortwave, mean_outgoing_longwave
    close (unit)
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
                                              u, v, log_surface_pressure, surface_temperature, deep_temperature)
    character(*), intent(in) :: case_directory
    integer, intent(in) :: year
    integer, intent(in) :: nlon(:)
    complex(real64), intent(in) :: zeta_spectral(0:, 0:, :), delta_spectral(0:, 0:, :)
    complex(real64), intent(in) :: temperature_spectral(0:, 0:, :)
    complex(real64), intent(in) :: log_surface_pressure_spectral(0:, 0:)
    real(real64), intent(in) :: zeta(:, :, :), delta(:, :, :), temperature(:, :, :)
    real(real64), intent(in) :: u(:, :, :), v(:, :, :), log_surface_pressure(:, :)
    real(real64), intent(in) :: surface_temperature(:, :), deep_temperature(:, :)
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
    call check_finite('yearly_deep_temperature', nlon, deep_temperature)
    call check_spectral_finite('yearly_log_surface_pressure_spectral', log_surface_pressure_spectral)
    call write_field(trim(case_directory)//'/yearly_log_surface_pressure_y'//year_text//'.bin', &
                     nlon, log_surface_pressure)
    call write_field(trim(case_directory)//'/yearly_surface_temperature_y'//year_text//'.bin', &
                     nlon, surface_temperature)
    call write_field(trim(case_directory)//'/yearly_deep_temperature_y'//year_text//'.bin', &
                     nlon, deep_temperature)
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

  subroutine write_dry_metadata(case_directory, initial_condition, truncation, time_step, duration, &
                                number_of_steps, snapshot_interval_steps, maximum_cfl, &
                                elapsed_wall_seconds, nlon, mu, &
                                pressure_half, delta_pressure, layer_l, alpha, reference_temperature, &
                                a_half, b_half, include_ground_temperatures)
    character(*), intent(in) :: case_directory, initial_condition
    integer, intent(in) :: truncation, number_of_steps, snapshot_interval_steps
    real(real64), intent(in) :: time_step, duration, maximum_cfl, elapsed_wall_seconds
    integer, intent(in) :: nlon(:)
    real(real64), intent(in) :: mu(:), pressure_half(0:), delta_pressure(:), layer_l(:), alpha(:)
    real(real64), intent(in) :: reference_temperature(:)
    real(real64), intent(in) :: a_half(0:), b_half(0:)
    logical, intent(in), optional :: include_ground_temperatures
    integer :: unit, number_of_snapshots
    logical :: write_ground_temperatures

    write_ground_temperatures = .false.
    if (present(include_ground_temperatures)) write_ground_temperatures = include_ground_temperatures
    number_of_snapshots = number_of_steps/snapshot_interval_steps + 1
    if (mod(number_of_steps, snapshot_interval_steps) /= 0) number_of_snapshots = number_of_snapshots + 1
    open (newunit=unit, file=trim(case_directory)//'/metadata.json', status='replace', action='write')
    write (unit, '(a)') '{'
    write (unit, '(a)') '  "schema_version": 1,'
    write (unit, '(a)') '  "equation": "dry_hydrostatic_atmosphere",'
    write (unit, '(a)') '  "initial_condition": "'//initial_condition//'",'
    write (unit, '(a,es24.16e3,a)') '  "duration_seconds": ', duration, ','
    write (unit, '(a,es24.16e3,a)') '  "time_step_seconds": ', time_step, ','
    write (unit, '(a,i0,a)') '  "number_of_steps": ', number_of_steps, ','
    write (unit, '(a,i0,a)') '  "snapshot_interval_steps": ', snapshot_interval_steps, ','
    write (unit, '(a,i0,a)') '  "number_of_snapshots": ', number_of_snapshots, ','
    write (unit, '(a,i0,a)') '  "spectral_truncation": ', truncation, ','
    write (unit, '(a,i0,a)') '  "number_of_levels": ', size(reference_temperature), ','
    write (unit, '(a)') '  "time_integrator": "semi-implicit RAW-filtered leapfrog",'
    write (unit, '(a)') '  "hyperdiffusion_e_folding_seconds": {"zeta":14400,"delta":3600,"temperature":14400},'
    write (unit, '(a,es24.16e3,a)') '  "maximum_advective_cfl": ', maximum_cfl, ','
    write (unit, '(a,es24.16e3,a)') '  "elapsed_wall_seconds": ', elapsed_wall_seconds, ','
    write (unit, '(a)', advance='no') '  "reference_half_level_pressure_pa": ['
    call write_inline_real_values(unit, pressure_half)
    write (unit, '(a)', advance='no') '  "reference_layer_pressure_thickness_pa": ['
    call write_inline_real_values(unit, delta_pressure)
    write (unit, '(a)', advance='no') '  "reference_layer_log_pressure_thickness": ['
    call write_inline_real_values(unit, layer_l)
    write (unit, '(a)', advance='no') '  "reference_alpha": ['
    call write_inline_real_values(unit, alpha)
    write (unit, '(a)', advance='no') '  "reference_temperature_k": ['
    call write_inline_real_values(unit, reference_temperature)
    write (unit, '(a)', advance='no') '  "hybrid_a_half_pa": ['
    call write_inline_real_values(unit, a_half)
    write (unit, '(a)', advance='no') '  "hybrid_b_half": ['
    call write_inline_real_values(unit, b_half)
    write (unit, '(a)') '  "grid": {'
    write (unit, '(a)') '    "type": "octahedral_gaussian",'
    call write_real_array(unit, 'mu', mu, .true.)
    call write_integer_array(unit, 'nlon', nlon, .false.)
    write (unit, '(a)') '  },'
    write (unit, '(a)') '  "output": {'
    write (unit, '(a)') '    "dtype": "float64 little-endian",'
    write (unit, '(a)') '    "layout": "flat ring-major",'
    write (unit, '(a)') '    "surface_pressure": "surface_pressure_{step:05d}.bin",'
    if (write_ground_temperatures) then
      write (unit, '(a)') '    "surface_temperature": "surface_temperature_{step:05d}.bin",'
      write (unit, '(a)') '    "deep_temperature": "deep_temperature_{step:05d}.bin",'
    end if
    write (unit, '(a)') '    "level_fields": "{name}_l{level:02d}_{step:05d}.bin"'
    write (unit, '(a)') '  }'
    write (unit, '(a)') '}'
    close (unit)
  end subroutine write_dry_metadata

  subroutine write_radiation_metadata(case_directory, truncation, time_step, duration, number_of_steps, &
                                      maximum_cfl, elapsed_wall_seconds, nlon, mu, pressure_half, &
                                      delta_pressure, layer_l, alpha, reference_temperature, a_half, b_half)
    character(*), intent(in) :: case_directory
    integer, intent(in) :: truncation, number_of_steps
    real(real64), intent(in) :: time_step, duration, maximum_cfl, elapsed_wall_seconds
    integer, intent(in) :: nlon(:)
    real(real64), intent(in) :: mu(:), pressure_half(0:), delta_pressure(:), layer_l(:), alpha(:)
    real(real64), intent(in) :: reference_temperature(:), a_half(0:), b_half(0:)
    integer :: unit

    open (newunit=unit, file=trim(case_directory)//'/metadata.json', status='replace', action='write')
    write (unit, '(a)') '{'
    write (unit, '(a)') '  "schema_version": 2,'
    write (unit, '(a)') '  "case_name": "dry_radiation",'
    write (unit, '(a)') '  "equation": "dry_hydrostatic_atmosphere",'
    write (unit, '(a)') '  "initial_condition": '// &
      '"Jablonowski-Williamson perturbation, flat terrain, two-layer ground",'
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
    write (unit, '(a,es24.16e3,a)') '    "solar_day_seconds": ', solar_day, ','
    write (unit, '(a,i0,a)') '    "days_per_month": ', days_per_month, ','
    write (unit, '(a,i0,a)') '    "months_per_year": ', months_per_year, ','
    write (unit, '(a,i0)') '    "days_per_year": ', days_per_year
    write (unit, '(a)') '  },'
    write (unit, '(a)') '  "radiation": {'
    write (unit, '(a,es24.16e3,a)') '    "solar_constant_w_m-2": ', solar_constant, ','
    write (unit, '(a,es24.16e3,a)') '    "surface_shortwave_albedo": ', surface_shortwave_albedo, ','
    write (unit, '(a,es24.16e3,a)') '    "axial_tilt_radians": ', axial_tilt, ','
    write (unit, '(a,es24.16e3,a)') '    "orbital_period_seconds": ', orbital_period, ','
    write (unit, '(a,es24.16e3,a)') '    "rotation_rate_rad_s": ', planetary_rotation_rate, ','
    write (unit, '(a,es24.16e3,a)') &
      '    "longwave_surface_optical_depth": ', longwave_surface_optical_depth, ','
    write (unit, '(a,es24.16e3,a)') &
      '    "stefan_boltzmann_constant_w_m-2_k-4": ', stefan_boltzmann_constant, ','
    write (unit, '(a,es24.16e3,a)') '    "dry_air_specific_heat_j_kg-1_k-1": ', dry_air_specific_heat, ','
    write (unit, '(a,es24.16e3,a)') '    "gravity_acceleration_m_s-2": ', dry_gravity_acceleration, ','
    write (unit, '(a)') '    "shortwave_atmosphere": "transparent",'
    write (unit, '(a)') '    "shortwave_forcing": "instantaneous local solar zenith angle",'
    write (unit, '(a)') '    "orbit_eccentricity": 0'
    write (unit, '(a)') '  },'
    write (unit, '(a)') '  "ground": {'
    write (unit, '(a,es24.16e3,a)') '    "surface_heat_capacity_j_m-2_k-1": ', surface_heat_capacity, ','
    write (unit, '(a,es24.16e3,a)') '    "deep_heat_capacity_j_m-2_k-1": ', deep_ground_heat_capacity, ','
    write (unit, '(a,es24.16e3,a)') '    "exchange_coefficient_w_m-2_k-1": ', ground_exchange_coefficient, ','
    write (unit, '(a,es24.16e3,a)') '    "surface_exchange_coefficient": ', surface_exchange_coefficient, ','
    write (unit, '(a,es24.16e3)') '    "gustiness_speed_m_s-1": ', gustiness_speed
    write (unit, '(a)') '  },'
    write (unit, '(a)') '  "numerics": {'
    write (unit, '(a,i0,a)') '    "spectral_truncation": ', truncation, ','
    write (unit, '(a)') '    "time_integrator": "semi-implicit RAW-filtered leapfrog",'
    write (unit, '(a,es24.16e3,a)') '    "maximum_advective_cfl": ', maximum_cfl, ','
    write (unit, '(a,es24.16e3)') '    "elapsed_wall_seconds": ', elapsed_wall_seconds
    write (unit, '(a)') '  },'
    write (unit, '(a)', advance='no') '  "reference_half_level_pressure_pa": ['
    call write_inline_real_values(unit, pressure_half)
    write (unit, '(a)', advance='no') '  "reference_layer_pressure_thickness_pa": ['
    call write_inline_real_values(unit, delta_pressure)
    write (unit, '(a)', advance='no') '  "reference_layer_log_pressure_thickness": ['
    call write_inline_real_values(unit, layer_l)
    write (unit, '(a)', advance='no') '  "reference_alpha": ['
    call write_inline_real_values(unit, alpha)
    write (unit, '(a)', advance='no') '  "reference_temperature_k": ['
    call write_inline_real_values(unit, reference_temperature)
    write (unit, '(a)', advance='no') '  "hybrid_a_half_pa": ['
    call write_inline_real_values(unit, a_half)
    write (unit, '(a)', advance='no') '  "hybrid_b_half": ['
    call write_inline_real_values(unit, b_half)
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
      '    "monthly_period_seconds": ', real(days_per_month, real64)*solar_day, ','
    write (unit, '(a)') '    "monthly_sampling": "online mean of every model step in each 30-solar-day interval",'
    write (unit, '(a)') '    "monthly_index_origin": "m0001 is 0001-04-01 through 0001-04-30",'
    write (unit, '(a,i0,a)') '    "monthly_count": ', &
      nint(duration/orbital_period)*months_per_year, ','
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

  subroutine write_inline_real_values(unit, values)
    integer, intent(in) :: unit
    real(real64), intent(in) :: values(:)
    integer :: i

    do i = 1, size(values)
      write (unit, '(es24.16e3)', advance='no') values(i)
      if (i < size(values)) write (unit, '(a)', advance='no') ','
    end do
    write (unit, '(a)') '],'
  end subroutine write_inline_real_values

  subroutine write_field(path, nlon, field)
    character(*), intent(in) :: path
    integer, intent(in) :: nlon(:)
    real(real64), intent(in) :: field(:, :)
    integer :: unit, j

    open (newunit=unit, file=path, status='replace', access='stream', &
          form='unformatted', action='write', convert='little_endian')
    do j = 1, size(nlon)
      write (unit) field(1:nlon(j), j)
    end do
    close (unit)
  end subroutine write_field

  subroutine write_rectangular_field(path, field)
    character(*), intent(in) :: path
    real(real64), intent(in) :: field(:, :)
    integer :: unit, k

    open (newunit=unit, file=path, status='replace', access='stream', &
          form='unformatted', action='write', convert='little_endian')
    do k = 1, size(field, 2)
      write (unit) field(:, k)
    end do
    close (unit)
  end subroutine write_rectangular_field

  subroutine write_spectral_field(path, field)
    character(*), intent(in) :: path
    complex(real64), intent(in) :: field(0:, 0:)
    integer :: unit, n, m

    open (newunit=unit, file=path, status='replace', access='stream', &
          form='unformatted', action='write', convert='little_endian')
    ! State arrays carry an n=T+1 work row for differentiated fields.  It is
    ! constrained to zero and is not part of the triangular T spectrum.
    do m = 0, ubound(field, 2)
      do n = 0, ubound(field, 1) - 1
        write (unit) real(field(n, m), real64), aimag(field(n, m))
      end do
    end do
    close (unit)
  end subroutine write_spectral_field

  subroutine write_run_metadata(case_directory, case_name, initial_condition_json, &
                                truncation, time_step, duration, number_of_steps, &
                                snapshot_interval_steps, maximum_cfl, elapsed_wall_seconds, &
                                nlon, mu)
    character(*), intent(in) :: case_directory, case_name, initial_condition_json
    integer, intent(in) :: truncation, number_of_steps, snapshot_interval_steps
    real(real64), intent(in) :: time_step, duration, maximum_cfl, elapsed_wall_seconds
    integer, intent(in) :: nlon(:)
    real(real64), intent(in) :: mu(:)
    integer :: unit, point_count, number_of_snapshots

    if (size(nlon) /= size(mu)) error stop 'metadata grid arrays have different sizes'
    if (snapshot_interval_steps <= 0) error stop 'snapshot interval must be positive'
    point_count = sum(nlon)
    number_of_snapshots = number_of_steps/snapshot_interval_steps + 1
    if (mod(number_of_steps, snapshot_interval_steps) /= 0) then
      number_of_snapshots = number_of_snapshots + 1
    end if
    open (newunit=unit, file=trim(case_directory)//'/metadata.json', &
          status='replace', action='write')
    write (unit, '(a)') '{'
    write (unit, '(a)') '  "schema_version": 1,'
    write (unit, '(3a)') '  "case_name": "', trim(case_name), '",'
    write (unit, '(3a)') '  "initial_condition": ', trim(initial_condition_json), ','
    write (unit, '(a)') '  "simulation": {'
    write (unit, '(a,es24.16e3,a)') '    "duration_seconds": ', duration, ','
    write (unit, '(a,es24.16e3,a)') '    "time_step_seconds": ', time_step, ','
    write (unit, '(a,i0,a)') '    "number_of_steps": ', number_of_steps, ','
    write (unit, '(a,i0,a)') '    "snapshot_interval_steps": ', snapshot_interval_steps, ','
    write (unit, '(a,es24.16e3,a)') &
      '    "snapshot_interval_seconds": ', time_step*real(snapshot_interval_steps, real64), ','
    write (unit, '(a,i0,a)') '    "number_of_snapshots": ', number_of_snapshots, ','
    write (unit, '(a)') '    "snapshot_time_seconds": "step * time_step_seconds"'
    write (unit, '(a)') '  },'
    write (unit, '(a)') '  "physical_constants": {'
    write (unit, '(a,es24.16e3,a)') '    "earth_radius_m": ', earth_radius, ','
    write (unit, '(a,es24.16e3)') '    "rotation_rate_rad_s": ', rotation_rate
    write (unit, '(a)') '  },'
    write (unit, '(a)') '  "numerics": {'
    write (unit, '(a,i0,a)') '    "spectral_truncation": ', truncation, ','
    write (unit, '(a)') '    "time_integrator": "RAW-filtered leapfrog",'
    write (unit, '(a)') '    "raw_filter": {'
    write (unit, '(a,es24.16e3,a)') '      "epsilon": ', raw_filter_epsilon, ','
    write (unit, '(a,es24.16e3)') '      "alpha": ', raw_filter_alpha
    write (unit, '(a)') '    },'
    write (unit, '(a)') '    "hyperdiffusion": {'
    write (unit, '(a,i0,a)') '      "order": ', hyperdiffusion_order, ','
    write (unit, '(a,es24.16e3)') &
      '      "e_folding_time_at_truncation_seconds": ', hyperdiffusion_timescale_seconds
    write (unit, '(a)') '    },'
    write (unit, '(a,es24.16e3)') '    "maximum_cfl": ', maximum_cfl
    write (unit, '(a)') '  },'
    write (unit, '(a)') '  "performance": {'
    write (unit, '(a,es24.16e3)') '    "elapsed_wall_seconds": ', elapsed_wall_seconds
    write (unit, '(a)') '  },'
    write (unit, '(a)') '  "grid": {'
    write (unit, '(a)') '    "type": "octahedral_gaussian",'
    write (unit, '(a)') '    "latitude_coordinate": "mu = sin(latitude_radians)",'
    write (unit, '(a)') &
      '    "longitude_radians": "2*pi*k/nlon[j], k=0,...,nlon[j]-1",'
    write (unit, '(a)') '    "ring_order": "south_to_north",'
    write (unit, '(a,i0,a)') '    "number_of_latitudes": ', size(nlon), ','
    write (unit, '(a,i0,a)') '    "point_count": ', point_count, ','
    call write_real_array(unit, 'mu', mu, .true.)
    call write_integer_array(unit, 'nlon', nlon, .true.)
    call write_ring_offsets(unit, nlon)
    write (unit, '(a)') '  },'
    write (unit, '(a)') '  "output": {'
    write (unit, '(a)') '    "dtype": "float64",'
    write (unit, '(a)') '    "byte_order": "little_endian",'
    write (unit, '(a)') '    "layout": "flat ring-major; longitude index varies fastest",'
    write (unit, '(a)') '    "file_patterns": {'
    write (unit, '(a)') '      "zeta": "zeta_{step:05d}.bin",'
    write (unit, '(a)') '      "u": "u_{step:05d}.bin",'
    write (unit, '(a)') '      "v": "v_{step:05d}.bin"'
    write (unit, '(a)') '    },'
    write (unit, '(a)') '    "units": {'
    write (unit, '(a)') '      "zeta": "s^-1",'
    write (unit, '(a)') '      "u": "m s^-1",'
    write (unit, '(a)') '      "v": "m s^-1"'
    write (unit, '(a)') '    }'
    write (unit, '(a)') '  }'
    write (unit, '(a)') '}'
    close (unit)
  end subroutine write_run_metadata

  subroutine write_shallow_water_metadata(case_directory, case_name, initial_condition_json, &
                                           truncation, time_step, duration, number_of_steps, &
                                           snapshot_interval_steps, maximum_cfl, elapsed_wall_seconds, &
                                           nlon, mu, gravity, depth, implicitness)
    character(*), intent(in) :: case_directory, case_name, initial_condition_json
    integer, intent(in) :: truncation, number_of_steps, snapshot_interval_steps
    real(real64), intent(in) :: time_step, duration, maximum_cfl, elapsed_wall_seconds
    real(real64), intent(in) :: gravity, depth, implicitness
    integer, intent(in) :: nlon(:)
    real(real64), intent(in) :: mu(:)
    integer :: unit, point_count, number_of_snapshots

    if (size(nlon) /= size(mu)) error stop 'metadata grid arrays have different sizes'
    if (snapshot_interval_steps <= 0) error stop 'snapshot interval must be positive'
    point_count = sum(nlon)
    number_of_snapshots = number_of_steps/snapshot_interval_steps + 1
    if (mod(number_of_steps, snapshot_interval_steps) /= 0) number_of_snapshots = number_of_snapshots + 1
    open (newunit=unit, file=trim(case_directory)//'/metadata.json', status='replace', action='write')
    write (unit, '(a)') '{'
    write (unit, '(a)') '  "schema_version": 1,'
    write (unit, '(a)') '  "equation": "shallow_water",'
    write (unit, '(3a)') '  "case_name": "', trim(case_name), '",'
    write (unit, '(3a)') '  "initial_condition": ', trim(initial_condition_json), ','
    write (unit, '(a)') '  "simulation": {'
    write (unit, '(a,es24.16e3,a)') '    "duration_seconds": ', duration, ','
    write (unit, '(a,es24.16e3,a)') '    "time_step_seconds": ', time_step, ','
    write (unit, '(a,i0,a)') '    "number_of_steps": ', number_of_steps, ','
    write (unit, '(a,i0,a)') '    "snapshot_interval_steps": ', snapshot_interval_steps, ','
    write (unit, '(a,es24.16e3,a)') &
      '    "snapshot_interval_seconds": ', time_step*real(snapshot_interval_steps, real64), ','
    write (unit, '(a,i0,a)') '    "number_of_snapshots": ', number_of_snapshots, ','
    write (unit, '(a)') '    "snapshot_time_seconds": "step * time_step_seconds"'
    write (unit, '(a)') '  },'
    write (unit, '(a)') '  "physical_constants": {'
    write (unit, '(a,es24.16e3,a)') '    "earth_radius_m": ', earth_radius, ','
    write (unit, '(a,es24.16e3,a)') '    "rotation_rate_rad_s": ', rotation_rate, ','
    write (unit, '(a,es24.16e3,a)') '    "gravity_acceleration_m_s-2": ', gravity, ','
    write (unit, '(a,es24.16e3)') '    "mean_depth_m": ', depth
    write (unit, '(a)') '  },'
    write (unit, '(a)') '  "numerics": {'
    write (unit, '(a,i0,a)') '    "spectral_truncation": ', truncation, ','
    write (unit, '(a)') '    "time_integrator": "semi-implicit RAW-filtered leapfrog",'
    write (unit, '(a,es24.16e3,a)') '    "gravity_wave_implicitness_beta": ', implicitness, ','
    write (unit, '(a)') '    "raw_filter": {'
    write (unit, '(a,es24.16e3,a)') '      "epsilon": ', raw_filter_epsilon, ','
    write (unit, '(a,es24.16e3)') '      "alpha": ', raw_filter_alpha
    write (unit, '(a)') '    },'
    write (unit, '(a)') '    "hyperdiffusion": {'
    write (unit, '(a)') '      "applied_to": ["zeta", "delta"],'
    write (unit, '(a,i0,a)') '      "order": ', hyperdiffusion_order, ','
    write (unit, '(a,es24.16e3)') &
      '      "e_folding_time_at_truncation_seconds": ', hyperdiffusion_timescale_seconds
    write (unit, '(a)') '    },'
    write (unit, '(a,es24.16e3,a)') '    "maximum_cfl": ', maximum_cfl, ','
    write (unit, '(a)') '    "cfl_kind": "advective"'
    write (unit, '(a)') '  },'
    write (unit, '(a)') '  "performance": {'
    write (unit, '(a,es24.16e3)') '    "elapsed_wall_seconds": ', elapsed_wall_seconds
    write (unit, '(a)') '  },'
    write (unit, '(a)') '  "grid": {'
    write (unit, '(a)') '    "type": "octahedral_gaussian",'
    write (unit, '(a)') '    "latitude_coordinate": "mu = sin(latitude_radians)",'
    write (unit, '(a)') '    "longitude_radians": "2*pi*k/nlon[j], k=0,...,nlon[j]-1",'
    write (unit, '(a)') '    "ring_order": "south_to_north",'
    write (unit, '(a,i0,a)') '    "number_of_latitudes": ', size(nlon), ','
    write (unit, '(a,i0,a)') '    "point_count": ', point_count, ','
    call write_real_array(unit, 'mu', mu, .true.)
    call write_integer_array(unit, 'nlon', nlon, .true.)
    call write_ring_offsets(unit, nlon)
    write (unit, '(a)') '  },'
    write (unit, '(a)') '  "output": {'
    write (unit, '(a)') '    "grid_fields": {'
    write (unit, '(a)') '      "dtype": "float64",'
    write (unit, '(a)') '      "byte_order": "little_endian",'
    write (unit, '(a)') '      "layout": "flat ring-major; longitude index varies fastest",'
    write (unit, '(a)') '      "file_patterns": {'
    write (unit, '(a)') '        "zeta": "zeta_{step:05d}.bin",'
    write (unit, '(a)') '        "delta": "delta_{step:05d}.bin",'
    write (unit, '(a)') '        "eta": "eta_{step:05d}.bin",'
    write (unit, '(a)') '        "u": "u_{step:05d}.bin",'
    write (unit, '(a)') '        "v": "v_{step:05d}.bin"'
    write (unit, '(a)') '      },'
    write (unit, '(a)') '      "units": {"zeta":"s^-1","delta":"s^-1","eta":"m","u":"m s^-1","v":"m s^-1"}'
    write (unit, '(a)') '    },'
    write (unit, '(a)') '    "spectral_fields": {'
    write (unit, '(a)') '      "dtype": "complex128 as interleaved float64 real,imag",'
    write (unit, '(a)') '      "byte_order": "little_endian",'
    write (unit, '(a)') '      "layout": "rectangular (n,m); n=0:T varies fastest, then m=0:T",'
    write (unit, '(a,i0,a,i0,a)') '      "shape": [', truncation + 1, ',', truncation + 1, '],'
    write (unit, '(a)') '      "file_patterns": {'
    write (unit, '(a)') '        "zeta": "zeta_spectral_{step:05d}.bin",'
    write (unit, '(a)') '        "delta": "delta_spectral_{step:05d}.bin",'
    write (unit, '(a)') '        "eta": "eta_spectral_{step:05d}.bin"'
    write (unit, '(a)') '      },'
    write (unit, '(a)') '      "units": {"zeta":"s^-1","delta":"s^-1","eta":"m"}'
    write (unit, '(a)') '    }'
    write (unit, '(a)') '  }'
    write (unit, '(a)') '}'
    close (unit)
  end subroutine write_shallow_water_metadata

  subroutine write_real_array(unit, name, values, trailing_comma)
    integer, intent(in) :: unit
    character(*), intent(in) :: name
    real(real64), intent(in) :: values(:)
    logical, intent(in) :: trailing_comma
    integer :: j

    write (unit, '(3a)', advance='no') '    "', trim(name), '": ['
    do j = 1, size(values)
      write (unit, '(es24.16e3)', advance='no') values(j)
      if (j < size(values)) write (unit, '(a)', advance='no') ','
    end do
    if (trailing_comma) then
      write (unit, '(a)') '],'
    else
      write (unit, '(a)') ']'
    end if
  end subroutine write_real_array

  subroutine write_integer_array(unit, name, values, trailing_comma)
    integer, intent(in) :: unit
    character(*), intent(in) :: name
    integer, intent(in) :: values(:)
    logical, intent(in) :: trailing_comma
    integer :: j

    write (unit, '(3a)', advance='no') '    "', trim(name), '": ['
    do j = 1, size(values)
      write (unit, '(i0)', advance='no') values(j)
      if (j < size(values)) write (unit, '(a)', advance='no') ','
    end do
    if (trailing_comma) then
      write (unit, '(a)') '],'
    else
      write (unit, '(a)') ']'
    end if
  end subroutine write_integer_array

  subroutine write_ring_offsets(unit, nlon)
    integer, intent(in) :: unit
    integer, intent(in) :: nlon(:)
    integer :: j, offset

    write (unit, '(a)', advance='no') '    "ring_offsets": [0'
    offset = 0
    do j = 1, size(nlon)
      offset = offset + nlon(j)
      write (unit, '(a,i0)', advance='no') ',', offset
    end do
    write (unit, '(a)') ']'
  end subroutine write_ring_offsets

  subroutine check_finite(name, nlon, field)
    character(*), intent(in) :: name
    integer, intent(in) :: nlon(:)
    real(real64), intent(in) :: field(:, :)
    integer :: j

    do j = 1, size(nlon)
      if (.not. all(ieee_is_finite(field(1:nlon(j), j)))) then
        error stop 'non-finite value in output field '//name
      end if
    end do
  end subroutine check_finite

  subroutine check_spectral_finite(name, field)
    character(*), intent(in) :: name
    complex(real64), intent(in) :: field(0:, 0:)

    if (.not. all(ieee_is_finite(real(field, real64))) .or. &
        .not. all(ieee_is_finite(aimag(field)))) then
      error stop 'non-finite value in output field '//name
    end if
  end subroutine check_spectral_finite

  subroutine check_rectangular_finite(name, field)
    character(*), intent(in) :: name
    real(real64), intent(in) :: field(:, :)

    if (.not. all(ieee_is_finite(field))) then
      error stop 'non-finite value in output field '//name
    end if
  end subroutine check_rectangular_finite

end module field_output
