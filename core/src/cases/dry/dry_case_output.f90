!> Output files and metadata schema of the dry-atmosphere cases.
module dry_case_output
  use iso_fortran_env, only: real64
  use field_binary_writer, only: write_field, check_finite
  use json_writer, only: write_inline_real_values, write_real_array, write_integer_array
  use numerics_config, only: dry_hyperdiffusion_config
  implicit none
  private
  public :: write_dry_snapshot, write_dry_metadata, write_dry_reference_atmosphere

contains

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

  subroutine write_dry_metadata(case_directory, initial_condition, truncation, time_step, duration, &
                                number_of_steps, snapshot_interval_steps, maximum_cfl, &
                                elapsed_wall_seconds, nlon, mu, hyperdiffusion, &
                                pressure_half, delta_pressure, layer_l, alpha, reference_temperature, &
                                a_half, b_half, include_ground_temperatures)
    character(*), intent(in) :: case_directory, initial_condition
    integer, intent(in) :: truncation, number_of_steps, snapshot_interval_steps
    real(real64), intent(in) :: time_step, duration, maximum_cfl, elapsed_wall_seconds
    integer, intent(in) :: nlon(:)
    real(real64), intent(in) :: mu(:)
    type(dry_hyperdiffusion_config), intent(in) :: hyperdiffusion
    real(real64), intent(in) :: pressure_half(0:), delta_pressure(:), layer_l(:), alpha(:)
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
    write (unit, '(a)') '  "hyperdiffusion_e_folding_seconds": {"zeta":'// &
      json_seconds(hyperdiffusion%vorticity_timescale_seconds)//',"delta":'// &
      json_seconds(hyperdiffusion%divergence_timescale_seconds)//',"temperature":'// &
      json_seconds(hyperdiffusion%temperature_timescale_seconds)//'},'
    write (unit, '(a,es24.16e3,a)') '  "maximum_advective_cfl": ', maximum_cfl, ','
    write (unit, '(a,es24.16e3,a)') '  "elapsed_wall_seconds": ', elapsed_wall_seconds, ','
    call write_dry_reference_atmosphere(unit, pressure_half, delta_pressure, layer_l, alpha, &
                                        reference_temperature, a_half, b_half)
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

  !> Reference-atmosphere and hybrid-coordinate entries shared by the dry and radiation schemas.
  subroutine write_dry_reference_atmosphere(unit, pressure_half, delta_pressure, layer_l, alpha, &
                                            reference_temperature, a_half, b_half)
    integer, intent(in) :: unit
    real(real64), intent(in) :: pressure_half(0:), delta_pressure(:), layer_l(:), alpha(:)
    real(real64), intent(in) :: reference_temperature(:), a_half(0:), b_half(0:)

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
  end subroutine write_dry_reference_atmosphere

  !> Whole seconds are written as JSON integers (the established schema); other values keep full precision.
  function json_seconds(seconds) result(text)
    real(real64), intent(in) :: seconds
    character(len=:), allocatable :: text
    character(len=32) :: buffer

    if (seconds == aint(seconds) .and. abs(seconds) < 1.0e9_real64) then
      write (buffer, '(i0)') nint(seconds)
    else
      write (buffer, '(es24.16e3)') seconds
    end if
    text = trim(adjustl(buffer))
  end function json_seconds

end module dry_case_output
