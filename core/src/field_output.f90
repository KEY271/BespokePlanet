module field_output
  use iso_fortran_env, only: real64
  use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
  use barotropic_vorticity, only: earth_radius, rotation_rate
  use raw_filter, only: raw_filter_epsilon, raw_filter_alpha
  use spectral_hyperdiffusion, only: hyperdiffusion_order, &
                                       hyperdiffusion_timescale_seconds
  implicit none
  private

  public :: make_directory
  public :: write_snapshot
  public :: write_run_metadata
  public :: write_shallow_water_snapshot
  public :: write_shallow_water_metadata

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

end module field_output
