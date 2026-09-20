!> Output files and metadata schema of the barotropic-vorticity cases.
module barotropic_case_output
  use iso_fortran_env, only: real64
  use field_binary_writer, only: write_field, check_finite
  use json_writer, only: write_real_array, write_integer_array, write_ring_offsets
  use numerics_config, only: model_numerics_config
  implicit none
  private
  public :: write_barotropic_snapshot, write_barotropic_metadata

contains

  subroutine write_barotropic_snapshot(case_directory, step, nlon, zeta, u, v)
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
  end subroutine write_barotropic_snapshot

  !> Metadata of one run.  All physical and numerical values are supplied by the
  !> case runner so that the file records the configuration actually integrated.
  subroutine write_barotropic_metadata(case_directory, case_name, initial_condition_json, &
                                       numerics, duration, number_of_steps, &
                                       snapshot_interval_steps, maximum_cfl, elapsed_wall_seconds, &
                                       nlon, mu, radius, rotation_rate)
    character(*), intent(in) :: case_directory, case_name, initial_condition_json
    type(model_numerics_config), intent(in) :: numerics
    integer, intent(in) :: number_of_steps, snapshot_interval_steps
    real(real64), intent(in) :: duration, maximum_cfl, elapsed_wall_seconds
    integer, intent(in) :: nlon(:)
    real(real64), intent(in) :: mu(:)
    real(real64), intent(in) :: radius, rotation_rate
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
    write (unit, '(a,es24.16e3,a)') '    "time_step_seconds": ', numerics%time_step, ','
    write (unit, '(a,i0,a)') '    "number_of_steps": ', number_of_steps, ','
    write (unit, '(a,i0,a)') '    "snapshot_interval_steps": ', snapshot_interval_steps, ','
    write (unit, '(a,es24.16e3,a)') &
      '    "snapshot_interval_seconds": ', numerics%time_step*real(snapshot_interval_steps, real64), ','
    write (unit, '(a,i0,a)') '    "number_of_snapshots": ', number_of_snapshots, ','
    write (unit, '(a)') '    "snapshot_time_seconds": "step * time_step_seconds"'
    write (unit, '(a)') '  },'
    write (unit, '(a)') '  "physical_constants": {'
    write (unit, '(a,es24.16e3,a)') '    "earth_radius_m": ', radius, ','
    write (unit, '(a,es24.16e3)') '    "rotation_rate_rad_s": ', rotation_rate
    write (unit, '(a)') '  },'
    write (unit, '(a)') '  "numerics": {'
    write (unit, '(a,i0,a)') '    "spectral_truncation": ', numerics%truncation, ','
    write (unit, '(a)') '    "time_integrator": "RAW-filtered leapfrog",'
    write (unit, '(a)') '    "raw_filter": {'
    write (unit, '(a,es24.16e3,a)') '      "epsilon": ', numerics%raw_filter%epsilon, ','
    write (unit, '(a,es24.16e3)') '      "alpha": ', numerics%raw_filter%alpha
    write (unit, '(a)') '    },'
    write (unit, '(a)') '    "hyperdiffusion": {'
    write (unit, '(a,i0,a)') '      "order": ', numerics%hyperdiffusion%order, ','
    write (unit, '(a,es24.16e3)') &
      '      "e_folding_time_at_truncation_seconds": ', numerics%hyperdiffusion%timescale_seconds
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
  end subroutine write_barotropic_metadata

end module barotropic_case_output
