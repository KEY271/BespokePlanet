module field_output
  use iso_fortran_env, only: real64
  use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
  use barotropic_vorticity, only: earth_radius, rotation_rate, &
                                  raw_filter_epsilon, raw_filter_alpha, &
                                  hyperdiffusion_order, &
                                  hyperdiffusion_timescale_seconds
  implicit none
  private

  public :: make_directory
  public :: write_snapshot
  public :: write_run_metadata

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

  subroutine write_run_metadata(case_directory, case_name, initial_condition_json, &
                                truncation, time_step, duration, number_of_steps, &
                                maximum_cfl, elapsed_wall_seconds, nlon, mu)
    character(*), intent(in) :: case_directory, case_name, initial_condition_json
    integer, intent(in) :: truncation, number_of_steps
    real(real64), intent(in) :: time_step, duration, maximum_cfl, elapsed_wall_seconds
    integer, intent(in) :: nlon(:)
    real(real64), intent(in) :: mu(:)
    integer :: unit, point_count

    if (size(nlon) /= size(mu)) error stop 'metadata grid arrays have different sizes'
    point_count = sum(nlon)
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
    write (unit, '(a,i0,a)') '    "number_of_snapshots": ', number_of_steps + 1, ','
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

end module field_output
