!> Unit checks for the T31 -> T63 restart of the land--sea cases
!> (docs/cases/land-sea.md#t63-の初期値, core/src/cases/radiation/land_sea_restart.f90).
program check_land_sea_restart
  use iso_fortran_env, only: real64
  use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
  use harmonics, only: harmonic_transform
  use topography, only: topography_config, topography_diagnostics, generate_topography
  use dry_physics_config, only: dry_model_physics_config, radiation_orbital_period
  use dry_case_initial_conditions, only: land_sea_case_physics, radiation_case_planet, set_land_sea_case_state
  use dry_atmosphere, only: dry_atmosphere_solver
  use dry_vertical_coordinate, only: dry_air_gas_constant
  use numerics_config, only: model_numerics_config
  use planet_parameters, only: planet_config, earth_gravity
  use field_binary_writer, only: write_field, write_spectral_field
  use filesystem, only: make_directory
  use land_sea_restart, only: land_sea_restart_summary, set_land_sea_state_from_snapshot, read_ring_field, &
                              read_spectral_snapshot_field, read_snapshot_time, pad_spectral_field, &
                              regrid_ring_field, terrain_log_surface_pressure_change
  implicit none

  real(real64), parameter :: pi = acos(-1.0_real64)
  character(len=*), parameter :: scratch = 'build/check_land_sea_restart'
  type(harmonic_transform) :: t31, t63

  call t31%init(31)
  call t63%init(63)
  call make_directory(scratch)
  call check_spectral_normalization()
  call check_regridding()
  call check_terrain_adjustment()
  call check_file_round_trip()
  call check_restart_from_snapshot()
  write (*, '(a)') 'check_land_sea_restart: all checks passed'

contains

  !> A band-limited field analysed at T31 and at T63 has the same coefficients, so
  !> zero padding reproduces the T31 field exactly on the T63 grid.
  subroutine check_spectral_normalization()
    real(real64), allocatable :: grid31(:, :), grid63(:, :), back(:, :)
    complex(real64), allocatable :: a31(:, :), a63(:, :), padded(:, :)

    call smooth_field(t31, grid31)
    call smooth_field(t63, grid63)
    call t31%grid_to_spectral(grid31, a31)
    call t63%grid_to_spectral(grid63, a63)
    call pad_spectral_field(a31, 63, padded)
    if (maxval(abs(padded(0:63, :) - a63(0:63, :))) > 1.0e-12_real64) then
      error stop 'spectral coefficients depend on the truncation'
    end if
    call t63%spectral_to_grid(padded, back)
    if (maxval(abs(back - grid63)) > 1.0e-12_real64) error stop 'padded T31 field differs on the T63 grid'
  end subroutine check_spectral_normalization

  subroutine check_regridding()
    real(real64), allocatable :: source(:, :), expected(:, :), regridded(:, :)
    logical, allocatable :: valid(:, :), needed(:, :)
    integer, allocatable :: nlon31(:), nlon63(:)
    integer :: i, j, fallback
    real(real64) :: longitude, error, polar_error

    nlon31 = t31%get_nlon()
    nlon63 = t63%get_nlon()
    ! A constant is reproduced exactly and a smooth field to second order.
    call t31%allocate_field(source)
    source = 287.5_real64
    call t63%allocate_field(regridded)
    call regrid_ring_field(t31%mu, nlon31, source, t63%mu, nlon63, regridded)
    do j = 1, size(nlon63)
      if (any(abs(regridded(1:nlon63(j), j) - 287.5_real64) > 1.0e-12_real64)) &
        error stop 'regridding does not preserve a constant'
    end do
    call smooth_field(t31, source)
    call smooth_field(t63, expected)
    call regrid_ring_field(t31%mu, nlon31, source, t63%mu, nlon63, regridded)
    ! Between the outermost T31 rings the error is second order (largest where the
    ! octahedral rings are short, 8 degrees at 69S); poleward of them
    ! the nearest ring is used, which is first order in the ring spacing.
    error = 0.0_real64
    polar_error = 0.0_real64
    do j = 1, size(nlon63)
      if (abs(t63%mu(j)) > abs(t31%mu(1))) then
        polar_error = max(polar_error, maxval(abs(regridded(1:nlon63(j), j) - expected(1:nlon63(j), j))))
      else
        error = max(error, maxval(abs(regridded(1:nlon63(j), j) - expected(1:nlon63(j), j))))
      end if
    end do
    if (error > 2.0e-3_real64 .or. polar_error > 2.5e-2_real64) error stop 'bilinear regridding error is too large'

    ! Masked: only the western hemisphere is valid.  Every value comes from valid
    ! points; eastern points far from the mask fall back to the nearest one.
    allocate (valid(size(source, 1), size(source, 2)))
    valid = .false.
    source = 1.0e6_real64
    do j = 1, size(nlon31)
      do i = 1, nlon31(j)
        longitude = 2.0_real64*pi*real(i - 1, real64)/real(nlon31(j), real64)
        if (longitude < pi) then
          valid(i, j) = .true.
          source(i, j) = 1.0_real64
        end if
      end do
    end do
    regridded = -1.0_real64
    call regrid_ring_field(t31%mu, nlon31, source, t63%mu, nlon63, regridded, valid, fallback_points=fallback)
    do j = 1, size(nlon63)
      if (any(regridded(1:nlon63(j), j) /= 1.0_real64)) error stop 'masked regridding used an invalid point'
    end do
    if (fallback == 0) error stop 'masked regridding never fell back to the nearest valid point'
    ! Points outside `needed` keep the unmasked interpolation instead.
    allocate (needed(size(regridded, 1), size(regridded, 2)))
    needed = .false.
    call regrid_ring_field(t31%mu, nlon31, source, t63%mu, nlon63, regridded, valid, needed, fallback)
    if (fallback /= 0) error stop 'masked regridding fell back outside the needed points'
    if (maxval(regridded(1:nlon63(64), 64)) /= 1.0e6_real64) error stop 'unneeded points were masked'
  end subroutine check_regridding

  !> Raising the terrain by 100 m lowers ln p_s by g*100/(R_d T_v).
  subroutine check_terrain_adjustment()
    complex(real64), allocatable :: phi31(:, :), phi63(:, :), temperature(:, :), humidity(:, :), change(:, :)
    real(real64), allocatable :: grid(:, :)
    integer, allocatable :: nlon63(:)
    real(real64) :: maximum_change, expected

    allocate (phi31(0:32, 0:31), phi63(0:64, 0:63), temperature(0:64, 0:63), humidity(0:64, 0:63))
    phi31 = 0.0_real64
    phi63 = 0.0_real64
    temperature = 0.0_real64
    humidity = 0.0_real64
    phi31(0, 0) = 500.0_real64*earth_gravity
    phi63(0, 0) = 600.0_real64*earth_gravity
    temperature(0, 0) = 280.0_real64
    nlon63 = t63%get_nlon()
    call terrain_log_surface_pressure_change(t63, nlon63, phi31, phi63, temperature, humidity, change, &
                                             maximum_change)
    expected = -100.0_real64*earth_gravity/(dry_air_gas_constant*280.0_real64)
    call t63%spectral_to_grid(change, grid)
    if (abs(maximum_change - abs(expected)) > 1.0e-12_real64 .or. &
        maxval(abs(grid(1:nlon63(1), 1) - expected)) > 1.0e-12_real64) then
      error stop 'terrain log surface pressure change is wrong'
    end if
    phi63(0:32, 0:31) = phi31
    call terrain_log_surface_pressure_change(t63, nlon63, phi31, phi63, temperature, humidity, change, &
                                             maximum_change)
    if (maximum_change /= 0.0_real64) error stop 'identical terrain changed ln p_s'
  end subroutine check_terrain_adjustment

  subroutine check_file_round_trip()
    real(real64), allocatable :: field(:, :), back(:, :)
    complex(real64), allocatable :: spectrum(:, :), spectral_back(:, :)
    integer, allocatable :: nlon31(:)
    integer :: unit, j

    nlon31 = t31%get_nlon()
    call smooth_field(t31, field)
    call write_field(scratch//'/ring.bin', nlon31, field)
    call t31%allocate_field(back)
    call read_ring_field(scratch//'/ring.bin', nlon31, back)
    do j = 1, size(nlon31)
      if (any(back(1:nlon31(j), j) /= field(1:nlon31(j), j))) error stop 'ring field round trip changed values'
    end do
    call t31%grid_to_spectral(field, spectrum)
    spectrum(32, :) = 0.0_real64
    call write_spectral_field(scratch//'/spectrum.bin', spectrum)
    call read_spectral_snapshot_field(scratch//'/spectrum.bin', 31, spectral_back)
    if (any(spectral_back /= spectrum)) error stop 'spectral field round trip changed values'
    open (newunit=unit, file=scratch//'/time.json', status='replace', action='write')
    write (unit, '(a)') '{'
    write (unit, '(a)') '  "time_seconds":  1.5552000000000000E+008,'
    write (unit, '(a)') '  "step": 129600'
    write (unit, '(a)') '}'
    close (unit)
    if (read_snapshot_time(scratch//'/time.json') /= 1.5552e8_real64) error stop 'snapshot time was misread'
  end subroutine check_file_round_trip

  !> Writes a T31 land--sea state (with sea ice added by hand) in the snapshot
  !> format, starts a T63 solver from it and integrates a few steps.
  subroutine check_restart_from_snapshot()
    type(topography_config) :: terrain
    type(topography_diagnostics) :: terrain_diagnostics
    type(dry_model_physics_config) :: physics
    type(planet_config) :: planet
    type(model_numerics_config) :: numerics
    type(dry_atmosphere_solver) :: solver31, solver63
    type(land_sea_restart_summary) :: summary
    real(real64), allocatable :: land31(:, :), analytic31(:, :), height31(:, :)
    real(real64), allocatable :: land63(:, :), analytic63(:, :), height63(:, :)
    real(real64), allocatable :: zeta_grid(:, :, :), delta_grid(:, :, :), temperature_grid(:, :, :)
    real(real64), allocatable :: surface_pressure(:, :), u(:, :, :), v(:, :, :), humidity_grid(:, :, :)
    real(real64), allocatable :: surface_temperature(:, :), deep_temperature(:, :), surface_water(:, :)
    real(real64), allocatable :: land_temperature(:, :), ocean_temperature(:, :), ice_fraction(:, :)
    real(real64), allocatable :: ice_volume(:, :), snow_water(:, :)
    complex(real64), allocatable :: phi31(:, :), phi63(:, :)
    complex(real64), allocatable :: zeta31(:, :, :), delta31(:, :, :), temperature31(:, :, :), log_ps31(:, :)
    complex(real64), allocatable :: humidity31(:, :, :)
    complex(real64), allocatable :: zeta63(:, :, :), delta63(:, :, :), temperature63(:, :, :), log_ps63(:, :)
    complex(real64), allocatable :: humidity63(:, :, :)
    integer, allocatable :: nlon31(:), nlon63(:)
    character(len=:), allocatable :: directory
    character(len=2) :: level_text
    integer :: k, unit, step, i, j
    real(real64) :: period

    directory = scratch//'/moist_land_sea_t31'
    call make_directory(directory)
    nlon31 = t31%get_nlon()
    nlon63 = t63%get_nlon()
    physics = land_sea_case_physics()
    planet = radiation_case_planet(physics)
    period = radiation_orbital_period(physics%radiation)
    call generate_topography(t31, terrain, land31, analytic31, phi31, height31, terrain_diagnostics)
    call generate_topography(t63, terrain, land63, analytic63, phi63, height63, terrain_diagnostics)

    numerics = model_numerics_config()
    numerics%truncation = 31
    numerics%time_step = 1200.0_real64
    call solver31%init_with_config(numerics)
    call set_land_sea_case_state(solver31, t31, physics, planet, phi31, land31)
    call solver31%get_spectral_state(zeta31, delta31, temperature31, log_ps31, specific_humidity=humidity31)
    call solver31%get_fields(zeta_grid, delta_grid, temperature_grid, surface_pressure, u, v, &
      surface_temperature=surface_temperature, deep_temperature=deep_temperature, surface_water=surface_water, &
      land_temperature=land_temperature, ocean_temperature=ocean_temperature, sea_ice_fraction=ice_fraction, &
      sea_ice_volume=ice_volume, snow_water=snow_water)
    ! Sea ice in equilibrium poleward of 60 degrees; snow on the land there.
    do j = 1, size(nlon31)
      if (abs(t31%mu(j)) < sin(pi/3.0_real64)) cycle
      do i = 1, nlon31(j)
        if (land31(i, j) < 1.0_real64) then
          ocean_temperature(i, j) = physics%sea_ice%freezing_temperature
          ice_fraction(i, j) = 0.8_real64
          ice_volume(i, j) = 1.2_real64
        end if
        if (land31(i, j) > 0.0_real64) snow_water(i, j) = 200.0_real64
      end do
    end do

    call write_field(directory//'/land_fraction.bin', nlon31, land31)
    do k = 1, size(zeta31, 3)
      write (level_text, '(i2.2)') k
      call write_spectral_field(directory//'/yearly_zeta_spectral_y0006_l'//level_text//'.bin', zeta31(:, :, k))
      call write_spectral_field(directory//'/yearly_delta_spectral_y0006_l'//level_text//'.bin', delta31(:, :, k))
      call write_spectral_field(directory//'/yearly_temperature_spectral_y0006_l'//level_text//'.bin', &
                                temperature31(:, :, k))
      call write_spectral_field(directory//'/yearly_specific_humidity_spectral_y0006_l'//level_text//'.bin', &
                                humidity31(:, :, k))
    end do
    call write_spectral_field(directory//'/yearly_log_surface_pressure_spectral_y0006.bin', log_ps31)
    call write_field(directory//'/yearly_surface_temperature_y0006.bin', nlon31, surface_temperature)
    call write_field(directory//'/yearly_deep_temperature_y0006.bin', nlon31, deep_temperature)
    call write_field(directory//'/yearly_surface_water_y0006.bin', nlon31, surface_water)
    call write_field(directory//'/yearly_land_temperature_y0006.bin', nlon31, land_temperature)
    call write_field(directory//'/yearly_ocean_temperature_y0006.bin', nlon31, ocean_temperature)
    call write_field(directory//'/yearly_sea_ice_fraction_y0006.bin', nlon31, ice_fraction)
    call write_field(directory//'/yearly_sea_ice_volume_y0006.bin', nlon31, ice_volume)
    call write_field(directory//'/yearly_snow_water_y0006.bin', nlon31, snow_water)
    open (newunit=unit, file=directory//'/yearly_time_y0006.json', status='replace', action='write')
    write (unit, '(a,es24.16e3,a)') '{"time_seconds": ', 5.0_real64*period, ', "step": 0}'
    close (unit)

    numerics%truncation = 63
    call solver63%init_with_config(numerics)
    call set_land_sea_state_from_snapshot(solver63, t63, physics, planet, phi63, land63, t31, phi31, land31, &
                                          directory, 6, period, summary)
    if (summary%source_truncation /= 31) error stop 'restart summary has the wrong source truncation'
    if (summary%ice_points == 0) error stop 'restart lost the sea ice'
    if (summary%maximum_log_surface_pressure_change <= 0.0_real64 .or. &
        summary%maximum_log_surface_pressure_change > 0.2_real64) then
      error stop 'restart terrain adjustment of ln p_s is implausible'
    end if

    ! The dynamics above the terrain adjustment are the T31 coefficients.
    call solver63%get_spectral_state(zeta63, delta63, temperature63, log_ps63, specific_humidity=humidity63)
    if (any(zeta63(0:31, 0:31, :) /= zeta31(0:31, 0:31, :)) .or. any(zeta63(32:, :, :) /= 0.0_real64) .or. &
        any(temperature63(0:31, 0:31, :) /= temperature31(0:31, 0:31, :)) .or. &
        any(humidity63(0:31, 0:31, :) /= humidity31(0:31, 0:31, :)) .or. &
        any(delta63(0:31, 0:31, :) /= delta31(0:31, 0:31, :))) then
      error stop 'restart did not zero-pad the spectral state'
    end if

    do step = 1, 3
      call solver63%advance()
    end do
    call solver63%get_fields(zeta_grid, delta_grid, temperature_grid, surface_pressure, u, v, &
      surface_temperature=surface_temperature, deep_temperature=deep_temperature, &
      specific_humidity=humidity_grid, surface_water=surface_water, land_temperature=land_temperature, &
      ocean_temperature=ocean_temperature, sea_ice_fraction=ice_fraction, sea_ice_volume=ice_volume, &
      snow_water=snow_water)
    do j = 1, size(nlon63)
      if (.not. all(ieee_is_finite(temperature_grid(1:nlon63(j), j, :))) .or. &
          .not. all(ieee_is_finite(surface_pressure(1:nlon63(j), j))) .or. &
          .not. all(ieee_is_finite(land_temperature(1:nlon63(j), j))) .or. &
          .not. all(ieee_is_finite(ocean_temperature(1:nlon63(j), j))) .or. &
          .not. all(ieee_is_finite(ice_volume(1:nlon63(j), j)))) then
        error stop 'short T63 integration from the T31 snapshot produced a non-finite field'
      end if
      if (any(surface_water(1:nlon63(j), j) < 0.0_real64) .or. &
          any(surface_water(1:nlon63(j), j) > physics%bucket%capacity) .or. &
          any(snow_water(1:nlon63(j), j) < 0.0_real64)) then
        error stop 'restarted land water is out of bounds'
      end if
      do i = 1, nlon63(j)
        if (land63(i, j) >= 1.0_real64 .and. ice_fraction(i, j) /= 0.0_real64) error stop 'restart put ice on land'
        if (land63(i, j) <= 0.0_real64 .and. (surface_water(i, j) /= 0.0_real64 .or. snow_water(i, j) /= 0.0_real64)) &
          error stop 'restart put land water on the open ocean'
      end do
    end do
    if (maxval(ice_fraction) <= 0.0_real64) error stop 'short T63 integration lost the restarted sea ice'
    if (maxval(snow_water) <= 0.0_real64) error stop 'short T63 integration lost the restarted snowpack'
  end subroutine check_restart_from_snapshot

  !> sin(phi) + cos(phi) cos(lambda) + 0.3 cos^2(phi) sin(2 lambda): degree <= 2.
  subroutine smooth_field(transform, field)
    type(harmonic_transform), intent(in) :: transform
    real(real64), allocatable, intent(out) :: field(:, :)
    integer, allocatable :: nlon(:)
    integer :: i, j
    real(real64) :: longitude, latitude

    call transform%allocate_field(field)
    field = 0.0_real64
    nlon = transform%get_nlon()
    do j = 1, size(nlon)
      latitude = asin(transform%mu(j))
      do i = 1, nlon(j)
        longitude = 2.0_real64*pi*real(i - 1, real64)/real(nlon(j), real64)
        field(i, j) = sin(latitude) + cos(latitude)*cos(longitude) + &
                      0.3_real64*cos(latitude)**2*sin(2.0_real64*longitude)
      end do
    end do
  end subroutine smooth_field

end program check_land_sea_restart
