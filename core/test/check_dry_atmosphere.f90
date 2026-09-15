program check_dry_atmosphere
  use iso_fortran_env, only: real64
  use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
  use harmonics, only: harmonic_transform
  use dry_vertical_coordinate, only: hybrid_sigma_coordinate, reference_surface_pressure, dry_air_kappa
  use dry_gravity_wave, only: dry_gravity_wave_solver, dry_gravity_wave_implicitness
  use dry_atmosphere, only: dry_atmosphere_solver
  use dry_held_suarez, only: held_suarez_forcing, held_suarez_initial_temperature, &
                              held_suarez_temperature_perturbation, held_suarez_sigma_boundary, &
                              held_suarez_upper_thermal_rate, held_suarez_lower_thermal_rate, &
                              held_suarez_friction_rate, held_suarez_minimum_equilibrium_temperature, &
                              held_suarez_equatorial_temperature, held_suarez_equator_to_pole_difference, &
                              held_suarez_vertical_difference
  implicit none

  integer, parameter :: truncation = 5
  real(real64), parameter :: time_step = 900.0_real64

  call check_reference_atmosphere()
  call check_precomputed_gravity_wave_inverse()
  call check_resting_atmosphere()
  call check_held_suarez_forcing()
  call check_held_suarez_state()
  call check_jablonowski_state()
  call check_jablonowski_steady_state()

contains

  subroutine check_reference_atmosphere()
    type(hybrid_sigma_coordinate) :: coordinate

    call coordinate%init_default()
    if (coordinate%number_of_levels /= 10) error stop 'default dry atmosphere does not have ten levels'
    if (abs(coordinate%reference_p_half(0) - 1000.0_real64) > 1.0e-12_real64) then
      error stop 'dry atmosphere has the wrong top pressure'
    end if
    if (abs(coordinate%reference_p_half(10) - reference_surface_pressure) > 1.0e-12_real64) then
      error stop 'dry atmosphere has the wrong reference surface pressure'
    end if
    if (any(coordinate%reference_delta_p <= 0.0_real64) .or. &
        any(coordinate%reference_temperature <= 0.0_real64)) then
      error stop 'dry reference atmosphere contains a nonphysical value'
    end if
    if (abs(coordinate%full_level_eta(1) - 0.03_real64) > 1.0e-14_real64 .or. &
        abs(coordinate%full_level_eta(10) - 0.94_real64) > 1.0e-14_real64) then
      error stop 'dry full-level eta values do not match the hybrid-sigma pressures'
    end if
  end subroutine check_reference_atmosphere

  subroutine check_precomputed_gravity_wave_inverse()
    type(hybrid_sigma_coordinate) :: coordinate
    type(dry_gravity_wave_solver) :: gravity_wave
    complex(real64), allocatable :: surface_pressure(:, :), delta(:, :, :), temperature(:, :, :)
    complex(real64), allocatable :: rhs_surface_pressure(:, :), rhs_delta(:, :, :), rhs_temperature(:, :, :)
    complex(real64), allocatable :: next_surface_pressure(:, :), next_delta(:, :, :), next_temperature(:, :, :)
    complex(real64), allocatable :: g_next_surface_pressure(:, :), g_next_delta(:, :, :), g_next_temperature(:, :, :)
    complex(real64), allocatable :: residual_surface_pressure(:, :), residual_delta(:, :, :)
    complex(real64), allocatable :: residual_temperature(:, :, :)
    real(real64) :: centered_interval

    call coordinate%init_default()
    call gravity_wave%init(truncation, time_step, coordinate)
    call allocate_zero_state(coordinate%number_of_levels, surface_pressure, delta, temperature)
    delta(3, 2, 4) = cmplx(2.0e-5_real64, -1.0e-5_real64, real64)
    temperature(3, 2, 7) = cmplx(0.7_real64, 0.2_real64, real64)
    surface_pressure(3, 2) = cmplx(1.0e-4_real64, -2.0e-4_real64, real64)
    call gravity_wave%apply(surface_pressure, delta, temperature, &
                            rhs_surface_pressure, rhs_delta, rhs_temperature)
    centered_interval = 2.0_real64*time_step
    call gravity_wave%solve(centered_interval, surface_pressure, delta, temperature, &
                            surface_pressure, delta, temperature, &
                            rhs_surface_pressure, rhs_delta, rhs_temperature, &
                            next_surface_pressure, next_delta, next_temperature)
    call gravity_wave%apply(next_surface_pressure, next_delta, next_temperature, &
                            g_next_surface_pressure, g_next_delta, g_next_temperature)
    residual_surface_pressure = next_surface_pressure - &
      centered_interval*dry_gravity_wave_implicitness*g_next_surface_pressure - &
      (surface_pressure + centered_interval*(1.0_real64 - dry_gravity_wave_implicitness)*rhs_surface_pressure)
    residual_delta = next_delta - centered_interval*dry_gravity_wave_implicitness*g_next_delta - &
      (delta + centered_interval*(1.0_real64 - dry_gravity_wave_implicitness)*rhs_delta)
    residual_temperature = next_temperature - centered_interval*dry_gravity_wave_implicitness*g_next_temperature - &
      (temperature + centered_interval*(1.0_real64 - dry_gravity_wave_implicitness)*rhs_temperature)
    if (maxval(abs(residual_surface_pressure)) > 1.0e-10_real64 .or. &
        maxval(abs(residual_delta)) > 1.0e-10_real64 .or. &
        maxval(abs(residual_temperature)) > 1.0e-10_real64) then
      error stop 'precomputed dry gravity-wave inverse does not solve its implicit system'
    end if
  end subroutine check_precomputed_gravity_wave_inverse

  subroutine check_resting_atmosphere()
    type(hybrid_sigma_coordinate) :: coordinate
    type(harmonic_transform) :: transform
    type(dry_atmosphere_solver) :: solver
    complex(real64), allocatable :: zeta(:, :, :), delta(:, :, :), temperature(:, :, :), log_ps(:, :)
    complex(real64), allocatable :: initial_temperature(:, :, :), initial_log_ps(:, :)
    complex(real64), allocatable :: level_spectral(:, :)
    real(real64), allocatable :: grid(:, :)
    integer :: k

    call coordinate%init_default()
    call transform%init(truncation)
    call allocate_zero_state(coordinate%number_of_levels, log_ps, delta, temperature)
    allocate (zeta, mold=delta)
    zeta = 0.0_real64
    call transform%allocate_field(grid)
    do k = 1, coordinate%number_of_levels
      grid = coordinate%reference_temperature(k)
      call transform%grid_to_spectral(grid, level_spectral)
      temperature(:, :, k) = level_spectral
    end do
    grid = log(reference_surface_pressure)
    call transform%grid_to_spectral(grid, log_ps)
    initial_temperature = temperature
    initial_log_ps = log_ps
    call solver%init(truncation, time_step)
    call solver%set_initial_state(zeta, delta, temperature, log_ps)
    call solver%advance()
    call solver%advance()
    call solver%get_spectral_state(zeta, delta, temperature, log_ps)
    if (maxval(abs(zeta)) > 1.0e-11_real64 .or. maxval(abs(delta)) > 1.0e-11_real64) then
      error stop 'resting dry reference atmosphere developed wind'
    end if
    if (maxval(abs(temperature - initial_temperature)) > 1.0e-8_real64 .or. &
        maxval(abs(log_ps - initial_log_ps)) > 1.0e-10_real64) then
      error stop 'resting dry reference atmosphere changed'
    end if
  end subroutine check_resting_atmosphere

  subroutine check_held_suarez_forcing()
    real(real64) :: forcing_u, forcing_v, temperature_tendency
    real(real64) :: sinphi, pressure, surface_pressure, temperature, u, v
    real(real64) :: pressure_ratio, cosphi_squared, sigma_weight, equilibrium_temperature
    real(real64) :: thermal_rate

    call held_suarez_forcing(0.0_real64, reference_surface_pressure, reference_surface_pressure, &
                             323.0_real64, 12.0_real64, -6.0_real64, &
                             forcing_u, forcing_v, temperature_tendency)
    if (abs(forcing_u + 12.0_real64*held_suarez_friction_rate) > 1.0e-15_real64 .or. &
        abs(forcing_v - 6.0_real64*held_suarez_friction_rate) > 1.0e-15_real64 .or. &
        abs(temperature_tendency + 8.0_real64*held_suarez_lower_thermal_rate) > 1.0e-15_real64) then
      error stop 'Held-Suarez lower-boundary forcing is incorrect'
    end if

    sinphi = 0.5_real64
    pressure = 0.5_real64*reference_surface_pressure
    surface_pressure = reference_surface_pressure
    temperature = 250.0_real64
    u = 12.0_real64
    v = -6.0_real64
    call held_suarez_forcing(sinphi, pressure, surface_pressure, temperature, u, v, &
                             forcing_u, forcing_v, temperature_tendency)
    pressure_ratio = pressure/reference_surface_pressure
    cosphi_squared = 1.0_real64 - sinphi**2
    sigma_weight = max(0.0_real64, (pressure/surface_pressure - held_suarez_sigma_boundary)/ &
                       (1.0_real64 - held_suarez_sigma_boundary))
    equilibrium_temperature = max(held_suarez_minimum_equilibrium_temperature, &
      (held_suarez_equatorial_temperature - held_suarez_equator_to_pole_difference*sinphi**2 - &
       held_suarez_vertical_difference*log(pressure_ratio)*cosphi_squared)*pressure_ratio**dry_air_kappa)
    thermal_rate = held_suarez_upper_thermal_rate + &
      (held_suarez_lower_thermal_rate - held_suarez_upper_thermal_rate)*sigma_weight*cosphi_squared**2
    if (abs(forcing_u) > 1.0e-15_real64 .or. abs(forcing_v) > 1.0e-15_real64 .or. &
        abs(temperature_tendency + thermal_rate*(temperature - equilibrium_temperature)) > 1.0e-15_real64) then
      error stop 'Held-Suarez free-atmosphere forcing is incorrect'
    end if
  end subroutine check_held_suarez_forcing

  subroutine check_held_suarez_state()
    type(harmonic_transform) :: transform
    type(dry_atmosphere_solver) :: forced_solver, unforced_solver
    complex(real64), allocatable :: zeta_spectral(:, :, :), delta_spectral(:, :, :)
    complex(real64), allocatable :: temperature_spectral(:, :, :), log_ps_spectral(:, :)
    complex(real64), allocatable :: unforced_temperature(:, :, :)
    real(real64), allocatable :: zeta(:, :, :), delta(:, :, :), temperature(:, :, :)
    real(real64), allocatable :: surface_pressure(:, :), u(:, :, :), v(:, :, :)
    integer, allocatable :: nlon(:)
    integer :: i, j, k
    real(real64) :: longitude, cosphi, expected_temperature

    call transform%init(truncation)
    nlon = transform%get_nlon()
    call forced_solver%init(truncation, time_step)
    call forced_solver%set_held_suarez_state()
    call forced_solver%get_fields(zeta, delta, temperature, surface_pressure, u, v)
    do k = 1, size(temperature, 3)
      do j = 1, size(nlon)
        cosphi = sqrt(max(0.0_real64, 1.0_real64 - transform%mu(j)**2))
        do i = 1, nlon(j)
          longitude = 2.0_real64*acos(-1.0_real64)*real(i - 1, real64)/real(nlon(j), real64)
          expected_temperature = held_suarez_initial_temperature + &
            held_suarez_temperature_perturbation*cosphi*sin(longitude)
          if (abs(temperature(i, j, k) - expected_temperature) > 1.0e-10_real64) then
            error stop 'Held-Suarez initial temperature is incorrect'
          end if
          if (abs(zeta(i, j, k)) > 1.0e-12_real64 .or. abs(delta(i, j, k)) > 1.0e-12_real64 .or. &
              abs(u(i, j, k)) > 1.0e-10_real64 .or. abs(v(i, j, k)) > 1.0e-10_real64 .or. &
              abs(surface_pressure(i, j) - reference_surface_pressure) > 1.0e-8_real64) then
            error stop 'Held-Suarez atmosphere does not start at rest with uniform surface pressure'
          end if
        end do
      end do
    end do

    call forced_solver%get_spectral_state(zeta_spectral, delta_spectral, temperature_spectral, log_ps_spectral)
    call unforced_solver%init(truncation, time_step)
    call unforced_solver%set_initial_state(zeta_spectral, delta_spectral, temperature_spectral, log_ps_spectral)
    call forced_solver%advance()
    call unforced_solver%advance()
    call forced_solver%get_spectral_state(zeta_spectral, delta_spectral, temperature_spectral, log_ps_spectral)
    call unforced_solver%get_spectral_state(zeta_spectral, delta_spectral, unforced_temperature, log_ps_spectral)
    if (maxval(abs(temperature_spectral - unforced_temperature)) <= 1.0e-8_real64) then
      error stop 'Held-Suarez thermal forcing was not enabled by its initial-state setup'
    end if
    call forced_solver%get_fields(zeta, delta, temperature, surface_pressure, u, v)
    if (.not. all(ieee_is_finite(temperature)) .or. .not. all(ieee_is_finite(surface_pressure)) .or. &
        .not. all(ieee_is_finite(u)) .or. .not. all(ieee_is_finite(v))) then
      error stop 'Held-Suarez state produced a non-finite value'
    end if
  end subroutine check_held_suarez_state

  subroutine check_jablonowski_state()
    type(dry_atmosphere_solver) :: solver
    real(real64), allocatable :: zeta(:, :, :), delta(:, :, :), temperature(:, :, :)
    real(real64), allocatable :: surface_pressure(:, :), u(:, :, :), v(:, :, :)
    real(real64) :: cfl
    integer :: step

    call solver%init(truncation, time_step)
    call solver%set_jablonowski_williamson_state(.true.)
    do step = 1, 4
      call solver%advance()
    end do
    call solver%get_fields(zeta, delta, temperature, surface_pressure, u, v, cfl)
    if (.not. all(ieee_is_finite(zeta)) .or. .not. all(ieee_is_finite(delta)) .or. &
        .not. all(ieee_is_finite(temperature)) .or. .not. all(ieee_is_finite(surface_pressure)) .or. &
        .not. all(ieee_is_finite(u)) .or. .not. all(ieee_is_finite(v)) .or. &
        .not. ieee_is_finite(cfl)) then
      error stop 'Jablonowski-Williamson dry state produced a non-finite value'
    end if
    if (maxval(abs(u)) <= 1.0_real64) error stop 'Jablonowski-Williamson jet is absent'
  end subroutine check_jablonowski_state

  subroutine check_jablonowski_steady_state()
    ! Without the wind perturbation the Jablonowski-Williamson base state is a steady
    ! solution.  It stays balanced only when eta is consistent with the hybrid-sigma
    ! pressures and the surface geopotential is supplied; otherwise |v| grows to
    ! several m/s and ps oscillates by several hPa within a day.
    integer, parameter :: steady_truncation = 21
    integer, parameter :: steady_steps = 96
    real(real64), parameter :: wind_tolerance = 1.0_real64
    real(real64), parameter :: surface_pressure_tolerance = 100.0_real64
    type(harmonic_transform) :: transform
    type(dry_atmosphere_solver) :: solver
    real(real64), allocatable :: zeta(:, :, :), delta(:, :, :), temperature(:, :, :)
    real(real64), allocatable :: surface_pressure(:, :), u(:, :, :), v(:, :, :), initial_u(:, :, :)
    integer, allocatable :: nlon(:)
    real(real64) :: maximum_v, maximum_u_change, maximum_surface_pressure_change
    integer :: step, j

    call transform%init(steady_truncation)
    nlon = transform%get_nlon()
    call solver%init(steady_truncation, time_step)
    call solver%set_jablonowski_williamson_state(.false.)
    call solver%get_fields(zeta, delta, temperature, surface_pressure, initial_u, v)
    do step = 1, steady_steps
      call solver%advance()
    end do
    call solver%get_fields(zeta, delta, temperature, surface_pressure, u, v)

    maximum_v = 0.0_real64
    maximum_u_change = 0.0_real64
    maximum_surface_pressure_change = 0.0_real64
    do j = 1, size(nlon)
      maximum_v = max(maximum_v, maxval(abs(v(1:nlon(j), j, :))))
      maximum_u_change = max(maximum_u_change, &
        maxval(abs(u(1:nlon(j), j, :) - initial_u(1:nlon(j), j, :))))
      maximum_surface_pressure_change = max(maximum_surface_pressure_change, &
        maxval(abs(surface_pressure(1:nlon(j), j) - reference_surface_pressure)))
    end do
    if (maximum_v > wind_tolerance .or. maximum_u_change > wind_tolerance .or. &
        maximum_surface_pressure_change > surface_pressure_tolerance) then
      write (*, '(a,3es12.4)') 'max |v|, max |du|, max |dps| = ', maximum_v, maximum_u_change, &
        maximum_surface_pressure_change
      error stop 'Jablonowski-Williamson base state is not steady'
    end if
  end subroutine check_jablonowski_steady_state

  subroutine allocate_zero_state(levels, surface_pressure, delta, temperature)
    integer, intent(in) :: levels
    complex(real64), allocatable, intent(out) :: surface_pressure(:, :), delta(:, :, :), temperature(:, :, :)

    allocate (surface_pressure(0:truncation + 1, 0:truncation))
    allocate (delta(0:truncation + 1, 0:truncation, levels))
    allocate (temperature(0:truncation + 1, 0:truncation, levels))
    surface_pressure = 0.0_real64
    delta = 0.0_real64
    temperature = 0.0_real64
  end subroutine allocate_zero_state

end program check_dry_atmosphere
