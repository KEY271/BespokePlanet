module dry_initial_conditions
  use iso_fortran_env, only: real64
  use harmonics, only: harmonic_transform
  use planet_parameters, only: earth_radius, earth_gravity, rotation_rate => earth_rotation_rate
  use spectral_vector_operators, only: flux_divergence, flux_curl
  use dry_vertical_coordinate, only: hybrid_sigma_coordinate, dry_air_gas_constant, &
                                     reference_surface_pressure, jablonowski_mean_temperature
  implicit none
  private

  real(real64), parameter, public :: jablonowski_maximum_wind = 35.0_real64
  real(real64), parameter, public :: jablonowski_jet_eta = 0.252_real64
  real(real64), parameter, public :: jablonowski_perturbation_wind = 1.0_real64

  public :: jablonowski_williamson_initial_state
  public :: jablonowski_williamson_topographic_initial_state
  public :: topographic_surface_pressure

contains

  subroutine jablonowski_williamson_initial_state(transform, truncation, coordinate, &
                                                  include_perturbation, zeta, delta, temperature, &
                                                  log_surface_pressure, surface_geopotential, &
                                                  planetary_rotation_rate, flat_terrain)
    type(harmonic_transform), intent(inout) :: transform
    integer, intent(in) :: truncation
    type(hybrid_sigma_coordinate), intent(in) :: coordinate
    logical, intent(in), optional :: include_perturbation
    complex(real64), allocatable, intent(out) :: zeta(:, :, :), delta(:, :, :), temperature(:, :, :)
    complex(real64), allocatable, intent(out) :: log_surface_pressure(:, :)
    complex(real64), allocatable, intent(out) :: surface_geopotential(:, :)
    real(real64), intent(in), optional :: planetary_rotation_rate
    ! When true, Phi_s = 0 and the base-state zonal wind is replaced by the wind in
    ! gradient-wind balance with the unchanged temperature over flat terrain.
    logical, intent(in), optional :: flat_terrain
    real(real64), allocatable :: u(:, :), v(:, :), temperature_grid(:, :), log_ps_grid(:, :)
    real(real64), allocatable :: surface_geopotential_grid(:, :)
    complex(real64), allocatable :: level_spectral(:, :)
    integer, allocatable :: nlon(:)
    integer :: i, j, k, number_of_levels
    logical :: perturb, flat
    real(real64) :: pi, longitude, eta, eta_v, vertical_shape, surface_shape
    real(real64) :: sinphi, cosphi, angular_cosine, distance, wind_perturbation
    real(real64) :: temperature_correction, active_rotation_rate, base_wind
    real(real64), parameter :: centre_longitude = 20.0_real64*acos(-1.0_real64)/180.0_real64
    real(real64), parameter :: centre_latitude = 40.0_real64*acos(-1.0_real64)/180.0_real64
    real(real64), parameter :: perturbation_radius = earth_radius/10.0_real64

    perturb = .true.
    if (present(include_perturbation)) perturb = include_perturbation
    active_rotation_rate = rotation_rate
    if (present(planetary_rotation_rate)) active_rotation_rate = planetary_rotation_rate
    flat = .false.
    if (present(flat_terrain)) flat = flat_terrain
    number_of_levels = coordinate%number_of_levels
    if (number_of_levels < 1) error stop 'Jablonowski-Williamson coordinate is not initialized'
    pi = acos(-1.0_real64)
    nlon = transform%get_nlon()
    call transform%allocate_field(u)
    call transform%allocate_field(v)
    call transform%allocate_field(temperature_grid)
    call transform%allocate_field(log_ps_grid)
    call transform%allocate_field(surface_geopotential_grid)
    allocate (zeta(0:truncation + 1, 0:truncation, number_of_levels))
    allocate (delta(0:truncation + 1, 0:truncation, number_of_levels))
    allocate (temperature(0:truncation + 1, 0:truncation, number_of_levels))
    zeta = 0.0_real64
    delta = 0.0_real64
    temperature = 0.0_real64

    ! The surface geopotential balances the nonzero surface wind under a uniform
    ! surface pressure.  It is the eta = 1 value of the balanced geopotential.
    ! With flat_terrain the geopotential is instead zero and the wind is adjusted below.
    surface_geopotential_grid = 0.0_real64
    surface_shape = cos(0.5_real64*pi*(1.0_real64 - jablonowski_jet_eta))**1.5_real64
    if (.not. flat) then
      do j = 1, size(nlon)
        sinphi = transform%mu(j)
        cosphi = sqrt(max(0.0_real64, 1.0_real64 - sinphi*sinphi))
        surface_geopotential_grid(1:nlon(j), j) = jablonowski_maximum_wind*surface_shape*( &
          jablonowski_maximum_wind*surface_shape*latitude_function_a(sinphi, cosphi) + &
          earth_radius*active_rotation_rate*latitude_function_b(sinphi, cosphi))
      end do
    end if
    call transform%grid_to_spectral(surface_geopotential_grid, surface_geopotential)

    do k = 1, number_of_levels
      u = 0.0_real64
      v = 0.0_real64
      temperature_grid = 0.0_real64
      eta = coordinate%full_level_eta(k)
      eta_v = 0.5_real64*pi*(eta - jablonowski_jet_eta)
      vertical_shape = max(0.0_real64, cos(eta_v))
      do j = 1, size(nlon)
        sinphi = transform%mu(j)
        cosphi = sqrt(max(0.0_real64, 1.0_real64 - sinphi*sinphi))
        if (flat) then
          base_wind = flat_terrain_balanced_wind(sinphi, cosphi, vertical_shape**1.5_real64, &
                                                 surface_shape, active_rotation_rate)
        else
          base_wind = jablonowski_maximum_wind*vertical_shape**1.5_real64*double_angle_sine_squared(sinphi, cosphi)
        end if
        temperature_correction = 0.75_real64*eta*pi*jablonowski_maximum_wind/ &
          dry_air_gas_constant*sin(eta_v)*sqrt(vertical_shape)*( &
          2.0_real64*jablonowski_maximum_wind*vertical_shape**1.5_real64* &
          latitude_function_a(sinphi, cosphi) + &
          earth_radius*active_rotation_rate*latitude_function_b(sinphi, cosphi))
        do i = 1, nlon(j)
          longitude = 2.0_real64*pi*real(i - 1, real64)/real(nlon(j), real64)
          wind_perturbation = 0.0_real64
          if (perturb) then
            angular_cosine = sin(centre_latitude)*sinphi + &
              cos(centre_latitude)*cosphi*cos(longitude - centre_longitude)
            distance = earth_radius*acos(max(-1.0_real64, min(1.0_real64, angular_cosine)))
            wind_perturbation = jablonowski_perturbation_wind* &
              exp(-(distance/perturbation_radius)**2)
          end if
          u(i, j) = base_wind + wind_perturbation
          temperature_grid(i, j) = coordinate%reference_temperature(k) + temperature_correction
        end do
      end do
      call flux_curl(transform, u, v, level_spectral)
      zeta(:, :, k) = level_spectral
      call flux_divergence(transform, u, v, level_spectral)
      delta(:, :, k) = level_spectral
      call transform%grid_to_spectral(temperature_grid, level_spectral)
      temperature(:, :, k) = level_spectral
      zeta(0, 0, k) = 0.0_real64
      delta(0, 0, k) = 0.0_real64
    end do

    log_ps_grid = log(reference_surface_pressure)
    call transform%grid_to_spectral(log_ps_grid, log_surface_pressure)
  end subroutine jablonowski_williamson_initial_state

  !> The flat-radiation basic state interpreted as a function of pressure and
  !> placed over an arbitrary fixed surface geopotential.  A zero field is sent
  !> through the legacy routine so flat cases remain bit-for-bit unchanged.
  subroutine jablonowski_williamson_topographic_initial_state(transform, truncation, coordinate, &
      surface_geopotential, planetary_rotation_rate, zeta, delta, temperature, log_surface_pressure)
    type(harmonic_transform), intent(inout) :: transform
    integer, intent(in) :: truncation
    type(hybrid_sigma_coordinate), intent(in) :: coordinate
    complex(real64), intent(in) :: surface_geopotential(0:, 0:)
    real(real64), intent(in) :: planetary_rotation_rate
    complex(real64), allocatable, intent(out) :: zeta(:, :, :), delta(:, :, :), temperature(:, :, :)
    complex(real64), allocatable, intent(out) :: log_surface_pressure(:, :)
    complex(real64), allocatable :: unused_geopotential(:, :), level_spectral(:, :)
    real(real64), allocatable :: surface_geopotential_grid(:, :), log_ps_grid(:, :)
    real(real64), allocatable :: u(:, :), v(:, :), temperature_grid(:, :)
    integer, allocatable :: nlon(:)
    integer :: i, j, k, levels
    real(real64) :: surface_pressure, eta, eta_v, vertical_shape, surface_shape
    real(real64) :: sinphi, cosphi, temperature_correction
    real(real64), parameter :: pi = acos(-1.0_real64)

    if (all(surface_geopotential == cmplx(0.0_real64, 0.0_real64, kind=real64))) then
      call jablonowski_williamson_initial_state(transform, truncation, coordinate, .false., &
        zeta, delta, temperature, log_surface_pressure, unused_geopotential, planetary_rotation_rate, &
        flat_terrain=.true.)
      return
    end if
    levels = coordinate%number_of_levels
    nlon = transform%get_nlon()
    call transform%spectral_to_grid(surface_geopotential, surface_geopotential_grid)
    call transform%allocate_field(log_ps_grid)
    call transform%allocate_field(u)
    call transform%allocate_field(v)
    call transform%allocate_field(temperature_grid)
    log_ps_grid = log(reference_surface_pressure)
    do j = 1, size(nlon)
      do i = 1, nlon(j)
        surface_pressure = topographic_surface_pressure(transform%mu(j), surface_geopotential_grid(i, j), &
                                                        planetary_rotation_rate)
        log_ps_grid(i, j) = log(surface_pressure)
      end do
    end do
    call transform%grid_to_spectral(log_ps_grid, log_surface_pressure)
    allocate (zeta(0:truncation + 1, 0:truncation, levels), delta(0:truncation + 1, 0:truncation, levels))
    allocate (temperature(0:truncation + 1, 0:truncation, levels))
    surface_shape = cos(0.5_real64*pi*(1.0_real64 - jablonowski_jet_eta))**1.5_real64
    do k = 1, levels
      u = 0.0_real64
      v = 0.0_real64
      temperature_grid = 0.0_real64
      do j = 1, size(nlon)
        sinphi = transform%mu(j)
        cosphi = sqrt(max(0.0_real64, 1.0_real64 - sinphi*sinphi))
        do i = 1, nlon(j)
          surface_pressure = exp(log_ps_grid(i, j))
          eta = 0.5_real64*(coordinate%a_half(k - 1) + coordinate%b_half(k - 1)*surface_pressure + &
            coordinate%a_half(k) + coordinate%b_half(k)*surface_pressure)/reference_surface_pressure
          eta_v = 0.5_real64*pi*(eta - jablonowski_jet_eta)
          vertical_shape = max(0.0_real64, cos(eta_v))
          u(i, j) = flat_terrain_balanced_wind(sinphi, cosphi, vertical_shape**1.5_real64, &
                                               surface_shape, planetary_rotation_rate)
          temperature_correction = 0.75_real64*eta*pi*jablonowski_maximum_wind/ &
            dry_air_gas_constant*sin(eta_v)*sqrt(vertical_shape)*( &
            2.0_real64*jablonowski_maximum_wind*vertical_shape**1.5_real64* &
            latitude_function_a(sinphi, cosphi) + &
            earth_radius*planetary_rotation_rate*latitude_function_b(sinphi, cosphi))
          temperature_grid(i, j) = jablonowski_mean_temperature(eta) + temperature_correction
        end do
      end do
      call flux_curl(transform, u, v, level_spectral)
      zeta(:, :, k) = level_spectral
      call flux_divergence(transform, u, v, level_spectral)
      delta(:, :, k) = level_spectral
      call transform%grid_to_spectral(temperature_grid, level_spectral)
      temperature(:, :, k) = level_spectral
      zeta(0, 0, k) = 0.0_real64
      delta(0, 0, k) = 0.0_real64
    end do
  end subroutine jablonowski_williamson_topographic_initial_state

  !> Surface pressure satisfying Phi*(phi,p_s)=Phi_s.  Model terrain is far
  !> below the JW tropopause, so the closed tropospheric mean geopotential is
  !> used in the Newton residual exactly as documented.
  real(real64) function topographic_surface_pressure(sinphi, surface_geopotential, &
                                                     active_rotation_rate) result(surface_pressure)
    real(real64), intent(in) :: sinphi, surface_geopotential, active_rotation_rate
    real(real64) :: eta, residual, temperature, cosphi
    integer :: iteration
    real(real64), parameter :: base_temperature = 288.0_real64
    real(real64), parameter :: lapse_rate = 0.005_real64

    if (surface_geopotential == 0.0_real64) then
      surface_pressure = reference_surface_pressure
      return
    end if
    surface_pressure = reference_surface_pressure*exp(-surface_geopotential/ &
      (dry_air_gas_constant*base_temperature))
    cosphi = sqrt(max(0.0_real64, 1.0_real64 - sinphi*sinphi))
    do iteration = 1, 20
      eta = surface_pressure/reference_surface_pressure
      temperature = jablonowski_temperature(sinphi, cosphi, eta, active_rotation_rate)
      residual = pressure_geopotential(sinphi, cosphi, eta, active_rotation_rate, &
                                       base_temperature, lapse_rate) - surface_geopotential
      if (abs(residual) < 1.0e-6_real64) exit
      surface_pressure = surface_pressure*exp(residual/(dry_air_gas_constant*temperature))
    end do
    eta = surface_pressure/reference_surface_pressure
    residual = pressure_geopotential(sinphi, cosphi, eta, active_rotation_rate, &
                                     base_temperature, lapse_rate) - surface_geopotential
    if (abs(residual) >= 1.0e-6_real64) error stop 'topographic surface pressure did not converge'
  end function topographic_surface_pressure

  pure real(real64) function pressure_geopotential(sinphi, cosphi, eta, active_rotation_rate, &
                                                   base_temperature, lapse_rate) result(value)
    real(real64), intent(in) :: sinphi, cosphi, eta, active_rotation_rate, base_temperature, lapse_rate
    real(real64) :: vertical_factor, surface_factor

    vertical_factor = max(0.0_real64, cos(0.5_real64*acos(-1.0_real64)* &
      (eta - jablonowski_jet_eta)))**1.5_real64
    surface_factor = cos(0.5_real64*acos(-1.0_real64)*(1.0_real64 - jablonowski_jet_eta))**1.5_real64
    value = earth_gravity*base_temperature/lapse_rate* &
      (1.0_real64 - eta**(dry_air_gas_constant*lapse_rate/earth_gravity)) + &
      jablonowski_maximum_wind*vertical_factor*(jablonowski_maximum_wind*vertical_factor* &
      latitude_function_a(sinphi, cosphi) + earth_radius*active_rotation_rate*latitude_function_b(sinphi, cosphi)) - &
      jablonowski_maximum_wind*surface_factor*(jablonowski_maximum_wind*surface_factor* &
      latitude_function_a(sinphi, cosphi) + earth_radius*active_rotation_rate*latitude_function_b(sinphi, cosphi))
  end function pressure_geopotential

  pure real(real64) function jablonowski_temperature(sinphi, cosphi, eta, active_rotation_rate) result(value)
    real(real64), intent(in) :: sinphi, cosphi, eta, active_rotation_rate
    real(real64) :: eta_v, vertical_shape

    eta_v = 0.5_real64*acos(-1.0_real64)*(eta - jablonowski_jet_eta)
    vertical_shape = max(0.0_real64, cos(eta_v))
    value = jablonowski_mean_temperature(eta) + 0.75_real64*eta*acos(-1.0_real64)* &
      jablonowski_maximum_wind/dry_air_gas_constant*sin(eta_v)*sqrt(vertical_shape)*( &
      2.0_real64*jablonowski_maximum_wind*vertical_shape**1.5_real64*latitude_function_a(sinphi, cosphi) + &
      earth_radius*active_rotation_rate*latitude_function_b(sinphi, cosphi))
  end function jablonowski_temperature

  pure real(real64) function double_angle_sine_squared(sinphi, cosphi) result(value)
    real(real64), intent(in) :: sinphi, cosphi

    value = (2.0_real64*sinphi*cosphi)**2
  end function double_angle_sine_squared

  ! Zonal wind in gradient-wind balance, (f + u tan(phi)/a) u = -(1/a) dPhi/dphi, with the
  ! Jablonowski-Williamson temperature over flat terrain (Phi_s = 0).  vertical_factor and
  ! surface_factor are cos^(3/2) eta_v at the level and at eta = 1.  The rationalised root
  ! avoids cancellation near the poles.
  pure real(real64) function flat_terrain_balanced_wind(sinphi, cosphi, vertical_factor, &
                                                       surface_factor, active_rotation_rate) result(value)
    real(real64), intent(in) :: sinphi, cosphi, vertical_factor, surface_factor, active_rotation_rate
    real(real64) :: jet_wind, surface_wind, planetary_speed, discriminant_excess

    jet_wind = jablonowski_maximum_wind*vertical_factor*double_angle_sine_squared(sinphi, cosphi)
    surface_wind = jablonowski_maximum_wind*surface_factor*double_angle_sine_squared(sinphi, cosphi)
    planetary_speed = earth_radius*active_rotation_rate*cosphi
    ! Nonnegative because cos^(3/2) eta_v is smallest at eta = 1.
    discriminant_excess = (jet_wind - surface_wind)*(jet_wind + surface_wind + 2.0_real64*planetary_speed)
    if (discriminant_excess <= 0.0_real64) then
      value = 0.0_real64
    else
      value = discriminant_excess/(planetary_speed + sqrt(planetary_speed**2 + discriminant_excess))
    end if
  end function flat_terrain_balanced_wind

  pure real(real64) function latitude_function_a(sinphi, cosphi) result(value)
    real(real64), intent(in) :: sinphi, cosphi

    value = -2.0_real64*sinphi**6*(cosphi**2 + 1.0_real64/3.0_real64) + 10.0_real64/63.0_real64
  end function latitude_function_a

  pure real(real64) function latitude_function_b(sinphi, cosphi) result(value)
    real(real64), intent(in) :: sinphi, cosphi

    value = 8.0_real64/5.0_real64*cosphi**3*(sinphi**2 + 2.0_real64/3.0_real64) - &
            acos(-1.0_real64)/4.0_real64
  end function latitude_function_b

end module dry_initial_conditions
