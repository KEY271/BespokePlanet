module dry_initial_conditions
  use iso_fortran_env, only: real64
  use harmonics, only: harmonic_transform
  use barotropic_vorticity, only: earth_radius, rotation_rate
  use shallow_water_nonlinear, only: flux_divergence, flux_curl
  use dry_vertical_coordinate, only: hybrid_sigma_coordinate, dry_air_gas_constant, &
                                     reference_surface_pressure
  implicit none
  private

  real(real64), parameter, public :: jablonowski_maximum_wind = 35.0_real64
  real(real64), parameter, public :: jablonowski_jet_eta = 0.252_real64
  real(real64), parameter, public :: jablonowski_perturbation_wind = 1.0_real64

  public :: jablonowski_williamson_initial_state

contains

  subroutine jablonowski_williamson_initial_state(transform, truncation, coordinate, &
                                                  include_perturbation, zeta, delta, temperature, &
                                                  log_surface_pressure, surface_geopotential)
    type(harmonic_transform), intent(inout) :: transform
    integer, intent(in) :: truncation
    type(hybrid_sigma_coordinate), intent(in) :: coordinate
    logical, intent(in), optional :: include_perturbation
    complex(real64), allocatable, intent(out) :: zeta(:, :, :), delta(:, :, :), temperature(:, :, :)
    complex(real64), allocatable, intent(out) :: log_surface_pressure(:, :)
    complex(real64), allocatable, intent(out) :: surface_geopotential(:, :)
    real(real64), allocatable :: u(:, :), v(:, :), temperature_grid(:, :), log_ps_grid(:, :)
    real(real64), allocatable :: surface_geopotential_grid(:, :)
    complex(real64), allocatable :: level_spectral(:, :)
    integer, allocatable :: nlon(:)
    integer :: i, j, k, number_of_levels
    logical :: perturb
    real(real64) :: pi, longitude, latitude, eta, eta_v, vertical_shape, surface_shape
    real(real64) :: sinphi, cosphi, angular_cosine, distance, wind_perturbation
    real(real64) :: temperature_correction
    real(real64), parameter :: centre_longitude = 20.0_real64*acos(-1.0_real64)/180.0_real64
    real(real64), parameter :: centre_latitude = 40.0_real64*acos(-1.0_real64)/180.0_real64
    real(real64), parameter :: perturbation_radius = earth_radius/10.0_real64

    perturb = .true.
    if (present(include_perturbation)) perturb = include_perturbation
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
    surface_geopotential_grid = 0.0_real64
    surface_shape = cos(0.5_real64*pi*(1.0_real64 - jablonowski_jet_eta))**1.5_real64
    do j = 1, size(nlon)
      sinphi = transform%mu(j)
      cosphi = sqrt(max(0.0_real64, 1.0_real64 - sinphi*sinphi))
      surface_geopotential_grid(1:nlon(j), j) = jablonowski_maximum_wind*surface_shape*( &
        jablonowski_maximum_wind*surface_shape*latitude_function_a(sinphi, cosphi) + &
        earth_radius*rotation_rate*latitude_function_b(sinphi, cosphi))
    end do
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
        latitude = asin(sinphi)
        temperature_correction = 0.75_real64*eta*pi*jablonowski_maximum_wind/ &
          dry_air_gas_constant*sin(eta_v)*sqrt(vertical_shape)*( &
          2.0_real64*jablonowski_maximum_wind*vertical_shape**1.5_real64* &
          latitude_function_a(sinphi, cosphi) + &
          earth_radius*rotation_rate*latitude_function_b(sinphi, cosphi))
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
          u(i, j) = jablonowski_maximum_wind*vertical_shape**1.5_real64* &
                    sin(2.0_real64*latitude)**2 + wind_perturbation
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
