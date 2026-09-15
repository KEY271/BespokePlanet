module dry_held_suarez
  use iso_fortran_env, only: real64
  use harmonics, only: harmonic_transform
  use dry_vertical_coordinate, only: hybrid_sigma_coordinate, dry_air_kappa, &
                                     reference_surface_pressure
  implicit none
  private

  real(real64), parameter, public :: held_suarez_initial_temperature = 264.0_real64
  real(real64), parameter, public :: held_suarez_temperature_perturbation = 0.1_real64
  real(real64), parameter, public :: held_suarez_sigma_boundary = 0.7_real64
  real(real64), parameter, public :: held_suarez_minimum_equilibrium_temperature = 200.0_real64
  real(real64), parameter, public :: held_suarez_equatorial_temperature = 315.0_real64
  real(real64), parameter, public :: held_suarez_equator_to_pole_difference = 60.0_real64
  real(real64), parameter, public :: held_suarez_vertical_difference = 10.0_real64
  real(real64), parameter, public :: held_suarez_day = 86400.0_real64
  real(real64), parameter, public :: held_suarez_friction_rate = 1.0_real64/held_suarez_day
  real(real64), parameter, public :: held_suarez_upper_thermal_rate = 1.0_real64/(40.0_real64*held_suarez_day)
  real(real64), parameter, public :: held_suarez_lower_thermal_rate = 1.0_real64/(4.0_real64*held_suarez_day)

  public :: held_suarez_initial_state
  public :: held_suarez_forcing

contains

  subroutine held_suarez_initial_state(transform, truncation, coordinate, zeta, delta, temperature, &
                                       log_surface_pressure, surface_geopotential)
    type(harmonic_transform), intent(inout) :: transform
    integer, intent(in) :: truncation
    type(hybrid_sigma_coordinate), intent(in) :: coordinate
    complex(real64), allocatable, intent(out) :: zeta(:, :, :), delta(:, :, :), temperature(:, :, :)
    complex(real64), allocatable, intent(out) :: log_surface_pressure(:, :), surface_geopotential(:, :)
    real(real64), allocatable :: temperature_grid(:, :), log_ps_grid(:, :)
    complex(real64), allocatable :: temperature_spectral(:, :)
    integer, allocatable :: nlon(:)
    integer :: i, j, k, number_of_levels
    real(real64) :: longitude, cosphi, pi

    number_of_levels = coordinate%number_of_levels
    if (number_of_levels < 1) error stop 'Held-Suarez coordinate is not initialized'
    pi = acos(-1.0_real64)
    nlon = transform%get_nlon()
    call transform%allocate_field(temperature_grid)
    call transform%allocate_field(log_ps_grid)
    temperature_grid = 0.0_real64
    do j = 1, size(nlon)
      cosphi = sqrt(max(0.0_real64, 1.0_real64 - transform%mu(j)**2))
      do i = 1, nlon(j)
        longitude = 2.0_real64*pi*real(i - 1, real64)/real(nlon(j), real64)
        temperature_grid(i, j) = held_suarez_initial_temperature + &
          held_suarez_temperature_perturbation*cosphi*sin(longitude)
      end do
    end do
    call transform%grid_to_spectral(temperature_grid, temperature_spectral)

    allocate (zeta(0:truncation + 1, 0:truncation, number_of_levels))
    allocate (delta(0:truncation + 1, 0:truncation, number_of_levels))
    allocate (temperature(0:truncation + 1, 0:truncation, number_of_levels))
    zeta = 0.0_real64
    delta = 0.0_real64
    do k = 1, number_of_levels
      temperature(:, :, k) = temperature_spectral
    end do

    log_ps_grid = log(reference_surface_pressure)
    call transform%grid_to_spectral(log_ps_grid, log_surface_pressure)
    allocate (surface_geopotential(0:truncation + 1, 0:truncation))
    surface_geopotential = 0.0_real64
  end subroutine held_suarez_initial_state

  pure subroutine held_suarez_forcing(sinphi, full_level_pressure, surface_pressure, temperature, &
                                      u, v, forcing_u, forcing_v, temperature_tendency)
    real(real64), intent(in) :: sinphi, full_level_pressure, surface_pressure, temperature, u, v
    real(real64), intent(out) :: forcing_u, forcing_v, temperature_tendency
    real(real64) :: cosphi_squared, pressure_ratio, sigma, boundary_weight
    real(real64) :: equilibrium_temperature, thermal_rate, friction_rate

    cosphi_squared = max(0.0_real64, 1.0_real64 - sinphi**2)
    pressure_ratio = full_level_pressure/reference_surface_pressure
    sigma = full_level_pressure/surface_pressure
    boundary_weight = max(0.0_real64, (sigma - held_suarez_sigma_boundary)/ &
                          (1.0_real64 - held_suarez_sigma_boundary))
    equilibrium_temperature = max(held_suarez_minimum_equilibrium_temperature, &
      (held_suarez_equatorial_temperature - held_suarez_equator_to_pole_difference*sinphi**2 - &
       held_suarez_vertical_difference*log(pressure_ratio)*cosphi_squared)*pressure_ratio**dry_air_kappa)
    thermal_rate = held_suarez_upper_thermal_rate + &
      (held_suarez_lower_thermal_rate - held_suarez_upper_thermal_rate)* &
      boundary_weight*cosphi_squared**2
    friction_rate = held_suarez_friction_rate*boundary_weight

    forcing_u = -friction_rate*u
    forcing_v = -friction_rate*v
    temperature_tendency = -thermal_rate*(temperature - equilibrium_temperature)
  end subroutine held_suarez_forcing

end module dry_held_suarez
