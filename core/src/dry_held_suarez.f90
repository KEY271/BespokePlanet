module dry_held_suarez
  use iso_fortran_env, only: real64
  use harmonics, only: harmonic_transform
  use dry_vertical_coordinate, only: hybrid_sigma_coordinate, dry_air_kappa, &
                                     reference_surface_pressure
  use dry_physics_config, only: held_suarez_config, surface_friction_config
  implicit none
  private

  !> Initial-condition values of the Held-Suarez case; the forcing coefficients
  !> live in held_suarez_config and surface_friction_config.
  real(real64), parameter, public :: held_suarez_initial_temperature = 264.0_real64
  real(real64), parameter, public :: held_suarez_temperature_perturbation = 0.1_real64

  public :: held_suarez_initial_state
  public :: held_suarez_forcing
  public :: held_suarez_friction
  public :: held_suarez_thermal_relaxation

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

  !> Fraction of the Rayleigh boundary layer at this level, zero above sigma_b.
  !> Shared by the boundary friction and the thermal relaxation rate.
  pure real(real64) function held_suarez_boundary_weight(sigma_boundary, full_level_pressure, &
                                                         surface_pressure) result(weight)
    real(real64), intent(in) :: sigma_boundary, full_level_pressure, surface_pressure
    real(real64) :: sigma

    sigma = full_level_pressure/surface_pressure
    weight = max(0.0_real64, (sigma - sigma_boundary)/(1.0_real64 - sigma_boundary))
  end function held_suarez_boundary_weight

  !> Boundary-layer Rayleigh drag on the horizontal wind.  This is the only part
  !> of the Held-Suarez forcing that the radiation case also applies.
  pure subroutine held_suarez_friction(config, full_level_pressure, surface_pressure, u, v, forcing_u, forcing_v)
    type(surface_friction_config), intent(in) :: config
    real(real64), intent(in) :: full_level_pressure, surface_pressure, u, v
    real(real64), intent(out) :: forcing_u, forcing_v
    real(real64) :: friction_rate

    friction_rate = config%friction_rate* &
                    held_suarez_boundary_weight(config%sigma_boundary, full_level_pressure, surface_pressure)
    forcing_u = -friction_rate*u
    forcing_v = -friction_rate*v
  end subroutine held_suarez_friction

  !> Newtonian relaxation of temperature towards the Held-Suarez equilibrium profile.
  pure subroutine held_suarez_thermal_relaxation(config, sinphi, full_level_pressure, surface_pressure, &
                                                 temperature, temperature_tendency)
    type(held_suarez_config), intent(in) :: config
    real(real64), intent(in) :: sinphi, full_level_pressure, surface_pressure, temperature
    real(real64), intent(out) :: temperature_tendency
    real(real64) :: cosphi_squared, pressure_ratio, boundary_weight
    real(real64) :: equilibrium_temperature, thermal_rate

    cosphi_squared = max(0.0_real64, 1.0_real64 - sinphi**2)
    pressure_ratio = full_level_pressure/reference_surface_pressure
    boundary_weight = held_suarez_boundary_weight(config%sigma_boundary, full_level_pressure, surface_pressure)
    equilibrium_temperature = max(config%minimum_equilibrium_temperature, &
      (config%equatorial_temperature - config%equator_to_pole_difference*sinphi**2 - &
       config%vertical_difference*log(pressure_ratio)*cosphi_squared)*pressure_ratio**dry_air_kappa)
    thermal_rate = config%upper_thermal_rate + &
      (config%lower_thermal_rate - config%upper_thermal_rate)* &
      boundary_weight*cosphi_squared**2
    temperature_tendency = -thermal_rate*(temperature - equilibrium_temperature)
  end subroutine held_suarez_thermal_relaxation

  !> Full Held-Suarez forcing; the composition of the two independent parts above.
  pure subroutine held_suarez_forcing(relaxation, friction, sinphi, full_level_pressure, surface_pressure, &
                                      temperature, u, v, forcing_u, forcing_v, temperature_tendency)
    type(held_suarez_config), intent(in) :: relaxation
    type(surface_friction_config), intent(in) :: friction
    real(real64), intent(in) :: sinphi, full_level_pressure, surface_pressure, temperature, u, v
    real(real64), intent(out) :: forcing_u, forcing_v, temperature_tendency

    call held_suarez_friction(friction, full_level_pressure, surface_pressure, u, v, forcing_u, forcing_v)
    call held_suarez_thermal_relaxation(relaxation, sinphi, full_level_pressure, surface_pressure, &
                                        temperature, temperature_tendency)
  end subroutine held_suarez_forcing

end module dry_held_suarez
