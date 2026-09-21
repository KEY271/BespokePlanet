!> Bulk surface fluxes between the lowest model level and the ground or slab
!> ocean (docs/tendency/ground.md, docs/tendency/slab-ocean.md,
!> docs/tendency/evaporation.md).
!>
!> The sensible heat and the evaporation share the lowest-level density
!> rho_N = p_N/(R T_N), the transfer coefficient and the gusty wind speed.  The
!> lowest level is compared with the surface after dry-adiabatic extrapolation
!> to the surface pressure, so that a neutral column exchanges no sensible heat.
module surface_exchange
  use iso_fortran_env, only: real64
  use dry_vertical_coordinate, only: dry_air_gas_constant, dry_air_kappa
  use dry_physics_config, only: radiation_config
  use moist_thermodynamics, only: saturation_specific_humidity
  implicit none
  private
  public :: lowest_full_level_pressure, surface_transfer_mass_flux
  public :: surface_sensible_heat_flux, surface_evaporation_flux

contains

  !> p_N = p_{N+1/2} exp(-alpha_N) of one column of half-level pressures.
  pure real(real64) function lowest_full_level_pressure(pressure_half) result(pressure)
    real(real64), intent(in) :: pressure_half(0:)
    integer :: levels
    real(real64) :: thickness, alpha

    levels = ubound(pressure_half, 1)
    thickness = pressure_half(levels) - pressure_half(levels - 1)
    alpha = 1.0_real64 - pressure_half(levels - 1)*log(pressure_half(levels)/pressure_half(levels - 1))/thickness
    pressure = pressure_half(levels)*exp(-alpha)
  end function lowest_full_level_pressure

  !> rho_N C_H sqrt(u^2 + v^2 + U_g^2) in kg m^-2 s^-1, shared by both fluxes.
  pure real(real64) function surface_transfer_mass_flux(config, lowest_pressure, lowest_temperature, &
                                                        lowest_u, lowest_v) result(mass_flux)
    type(radiation_config), intent(in) :: config
    real(real64), intent(in) :: lowest_pressure, lowest_temperature, lowest_u, lowest_v

    mass_flux = lowest_pressure/(dry_air_gas_constant*lowest_temperature)*config%surface_exchange_coefficient* &
      sqrt(lowest_u**2 + lowest_v**2 + config%gustiness_speed**2)
  end function surface_transfer_mass_flux

  !> Upward sensible heat flux H = rho_N c_p C_H |U| [T_s - T_N (p_s/p_N)^kappa] in W m^-2.
  pure real(real64) function surface_sensible_heat_flux(config, pressure_half, lowest_temperature, &
                                                        surface_temperature, lowest_u, lowest_v) result(flux)
    type(radiation_config), intent(in) :: config
    real(real64), intent(in) :: pressure_half(0:), lowest_temperature, surface_temperature, lowest_u, lowest_v
    real(real64) :: lowest_pressure, surface_pressure

    lowest_pressure = lowest_full_level_pressure(pressure_half)
    surface_pressure = pressure_half(ubound(pressure_half, 1))
    flux = surface_transfer_mass_flux(config, lowest_pressure, lowest_temperature, lowest_u, lowest_v)* &
      config%dry_air_specific_heat* &
      (surface_temperature - lowest_temperature*(surface_pressure/lowest_pressure)**dry_air_kappa)
  end function surface_sensible_heat_flux

  !> Upward water vapour flux E = beta rho_N C_E |U| [q_s(T_s, p_s) - max(q_N, 0)] in kg m^-2 s^-1.
  !> Negative values (dew) are kept.
  pure real(real64) function surface_evaporation_flux(config, surface_wetness, pressure_half, lowest_temperature, &
                                                      lowest_humidity, surface_temperature, lowest_u, lowest_v) &
    result(flux)
    type(radiation_config), intent(in) :: config
    real(real64), intent(in) :: surface_wetness, pressure_half(0:), lowest_temperature, lowest_humidity
    real(real64), intent(in) :: surface_temperature, lowest_u, lowest_v
    real(real64) :: lowest_pressure, surface_pressure

    lowest_pressure = lowest_full_level_pressure(pressure_half)
    surface_pressure = pressure_half(ubound(pressure_half, 1))
    flux = surface_wetness*surface_transfer_mass_flux(config, lowest_pressure, lowest_temperature, lowest_u, lowest_v)* &
      (saturation_specific_humidity(surface_temperature, surface_pressure) - max(lowest_humidity, 0.0_real64))
  end function surface_evaporation_flux

end module surface_exchange
