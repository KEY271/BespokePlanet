!> Constants of water vapour and the saturation relations shared by evaporation,
!> large-scale condensation and the moist convective adjustment
!> (docs/tendency/saturation-specific-humidity.md).
!>
!> The latent heat is constant, there is no ice phase, and c_p and kappa keep
!> their dry-air values.
module moist_thermodynamics
  use iso_fortran_env, only: real64
  use dry_vertical_coordinate, only: dry_air_gas_constant, dry_air_kappa, reference_surface_pressure
  implicit none
  private

  real(real64), parameter, public :: dry_air_specific_heat = 1004.0_real64
  real(real64), parameter, public :: water_vapor_gas_constant = 461.5_real64
  !> epsilon = R/R_v
  real(real64), parameter, public :: gas_constant_ratio = dry_air_gas_constant/water_vapor_gas_constant
  !> delta_v = 1/epsilon - 1, so that T_v = (1 + delta_v q) T.
  real(real64), parameter, public :: virtual_temperature_coefficient = 1.0_real64/gas_constant_ratio - 1.0_real64
  real(real64), parameter, public :: latent_heat_of_condensation = 2.5e6_real64
  !> Saturation vapour pressure e_0 at the reference temperature T_0.
  real(real64), parameter, public :: reference_saturation_vapor_pressure = 610.78_real64
  real(real64), parameter, public :: saturation_reference_temperature = 273.16_real64
  !> Newton iteration bounds of the lifting-condensation-level solve.
  integer, parameter, public :: lifting_condensation_maximum_iterations = 20
  real(real64), parameter, public :: lifting_condensation_temperature_tolerance = 1.0e-3_real64

  public :: saturation_vapor_pressure, saturation_specific_humidity
  public :: saturation_specific_humidity_derivative
  public :: virtual_temperature, equivalent_potential_temperature
  public :: saturated_log_equivalent_potential_temperature_derivative
  public :: lifting_condensation_level
  public :: full_level_pressures

contains

  !> Clausius-Clapeyron with constant L, integrated from (T_0, e_0).
  pure real(real64) function saturation_vapor_pressure(temperature) result(pressure)
    real(real64), intent(in) :: temperature
    pressure = reference_saturation_vapor_pressure*exp(latent_heat_of_condensation/water_vapor_gas_constant* &
      (1.0_real64/saturation_reference_temperature - 1.0_real64/temperature))
  end function saturation_vapor_pressure

  !> q_s = eps e_s/(p - (1 - eps) e_s); one where p <= e_s, so that such a layer is never saturated.
  pure real(real64) function saturation_specific_humidity(temperature, pressure) result(humidity)
    real(real64), intent(in) :: temperature, pressure
    real(real64) :: vapor_pressure

    vapor_pressure = saturation_vapor_pressure(temperature)
    if (pressure <= vapor_pressure) then
      humidity = 1.0_real64
    else
      humidity = gas_constant_ratio*vapor_pressure/(pressure - (1.0_real64 - gas_constant_ratio)*vapor_pressure)
    end if
  end function saturation_specific_humidity

  !> dq_s/dT = q_s L/(R_v T^2) p/(p - (1 - eps) e_s); zero where p <= e_s.
  pure real(real64) function saturation_specific_humidity_derivative(temperature, pressure) result(derivative)
    real(real64), intent(in) :: temperature, pressure
    real(real64) :: vapor_pressure, humidity

    vapor_pressure = saturation_vapor_pressure(temperature)
    if (pressure <= vapor_pressure) then
      derivative = 0.0_real64
    else
      humidity = gas_constant_ratio*vapor_pressure/(pressure - (1.0_real64 - gas_constant_ratio)*vapor_pressure)
      derivative = humidity*latent_heat_of_condensation/(water_vapor_gas_constant*temperature**2)* &
        pressure/(pressure - (1.0_real64 - gas_constant_ratio)*vapor_pressure)
    end if
  end function saturation_specific_humidity_derivative

  !> T_v = (1 + delta_v max(q, 0)) T.  Negative spectral-truncation humidity is clipped here.
  pure real(real64) function virtual_temperature(temperature, humidity) result(value)
    real(real64), intent(in) :: temperature, humidity
    value = (1.0_real64 + virtual_temperature_coefficient*max(humidity, 0.0_real64))*temperature
  end function virtual_temperature

  !> theta_e = T (p_0/p)^kappa exp(L q/(c_p T)), the simplified pseudo-adiabatic invariant.
  pure real(real64) function equivalent_potential_temperature(temperature, humidity, pressure) result(theta_e)
    real(real64), intent(in) :: temperature, humidity, pressure
    theta_e = temperature*(reference_surface_pressure/pressure)**dry_air_kappa* &
      exp(latent_heat_of_condensation*humidity/(dry_air_specific_heat*temperature))
  end function equivalent_potential_temperature

  !> d ln theta_e/dT of a saturated parcel at fixed pressure:
  !>   1/T + L/(c_p T) (dq_s/dT - q_s/T), which is always positive.
  pure real(real64) function saturated_log_equivalent_potential_temperature_derivative(temperature, pressure) &
    result(derivative)
    real(real64), intent(in) :: temperature, pressure
    derivative = 1.0_real64/temperature + latent_heat_of_condensation/(dry_air_specific_heat*temperature)* &
      (saturation_specific_humidity_derivative(temperature, pressure) - &
       saturation_specific_humidity(temperature, pressure)/temperature)
  end function saturated_log_equivalent_potential_temperature_derivative

  !> Pressure and temperature at which a parcel (T_N, q_N > 0, p_N) lifted dry
  !> adiabatically just saturates.  A parcel that is already saturated stays where it is.
  pure subroutine lifting_condensation_level(parcel_temperature, parcel_humidity, parcel_pressure, &
                                             level_temperature, level_pressure)
    real(real64), intent(in) :: parcel_temperature, parcel_humidity, parcel_pressure
    real(real64), intent(out) :: level_temperature, level_pressure
    real(real64) :: temperature, pressure, vapor_pressure, residual, slope, change
    integer :: iteration

    if (parcel_humidity >= saturation_specific_humidity(parcel_temperature, parcel_pressure)) then
      level_temperature = parcel_temperature
      level_pressure = parcel_pressure
      return
    end if
    temperature = parcel_temperature
    pressure = parcel_pressure
    do iteration = 1, lifting_condensation_maximum_iterations
      pressure = parcel_pressure*(temperature/parcel_temperature)**(1.0_real64/dry_air_kappa)
      residual = log(saturation_specific_humidity(temperature, pressure)) - log(parcel_humidity)
      vapor_pressure = saturation_vapor_pressure(temperature)
      slope = pressure/(pressure - (1.0_real64 - gas_constant_ratio)*vapor_pressure)/temperature* &
        (latent_heat_of_condensation/(water_vapor_gas_constant*temperature) - 1.0_real64/dry_air_kappa)
      change = -residual/slope
      temperature = temperature + change
      if (abs(change) < lifting_condensation_temperature_tolerance) exit
    end do
    level_temperature = temperature
    level_pressure = parcel_pressure*(temperature/parcel_temperature)**(1.0_real64/dry_air_kappa)
  end subroutine lifting_condensation_level

  !> Full-level pressures p_k = p_{k+1/2} exp(-alpha_k) of one column and the layer thicknesses.
  pure subroutine full_level_pressures(pressure_half, full_level_pressure, delta_pressure)
    real(real64), intent(in) :: pressure_half(0:)
    real(real64), intent(out) :: full_level_pressure(:), delta_pressure(:)
    real(real64) :: alpha
    integer :: k

    do k = 1, size(full_level_pressure)
      delta_pressure(k) = pressure_half(k) - pressure_half(k - 1)
      alpha = 1.0_real64 - pressure_half(k - 1)*log(pressure_half(k)/pressure_half(k - 1))/delta_pressure(k)
      full_level_pressure(k) = pressure_half(k)*exp(-alpha)
    end do
  end subroutine full_level_pressures

end module moist_thermodynamics
