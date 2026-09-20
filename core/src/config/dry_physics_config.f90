!> Configuration of the dry-atmosphere physical processes.
!>
!> Each process has one type that records whether it is active and the
!> coefficients it uses.  The default initializers are the only place where the
!> default coefficients are written; the physics modules receive the values
!> through these types and carry no constants of their own.
module dry_physics_config
  use iso_fortran_env, only: real64
  use planet_parameters, only: earth_gravity
  implicit none
  private

  real(real64), parameter, public :: dry_reference_temperature = 300.0_real64
  real(real64), parameter :: pi = acos(-1.0_real64)
  real(real64), parameter :: day_seconds = 86400.0_real64

  !> Held-Suarez Newtonian relaxation of temperature.
  type, public :: held_suarez_config
    logical :: enabled = .false.
    real(real64) :: sigma_boundary = 0.7_real64
    real(real64) :: minimum_equilibrium_temperature = 200.0_real64
    real(real64) :: equatorial_temperature = 315.0_real64
    real(real64) :: equator_to_pole_difference = 60.0_real64
    real(real64) :: vertical_difference = 10.0_real64
    real(real64) :: upper_thermal_rate = 1.0_real64/(40.0_real64*day_seconds)
    real(real64) :: lower_thermal_rate = 1.0_real64/(4.0_real64*day_seconds)
  end type held_suarez_config

  !> Boundary-layer Rayleigh drag of the Held-Suarez form.  It is a process of
  !> its own because the radiation case applies it without the thermal relaxation.
  type, public :: surface_friction_config
    logical :: enabled = .false.
    real(real64) :: sigma_boundary = 0.7_real64
    real(real64) :: friction_rate = 1.0_real64/day_seconds
  end type surface_friction_config

  !> Rayleigh damping of the wind in the topmost levels (sponge layer).
  type, public :: rayleigh_friction_config
    logical :: enabled = .false.
    integer :: top_levels = 2
    real(real64) :: rate = 1.0_real64/day_seconds
  end type rayleigh_friction_config

  !> Grey longwave radiation, prescribed ozone shortwave absorption, and the
  !> two-layer ground energy budget, together with the calendar and orbit that
  !> drive the solar forcing.
  type, public :: radiation_config
    logical :: enabled = .false.
    real(real64) :: solar_day = day_seconds
    integer :: days_per_month = 30
    integer :: months_per_year = 12
    real(real64) :: axial_tilt = 23.4_real64*pi/180.0_real64
    real(real64) :: solar_constant = 1361.0_real64
    real(real64) :: surface_shortwave_albedo = 0.3_real64
    real(real64) :: longwave_surface_optical_depth = 1.0_real64
    real(real64) :: stefan_boltzmann_constant = 5.670374419e-8_real64
    real(real64) :: dry_air_specific_heat = 1004.0_real64
    real(real64) :: gravity_acceleration = earth_gravity
    real(real64) :: ultraviolet_shortwave_fraction = 0.02_real64
    real(real64) :: ozone_shortwave_optical_depth = 1.5_real64
    real(real64) :: ozone_longwave_optical_depth = 0.005_real64
    real(real64) :: ozone_pressure_lower_bound = 1.0e2_real64
    real(real64) :: ozone_pressure_upper_bound = 1.0e4_real64
    real(real64) :: ozone_peak_pressure = 1.0e3_real64
    real(real64) :: ozone_log_pressure_width = log(3.0_real64)
    real(real64) :: surface_heat_capacity = 2.0e6_real64
    real(real64) :: deep_ground_heat_capacity = 2.0e7_real64
    real(real64) :: ground_exchange_coefficient = 2.0_real64
    real(real64) :: surface_exchange_coefficient = 1.0e-3_real64
    real(real64) :: gustiness_speed = 1.0_real64
  end type radiation_config

  !> Enthalpy-conserving dry convective adjustment.
  type, public :: convection_config
    logical :: enabled = .false.
    real(real64) :: adjustment_time = 4.0_real64*3600.0_real64
  end type convection_config

  type, public :: dry_model_physics_config
    type(held_suarez_config) :: held_suarez
    type(surface_friction_config) :: surface_friction
    type(rayleigh_friction_config) :: rayleigh_friction
    type(radiation_config) :: radiation
    type(convection_config) :: convection
  end type dry_model_physics_config

  public :: radiation_days_per_year, radiation_orbital_period, radiation_planet_rotation_rate

contains

  pure integer function radiation_days_per_year(config) result(days)
    type(radiation_config), intent(in) :: config
    days = config%days_per_month*config%months_per_year
  end function radiation_days_per_year

  pure real(real64) function radiation_orbital_period(config) result(period)
    type(radiation_config), intent(in) :: config
    period = real(radiation_days_per_year(config), real64)*config%solar_day
  end function radiation_orbital_period

  !> Sidereal rotation rate of a planet whose solar day and year are those of the calendar.
  pure real(real64) function radiation_planet_rotation_rate(config) result(rate)
    type(radiation_config), intent(in) :: config
    rate = 2.0_real64*pi*(1.0_real64/config%solar_day + 1.0_real64/radiation_orbital_period(config))
  end function radiation_planet_rotation_rate

end module dry_physics_config
