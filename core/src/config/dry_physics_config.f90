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
  !> selectable two-layer-ground or slab-ocean surface energy budget, together
  !> with the calendar and orbit that drive the solar forcing.
  type, public :: radiation_config
    logical :: enabled = .false.
    real(real64) :: solar_day = day_seconds
    integer :: days_per_month = 30
    integer :: months_per_year = 12
    real(real64) :: axial_tilt = 23.4_real64*pi/180.0_real64
    real(real64) :: solar_constant = 1361.0_real64
    !> Albedo of the legacy single-surface path.  Without diagnosed clouds it
    !> is the planetary value with the cloud reflection folded in; the moist
    !> cases replace it by the open-ocean value (docs/tendency/shortwave-radiation.md).
    real(real64) :: surface_shortwave_albedo = 0.3_real64
    !> Reflectance of an overcast column; the downward shortwave reaching the
    !> troposphere is reflected once by C times this value.
    real(real64) :: cloud_shortwave_albedo = 0.43_real64
    !> When enabled, land and ocean temperatures are advanced separately and
    !> their fluxes are weighted by the fixed land fraction. The legacy ground/slab switch remains
    !> the exact path used by existing cases.  The land and ocean albedos are
    !> surface values without clouds.
    logical :: land_sea_mixing_enabled = .false.
    real(real64) :: land_shortwave_albedo = 0.2_real64
    real(real64) :: ocean_shortwave_albedo = 0.06_real64
    !> Grey longwave optical depth d tau/dp = (a mu + b q)/p_0: a well-mixed
    !> absorber (a mu) plus water vapour (b q), with the Byrne & O'Gorman (2013)
    !> coefficients as implemented in Isca (docs/tendency/longwave-radiation.md).
    real(real64) :: longwave_well_mixed_optical_depth = 0.8678_real64
    real(real64) :: longwave_water_vapor_optical_depth = 1997.9_real64
    real(real64) :: longwave_well_mixed_scaling = 1.0_real64
    real(real64) :: longwave_reference_pressure = 1.0e5_real64
    !> Surface value q_0 of the fixed reference humidity q_ref = q_0 (p/p_s)^3
    !> that the cases without prognostic water vapour give the longwave radiation.
    real(real64) :: longwave_reference_surface_humidity = 0.010_real64
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
    logical :: slab_ocean_enabled = .false.
    real(real64) :: slab_ocean_depth = 30.0_real64
    real(real64) :: seawater_density = 1000.0_real64
    real(real64) :: seawater_specific_heat = 4186.0_real64
    real(real64) :: surface_exchange_coefficient = 1.0e-3_real64
    real(real64) :: gustiness_speed = 1.0_real64
  end type radiation_config

  !> Enthalpy-conserving dry convective adjustment.
  type, public :: convection_config
    logical :: enabled = .false.
    real(real64) :: adjustment_time = 4.0_real64*3600.0_real64
  end type convection_config

  !> Water vapour as a prognostic variable: specific humidity is advected,
  !> diffused and filtered like temperature, and the virtual temperature enters
  !> the hydrostatic relation, the pressure-gradient term and the adiabatic
  !> heating.  The other moist processes require this to be enabled.
  type, public :: moisture_config
    logical :: enabled = .false.
    !> Initial state q = RH q_s(T, p) on the levels whose full-level pressure is
    !> at least initial_humidity_top_pressure; the levels above start dry.
    real(real64) :: initial_relative_humidity = 0.7_real64
    real(real64) :: initial_humidity_top_pressure = 2.0e4_real64
  end type moisture_config

  !> Bulk evaporation from the surface into the lowest model level, with the
  !> latent heat taken from the surface energy budget.  The bulk coefficient and
  !> gustiness are shared with the sensible heat flux of radiation_config.
  type, public :: evaporation_config
    logical :: enabled = .false.
    !> Surface wetness beta in [0, 1]; a slab ocean is saturated (1).
    real(real64) :: surface_wetness = 1.0_real64
    real(real64) :: land_surface_wetness = 0.5_real64
    real(real64) :: ocean_surface_wetness = 1.0_real64
  end type evaporation_config

  !> One-layer finite land-water reservoir (docs/tendency/bucket.md).
  type, public :: bucket_config
    logical :: enabled = .false.
    real(real64) :: capacity = 150.0_real64
    real(real64) :: initial_water = 75.0_real64
    real(real64) :: dry_threshold_fraction = 0.1_real64
  end type bucket_config

  !> Simplified Betts-Miller moist convective adjustment (Frierson 2007).
  type, public :: moist_convection_config
    logical :: enabled = .false.
    real(real64) :: adjustment_time = 2.0_real64*3600.0_real64
    real(real64) :: reference_relative_humidity = 0.7_real64
  end type moist_convection_config

  !> Grid-scale condensation of supersaturated water vapour; the condensate falls
  !> out immediately.  The condensation completes within one leapfrog step.
  type, public :: condensation_config
    logical :: enabled = .false.
    !> Relative tolerance of the saturation solve, |q - q_s| <= tolerance*q_s.
    real(real64) :: saturation_tolerance = 1.0e-4_real64
    integer :: maximum_iterations = 10
  end type condensation_config

  !> Diagnostic column cloud cover from the column-maximum relative humidity and
  !> the convective precipitation (Slingo 1987; docs/tendency/cloud.md).  Used by
  !> the shortwave reflection only.
  type, public :: cloud_config
    logical :: enabled = .false.
    real(real64) :: critical_relative_humidity = 0.8_real64
    real(real64) :: convective_intercept = 0.245_real64
    real(real64) :: convective_slope = 0.125_real64
    !> P_0 = 1 mm day^-1 in kg m^-2 s^-1.
    real(real64) :: convective_reference_precipitation = 1.0_real64/day_seconds
    real(real64) :: convective_maximum_cover = 0.8_real64
  end type cloud_config

  !> Snow-free, zero-heat-capacity sea ice; volume is per unit ocean area.
  type, public :: sea_ice_config
    logical :: enabled = .false.
    real(real64) :: freezing_temperature = 271.35_real64
    real(real64) :: melting_temperature = 273.15_real64
    real(real64) :: new_ice_thickness = 0.5_real64
    real(real64) :: conductivity = 2.0_real64
    real(real64) :: density = 917.0_real64
    real(real64) :: latent_heat = 3.34e5_real64
    real(real64) :: albedo = 0.60_real64
    real(real64) :: temperature_tolerance = 1.0e-7_real64
    real(real64) :: flux_tolerance = 1.0e-5_real64
    integer :: maximum_iterations = 100
  end type sea_ice_config

  type, public :: dry_model_physics_config
    type(sea_ice_config) :: sea_ice
    type(held_suarez_config) :: held_suarez
    type(surface_friction_config) :: surface_friction
    type(rayleigh_friction_config) :: rayleigh_friction
    type(radiation_config) :: radiation
    type(convection_config) :: convection
    type(moisture_config) :: moisture
    type(evaporation_config) :: evaporation
    type(bucket_config) :: bucket
    type(moist_convection_config) :: moist_convection
    type(condensation_config) :: condensation
    type(cloud_config) :: cloud
  end type dry_model_physics_config

  public :: radiation_days_per_year, radiation_orbital_period, radiation_planet_rotation_rate
  public :: radiation_surface_heat_capacity
  public :: mixed_surface_properties

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

  !> Heat capacity of the active surface layer.  The slab-ocean case replaces
  !> the two-layer ground surface with one well-mixed water column.
  pure real(real64) function radiation_surface_heat_capacity(config) result(heat_capacity)
    type(radiation_config), intent(in) :: config

    if (config%slab_ocean_enabled) then
      heat_capacity = config%seawater_density*config%seawater_specific_heat*config%slab_ocean_depth
    else
      heat_capacity = config%surface_heat_capacity
    end if
  end function radiation_surface_heat_capacity

  !> Time-independent surface coefficients at one mixed land--ocean grid point.
  !> Wetness is retained for legacy callers; bucket-enabled evaporation diagnoses
  !> its land value from the prognostic surface water instead.
  pure subroutine mixed_surface_properties(radiation, evaporation, land_fraction, heat_capacity, albedo, &
                                           wetness, ground_exchange)
    type(radiation_config), intent(in) :: radiation
    type(evaporation_config), intent(in) :: evaporation
    real(real64), intent(in) :: land_fraction
    real(real64), intent(out) :: heat_capacity, albedo, wetness, ground_exchange
    real(real64) :: ocean_heat_capacity

    ocean_heat_capacity = radiation%seawater_density*radiation%seawater_specific_heat*radiation%slab_ocean_depth
    heat_capacity = land_fraction*radiation%surface_heat_capacity + &
      (1.0_real64 - land_fraction)*ocean_heat_capacity
    albedo = land_fraction*radiation%land_shortwave_albedo + &
      (1.0_real64 - land_fraction)*radiation%ocean_shortwave_albedo
    wetness = land_fraction*evaporation%land_surface_wetness + &
      (1.0_real64 - land_fraction)*evaporation%ocean_surface_wetness
    ground_exchange = land_fraction*radiation%ground_exchange_coefficient
  end subroutine mixed_surface_properties

end module dry_physics_config
