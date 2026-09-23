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

  !> Radiation schemes: the grey longwave with UV-only ozone shortwave absorption
  !> (docs/tendency/longwave-radiation.md, shortwave-radiation.md) and the band
  !> scheme (docs/tendency/band-radiation.md).
  integer, parameter, public :: radiation_scheme_gray = 1
  integer, parameter, public :: radiation_scheme_band = 2

  !> Coefficients of the band radiation (docs/tendency/band-radiation.md).  The
  !> five longwave sub-bands are 1 window, 2 CO2 centre, 3 CO2 wings, 4 weak and
  !> 5 strong water vapour; the four shortwave bands are 1 UV, 2 visible, 3 weak
  !> and 4 strong near-infrared water vapour, split into twelve g points.
  !> Mass absorption coefficients are vertical and in m^2 kg^-1.
  type, public :: band_radiation_config
    real(real64) :: diffusivity = 1.66_real64
    !> Volume mixing ratio chi of CO2 and the reference chi_ref of the fitted coefficients.
    real(real64) :: co2_volume_mixing_ratio = 400.0e-6_real64
    real(real64) :: co2_reference_volume_mixing_ratio = 400.0e-6_real64
    !> M_CO2/M_air, converting the volume to the mass mixing ratio.
    real(real64) :: co2_molar_mass_ratio = 44.0095_real64/28.9647_real64
    !> Ozone column 300 DU in kg m^-2 (1 DU = 2.1415e-5 kg m^-2).
    real(real64) :: ozone_column = 300.0_real64*2.1415e-5_real64
    !> epsilon of the vapour pressure e = p q/epsilon of the self continuum.
    real(real64) :: water_vapor_molar_mass_ratio = 0.622_real64
    !> Temperature factor exp[T_*(1/T - 1/T_ref)] of the self continuum.
    real(real64) :: self_continuum_temperature = 1800.0_real64
    real(real64) :: self_continuum_reference_temperature = 296.0_real64
    real(real64) :: reference_pressure = 1.0e5_real64
    !> Second radiation constant c_2 = h c/k_B in cm K; the band edges are in cm^-1.
    real(real64) :: second_radiation_constant = 1.4387769_real64
    real(real64) :: longwave_edges(8) = [350.0_real64, 500.0_real64, 630.0_real64, 700.0_real64, &
                                         820.0_real64, 1180.0_real64, 1390.0_real64, 1800.0_real64]
    real(real64) :: longwave_line(5) = [0.0_real64, 1.0_real64, 0.067_real64, 0.39_real64, 8.5_real64]
    real(real64) :: longwave_line_pressure_exponent(5) = [1.0_real64, 1.0_real64, 1.0_real64, 1.0_real64, 0.0_real64]
    real(real64) :: longwave_self_continuum(5) = [0.69_real64, 0.0_real64, 0.0_real64, 0.0_real64, 0.0_real64]
    real(real64) :: longwave_co2(5) = [0.023_real64, 5.1_real64, 0.21_real64, 0.0_real64, 0.0_real64]
    real(real64) :: longwave_co2_pressure_exponent(5) = [1.0_real64, 0.0_real64, 0.7_real64, 1.0_real64, 1.0_real64]
    real(real64) :: longwave_co2_concentration_exponent(5) = [0.45_real64, 0.0_real64, 0.18_real64, 1.0_real64, &
                                                               1.0_real64]
    real(real64) :: longwave_ozone(5) = [7.6_real64, 0.0_real64, 0.0_real64, 0.0_real64, 0.0_real64]
    !> Fraction phi_b of the top-of-atmosphere solar flux in each shortwave band.
    real(real64) :: shortwave_band_fraction(4) = [0.0358_real64, 0.5074_real64, 0.2534_real64, 0.2034_real64]
    integer :: shortwave_gpoint_band(12) = [1, 1, 1, 2, 2, 3, 3, 3, 4, 4, 4, 4]
    real(real64) :: shortwave_gpoint_weight(12) = [0.27_real64, 0.25_real64, 0.48_real64, &
                                                   0.945_real64, 0.055_real64, &
                                                   0.104_real64, 0.314_real64, 0.582_real64, &
                                                   0.117_real64, 0.218_real64, 0.278_real64, 0.387_real64]
    real(real64) :: shortwave_gpoint_water_vapor(12) = [0.0_real64, 0.0_real64, 0.0_real64, &
                                                        0.0_real64, 0.0092_real64, &
                                                        0.21_real64, 0.0123_real64, 0.0_real64, &
                                                        8.8_real64, 0.49_real64, 0.0227_real64, 0.0_real64]
    real(real64) :: shortwave_gpoint_ozone(12) = [2300.0_real64, 170.0_real64, 6.8_real64, &
                                                  1.84_real64, 1.84_real64, &
                                                  0.0_real64, 0.0_real64, 0.0_real64, &
                                                  0.0_real64, 0.0_real64, 0.0_real64, 0.0_real64]
    real(real64) :: shortwave_water_vapor_pressure_exponent = 0.6_real64
    !> Rayleigh layer above the surface: R(mu_0) = a/(1 + b mu_0) for the direct
    !> beam and R* for the diffuse light from below.
    real(real64) :: rayleigh_direct_amplitude(4) = [0.81_real64, 0.63_real64, 0.0_real64, 0.0_real64]
    real(real64) :: rayleigh_direct_slope(4) = [1.85_real64, 8.2_real64, 0.0_real64, 0.0_real64]
    real(real64) :: rayleigh_diffuse_reflectance(4) = [0.38_real64, 0.096_real64, 0.0_real64, 0.0_real64]
    !> Magnification M(mu_0) = a/sqrt(b mu_0^2 + 1) of the direct-beam path (Lacis and Hansen 1974).
    real(real64) :: magnification_numerator = 35.0_real64
    real(real64) :: magnification_quadratic = 1224.0_real64
    !> exp(-x) is taken as zero beyond this optical path.
    real(real64) :: maximum_optical_path = 200.0_real64
    !> Radiative properties of the large-scale and the convective cloud (docs/tendency/band-radiation.md, 6.3).
    !> In the shortwave a cloud is a conservative delta-Eddington layer of visible
    !> optical depth tau and asymmetry factor g: tau = 3.9 is the ISCCP D-series
    !> global mean (Rossow and Schiffer 1999), g = 0.85 that of water droplets of
    !> r_e = 10 um (Hansen and Travis 1974).
    real(real64) :: large_scale_cloud_optical_depth = 3.9_real64
    real(real64) :: convective_cloud_optical_depth = 3.9_real64
    real(real64) :: cloud_asymmetry_factor = 0.85_real64
    real(real64) :: large_scale_cloud_longwave_emissivity = 1.0_real64
    real(real64) :: convective_cloud_longwave_emissivity = 1.0_real64
  end type band_radiation_config

  !> Grey longwave radiation, prescribed ozone shortwave absorption, and the
  !> selectable two-layer-ground or slab-ocean surface energy budget, together
  !> with the calendar and orbit that drive the solar forcing.  With the band
  !> scheme the radiative transfer follows `band` instead of the grey
  !> coefficients; the surface budget and the calendar are shared.
  type, public :: radiation_config
    logical :: enabled = .false.
    integer :: scheme = radiation_scheme_gray
    type(band_radiation_config) :: band
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
    !> Initial mixed-layer temperature T_p + (T_e - T_p) cos^2(phi) of the cases that
    !> advance the ocean as its own tile, as in the SpeedyWeather AquaPlanet ocean
    !> (docs/cases/moist.md).  The legacy single-surface slab ocean keeps T_s = T_N.
    real(real64) :: initial_ocean_equator_temperature = 302.0_real64
    real(real64) :: initial_ocean_pole_temperature = 273.0_real64
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

  !> Large-scale snowfall, its melting in warm layers, and the land snowpack
  !> (docs/tendency/snow.md).  Every temperature is in K.  The snowpack is a
  !> water-equivalent mass per unit land area; no snow collects on sea ice.
  type, public :: snow_config
    logical :: enabled = .false.
    !> Condensate of a layer colder than this (-10 C) is snow.
    real(real64) :: formation_temperature = 263.15_real64
    !> Falling snow melts in a layer warmer than this (5 C), using the excess sensible heat.
    real(real64) :: atmospheric_melting_temperature = 278.15_real64
    !> The snowpack melts with the heat of the land surface above this temperature.
    real(real64) :: surface_melting_temperature = 275.0_real64
    !> S_0 of the snow cover f = S/(S + S_0): 0.05 m water equivalent.
    real(real64) :: masking_water_equivalent = 50.0_real64
    !> Land albedo increase at full snow cover, alpha_L + f*albedo_increase.
    real(real64) :: albedo_increase = 0.4_real64
    real(real64) :: latent_heat_of_fusion = 3.34e5_real64
    !> Initial snowpack on land, kg m^-2 of land area.
    real(real64) :: initial_water_equivalent = 0.0_real64
  end type snow_config

  !> Prescribed, time-independent ocean heat convergence standing in for the
  !> ocean heat transport (docs/tendency/q-flux.md).  The northward transport
  !> (3 sqrt(3)/2) maximum_transport sin(phi) cos(phi)^2 peaks at sin(phi) = 1/sqrt(3).
  type, public :: q_flux_config
    logical :: enabled = .false.
    !> Peak northward ocean heat transport, W.
    real(real64) :: maximum_transport = 1.5e15_real64
  end type q_flux_config

  type, public :: dry_model_physics_config
    type(sea_ice_config) :: sea_ice
    type(q_flux_config) :: q_flux
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
    type(snow_config) :: snow
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
