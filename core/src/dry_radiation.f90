module dry_radiation
  use iso_fortran_env, only: real64
  use sea_ice, only: sea_ice_checks
  use dry_physics_config, only: radiation_config, radiation_orbital_period, radiation_planet_rotation_rate, &
                                radiation_surface_heat_capacity
  use surface_exchange, only: surface_sensible_heat_flux
  use dry_physics_config, only: radiation_scheme_band
  use cloud_diagnostics, only: cloud_layers
  use band_radiation, only: band_column_fluxes, band_longwave_subbands, band_longwave_optics, &
                            band_longwave_downward, band_longwave_upward, band_shortwave, band_surface_emission, &
                            validate_cloud_layers
  implicit none
  private

  real(real64), parameter :: pi = acos(-1.0_real64)

  !> Band-radiation part of a diagnostic sample (docs/tendency/band-radiation.md,
  !> 10.5): global means in W m^-2 and the cloud sub-column fractions, and the grid
  !> fields of the monthly means.  The cloud pressures are the full-level pressure
  !> of each cloud times its fraction, so that their sums over a month divided by
  !> the summed fraction give the mean cloud pressure.  Unallocated and zero with
  !> the grey scheme.
  type, public :: band_radiation_sample
    logical :: enabled = .false.
    real(real64) :: mean_clear_reflected_shortwave = 0.0_real64
    real(real64) :: mean_clear_outgoing_longwave = 0.0_real64
    real(real64) :: mean_window_outgoing_longwave = 0.0_real64
    real(real64) :: mean_surface_incident_shortwave = 0.0_real64
    real(real64) :: mean_atmospheric_shortwave_absorption = 0.0_real64
    real(real64) :: mean_surface_downward_longwave = 0.0_real64
    real(real64) :: mean_surface_upward_longwave = 0.0_real64
    real(real64) :: mean_large_scale_cloud_fraction = 0.0_real64
    real(real64) :: mean_convective_cloud_fraction = 0.0_real64
    real(real64), allocatable :: outgoing_longwave(:, :), reflected_shortwave(:, :)
    real(real64), allocatable :: clear_outgoing_longwave(:, :), clear_reflected_shortwave(:, :)
    real(real64), allocatable :: large_scale_cloud_fraction(:, :), convective_cloud_fraction(:, :)
    real(real64), allocatable :: large_scale_cloud_pressure(:, :), convective_cloud_pressure(:, :)
  end type band_radiation_sample

  !> One instantaneous sample of the global and zonal diagnostics that the
  !> radiation, slab-ocean, moist and land--sea cases aggregate.  The moist fields are zero
  !> and unallocated unless the moist processes are active.
  type, public :: radiation_diagnostics
    type(sea_ice_checks) :: ice_checks
    type(band_radiation_sample) :: band
    real(real64), allocatable :: land_temperature(:, :)
    real(real64), allocatable :: ocean_temperature(:, :)
    real(real64), allocatable :: sea_ice_fraction(:, :)
    real(real64), allocatable :: sea_ice_volume(:, :)
    real(real64), allocatable :: sea_ice_temperature(:, :)
    real(real64), allocatable :: sea_ice_thickness(:, :)
    !> Near-surface air temperature T_N (p_s/p_N)^kappa (K), allocated with the surface tiles.
    real(real64), allocatable :: surface_air_temperature(:, :)
    real(real64) :: sea_ice_area = 0.0_real64
    real(real64) :: sea_ice_total_volume = 0.0_real64
    real(real64) :: mean_sea_ice_thickness = 0.0_real64
    real(real64) :: time_seconds = 0.0_real64
    real(real64) :: mean_atmospheric_temperature = 0.0_real64
    real(real64) :: mean_surface_temperature = 0.0_real64
    real(real64) :: mean_deep_temperature = 0.0_real64
    real(real64) :: mean_land_surface_temperature = 0.0_real64
    real(real64) :: mean_ocean_surface_temperature = 0.0_real64
    real(real64) :: mean_land_precipitation = 0.0_real64
    real(real64) :: mean_ocean_precipitation = 0.0_real64
    real(real64) :: mean_land_evaporation = 0.0_real64
    real(real64) :: mean_ocean_evaporation = 0.0_real64
    real(real64) :: mean_surface_water = 0.0_real64
    real(real64) :: mean_surface_wetness = 0.0_real64
    real(real64) :: mean_runoff = 0.0_real64
    real(real64) :: dry_land_fraction = 0.0_real64
    real(real64) :: mean_water_budget_residual = 0.0_real64
    real(real64) :: maximum_water_budget_residual = 0.0_real64
    !> Snow (docs/tendency/snow.md): global means of the large-scale snowfall reaching the
    !> surface and of the snow melted in the air (kg m^-2 s^-1), and land-area means of the
    !> snowfall, the snowpack (kg m^-2), its cover fraction and its melt (kg m^-2 s^-1).
    real(real64) :: mean_snowfall = 0.0_real64
    real(real64) :: mean_atmospheric_snow_melt = 0.0_real64
    real(real64) :: mean_land_snowfall = 0.0_real64
    real(real64) :: mean_snow_water = 0.0_real64
    real(real64) :: mean_snow_fraction = 0.0_real64
    real(real64) :: mean_snow_melt = 0.0_real64
    real(real64) :: maximum_snow_budget_residual = 0.0_real64
    real(real64) :: mean_kinetic_energy = 0.0_real64
    real(real64) :: mean_surface_pressure = 0.0_real64
    real(real64) :: mean_incoming_shortwave = 0.0_real64
    real(real64) :: mean_reflected_shortwave = 0.0_real64
    real(real64) :: mean_outgoing_longwave = 0.0_real64
    !> Water budget in kg m^-2 s^-1 (precipitation, evaporation) and kg m^-2 (column water).
    real(real64) :: mean_convective_precipitation = 0.0_real64
    real(real64) :: mean_large_scale_precipitation = 0.0_real64
    real(real64) :: mean_evaporation = 0.0_real64
    real(real64) :: mean_latent_heat_flux = 0.0_real64
    real(real64) :: mean_precipitable_water = 0.0_real64
    real(real64) :: mean_signed_column_water = 0.0_real64
    real(real64) :: mean_negative_column_water = 0.0_real64
    !> Effective column cloud cover (docs/tendency/cloud.md), dimensionless.
    real(real64) :: mean_cloud_cover = 0.0_real64
    !> Largest wind speed of the sample over all grid points and levels, and where it occurs.
    real(real64) :: maximum_wind_speed = 0.0_real64
    real(real64) :: maximum_wind_longitude_degrees = 0.0_real64
    real(real64) :: maximum_wind_latitude_degrees = 0.0_real64
    integer :: maximum_wind_level = 0
    real(real64) :: maximum_wind_eta = 0.0_real64
    real(real64), allocatable :: surface_temperature(:, :)
    real(real64), allocatable :: deep_temperature(:, :)
    real(real64), allocatable :: surface_pressure(:, :)
    real(real64), allocatable :: zonal_temperature(:, :)
    real(real64), allocatable :: zonal_u(:, :)
    real(real64), allocatable :: zonal_v(:, :)
    real(real64), allocatable :: zonal_uv(:, :)
    real(real64), allocatable :: zonal_vt(:, :)
    !> Moist grid fields: total precipitation and evaporation (kg m^-2 s^-1), precipitable water (kg m^-2).
    real(real64), allocatable :: precipitation(:, :)
    real(real64), allocatable :: evaporation(:, :)
    real(real64), allocatable :: precipitable_water(:, :)
    real(real64), allocatable :: zonal_humidity(:, :)
    real(real64), allocatable :: zonal_vq(:, :)
    real(real64), allocatable :: cloud_cover(:, :)
    real(real64), allocatable :: surface_water(:, :), surface_wetness(:, :), runoff(:, :)
    !> Snow grid fields, allocated only with snow: snowpack (kg m^-2 of land), cover
    !> fraction, surface snowfall (kg m^-2 s^-1) and land snowmelt (kg m^-2 s^-1 of land).
    real(real64), allocatable :: snow_water(:, :), snow_fraction(:, :), snowfall(:, :), snow_melt(:, :)
  end type radiation_diagnostics

  public :: radiation_downward_column, radiation_upward_column
  public :: radiation_band_downward_column, radiation_band_upward_column
  public :: shortwave_downward_flux, solar_zenith_cosine
  public :: ozone_layer_optical_depth, ozone_layer_fraction, ozone_layer_fractions
  public :: ozone_longwave_layer_optical_depth
  public :: reference_layer_humidity
  public :: gas_longwave_layer_optical_depth
  public :: radiation_tendency
  public :: radiation_calendar_date
  public :: move_radiation_diagnostics

contains

  !> Transfers a diagnostic sample without copying its grid arrays; source arrays become unallocated.
  subroutine move_radiation_diagnostics(source, destination)
    type(radiation_diagnostics), intent(inout) :: source
    type(radiation_diagnostics), intent(out) :: destination

    destination%time_seconds = source%time_seconds
    destination%ice_checks = source%ice_checks
    destination%band = source%band
    destination%mean_atmospheric_temperature = source%mean_atmospheric_temperature
    destination%mean_surface_temperature = source%mean_surface_temperature
    destination%mean_deep_temperature = source%mean_deep_temperature
    destination%mean_land_surface_temperature = source%mean_land_surface_temperature
    destination%mean_ocean_surface_temperature = source%mean_ocean_surface_temperature
    destination%mean_land_precipitation = source%mean_land_precipitation
    destination%mean_ocean_precipitation = source%mean_ocean_precipitation
    destination%mean_land_evaporation = source%mean_land_evaporation
    destination%mean_ocean_evaporation = source%mean_ocean_evaporation
    destination%sea_ice_area = source%sea_ice_area
    destination%sea_ice_total_volume = source%sea_ice_total_volume
    destination%mean_sea_ice_thickness = source%mean_sea_ice_thickness
    if (allocated(source%land_temperature)) call move_alloc(source%land_temperature, destination%land_temperature)
    if (allocated(source%ocean_temperature)) call move_alloc(source%ocean_temperature, destination%ocean_temperature)
    if (allocated(source%sea_ice_fraction)) call move_alloc(source%sea_ice_fraction, destination%sea_ice_fraction)
    if (allocated(source%sea_ice_volume)) call move_alloc(source%sea_ice_volume, destination%sea_ice_volume)
    if (allocated(source%sea_ice_temperature)) call move_alloc(source%sea_ice_temperature, destination%sea_ice_temperature)
    if (allocated(source%sea_ice_thickness)) call move_alloc(source%sea_ice_thickness, destination%sea_ice_thickness)
    if (allocated(source%surface_air_temperature)) &
      call move_alloc(source%surface_air_temperature, destination%surface_air_temperature)
    destination%mean_surface_water = source%mean_surface_water
    destination%mean_surface_wetness = source%mean_surface_wetness
    destination%mean_runoff = source%mean_runoff
    destination%dry_land_fraction = source%dry_land_fraction
    destination%mean_water_budget_residual = source%mean_water_budget_residual
    destination%maximum_water_budget_residual = source%maximum_water_budget_residual
    destination%mean_snowfall = source%mean_snowfall
    destination%mean_atmospheric_snow_melt = source%mean_atmospheric_snow_melt
    destination%mean_land_snowfall = source%mean_land_snowfall
    destination%mean_snow_water = source%mean_snow_water
    destination%mean_snow_fraction = source%mean_snow_fraction
    destination%mean_snow_melt = source%mean_snow_melt
    destination%maximum_snow_budget_residual = source%maximum_snow_budget_residual
    destination%mean_kinetic_energy = source%mean_kinetic_energy
    destination%mean_surface_pressure = source%mean_surface_pressure
    destination%mean_incoming_shortwave = source%mean_incoming_shortwave
    destination%mean_reflected_shortwave = source%mean_reflected_shortwave
    destination%mean_outgoing_longwave = source%mean_outgoing_longwave
    destination%mean_convective_precipitation = source%mean_convective_precipitation
    destination%mean_large_scale_precipitation = source%mean_large_scale_precipitation
    destination%mean_evaporation = source%mean_evaporation
    destination%mean_latent_heat_flux = source%mean_latent_heat_flux
    destination%mean_precipitable_water = source%mean_precipitable_water
    destination%mean_signed_column_water = source%mean_signed_column_water
    destination%mean_negative_column_water = source%mean_negative_column_water
    destination%mean_cloud_cover = source%mean_cloud_cover
    destination%maximum_wind_speed = source%maximum_wind_speed
    destination%maximum_wind_longitude_degrees = source%maximum_wind_longitude_degrees
    destination%maximum_wind_latitude_degrees = source%maximum_wind_latitude_degrees
    destination%maximum_wind_level = source%maximum_wind_level
    destination%maximum_wind_eta = source%maximum_wind_eta
    if (allocated(source%precipitation)) call move_alloc(source%precipitation, destination%precipitation)
    if (allocated(source%evaporation)) call move_alloc(source%evaporation, destination%evaporation)
    if (allocated(source%precipitable_water)) call move_alloc(source%precipitable_water, destination%precipitable_water)
    if (allocated(source%zonal_humidity)) call move_alloc(source%zonal_humidity, destination%zonal_humidity)
    if (allocated(source%zonal_vq)) call move_alloc(source%zonal_vq, destination%zonal_vq)
    if (allocated(source%cloud_cover)) call move_alloc(source%cloud_cover, destination%cloud_cover)
    if (allocated(source%surface_water)) call move_alloc(source%surface_water, destination%surface_water)
    if (allocated(source%surface_wetness)) call move_alloc(source%surface_wetness, destination%surface_wetness)
    if (allocated(source%runoff)) call move_alloc(source%runoff, destination%runoff)
    if (allocated(source%snow_water)) call move_alloc(source%snow_water, destination%snow_water)
    if (allocated(source%snow_fraction)) call move_alloc(source%snow_fraction, destination%snow_fraction)
    if (allocated(source%snowfall)) call move_alloc(source%snowfall, destination%snowfall)
    if (allocated(source%snow_melt)) call move_alloc(source%snow_melt, destination%snow_melt)
    if (allocated(source%surface_temperature)) call move_alloc(source%surface_temperature, destination%surface_temperature)
    if (allocated(source%deep_temperature)) call move_alloc(source%deep_temperature, destination%deep_temperature)
    if (allocated(source%surface_pressure)) call move_alloc(source%surface_pressure, destination%surface_pressure)
    if (allocated(source%zonal_temperature)) call move_alloc(source%zonal_temperature, destination%zonal_temperature)
    if (allocated(source%zonal_u)) call move_alloc(source%zonal_u, destination%zonal_u)
    if (allocated(source%zonal_v)) call move_alloc(source%zonal_v, destination%zonal_v)
    if (allocated(source%zonal_uv)) call move_alloc(source%zonal_uv, destination%zonal_uv)
    if (allocated(source%zonal_vt)) call move_alloc(source%zonal_vt, destination%zonal_vt)
  end subroutine move_radiation_diagnostics

  real(real64) function shortwave_downward_flux(config, sin_latitude, longitude, time_seconds) result(flux)
    type(radiation_config), intent(in) :: config
    real(real64), intent(in) :: sin_latitude, longitude, time_seconds

    flux = config%solar_constant*max(0.0_real64, solar_zenith_cosine(config, sin_latitude, longitude, time_seconds))
  end function shortwave_downward_flux

  !> Cosine of the solar zenith angle mu_0, negative on the night side.
  real(real64) function solar_zenith_cosine(config, sin_latitude, longitude, time_seconds) result(cos_zenith)
    type(radiation_config), intent(in) :: config
    real(real64), intent(in) :: sin_latitude, longitude, time_seconds
    real(real64) :: orbital_longitude, solar_right_ascension, orbital_period, axial_tilt
    real(real64) :: sin_declination, cos_declination, cos_latitude, hour_angle

    if (abs(sin_latitude) > 1.0_real64) error stop 'radiation latitude sine is outside [-1,1]'
    orbital_period = radiation_orbital_period(config)
    axial_tilt = config%axial_tilt
    orbital_longitude = modulo(2.0_real64*pi*time_seconds/orbital_period, 2.0_real64*pi)
    solar_right_ascension = atan2(cos(axial_tilt)*sin(orbital_longitude), cos(orbital_longitude))
    sin_declination = sin(axial_tilt)*sin(orbital_longitude)
    cos_declination = sqrt(max(0.0_real64, 1.0_real64 - sin_declination**2))
    cos_latitude = sqrt(max(0.0_real64, 1.0_real64 - sin_latitude**2))
    hour_angle = radiation_planet_rotation_rate(config)*time_seconds + longitude - solar_right_ascension
    cos_zenith = sin_latitude*sin_declination + cos_latitude*cos_declination*cos(hour_angle)
  end function solar_zenith_cosine

  !> 1/(erf(x_upper) - erf(x_lower)) of the whole prescribed ozone column, the
  !> normalization every layer fraction shares.
  pure real(real64) function ozone_profile_normalization(config) result(normalization)
    type(radiation_config), intent(in) :: config
    real(real64) :: x_profile_lower, x_profile_upper

    x_profile_lower = log(config%ozone_pressure_lower_bound/config%ozone_peak_pressure)/ &
      (sqrt(2.0_real64)*config%ozone_log_pressure_width)
    x_profile_upper = log(config%ozone_pressure_upper_bound/config%ozone_peak_pressure)/ &
      (sqrt(2.0_real64)*config%ozone_log_pressure_width)
    normalization = 1.0_real64/(erf(x_profile_upper) - erf(x_profile_lower))
  end function ozone_profile_normalization

  !> Fraction of the prescribed ozone column within one pressure layer, given the
  !> column normalization above.
  pure real(real64) function ozone_layer_fraction_normalized(config, pressure_top, pressure_bottom, &
                                                             normalization) result(fraction)
    type(radiation_config), intent(in) :: config
    real(real64), intent(in) :: pressure_top, pressure_bottom, normalization
    real(real64) :: pressure_lower, pressure_upper, x_lower, x_upper

    pressure_lower = max(pressure_top, config%ozone_pressure_lower_bound)
    pressure_upper = min(pressure_bottom, config%ozone_pressure_upper_bound)
    if (pressure_upper <= pressure_lower) then
      fraction = 0.0_real64
      return
    end if
    x_lower = log(pressure_lower/config%ozone_peak_pressure)/(sqrt(2.0_real64)*config%ozone_log_pressure_width)
    x_upper = log(pressure_upper/config%ozone_peak_pressure)/(sqrt(2.0_real64)*config%ozone_log_pressure_width)
    fraction = (erf(x_upper) - erf(x_lower))*normalization
  end function ozone_layer_fraction_normalized

  !> Fractions of the prescribed ozone column in every layer of a column, equal
  !> bit for bit to ozone_layer_fraction_normalized layer by layer.  Each
  !> interface shares its erf between the two layers it bounds, and an interface
  !> outside the ozone layer takes the erf of the bound it is clipped to, which
  !> the normalization needs anyway; only the interfaces within the ozone layer
  !> call erf.
  pure subroutine ozone_layer_fractions(config, pressure_half, fraction)
    type(radiation_config), intent(in) :: config
    real(real64), intent(in) :: pressure_half(0:)
    real(real64), intent(out) :: fraction(:)
    real(real64) :: interface_erf(0:size(fraction)), normalization, lower_erf, upper_erf
    integer :: k

    lower_erf = erf(log(config%ozone_pressure_lower_bound/config%ozone_peak_pressure)/ &
                    (sqrt(2.0_real64)*config%ozone_log_pressure_width))
    upper_erf = erf(log(config%ozone_pressure_upper_bound/config%ozone_peak_pressure)/ &
                    (sqrt(2.0_real64)*config%ozone_log_pressure_width))
    normalization = 1.0_real64/(upper_erf - lower_erf)
    do k = 0, size(fraction)
      if (pressure_half(k) <= config%ozone_pressure_lower_bound) then
        interface_erf(k) = lower_erf
      else if (pressure_half(k) >= config%ozone_pressure_upper_bound) then
        interface_erf(k) = upper_erf
      else
        interface_erf(k) = erf(log(pressure_half(k)/config%ozone_peak_pressure)/ &
                               (sqrt(2.0_real64)*config%ozone_log_pressure_width))
      end if
    end do
    do k = 1, size(fraction)
      if (min(pressure_half(k), config%ozone_pressure_upper_bound) <= &
          max(pressure_half(k - 1), config%ozone_pressure_lower_bound)) then
        fraction(k) = 0.0_real64
      else
        fraction(k) = (interface_erf(k) - interface_erf(k - 1))*normalization
      end if
    end do
  end subroutine ozone_layer_fractions

  !> Fraction of the prescribed ozone column within one pressure layer.
  pure real(real64) function ozone_layer_fraction(config, pressure_top, pressure_bottom) result(fraction)
    type(radiation_config), intent(in) :: config
    real(real64), intent(in) :: pressure_top, pressure_bottom

    fraction = ozone_layer_fraction_normalized(config, pressure_top, pressure_bottom, &
                                               ozone_profile_normalization(config))
  end function ozone_layer_fraction

  !> Vertical UV optical depth of the prescribed ozone profile within one pressure layer.
  pure real(real64) function ozone_layer_optical_depth(config, pressure_top, pressure_bottom) result(optical_depth)
    type(radiation_config), intent(in) :: config
    real(real64), intent(in) :: pressure_top, pressure_bottom

    optical_depth = config%ozone_shortwave_optical_depth*ozone_layer_fraction(config, pressure_top, pressure_bottom)
  end function ozone_layer_optical_depth

  !> Vertical longwave optical depth added by ozone within one pressure layer.
  pure real(real64) function ozone_longwave_layer_optical_depth(config, pressure_top, pressure_bottom) &
    result(optical_depth)
    type(radiation_config), intent(in) :: config
    real(real64), intent(in) :: pressure_top, pressure_bottom

    optical_depth = config%ozone_longwave_optical_depth*ozone_layer_fraction(config, pressure_top, pressure_bottom)
  end function ozone_longwave_layer_optical_depth

  !> Mass-weighted mean of the reference humidity q_ref = q_0 (p/p_s)^3 over one
  !> layer, q_0 (p_bottom^4 - p_top^4)/(4 p_s^3 dp).  The cases without prognostic
  !> water vapour give this to the longwave radiation in place of q.
  pure real(real64) function reference_layer_humidity(config, pressure_top, pressure_bottom, surface_pressure) &
    result(humidity)
    type(radiation_config), intent(in) :: config
    real(real64), intent(in) :: pressure_top, pressure_bottom, surface_pressure

    humidity = config%longwave_reference_surface_humidity*(pressure_bottom**4 - pressure_top**4)/ &
      (4.0_real64*surface_pressure**3*(pressure_bottom - pressure_top))
  end function reference_layer_humidity

  !> Longwave optical depth (a mu + b q^+) dp/p_0 of the gases other than ozone
  !> within one layer; negative humidity counts as zero.
  pure real(real64) function gas_longwave_layer_optical_depth(config, pressure_top, pressure_bottom, humidity) &
    result(optical_depth)
    type(radiation_config), intent(in) :: config
    real(real64), intent(in) :: pressure_top, pressure_bottom, humidity

    optical_depth = (config%longwave_well_mixed_optical_depth*config%longwave_well_mixed_scaling + &
      config%longwave_water_vapor_optical_depth*max(humidity, 0.0_real64))* &
      (pressure_bottom - pressure_top)/config%longwave_reference_pressure
  end function gas_longwave_layer_optical_depth

  subroutine radiation_calendar_date(config, time_seconds, year, month, day, seconds_of_day)
    type(radiation_config), intent(in) :: config
    real(real64), intent(in) :: time_seconds
    integer, intent(out) :: year, month, day
    real(real64), intent(out) :: seconds_of_day
    integer :: elapsed_days, zero_based_month, days_per_month, months_per_year
    real(real64) :: solar_day

    solar_day = config%solar_day
    days_per_month = config%days_per_month
    months_per_year = config%months_per_year

    if (time_seconds < 0.0_real64) error stop 'radiation calendar time must be nonnegative'
    elapsed_days = floor(time_seconds/solar_day)
    zero_based_month = 3 + elapsed_days/days_per_month
    year = 1 + zero_based_month/months_per_year
    month = modulo(zero_based_month, months_per_year) + 1
    day = modulo(elapsed_days, days_per_month) + 1
    seconds_of_day = modulo(time_seconds, solar_day)
  end subroutine radiation_calendar_date

  !> Column radiation, surface fluxes and the selected surface energy budget.
  !>
  !> The caller supplies one internally consistent RAW-filtered previous-time column for every
  !> state-dependent term.  Solar geometry is evaluated at time_seconds because it is prescribed
  !> rather than prognostic.  The sensible heat flux compares the surface with the lowest level
  !> extrapolated dry-adiabatically to the surface pressure, so a neutral column exchanges none.
  !> An optional latent heat flux L E (W m^-2, upward) is taken from the surface budget; the
  !> matching water vapour source is added to the atmosphere by the evaporation tendency.
  !> The grey longwave optical depth of each layer is (a mu + b q^+) dp/p_0 plus ozone, with q
  !> the optional specific humidity column; without it the fixed reference humidity
  !> q_0 (p/p_s)^3 stands in for the water vapour.
  !> The optional cloud cover C reflects C times the cloud albedo of the downward shortwave
  !> below the ozone layer once (docs/tendency/shortwave-radiation.md); the rest reaches the
  !> surface.  Clouds absorb nothing and do not enter the longwave.  reflected_shortwave is
  !> the cloud plus surface reflection, which leaves at the top of the atmosphere unabsorbed.
  !> Surface-independent optical properties, downward LW and SW absorption.
  subroutine radiation_downward_column(config, pressure_half, temperature, sin_latitude, longitude, time_seconds, &
                                       transmission, emission, downward_longwave, shortwave_downward, &
                                       incoming_shortwave, shortwave_heating, specific_humidity)
    type(radiation_config), intent(in) :: config
    real(real64), intent(in) :: pressure_half(0:), temperature(:), sin_latitude, longitude, time_seconds
    real(real64), intent(in), optional :: specific_humidity(:)
    real(real64), intent(out) :: transmission(:), emission(:), downward_longwave(0:)
    real(real64), intent(out) :: shortwave_downward, incoming_shortwave, shortwave_heating(:)
    real(real64) :: ozone_fraction(size(temperature))
    real(real64) :: pressure_thickness, layer_humidity, layer_longwave_optical_depth
    real(real64) :: layer_shortwave_optical_depth, shortwave_transmission, shortwave_absorbed
    real(real64) :: ultraviolet_downward, non_ultraviolet_downward
    real(real64) :: stefan_boltzmann_constant, dry_gravity_acceleration, dry_air_specific_heat
    integer :: number_of_levels, k
    number_of_levels = size(temperature)
    if (number_of_levels < 1 .or. size(pressure_half) /= number_of_levels + 1) &
      error stop 'radiation column has inconsistent vertical dimensions'
    if (any(temperature <= 0.0_real64)) error stop 'radiation temperature must be positive'
    if (present(specific_humidity)) then
      if (size(specific_humidity) /= number_of_levels) error stop 'radiation humidity shape mismatch'
    end if
    stefan_boltzmann_constant = config%stefan_boltzmann_constant
    dry_gravity_acceleration = config%gravity_acceleration
    dry_air_specific_heat = config%dry_air_specific_heat
    shortwave_heating = 0.0_real64
    ! The ozone fraction of a layer is shared by its longwave and shortwave optical depths.
    call ozone_layer_fractions(config, pressure_half, ozone_fraction)
    do k = 1, number_of_levels
      pressure_thickness = pressure_half(k) - pressure_half(k - 1)
      if (pressure_thickness <= 0.0_real64) then
        error stop 'radiation pressures must increase downward'
      end if
      if (present(specific_humidity)) then
        layer_humidity = specific_humidity(k)
      else
        layer_humidity = reference_layer_humidity(config, pressure_half(k - 1), pressure_half(k), &
                                                  pressure_half(number_of_levels))
      end if
      layer_longwave_optical_depth = &
        gas_longwave_layer_optical_depth(config, pressure_half(k - 1), pressure_half(k), layer_humidity) + &
        config%ozone_longwave_optical_depth*ozone_fraction(k)
      transmission(k) = exp(-layer_longwave_optical_depth)
      emission(k) = (1.0_real64 - transmission(k))*stefan_boltzmann_constant*temperature(k)**4
    end do

    downward_longwave(0) = 0.0_real64
    do k = 1, number_of_levels
      downward_longwave(k) = transmission(k)*downward_longwave(k - 1) + emission(k)
    end do
    incoming_shortwave = shortwave_downward_flux(config, sin_latitude, longitude, time_seconds)
    shortwave_downward = incoming_shortwave
    if (incoming_shortwave > 0.0_real64) then
      non_ultraviolet_downward = (1.0_real64 - config%ultraviolet_shortwave_fraction)*incoming_shortwave
      ultraviolet_downward = config%ultraviolet_shortwave_fraction*incoming_shortwave
      do k = 1, number_of_levels
        layer_shortwave_optical_depth = config%ozone_shortwave_optical_depth*ozone_fraction(k)
        shortwave_transmission = exp(-layer_shortwave_optical_depth)
        shortwave_absorbed = ultraviolet_downward*(1.0_real64 - shortwave_transmission)
        pressure_thickness = pressure_half(k) - pressure_half(k - 1)
        shortwave_heating(k) = shortwave_heating(k) + &
          dry_gravity_acceleration/(dry_air_specific_heat*pressure_thickness)*shortwave_absorbed
        ultraviolet_downward = shortwave_transmission*ultraviolet_downward
      end do
      shortwave_downward = non_ultraviolet_downward + ultraviolet_downward
    end if

  end subroutine radiation_downward_column

  !> The same area-averaged upward LW boundary is used by the air and surface budgets.
  subroutine radiation_upward_column(config, pressure_half, transmission, emission, downward, surface_upward, &
                                     heating, outgoing)
    type(radiation_config), intent(in) :: config
    real(real64), intent(in) :: pressure_half(0:), transmission(:), emission(:), downward(0:), surface_upward
    real(real64), intent(out) :: heating(:), outgoing
    real(real64) :: upward(0:size(emission)), net(0:size(emission))
    integer :: k, levels
    levels = size(emission)
    upward(levels) = surface_upward
    do k = levels, 1, -1
      upward(k - 1) = transmission(k)*upward(k) + emission(k)
    end do
    net = upward - downward
    do k = 1, levels
      heating(k) = config%gravity_acceleration/(config%dry_air_specific_heat* &
        (pressure_half(k) - pressure_half(k - 1)))*(net(k) - net(k - 1))
    end do
    outgoing = upward(0)
  end subroutine radiation_upward_column

  !> Downward stage of the band radiation (docs/tendency/band-radiation.md): the
  !> longwave transmittances and sources of the clear sub-column, the area-mean
  !> downward longwave, and the whole shortwave for a surface of area-mean
  !> albedo `surface_albedo`.  Without `specific_humidity` the fixed reference
  !> humidity stands in for the water vapour, as in the grey scheme.
  subroutine radiation_band_downward_column(config, pressure_half, temperature, sin_latitude, longitude, &
                                            time_seconds, surface_albedo, clouds, transmission, source, &
                                            downward_longwave, shortwave_heating, fluxes, specific_humidity)
    type(radiation_config), intent(in) :: config
    real(real64), intent(in) :: pressure_half(0:), temperature(:), sin_latitude, longitude, time_seconds
    real(real64), intent(in) :: surface_albedo
    type(cloud_layers), intent(in) :: clouds
    real(real64), intent(out) :: transmission(:, :), source(:, :), downward_longwave(0:), shortwave_heating(:)
    type(band_column_fluxes), intent(out) :: fluxes
    real(real64), intent(in), optional :: specific_humidity(:)
    real(real64) :: humidity(size(temperature)), ozone_fraction(size(temperature))
    real(real64) :: cos_zenith, incoming
    integer :: levels, k

    levels = size(temperature)
    if (levels < 1 .or. size(pressure_half) /= levels + 1) &
      error stop 'band radiation column has inconsistent vertical dimensions'
    if (size(transmission, 1) /= band_longwave_subbands .or. size(transmission, 2) /= levels .or. &
        any(shape(source) /= shape(transmission)) .or. size(shortwave_heating) /= levels) &
      error stop 'band radiation work arrays have inconsistent shapes'
    if (any(temperature <= 0.0_real64)) error stop 'radiation temperature must be positive'
    if (surface_albedo < 0.0_real64 .or. surface_albedo > 1.0_real64) error stop 'surface albedo is outside [0,1]'
    call validate_cloud_layers(clouds, levels)
    do k = 1, levels
      if (pressure_half(k) <= pressure_half(k - 1)) error stop 'radiation pressures must increase downward'
    end do
    call ozone_layer_fractions(config, pressure_half, ozone_fraction)
    if (present(specific_humidity)) then
      if (size(specific_humidity) /= levels) error stop 'radiation humidity shape mismatch'
      humidity = max(specific_humidity, 0.0_real64)
    else
      do k = 1, levels
        humidity(k) = reference_layer_humidity(config, pressure_half(k - 1), pressure_half(k), pressure_half(levels))
      end do
    end if
    call band_longwave_optics(config%band, config%stefan_boltzmann_constant, config%gravity_acceleration, &
                              pressure_half, temperature, humidity, ozone_fraction, transmission, source)
    call band_longwave_downward(config%band, clouds, transmission, source, downward_longwave)
    cos_zenith = solar_zenith_cosine(config, sin_latitude, longitude, time_seconds)
    incoming = config%solar_constant*max(0.0_real64, cos_zenith)
    call band_shortwave(config%band, config%gravity_acceleration, config%dry_air_specific_heat, pressure_half, &
                        humidity, ozone_fraction, incoming, cos_zenith, surface_albedo, clouds, &
                        shortwave_heating, fluxes)
    fluxes%surface_downward_longwave = downward_longwave(levels)
  end subroutine radiation_band_downward_column

  !> Upward stage of the band radiation: `surface_emission` is the area-weighted
  !> black-body emission of the surface in each longwave sub-band.  Returns the
  !> longwave heating and completes the longwave entries of `fluxes`.
  subroutine radiation_band_upward_column(config, pressure_half, clouds, transmission, source, downward_longwave, &
                                          surface_emission, heating, fluxes)
    type(radiation_config), intent(in) :: config
    real(real64), intent(in) :: pressure_half(0:), transmission(:, :), source(:, :), downward_longwave(0:)
    type(cloud_layers), intent(in) :: clouds
    real(real64), intent(in) :: surface_emission(:)
    real(real64), intent(out) :: heating(:)
    type(band_column_fluxes), intent(inout) :: fluxes

    if (size(surface_emission) /= band_longwave_subbands) error stop 'band surface emission has the wrong length'
    call band_longwave_upward(config%band, clouds, config%gravity_acceleration, config%dry_air_specific_heat, &
                              pressure_half, transmission, source, downward_longwave, surface_emission, heating, fluxes)
  end subroutine radiation_band_upward_column

  subroutine radiation_tendency(config, pressure_half, temperature, surface_temperature, &
                                deep_temperature, lowest_u, lowest_v, sin_latitude, &
                                longitude, time_seconds, temperature_tendency, &
                                surface_temperature_tendency, deep_temperature_tendency, &
                                incoming_shortwave, reflected_shortwave, outgoing_longwave, &
                                latent_heat_flux, specific_humidity, cloud_cover, land_fraction, clouds, &
                                band_fluxes)
    type(radiation_config), intent(in) :: config
    real(real64), intent(in) :: pressure_half(0:), temperature(:)
    real(real64), intent(in) :: surface_temperature, deep_temperature
    real(real64), intent(in) :: lowest_u, lowest_v
    real(real64), intent(in) :: sin_latitude, longitude, time_seconds
    real(real64), intent(out) :: temperature_tendency(:)
    real(real64), intent(out) :: surface_temperature_tendency, deep_temperature_tendency
    real(real64), intent(out) :: incoming_shortwave, reflected_shortwave, outgoing_longwave
    real(real64), intent(in), optional :: latent_heat_flux
    real(real64), intent(in), optional :: specific_humidity(:)
    real(real64), intent(in), optional :: cloud_cover
    !> Fraction of the grid cell covered by land.  It is required by the mixed
    !> surface path and absent on the legacy ground/slab paths.
    real(real64), intent(in), optional :: land_fraction
    !> Cloud sub-columns of the band scheme, which ignores `cloud_cover`, and the
    !> band column fluxes for the diagnostics.
    type(cloud_layers), intent(in), optional :: clouds
    type(band_column_fluxes), intent(out), optional :: band_fluxes
    type(cloud_layers) :: column_clouds
    type(band_column_fluxes) :: fluxes
    real(real64) :: band_transmission(band_longwave_subbands, size(temperature))
    real(real64) :: band_source(band_longwave_subbands, size(temperature))
    real(real64) :: upward_longwave(0:size(temperature)), downward_longwave(0:size(temperature))
    real(real64) :: transmission(size(temperature)), emission(size(temperature))
    real(real64) :: sw_heating(size(temperature))
    real(real64) :: pressure_thickness, layer_longwave_optical_depth, layer_humidity, sensible_heat, surface_deep_heat
    real(real64) :: layer_shortwave_optical_depth, shortwave_transmission
    real(real64) :: shortwave_downward, ultraviolet_downward, non_ultraviolet_downward, shortwave_absorbed
    real(real64) :: cloud_reflected, surface_downward, surface_reflected, cloud
    real(real64) :: stefan_boltzmann_constant, dry_gravity_acceleration, dry_air_specific_heat
    real(real64) :: surface_heat_capacity, surface_latent_heat, surface_albedo
    real(real64) :: ocean_heat_capacity, ground_exchange, land
    real(real64) :: ozone_fraction(size(temperature)), ozone_normalization
    integer :: k, number_of_levels

    surface_latent_heat = 0.0_real64
    if (present(latent_heat_flux)) surface_latent_heat = latent_heat_flux
    cloud = 0.0_real64
    if (present(cloud_cover)) then
      if (cloud_cover < 0.0_real64 .or. cloud_cover > 1.0_real64) error stop 'radiation cloud cover is outside [0,1]'
      cloud = cloud_cover
    end if
    stefan_boltzmann_constant = config%stefan_boltzmann_constant
    dry_gravity_acceleration = config%gravity_acceleration
    dry_air_specific_heat = config%dry_air_specific_heat
    if (config%land_sea_mixing_enabled) then
      if (.not. present(land_fraction)) error stop 'mixed land-sea radiation requires land fraction'
      if (land_fraction < 0.0_real64 .or. land_fraction > 1.0_real64) then
        error stop 'radiation land fraction is outside [0,1]'
      end if
      land = land_fraction
      ocean_heat_capacity = config%seawater_density*config%seawater_specific_heat*config%slab_ocean_depth
      surface_heat_capacity = land*config%surface_heat_capacity + (1.0_real64 - land)*ocean_heat_capacity
      surface_albedo = land*config%land_shortwave_albedo + &
        (1.0_real64 - land)*config%ocean_shortwave_albedo
      ground_exchange = land*config%ground_exchange_coefficient
    else
      surface_heat_capacity = radiation_surface_heat_capacity(config)
      surface_albedo = config%surface_shortwave_albedo
      if (config%slab_ocean_enabled) then
        ground_exchange = 0.0_real64
      else
        ground_exchange = config%ground_exchange_coefficient
      end if
    end if
    number_of_levels = size(temperature)
    if (number_of_levels < 1 .or. size(pressure_half) /= number_of_levels + 1 .or. &
        size(temperature_tendency) /= number_of_levels) then
      error stop 'radiation column has inconsistent vertical dimensions'
    end if
    if (pressure_half(number_of_levels) <= pressure_half(0) .or. any(temperature <= 0.0_real64) .or. &
        surface_temperature <= 0.0_real64 .or. &
        ((config%land_sea_mixing_enabled .or. .not. config%slab_ocean_enabled) .and. &
         deep_temperature <= 0.0_real64)) then
      error stop 'radiation column contains a nonphysical state'
    end if
    if (surface_heat_capacity <= 0.0_real64) then
      error stop 'radiation surface heat capacity must be positive'
    end if
    if (config%longwave_reference_pressure <= 0.0_real64) then
      error stop 'radiation longwave reference pressure must be positive'
    end if
    if (present(specific_humidity)) then
      if (size(specific_humidity) /= number_of_levels) then
        error stop 'radiation humidity column has an inconsistent length'
      end if
    end if

    if (config%scheme == radiation_scheme_band) then
      column_clouds = cloud_layers()
      if (present(clouds)) then
        column_clouds = clouds
      else if (cloud > 0.0_real64) then
        error stop 'band radiation takes the cloud sub-columns, not the effective cloud cover'
      end if
      call radiation_band_downward_column(config, pressure_half, temperature, sin_latitude, longitude, time_seconds, &
        surface_albedo, column_clouds, band_transmission, band_source, downward_longwave, sw_heating, fluxes, &
        specific_humidity)
      call radiation_band_upward_column(config, pressure_half, column_clouds, band_transmission, band_source, &
        downward_longwave, band_surface_emission(config%band, stefan_boltzmann_constant, surface_temperature), &
        temperature_tendency, fluxes)
      temperature_tendency = temperature_tendency + sw_heating
      sensible_heat = surface_sensible_heat_flux(config, pressure_half, temperature(number_of_levels), &
                                                 surface_temperature, lowest_u, lowest_v)
      pressure_thickness = pressure_half(number_of_levels) - pressure_half(number_of_levels - 1)
      temperature_tendency(number_of_levels) = temperature_tendency(number_of_levels) + &
        dry_gravity_acceleration/(dry_air_specific_heat*pressure_thickness)*sensible_heat
      incoming_shortwave = fluxes%incoming_shortwave
      reflected_shortwave = fluxes%reflected_shortwave
      outgoing_longwave = fluxes%outgoing_longwave
      surface_deep_heat = ground_exchange*(surface_temperature - deep_temperature)
      ! The atmosphere sees the sub-band sum of sigma T_s^4, equal to it within rounding.
      surface_temperature_tendency = ((1.0_real64 - surface_albedo)*fluxes%surface_incident_shortwave + &
        downward_longwave(number_of_levels) - stefan_boltzmann_constant*surface_temperature**4 - &
        surface_deep_heat - sensible_heat - surface_latent_heat)/surface_heat_capacity
      if (.not. config%land_sea_mixing_enabled .and. config%slab_ocean_enabled) then
        deep_temperature_tendency = 0.0_real64
      else
        deep_temperature_tendency = surface_deep_heat/config%deep_ground_heat_capacity
      end if
      if (present(band_fluxes)) band_fluxes = fluxes
      return
    end if

    call radiation_downward_column(config, pressure_half, temperature, sin_latitude, longitude, time_seconds, &
      transmission, emission, downward_longwave, shortwave_downward, incoming_shortwave, sw_heating, specific_humidity)
    upward_longwave(number_of_levels) = stefan_boltzmann_constant*surface_temperature**4
    call radiation_upward_column(config, pressure_half, transmission, emission, downward_longwave, &
      upward_longwave(number_of_levels), temperature_tendency, outgoing_longwave)
    upward_longwave(0) = outgoing_longwave
    temperature_tendency = temperature_tendency + sw_heating

    sensible_heat = surface_sensible_heat_flux(config, pressure_half, temperature(number_of_levels), &
                                               surface_temperature, lowest_u, lowest_v)
    pressure_thickness = pressure_half(number_of_levels) - pressure_half(number_of_levels - 1)
    temperature_tendency(number_of_levels) = temperature_tendency(number_of_levels) + &
      dry_gravity_acceleration/(dry_air_specific_heat*pressure_thickness)*sensible_heat

    ! The cloud reflects once below the ozone layer; with C = 0 these reduce bit for bit to the
    ! cloud-free expressions.
    cloud_reflected = cloud*config%cloud_shortwave_albedo*shortwave_downward
    surface_downward = (1.0_real64 - cloud*config%cloud_shortwave_albedo)*shortwave_downward
    surface_reflected = surface_albedo*surface_downward
    reflected_shortwave = cloud_reflected + surface_reflected
    outgoing_longwave = upward_longwave(0)
    surface_deep_heat = ground_exchange*(surface_temperature - deep_temperature)
    ! The surface loses exactly the upward longwave the lowest layer sees, so the column budget closes.
    surface_temperature_tendency = (surface_downward - surface_reflected + &
      downward_longwave(number_of_levels) - upward_longwave(number_of_levels) - &
      surface_deep_heat - sensible_heat - surface_latent_heat)/surface_heat_capacity
    if (.not. config%land_sea_mixing_enabled .and. config%slab_ocean_enabled) then
      deep_temperature_tendency = 0.0_real64
    else
      deep_temperature_tendency = surface_deep_heat/config%deep_ground_heat_capacity
    end if
  end subroutine radiation_tendency

end module dry_radiation
