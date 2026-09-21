module dry_radiation
  use iso_fortran_env, only: real64
  use dry_physics_config, only: radiation_config, radiation_orbital_period, radiation_planet_rotation_rate, &
                                radiation_surface_heat_capacity
  use surface_exchange, only: surface_sensible_heat_flux
  implicit none
  private

  real(real64), parameter :: pi = acos(-1.0_real64)

  !> One instantaneous sample of the global and zonal diagnostics that the
  !> radiation, slab-ocean and moist cases aggregate.  The moist fields are zero
  !> and unallocated unless the moist processes are active.
  type, public :: radiation_diagnostics
    real(real64) :: time_seconds = 0.0_real64
    real(real64) :: mean_atmospheric_temperature = 0.0_real64
    real(real64) :: mean_surface_temperature = 0.0_real64
    real(real64) :: mean_deep_temperature = 0.0_real64
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
    !> Largest wind speed of the sample over all grid points and levels, and where it occurs.
    real(real64) :: maximum_wind_speed = 0.0_real64
    real(real64) :: maximum_wind_longitude_degrees = 0.0_real64
    real(real64) :: maximum_wind_latitude_degrees = 0.0_real64
    integer :: maximum_wind_level = 0
    real(real64) :: maximum_wind_eta = 0.0_real64
    real(real64), allocatable :: surface_temperature(:, :)
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
  end type radiation_diagnostics

  public :: shortwave_downward_flux
  public :: ozone_layer_optical_depth
  public :: ozone_longwave_layer_optical_depth
  public :: radiation_tendency
  public :: radiation_calendar_date
  public :: move_radiation_diagnostics

contains

  !> Transfers a diagnostic sample without copying its grid arrays; source arrays become unallocated.
  subroutine move_radiation_diagnostics(source, destination)
    type(radiation_diagnostics), intent(inout) :: source
    type(radiation_diagnostics), intent(out) :: destination

    destination%time_seconds = source%time_seconds
    destination%mean_atmospheric_temperature = source%mean_atmospheric_temperature
    destination%mean_surface_temperature = source%mean_surface_temperature
    destination%mean_deep_temperature = source%mean_deep_temperature
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
    if (allocated(source%surface_temperature)) call move_alloc(source%surface_temperature, destination%surface_temperature)
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
    real(real64) :: orbital_longitude, solar_right_ascension, orbital_period, axial_tilt
    real(real64) :: sin_declination, cos_declination, cos_latitude, hour_angle, cos_zenith

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
    flux = config%solar_constant*max(0.0_real64, cos_zenith)
  end function shortwave_downward_flux

  !> Fraction of the prescribed ozone column within one pressure layer.
  pure real(real64) function ozone_layer_fraction(config, pressure_top, pressure_bottom) result(fraction)
    type(radiation_config), intent(in) :: config
    real(real64), intent(in) :: pressure_top, pressure_bottom
    real(real64) :: pressure_lower, pressure_upper, x_lower, x_upper
    real(real64) :: x_profile_lower, x_profile_upper
    real(real64) :: ozone_pressure_lower_bound, ozone_pressure_upper_bound
    real(real64) :: ozone_peak_pressure, ozone_log_pressure_width

    ozone_pressure_lower_bound = config%ozone_pressure_lower_bound
    ozone_pressure_upper_bound = config%ozone_pressure_upper_bound
    ozone_peak_pressure = config%ozone_peak_pressure
    ozone_log_pressure_width = config%ozone_log_pressure_width
    pressure_lower = max(pressure_top, ozone_pressure_lower_bound)
    pressure_upper = min(pressure_bottom, ozone_pressure_upper_bound)
    if (pressure_upper <= pressure_lower) then
      fraction = 0.0_real64
      return
    end if

    x_lower = log(pressure_lower/ozone_peak_pressure)/(sqrt(2.0_real64)*ozone_log_pressure_width)
    x_upper = log(pressure_upper/ozone_peak_pressure)/(sqrt(2.0_real64)*ozone_log_pressure_width)
    x_profile_lower = log(ozone_pressure_lower_bound/ozone_peak_pressure)/ &
      (sqrt(2.0_real64)*ozone_log_pressure_width)
    x_profile_upper = log(ozone_pressure_upper_bound/ozone_peak_pressure)/ &
      (sqrt(2.0_real64)*ozone_log_pressure_width)
    fraction = (erf(x_upper) - erf(x_lower))/ &
      (erf(x_profile_upper) - erf(x_profile_lower))
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
  subroutine radiation_tendency(config, pressure_half, temperature, surface_temperature, &
                                deep_temperature, lowest_u, lowest_v, sin_latitude, &
                                longitude, time_seconds, temperature_tendency, &
                                surface_temperature_tendency, deep_temperature_tendency, &
                                incoming_shortwave, reflected_shortwave, outgoing_longwave, &
                                latent_heat_flux)
    type(radiation_config), intent(in) :: config
    real(real64), intent(in) :: pressure_half(0:), temperature(:)
    real(real64), intent(in) :: surface_temperature, deep_temperature
    real(real64), intent(in) :: lowest_u, lowest_v
    real(real64), intent(in) :: sin_latitude, longitude, time_seconds
    real(real64), intent(out) :: temperature_tendency(:)
    real(real64), intent(out) :: surface_temperature_tendency, deep_temperature_tendency
    real(real64), intent(out) :: incoming_shortwave, reflected_shortwave, outgoing_longwave
    real(real64), intent(in), optional :: latent_heat_flux
    real(real64) :: upward_longwave(0:size(temperature)), downward_longwave(0:size(temperature))
    real(real64) :: transmission(size(temperature)), emission(size(temperature))
    real(real64) :: net_longwave(0:size(temperature))
    real(real64) :: pressure_thickness, layer_longwave_optical_depth, sensible_heat, surface_deep_heat
    real(real64) :: layer_shortwave_optical_depth, shortwave_transmission
    real(real64) :: shortwave_downward, ultraviolet_downward, non_ultraviolet_downward, shortwave_absorbed
    real(real64) :: stefan_boltzmann_constant, dry_gravity_acceleration, dry_air_specific_heat
    real(real64) :: surface_heat_capacity, surface_latent_heat
    integer :: k, number_of_levels

    surface_latent_heat = 0.0_real64
    if (present(latent_heat_flux)) surface_latent_heat = latent_heat_flux
    stefan_boltzmann_constant = config%stefan_boltzmann_constant
    dry_gravity_acceleration = config%gravity_acceleration
    dry_air_specific_heat = config%dry_air_specific_heat
    surface_heat_capacity = radiation_surface_heat_capacity(config)
    number_of_levels = size(temperature)
    if (number_of_levels < 1 .or. size(pressure_half) /= number_of_levels + 1 .or. &
        size(temperature_tendency) /= number_of_levels) then
      error stop 'radiation column has inconsistent vertical dimensions'
    end if
    if (pressure_half(number_of_levels) <= pressure_half(0) .or. any(temperature <= 0.0_real64) .or. &
        surface_temperature <= 0.0_real64 .or. &
        (.not. config%slab_ocean_enabled .and. deep_temperature <= 0.0_real64)) then
      error stop 'radiation column contains a nonphysical state'
    end if
    if (surface_heat_capacity <= 0.0_real64) then
      error stop 'radiation surface heat capacity must be positive'
    end if

    do k = 1, number_of_levels
      pressure_thickness = pressure_half(k) - pressure_half(k - 1)
      if (pressure_thickness <= 0.0_real64) then
        error stop 'radiation pressures must increase downward'
      end if
      layer_longwave_optical_depth = config%longwave_surface_optical_depth*pressure_thickness/ &
        (pressure_half(number_of_levels) - pressure_half(0)) + &
        ozone_longwave_layer_optical_depth(config, pressure_half(k - 1), pressure_half(k))
      transmission(k) = exp(-layer_longwave_optical_depth)
      emission(k) = (1.0_real64 - transmission(k))*stefan_boltzmann_constant*temperature(k)**4
    end do

    downward_longwave(0) = 0.0_real64
    do k = 1, number_of_levels
      downward_longwave(k) = transmission(k)*downward_longwave(k - 1) + emission(k)
    end do
    upward_longwave(number_of_levels) = stefan_boltzmann_constant*surface_temperature**4
    do k = number_of_levels, 1, -1
      upward_longwave(k - 1) = transmission(k)*upward_longwave(k) + emission(k)
    end do
    net_longwave = upward_longwave - downward_longwave

    do k = 1, number_of_levels
      pressure_thickness = pressure_half(k) - pressure_half(k - 1)
      temperature_tendency(k) = dry_gravity_acceleration/(dry_air_specific_heat*pressure_thickness)* &
        (net_longwave(k) - net_longwave(k - 1))
    end do

    incoming_shortwave = shortwave_downward_flux(config, sin_latitude, longitude, time_seconds)
    shortwave_downward = incoming_shortwave
    if (incoming_shortwave > 0.0_real64) then
      non_ultraviolet_downward = (1.0_real64 - config%ultraviolet_shortwave_fraction)*incoming_shortwave
      ultraviolet_downward = config%ultraviolet_shortwave_fraction*incoming_shortwave
      do k = 1, number_of_levels
        layer_shortwave_optical_depth = ozone_layer_optical_depth(config, pressure_half(k - 1), pressure_half(k))
        shortwave_transmission = exp(-layer_shortwave_optical_depth)
        shortwave_absorbed = ultraviolet_downward*(1.0_real64 - shortwave_transmission)
        pressure_thickness = pressure_half(k) - pressure_half(k - 1)
        temperature_tendency(k) = temperature_tendency(k) + &
          dry_gravity_acceleration/(dry_air_specific_heat*pressure_thickness)*shortwave_absorbed
        ultraviolet_downward = shortwave_transmission*ultraviolet_downward
      end do
      shortwave_downward = non_ultraviolet_downward + ultraviolet_downward
    end if

    sensible_heat = surface_sensible_heat_flux(config, pressure_half, temperature(number_of_levels), &
                                               surface_temperature, lowest_u, lowest_v)
    pressure_thickness = pressure_half(number_of_levels) - pressure_half(number_of_levels - 1)
    temperature_tendency(number_of_levels) = temperature_tendency(number_of_levels) + &
      dry_gravity_acceleration/(dry_air_specific_heat*pressure_thickness)*sensible_heat

    reflected_shortwave = config%surface_shortwave_albedo*shortwave_downward
    outgoing_longwave = upward_longwave(0)
    if (config%slab_ocean_enabled) then
      ! A slab ocean is one isolated, well-mixed water column.  It has no
      ! bottom heat exchange and its otherwise-unused deep state stays fixed.
      surface_deep_heat = 0.0_real64
    else
      surface_deep_heat = config%ground_exchange_coefficient*(surface_temperature - deep_temperature)
    end if
    ! The surface loses exactly the upward longwave the lowest layer sees, so the column budget closes.
    surface_temperature_tendency = (shortwave_downward - reflected_shortwave + &
      downward_longwave(number_of_levels) - upward_longwave(number_of_levels) - &
      surface_deep_heat - sensible_heat - surface_latent_heat)/surface_heat_capacity
    if (config%slab_ocean_enabled) then
      deep_temperature_tendency = 0.0_real64
    else
      deep_temperature_tendency = surface_deep_heat/config%deep_ground_heat_capacity
    end if
  end subroutine radiation_tendency

end module dry_radiation
