module dry_radiation
  use iso_fortran_env, only: real64
  use dry_vertical_coordinate, only: dry_air_gas_constant
  implicit none
  private

  real(real64), parameter, public :: pi = acos(-1.0_real64)
  real(real64), parameter, public :: stefan_boltzmann_constant = 5.670374419e-8_real64
  real(real64), parameter, public :: longwave_surface_optical_depth = 1.0_real64
  real(real64), parameter, public :: dry_air_specific_heat = 1004.0_real64
  real(real64), parameter, public :: dry_gravity_acceleration = 9.80616_real64
  real(real64), parameter, public :: solar_constant = 1361.0_real64
  real(real64), parameter, public :: surface_shortwave_albedo = 0.3_real64
  real(real64), parameter, public :: ultraviolet_shortwave_fraction = 0.02_real64
  real(real64), parameter, public :: ozone_shortwave_optical_depth = 1.5_real64
  real(real64), parameter, public :: ozone_longwave_optical_depth = 0.005_real64
  real(real64), parameter, public :: ozone_pressure_lower_bound = 1.0e2_real64
  real(real64), parameter, public :: ozone_pressure_upper_bound = 1.0e4_real64
  real(real64), parameter, public :: ozone_peak_pressure = 1.0e3_real64
  real(real64), parameter, public :: ozone_log_pressure_width = log(3.0_real64)
  real(real64), parameter, public :: surface_heat_capacity = 2.0e6_real64
  real(real64), parameter, public :: deep_ground_heat_capacity = 2.0e7_real64
  real(real64), parameter, public :: ground_exchange_coefficient = 2.0_real64
  real(real64), parameter, public :: surface_exchange_coefficient = 1.0e-3_real64
  real(real64), parameter, public :: gustiness_speed = 1.0_real64
  real(real64), parameter, public :: axial_tilt = 23.4_real64*pi/180.0_real64
  real(real64), parameter, public :: solar_day = 86400.0_real64
  integer, parameter, public :: radiation_top_rayleigh_levels = 2
  real(real64), parameter, public :: radiation_top_rayleigh_rate = 1.0_real64/solar_day
  integer, parameter, public :: days_per_month = 30
  integer, parameter, public :: months_per_year = 12
  integer, parameter, public :: days_per_year = days_per_month*months_per_year
  real(real64), parameter, public :: orbital_period = real(days_per_year, real64)*solar_day
  real(real64), parameter, public :: planetary_rotation_rate = &
    2.0_real64*pi*(1.0_real64/solar_day + 1.0_real64/orbital_period)

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
    real(real64), allocatable :: surface_temperature(:, :)
    real(real64), allocatable :: surface_pressure(:, :)
    real(real64), allocatable :: zonal_temperature(:, :)
    real(real64), allocatable :: zonal_u(:, :)
    real(real64), allocatable :: zonal_v(:, :)
    real(real64), allocatable :: zonal_uv(:, :)
    real(real64), allocatable :: zonal_vt(:, :)
  end type radiation_diagnostics

  type, public :: radiation_monthly_accumulator
    private
    integer :: sample_count = 0
    real(real64), allocatable :: surface_temperature_sum(:, :)
    real(real64), allocatable :: surface_pressure_sum(:, :)
    real(real64), allocatable :: zonal_temperature_sum(:, :)
    real(real64), allocatable :: zonal_u_sum(:, :)
    real(real64), allocatable :: zonal_v_sum(:, :)
    real(real64), allocatable :: zonal_uv_sum(:, :)
    real(real64), allocatable :: zonal_vt_sum(:, :)
  contains
    procedure, public :: add => add_monthly_sample
    procedure, public :: take => take_monthly_means
    procedure, public :: reset => reset_monthly_accumulator
    procedure, public :: count => monthly_sample_count
  end type radiation_monthly_accumulator

  !> Online equal-weight mean of the global scalar diagnostics over one output interval.
  type, public :: radiation_daily_accumulator
    private
    integer :: sample_count = 0
    real(real64) :: first_time_seconds = 0.0_real64
    real(real64) :: atmospheric_temperature_sum = 0.0_real64
    real(real64) :: surface_temperature_sum = 0.0_real64
    real(real64) :: deep_temperature_sum = 0.0_real64
    real(real64) :: kinetic_energy_sum = 0.0_real64
    real(real64) :: surface_pressure_sum = 0.0_real64
    real(real64) :: incoming_shortwave_sum = 0.0_real64
    real(real64) :: reflected_shortwave_sum = 0.0_real64
    real(real64) :: outgoing_longwave_sum = 0.0_real64
  contains
    procedure, public :: add => add_daily_sample
    procedure, public :: take => take_daily_means
    procedure, public :: reset => reset_daily_accumulator
    procedure, public :: count => daily_sample_count
  end type radiation_daily_accumulator

  public :: shortwave_downward_flux
  public :: ozone_layer_optical_depth
  public :: ozone_longwave_layer_optical_depth
  public :: radiation_tendency
  public :: radiation_calendar_date
  public :: radiation_rayleigh_rate

contains

  pure real(real64) function radiation_rayleigh_rate(level) result(rate)
    integer, intent(in) :: level

    rate = 0.0_real64
    if (level >= 1 .and. level <= radiation_top_rayleigh_levels) then
      rate = radiation_top_rayleigh_rate
    end if
  end function radiation_rayleigh_rate

  real(real64) function shortwave_downward_flux(sin_latitude, longitude, time_seconds) result(flux)
    real(real64), intent(in) :: sin_latitude, longitude, time_seconds
    real(real64) :: orbital_longitude, solar_right_ascension
    real(real64) :: sin_declination, cos_declination, cos_latitude, hour_angle, cos_zenith

    if (abs(sin_latitude) > 1.0_real64) error stop 'radiation latitude sine is outside [-1,1]'
    orbital_longitude = modulo(2.0_real64*pi*time_seconds/orbital_period, 2.0_real64*pi)
    solar_right_ascension = atan2(cos(axial_tilt)*sin(orbital_longitude), cos(orbital_longitude))
    sin_declination = sin(axial_tilt)*sin(orbital_longitude)
    cos_declination = sqrt(max(0.0_real64, 1.0_real64 - sin_declination**2))
    cos_latitude = sqrt(max(0.0_real64, 1.0_real64 - sin_latitude**2))
    hour_angle = planetary_rotation_rate*time_seconds + longitude - solar_right_ascension
    cos_zenith = sin_latitude*sin_declination + cos_latitude*cos_declination*cos(hour_angle)
    flux = solar_constant*max(0.0_real64, cos_zenith)
  end function shortwave_downward_flux

  !> Fraction of the prescribed ozone column within one pressure layer.
  pure real(real64) function ozone_layer_fraction(pressure_top, pressure_bottom) result(fraction)
    real(real64), intent(in) :: pressure_top, pressure_bottom
    real(real64) :: pressure_lower, pressure_upper, x_lower, x_upper
    real(real64) :: x_profile_lower, x_profile_upper

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
  pure real(real64) function ozone_layer_optical_depth(pressure_top, pressure_bottom) result(optical_depth)
    real(real64), intent(in) :: pressure_top, pressure_bottom

    optical_depth = ozone_shortwave_optical_depth*ozone_layer_fraction(pressure_top, pressure_bottom)
  end function ozone_layer_optical_depth

  !> Vertical longwave optical depth added by ozone within one pressure layer.
  pure real(real64) function ozone_longwave_layer_optical_depth(pressure_top, pressure_bottom) result(optical_depth)
    real(real64), intent(in) :: pressure_top, pressure_bottom

    optical_depth = ozone_longwave_optical_depth*ozone_layer_fraction(pressure_top, pressure_bottom)
  end function ozone_longwave_layer_optical_depth

  subroutine radiation_calendar_date(time_seconds, year, month, day, seconds_of_day)
    real(real64), intent(in) :: time_seconds
    integer, intent(out) :: year, month, day
    real(real64), intent(out) :: seconds_of_day
    integer :: elapsed_days, zero_based_month

    if (time_seconds < 0.0_real64) error stop 'radiation calendar time must be nonnegative'
    elapsed_days = floor(time_seconds/solar_day)
    zero_based_month = 3 + elapsed_days/days_per_month
    year = 1 + zero_based_month/months_per_year
    month = modulo(zero_based_month, months_per_year) + 1
    day = modulo(elapsed_days, days_per_month) + 1
    seconds_of_day = modulo(time_seconds, solar_day)
  end subroutine radiation_calendar_date

  !> Column radiation, surface fluxes and ground temperatures.
  !>
  !> The caller supplies one internally consistent RAW-filtered previous-time column for every
  !> state-dependent term.  Solar geometry is evaluated at time_seconds because it is prescribed
  !> rather than prognostic.
  subroutine radiation_tendency(pressure_half, temperature, surface_temperature, &
                                deep_temperature, lowest_u, lowest_v, sin_latitude, &
                                longitude, time_seconds, temperature_tendency, &
                                surface_temperature_tendency, deep_temperature_tendency, &
                                incoming_shortwave, reflected_shortwave, outgoing_longwave)
    real(real64), intent(in) :: pressure_half(0:), temperature(:)
    real(real64), intent(in) :: surface_temperature, deep_temperature
    real(real64), intent(in) :: lowest_u, lowest_v
    real(real64), intent(in) :: sin_latitude, longitude, time_seconds
    real(real64), intent(out) :: temperature_tendency(:)
    real(real64), intent(out) :: surface_temperature_tendency, deep_temperature_tendency
    real(real64), intent(out) :: incoming_shortwave, reflected_shortwave, outgoing_longwave
    real(real64) :: upward_longwave(0:size(temperature)), downward_longwave(0:size(temperature))
    real(real64) :: transmission(size(temperature)), emission(size(temperature))
    real(real64) :: net_longwave(0:size(temperature))
    real(real64) :: pressure_thickness, layer_longwave_optical_depth, sensible_heat, surface_deep_heat
    real(real64) :: layer_shortwave_optical_depth, shortwave_transmission
    real(real64) :: shortwave_downward, ultraviolet_downward, non_ultraviolet_downward, shortwave_absorbed
    integer :: k, number_of_levels

    number_of_levels = size(temperature)
    if (number_of_levels < 1 .or. size(pressure_half) /= number_of_levels + 1 .or. &
        size(temperature_tendency) /= number_of_levels) then
      error stop 'radiation column has inconsistent vertical dimensions'
    end if
    if (pressure_half(number_of_levels) <= pressure_half(0) .or. any(temperature <= 0.0_real64) .or. &
        surface_temperature <= 0.0_real64 .or. deep_temperature <= 0.0_real64) then
      error stop 'radiation column contains a nonphysical state'
    end if

    do k = 1, number_of_levels
      pressure_thickness = pressure_half(k) - pressure_half(k - 1)
      if (pressure_thickness <= 0.0_real64) then
        error stop 'radiation pressures must increase downward'
      end if
      layer_longwave_optical_depth = longwave_surface_optical_depth*pressure_thickness/ &
        (pressure_half(number_of_levels) - pressure_half(0)) + &
        ozone_longwave_layer_optical_depth(pressure_half(k - 1), pressure_half(k))
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

    incoming_shortwave = shortwave_downward_flux(sin_latitude, longitude, time_seconds)
    shortwave_downward = incoming_shortwave
    if (incoming_shortwave > 0.0_real64) then
      non_ultraviolet_downward = (1.0_real64 - ultraviolet_shortwave_fraction)*incoming_shortwave
      ultraviolet_downward = ultraviolet_shortwave_fraction*incoming_shortwave
      do k = 1, number_of_levels
        layer_shortwave_optical_depth = ozone_layer_optical_depth(pressure_half(k - 1), pressure_half(k))
        shortwave_transmission = exp(-layer_shortwave_optical_depth)
        shortwave_absorbed = ultraviolet_downward*(1.0_real64 - shortwave_transmission)
        pressure_thickness = pressure_half(k) - pressure_half(k - 1)
        temperature_tendency(k) = temperature_tendency(k) + &
          dry_gravity_acceleration/(dry_air_specific_heat*pressure_thickness)*shortwave_absorbed
        ultraviolet_downward = shortwave_transmission*ultraviolet_downward
      end do
      shortwave_downward = non_ultraviolet_downward + ultraviolet_downward
    end if

    sensible_heat = pressure_half(number_of_levels)/(dry_air_gas_constant*temperature(number_of_levels))* &
      dry_air_specific_heat*surface_exchange_coefficient* &
      sqrt(lowest_u**2 + lowest_v**2 + gustiness_speed**2)* &
      (surface_temperature - temperature(number_of_levels))
    pressure_thickness = pressure_half(number_of_levels) - pressure_half(number_of_levels - 1)
    temperature_tendency(number_of_levels) = temperature_tendency(number_of_levels) + &
      dry_gravity_acceleration/(dry_air_specific_heat*pressure_thickness)*sensible_heat

    reflected_shortwave = surface_shortwave_albedo*shortwave_downward
    outgoing_longwave = upward_longwave(0)
    surface_deep_heat = ground_exchange_coefficient*(surface_temperature - deep_temperature)
    ! The surface loses exactly the upward longwave the lowest layer sees, so the column budget closes.
    surface_temperature_tendency = (shortwave_downward - reflected_shortwave + &
      downward_longwave(number_of_levels) - upward_longwave(number_of_levels) - &
      surface_deep_heat - sensible_heat)/surface_heat_capacity
    deep_temperature_tendency = surface_deep_heat/deep_ground_heat_capacity
  end subroutine radiation_tendency

  subroutine add_monthly_sample(this, sample)
    class(radiation_monthly_accumulator), intent(inout) :: this
    type(radiation_diagnostics), intent(in) :: sample

    if (.not. allocated(sample%surface_temperature) .or. .not. allocated(sample%surface_pressure) .or. &
        .not. allocated(sample%zonal_temperature) .or. .not. allocated(sample%zonal_u) .or. &
        .not. allocated(sample%zonal_v) .or. .not. allocated(sample%zonal_uv) .or. &
        .not. allocated(sample%zonal_vt)) error stop 'incomplete radiation diagnostic sample'
    if (.not. allocated(this%surface_temperature_sum)) then
      allocate (this%surface_temperature_sum, mold=sample%surface_temperature)
      allocate (this%surface_pressure_sum, mold=sample%surface_pressure)
      allocate (this%zonal_temperature_sum, mold=sample%zonal_temperature)
      allocate (this%zonal_u_sum, mold=sample%zonal_u)
      allocate (this%zonal_v_sum, mold=sample%zonal_v)
      allocate (this%zonal_uv_sum, mold=sample%zonal_uv)
      allocate (this%zonal_vt_sum, mold=sample%zonal_vt)
      call this%reset()
    end if
    this%surface_temperature_sum = this%surface_temperature_sum + sample%surface_temperature
    this%surface_pressure_sum = this%surface_pressure_sum + sample%surface_pressure
    this%zonal_temperature_sum = this%zonal_temperature_sum + sample%zonal_temperature
    this%zonal_u_sum = this%zonal_u_sum + sample%zonal_u
    this%zonal_v_sum = this%zonal_v_sum + sample%zonal_v
    this%zonal_uv_sum = this%zonal_uv_sum + sample%zonal_uv
    this%zonal_vt_sum = this%zonal_vt_sum + sample%zonal_vt
    this%sample_count = this%sample_count + 1
  end subroutine add_monthly_sample

  subroutine take_monthly_means(this, surface_temperature, surface_pressure, zonal_temperature, &
                                zonal_u, zonal_v, eddy_uv, eddy_vt)
    class(radiation_monthly_accumulator), intent(inout) :: this
    real(real64), allocatable, intent(out) :: surface_temperature(:, :), surface_pressure(:, :)
    real(real64), allocatable, intent(out) :: zonal_temperature(:, :), zonal_u(:, :), zonal_v(:, :)
    real(real64), allocatable, intent(out) :: eddy_uv(:, :), eddy_vt(:, :)
    real(real64) :: inverse_count

    if (this%sample_count <= 0) error stop 'radiation monthly accumulator is empty'
    inverse_count = 1.0_real64/real(this%sample_count, real64)
    surface_temperature = this%surface_temperature_sum*inverse_count
    surface_pressure = this%surface_pressure_sum*inverse_count
    zonal_temperature = this%zonal_temperature_sum*inverse_count
    zonal_u = this%zonal_u_sum*inverse_count
    zonal_v = this%zonal_v_sum*inverse_count
    eddy_uv = this%zonal_uv_sum*inverse_count - zonal_u*zonal_v
    eddy_vt = this%zonal_vt_sum*inverse_count - zonal_v*zonal_temperature
    call this%reset()
  end subroutine take_monthly_means

  subroutine reset_monthly_accumulator(this)
    class(radiation_monthly_accumulator), intent(inout) :: this

    this%sample_count = 0
    if (allocated(this%surface_temperature_sum)) this%surface_temperature_sum = 0.0_real64
    if (allocated(this%surface_pressure_sum)) this%surface_pressure_sum = 0.0_real64
    if (allocated(this%zonal_temperature_sum)) this%zonal_temperature_sum = 0.0_real64
    if (allocated(this%zonal_u_sum)) this%zonal_u_sum = 0.0_real64
    if (allocated(this%zonal_v_sum)) this%zonal_v_sum = 0.0_real64
    if (allocated(this%zonal_uv_sum)) this%zonal_uv_sum = 0.0_real64
    if (allocated(this%zonal_vt_sum)) this%zonal_vt_sum = 0.0_real64
  end subroutine reset_monthly_accumulator

  integer function monthly_sample_count(this) result(count)
    class(radiation_monthly_accumulator), intent(in) :: this

    count = this%sample_count
  end function monthly_sample_count

  subroutine add_daily_sample(this, sample)
    class(radiation_daily_accumulator), intent(inout) :: this
    type(radiation_diagnostics), intent(in) :: sample

    if (this%sample_count == 0) this%first_time_seconds = sample%time_seconds
    this%atmospheric_temperature_sum = this%atmospheric_temperature_sum + sample%mean_atmospheric_temperature
    this%surface_temperature_sum = this%surface_temperature_sum + sample%mean_surface_temperature
    this%deep_temperature_sum = this%deep_temperature_sum + sample%mean_deep_temperature
    this%kinetic_energy_sum = this%kinetic_energy_sum + sample%mean_kinetic_energy
    this%surface_pressure_sum = this%surface_pressure_sum + sample%mean_surface_pressure
    this%incoming_shortwave_sum = this%incoming_shortwave_sum + sample%mean_incoming_shortwave
    this%reflected_shortwave_sum = this%reflected_shortwave_sum + sample%mean_reflected_shortwave
    this%outgoing_longwave_sum = this%outgoing_longwave_sum + sample%mean_outgoing_longwave
    this%sample_count = this%sample_count + 1
  end subroutine add_daily_sample

  !> Returns the interval means (time_seconds is the first sample time) and resets the sums.
  subroutine take_daily_means(this, means)
    class(radiation_daily_accumulator), intent(inout) :: this
    type(radiation_diagnostics), intent(out) :: means
    real(real64) :: inverse_count

    if (this%sample_count <= 0) error stop 'radiation daily accumulator is empty'
    inverse_count = 1.0_real64/real(this%sample_count, real64)
    means%time_seconds = this%first_time_seconds
    means%mean_atmospheric_temperature = this%atmospheric_temperature_sum*inverse_count
    means%mean_surface_temperature = this%surface_temperature_sum*inverse_count
    means%mean_deep_temperature = this%deep_temperature_sum*inverse_count
    means%mean_kinetic_energy = this%kinetic_energy_sum*inverse_count
    means%mean_surface_pressure = this%surface_pressure_sum*inverse_count
    means%mean_incoming_shortwave = this%incoming_shortwave_sum*inverse_count
    means%mean_reflected_shortwave = this%reflected_shortwave_sum*inverse_count
    means%mean_outgoing_longwave = this%outgoing_longwave_sum*inverse_count
    call this%reset()
  end subroutine take_daily_means

  subroutine reset_daily_accumulator(this)
    class(radiation_daily_accumulator), intent(inout) :: this

    this%sample_count = 0
    this%first_time_seconds = 0.0_real64
    this%atmospheric_temperature_sum = 0.0_real64
    this%surface_temperature_sum = 0.0_real64
    this%deep_temperature_sum = 0.0_real64
    this%kinetic_energy_sum = 0.0_real64
    this%surface_pressure_sum = 0.0_real64
    this%incoming_shortwave_sum = 0.0_real64
    this%reflected_shortwave_sum = 0.0_real64
    this%outgoing_longwave_sum = 0.0_real64
  end subroutine reset_daily_accumulator

  integer function daily_sample_count(this) result(count)
    class(radiation_daily_accumulator), intent(in) :: this

    count = this%sample_count
  end function daily_sample_count

end module dry_radiation
