!> Case-specific time aggregation of the radiation diagnostic samples.
!>
!> The physics module produces one instantaneous sample per step and knows
!> nothing about output calendars.  Deciding what a "day" or a "month" means,
!> and accumulating the samples over those intervals, belongs to the radiation
!> case and lives here.  The moist fields are accumulated only when the samples
!> carry them.
module radiation_diagnostics_collector
  use iso_fortran_env, only: real64
  use sea_ice, only: sea_ice_checks
  use dry_radiation, only: radiation_diagnostics
  implicit none
  private

  !> Monthly means of the grid and zonal fields, returned in one bundle.
  type, public :: radiation_monthly_means
    real(real64), allocatable :: land_temperature(:, :)
    real(real64), allocatable :: ocean_temperature(:, :)
    real(real64), allocatable :: sea_ice_fraction(:, :)
    real(real64), allocatable :: sea_ice_volume(:, :)
    real(real64), allocatable :: sea_ice_temperature(:, :)
    real(real64), allocatable :: sea_ice_thickness(:, :)
    real(real64), allocatable :: surface_temperature(:, :), deep_temperature(:, :), surface_pressure(:, :)
    real(real64), allocatable :: zonal_temperature(:, :), zonal_u(:, :), zonal_v(:, :)
    real(real64), allocatable :: eddy_uv(:, :), eddy_vt(:, :)
    !> Moist fields; unallocated in a dry run.
    real(real64), allocatable :: precipitation(:, :), evaporation(:, :), precipitable_water(:, :)
    real(real64), allocatable :: cloud_cover(:, :)
    real(real64), allocatable :: surface_water(:, :), surface_wetness(:, :), runoff(:, :)
    real(real64), allocatable :: zonal_humidity(:, :), eddy_vq(:, :)
  end type radiation_monthly_means

  !> Online equal-weight mean of the zonal and surface fields over one output interval.
  type, public :: radiation_monthly_accumulator
    private
    integer :: sample_count = 0
    logical :: moist = .false.
    logical :: tiles = .false.
    real(real64), allocatable :: land_temperature_sum(:, :)
    real(real64), allocatable :: ocean_temperature_sum(:, :)
    real(real64), allocatable :: sea_ice_fraction_sum(:, :)
    real(real64), allocatable :: sea_ice_volume_sum(:, :)
    real(real64), allocatable :: sea_ice_temperature_sum(:, :)
    real(real64), allocatable :: surface_temperature_sum(:, :)
    real(real64), allocatable :: deep_temperature_sum(:, :)
    real(real64), allocatable :: surface_pressure_sum(:, :)
    real(real64), allocatable :: zonal_temperature_sum(:, :)
    real(real64), allocatable :: zonal_u_sum(:, :)
    real(real64), allocatable :: zonal_v_sum(:, :)
    real(real64), allocatable :: zonal_uv_sum(:, :)
    real(real64), allocatable :: zonal_vt_sum(:, :)
    real(real64), allocatable :: precipitation_sum(:, :)
    real(real64), allocatable :: evaporation_sum(:, :)
    real(real64), allocatable :: precipitable_water_sum(:, :)
    real(real64), allocatable :: cloud_cover_sum(:, :)
    real(real64), allocatable :: surface_water_sum(:, :), surface_wetness_sum(:, :), runoff_sum(:, :)
    real(real64), allocatable :: zonal_humidity_sum(:, :)
    real(real64), allocatable :: zonal_vq_sum(:, :)
  contains
    procedure, public :: add => add_monthly_sample
    procedure, public :: take => take_monthly_means
    procedure, public :: reset => reset_monthly_accumulator
    procedure, public :: count => monthly_sample_count
  end type radiation_monthly_accumulator

  !> Online equal-weight mean of the global scalar diagnostics over one output
  !> interval, together with the interval maximum of the wind speed.
  type, public :: radiation_daily_accumulator
    private
    type(sea_ice_checks) :: ice_checks
    integer :: sample_count = 0
    real(real64) :: sea_ice_area_sum = 0.0_real64
    real(real64) :: sea_ice_total_volume_sum = 0.0_real64
    real(real64) :: first_time_seconds = 0.0_real64
    real(real64) :: atmospheric_temperature_sum = 0.0_real64
    real(real64) :: surface_temperature_sum = 0.0_real64
    real(real64) :: deep_temperature_sum = 0.0_real64
    real(real64) :: land_surface_temperature_sum = 0.0_real64
    real(real64) :: ocean_surface_temperature_sum = 0.0_real64
    real(real64) :: land_precipitation_sum = 0.0_real64
    real(real64) :: ocean_precipitation_sum = 0.0_real64
    real(real64) :: land_evaporation_sum = 0.0_real64
    real(real64) :: ocean_evaporation_sum = 0.0_real64
    real(real64) :: surface_water_sum = 0.0_real64
    real(real64) :: surface_wetness_sum = 0.0_real64
    real(real64) :: runoff_sum = 0.0_real64
    real(real64) :: dry_land_fraction_sum = 0.0_real64
    real(real64) :: water_budget_residual_sum = 0.0_real64
    real(real64) :: maximum_water_budget_residual = 0.0_real64
    real(real64) :: kinetic_energy_sum = 0.0_real64
    real(real64) :: surface_pressure_sum = 0.0_real64
    real(real64) :: incoming_shortwave_sum = 0.0_real64
    real(real64) :: reflected_shortwave_sum = 0.0_real64
    real(real64) :: outgoing_longwave_sum = 0.0_real64
    real(real64) :: convective_precipitation_sum = 0.0_real64
    real(real64) :: large_scale_precipitation_sum = 0.0_real64
    real(real64) :: evaporation_sum = 0.0_real64
    real(real64) :: latent_heat_flux_sum = 0.0_real64
    real(real64) :: precipitable_water_sum = 0.0_real64
    real(real64) :: signed_column_water_sum = 0.0_real64
    real(real64) :: negative_column_water_sum = 0.0_real64
    real(real64) :: cloud_cover_sum = 0.0_real64
    real(real64) :: maximum_wind_speed = -1.0_real64
    real(real64) :: maximum_wind_longitude_degrees = 0.0_real64
    real(real64) :: maximum_wind_latitude_degrees = 0.0_real64
    integer :: maximum_wind_level = 0
    real(real64) :: maximum_wind_eta = 0.0_real64
  contains
    procedure, public :: add => add_daily_sample
    procedure, public :: take => take_daily_means
    procedure, public :: reset => reset_daily_accumulator
    procedure, public :: count => daily_sample_count
  end type radiation_daily_accumulator

  !> The collector the radiation case runner drives: one sample in, daily and
  !> monthly means out at the boundaries the runner chooses.
  type, public :: radiation_case_diagnostics
    private
    type(radiation_daily_accumulator) :: daily
    type(radiation_monthly_accumulator) :: monthly
  contains
    procedure, public :: add => add_sample
    procedure, public :: take_daily
    procedure, public :: take_monthly
    procedure, public :: reset
  end type radiation_case_diagnostics

contains

  subroutine add_sample(this, sample)
    class(radiation_case_diagnostics), intent(inout) :: this
    type(radiation_diagnostics), intent(in) :: sample
    call this%daily%add(sample)
    call this%monthly%add(sample)
  end subroutine add_sample

  subroutine take_daily(this, means)
    class(radiation_case_diagnostics), intent(inout) :: this
    type(radiation_diagnostics), intent(out) :: means
    call this%daily%take(means)
  end subroutine take_daily

  subroutine take_monthly(this, means)
    class(radiation_case_diagnostics), intent(inout) :: this
    type(radiation_monthly_means), intent(out) :: means
    call this%monthly%take(means)
  end subroutine take_monthly

  subroutine reset(this)
    class(radiation_case_diagnostics), intent(inout) :: this
    call this%daily%reset()
    call this%monthly%reset()
  end subroutine reset

  subroutine add_monthly_sample(this, sample)
    class(radiation_monthly_accumulator), intent(inout) :: this
    type(radiation_diagnostics), intent(in) :: sample

    if (.not. allocated(sample%surface_temperature) .or. .not. allocated(sample%deep_temperature) .or. &
        .not. allocated(sample%surface_pressure) .or. &
        .not. allocated(sample%zonal_temperature) .or. .not. allocated(sample%zonal_u) .or. &
        .not. allocated(sample%zonal_v) .or. .not. allocated(sample%zonal_uv) .or. &
        .not. allocated(sample%zonal_vt)) error stop 'incomplete radiation diagnostic sample'
    if (.not. allocated(this%surface_temperature_sum)) then
      allocate (this%surface_temperature_sum, mold=sample%surface_temperature)
      allocate (this%deep_temperature_sum, mold=sample%deep_temperature)
      allocate (this%surface_pressure_sum, mold=sample%surface_pressure)
      allocate (this%zonal_temperature_sum, mold=sample%zonal_temperature)
      allocate (this%zonal_u_sum, mold=sample%zonal_u)
      allocate (this%zonal_v_sum, mold=sample%zonal_v)
      allocate (this%zonal_uv_sum, mold=sample%zonal_uv)
      allocate (this%zonal_vt_sum, mold=sample%zonal_vt)
      allocate (this%surface_water_sum, mold=sample%surface_water)
      allocate (this%surface_wetness_sum, mold=sample%surface_wetness)
      allocate (this%runoff_sum, mold=sample%runoff)
      this%tiles = allocated(sample%sea_ice_fraction)
      if (this%tiles) then
        allocate (this%land_temperature_sum, mold=sample%land_temperature)
        allocate (this%ocean_temperature_sum, mold=sample%ocean_temperature)
        allocate (this%sea_ice_fraction_sum, mold=sample%sea_ice_fraction)
        allocate (this%sea_ice_volume_sum, mold=sample%sea_ice_volume)
        allocate (this%sea_ice_temperature_sum, mold=sample%sea_ice_temperature)
      end if
      this%moist = allocated(sample%precipitation)
      if (this%moist) then
        if (.not. allocated(sample%evaporation) .or. .not. allocated(sample%precipitable_water) .or. &
            .not. allocated(sample%zonal_humidity) .or. .not. allocated(sample%zonal_vq) .or. &
            .not. allocated(sample%cloud_cover)) then
          error stop 'incomplete moist diagnostic sample'
        end if
        allocate (this%precipitation_sum, mold=sample%precipitation)
        allocate (this%evaporation_sum, mold=sample%evaporation)
        allocate (this%precipitable_water_sum, mold=sample%precipitable_water)
        allocate (this%cloud_cover_sum, mold=sample%cloud_cover)
        allocate (this%zonal_humidity_sum, mold=sample%zonal_humidity)
        allocate (this%zonal_vq_sum, mold=sample%zonal_vq)
      end if
      call this%reset()
    end if
    if (this%moist .neqv. allocated(sample%precipitation)) then
      error stop 'radiation diagnostic samples must all carry, or all lack, the moist fields'
    end if
    if (this%tiles .neqv. allocated(sample%sea_ice_fraction)) error stop 'inconsistent tile diagnostics'
    if (this%tiles) then
      this%land_temperature_sum = this%land_temperature_sum + sample%land_temperature
      this%ocean_temperature_sum = this%ocean_temperature_sum + sample%ocean_temperature
      this%sea_ice_fraction_sum = this%sea_ice_fraction_sum + sample%sea_ice_fraction
      this%sea_ice_volume_sum = this%sea_ice_volume_sum + sample%sea_ice_volume
      this%sea_ice_temperature_sum = this%sea_ice_temperature_sum + sample%sea_ice_fraction*sample%sea_ice_temperature
    end if
    this%surface_temperature_sum = this%surface_temperature_sum + sample%surface_temperature
    this%deep_temperature_sum = this%deep_temperature_sum + sample%deep_temperature
    this%surface_pressure_sum = this%surface_pressure_sum + sample%surface_pressure
    this%zonal_temperature_sum = this%zonal_temperature_sum + sample%zonal_temperature
    this%zonal_u_sum = this%zonal_u_sum + sample%zonal_u
    this%zonal_v_sum = this%zonal_v_sum + sample%zonal_v
    this%zonal_uv_sum = this%zonal_uv_sum + sample%zonal_uv
    this%zonal_vt_sum = this%zonal_vt_sum + sample%zonal_vt
    this%surface_water_sum = this%surface_water_sum + sample%surface_water
    this%surface_wetness_sum = this%surface_wetness_sum + sample%surface_wetness
    this%runoff_sum = this%runoff_sum + sample%runoff
    if (this%moist) then
      this%precipitation_sum = this%precipitation_sum + sample%precipitation
      this%evaporation_sum = this%evaporation_sum + sample%evaporation
      this%precipitable_water_sum = this%precipitable_water_sum + sample%precipitable_water
      this%cloud_cover_sum = this%cloud_cover_sum + sample%cloud_cover
      this%zonal_humidity_sum = this%zonal_humidity_sum + sample%zonal_humidity
      this%zonal_vq_sum = this%zonal_vq_sum + sample%zonal_vq
    end if
    this%sample_count = this%sample_count + 1
  end subroutine add_monthly_sample

  subroutine take_monthly_means(this, means)
    class(radiation_monthly_accumulator), intent(inout) :: this
    type(radiation_monthly_means), intent(out) :: means
    real(real64) :: inverse_count

    if (this%sample_count <= 0) error stop 'radiation monthly accumulator is empty'
    inverse_count = 1.0_real64/real(this%sample_count, real64)
    if (this%tiles) then
      means%land_temperature = this%land_temperature_sum*inverse_count
      means%ocean_temperature = this%ocean_temperature_sum*inverse_count
      means%sea_ice_fraction = this%sea_ice_fraction_sum*inverse_count
      means%sea_ice_volume = this%sea_ice_volume_sum*inverse_count
      allocate (means%sea_ice_temperature, mold=this%sea_ice_fraction_sum)
      allocate (means%sea_ice_thickness, mold=this%sea_ice_fraction_sum)
      means%sea_ice_temperature = 0.0_real64
      means%sea_ice_thickness = 0.0_real64
      where (this%sea_ice_fraction_sum > 0.0_real64)
        means%sea_ice_temperature = this%sea_ice_temperature_sum/this%sea_ice_fraction_sum
        means%sea_ice_thickness = this%sea_ice_volume_sum/this%sea_ice_fraction_sum
      end where
    end if
    means%surface_temperature = this%surface_temperature_sum*inverse_count
    means%deep_temperature = this%deep_temperature_sum*inverse_count
    means%surface_pressure = this%surface_pressure_sum*inverse_count
    means%zonal_temperature = this%zonal_temperature_sum*inverse_count
    means%zonal_u = this%zonal_u_sum*inverse_count
    means%zonal_v = this%zonal_v_sum*inverse_count
    means%eddy_uv = this%zonal_uv_sum*inverse_count - means%zonal_u*means%zonal_v
    means%eddy_vt = this%zonal_vt_sum*inverse_count - means%zonal_v*means%zonal_temperature
    means%surface_water = this%surface_water_sum*inverse_count
    means%surface_wetness = this%surface_wetness_sum*inverse_count
    means%runoff = this%runoff_sum*inverse_count
    if (this%moist) then
      means%precipitation = this%precipitation_sum*inverse_count
      means%evaporation = this%evaporation_sum*inverse_count
      means%precipitable_water = this%precipitable_water_sum*inverse_count
      means%cloud_cover = this%cloud_cover_sum*inverse_count
      means%zonal_humidity = this%zonal_humidity_sum*inverse_count
      means%eddy_vq = this%zonal_vq_sum*inverse_count - means%zonal_v*means%zonal_humidity
    end if
    call this%reset()
  end subroutine take_monthly_means

  subroutine reset_monthly_accumulator(this)
    class(radiation_monthly_accumulator), intent(inout) :: this

    this%sample_count = 0
    if (allocated(this%land_temperature_sum)) this%land_temperature_sum = 0.0_real64
    if (allocated(this%ocean_temperature_sum)) this%ocean_temperature_sum = 0.0_real64
    if (allocated(this%sea_ice_fraction_sum)) this%sea_ice_fraction_sum = 0.0_real64
    if (allocated(this%sea_ice_volume_sum)) this%sea_ice_volume_sum = 0.0_real64
    if (allocated(this%sea_ice_temperature_sum)) this%sea_ice_temperature_sum = 0.0_real64
    if (allocated(this%surface_temperature_sum)) this%surface_temperature_sum = 0.0_real64
    if (allocated(this%deep_temperature_sum)) this%deep_temperature_sum = 0.0_real64
    if (allocated(this%surface_pressure_sum)) this%surface_pressure_sum = 0.0_real64
    if (allocated(this%zonal_temperature_sum)) this%zonal_temperature_sum = 0.0_real64
    if (allocated(this%zonal_u_sum)) this%zonal_u_sum = 0.0_real64
    if (allocated(this%zonal_v_sum)) this%zonal_v_sum = 0.0_real64
    if (allocated(this%zonal_uv_sum)) this%zonal_uv_sum = 0.0_real64
    if (allocated(this%zonal_vt_sum)) this%zonal_vt_sum = 0.0_real64
    if (allocated(this%surface_water_sum)) this%surface_water_sum = 0.0_real64
    if (allocated(this%surface_wetness_sum)) this%surface_wetness_sum = 0.0_real64
    if (allocated(this%runoff_sum)) this%runoff_sum = 0.0_real64
    if (allocated(this%precipitation_sum)) this%precipitation_sum = 0.0_real64
    if (allocated(this%evaporation_sum)) this%evaporation_sum = 0.0_real64
    if (allocated(this%precipitable_water_sum)) this%precipitable_water_sum = 0.0_real64
    if (allocated(this%cloud_cover_sum)) this%cloud_cover_sum = 0.0_real64
    if (allocated(this%zonal_humidity_sum)) this%zonal_humidity_sum = 0.0_real64
    if (allocated(this%zonal_vq_sum)) this%zonal_vq_sum = 0.0_real64
  end subroutine reset_monthly_accumulator

  integer function monthly_sample_count(this) result(count)
    class(radiation_monthly_accumulator), intent(in) :: this

    count = this%sample_count
  end function monthly_sample_count

  subroutine add_daily_sample(this, sample)
    class(radiation_daily_accumulator), intent(inout) :: this
    type(radiation_diagnostics), intent(in) :: sample

    if (this%sample_count == 0) this%first_time_seconds = sample%time_seconds
    call this%ice_checks%merge(sample%ice_checks)
    this%sea_ice_area_sum = this%sea_ice_area_sum + sample%sea_ice_area
    this%sea_ice_total_volume_sum = this%sea_ice_total_volume_sum + sample%sea_ice_total_volume
    this%atmospheric_temperature_sum = this%atmospheric_temperature_sum + sample%mean_atmospheric_temperature
    this%surface_temperature_sum = this%surface_temperature_sum + sample%mean_surface_temperature
    this%deep_temperature_sum = this%deep_temperature_sum + sample%mean_deep_temperature
    this%land_surface_temperature_sum = this%land_surface_temperature_sum + sample%mean_land_surface_temperature
    this%ocean_surface_temperature_sum = this%ocean_surface_temperature_sum + sample%mean_ocean_surface_temperature
    this%land_precipitation_sum = this%land_precipitation_sum + sample%mean_land_precipitation
    this%ocean_precipitation_sum = this%ocean_precipitation_sum + sample%mean_ocean_precipitation
    this%land_evaporation_sum = this%land_evaporation_sum + sample%mean_land_evaporation
    this%ocean_evaporation_sum = this%ocean_evaporation_sum + sample%mean_ocean_evaporation
    this%surface_water_sum = this%surface_water_sum + sample%mean_surface_water
    this%surface_wetness_sum = this%surface_wetness_sum + sample%mean_surface_wetness
    this%runoff_sum = this%runoff_sum + sample%mean_runoff
    this%dry_land_fraction_sum = this%dry_land_fraction_sum + sample%dry_land_fraction
    this%water_budget_residual_sum = this%water_budget_residual_sum + sample%mean_water_budget_residual
    this%maximum_water_budget_residual = max(this%maximum_water_budget_residual, &
      sample%maximum_water_budget_residual)
    this%kinetic_energy_sum = this%kinetic_energy_sum + sample%mean_kinetic_energy
    this%surface_pressure_sum = this%surface_pressure_sum + sample%mean_surface_pressure
    this%incoming_shortwave_sum = this%incoming_shortwave_sum + sample%mean_incoming_shortwave
    this%reflected_shortwave_sum = this%reflected_shortwave_sum + sample%mean_reflected_shortwave
    this%outgoing_longwave_sum = this%outgoing_longwave_sum + sample%mean_outgoing_longwave
    this%convective_precipitation_sum = this%convective_precipitation_sum + sample%mean_convective_precipitation
    this%large_scale_precipitation_sum = this%large_scale_precipitation_sum + sample%mean_large_scale_precipitation
    this%evaporation_sum = this%evaporation_sum + sample%mean_evaporation
    this%latent_heat_flux_sum = this%latent_heat_flux_sum + sample%mean_latent_heat_flux
    this%precipitable_water_sum = this%precipitable_water_sum + sample%mean_precipitable_water
    this%signed_column_water_sum = this%signed_column_water_sum + sample%mean_signed_column_water
    this%negative_column_water_sum = this%negative_column_water_sum + sample%mean_negative_column_water
    this%cloud_cover_sum = this%cloud_cover_sum + sample%mean_cloud_cover
    ! The interval maximum keeps the first sample that reaches it.
    if (sample%maximum_wind_speed > this%maximum_wind_speed) then
      this%maximum_wind_speed = sample%maximum_wind_speed
      this%maximum_wind_longitude_degrees = sample%maximum_wind_longitude_degrees
      this%maximum_wind_latitude_degrees = sample%maximum_wind_latitude_degrees
      this%maximum_wind_level = sample%maximum_wind_level
      this%maximum_wind_eta = sample%maximum_wind_eta
    end if
    this%sample_count = this%sample_count + 1
  end subroutine add_daily_sample

  !> Returns the interval means (time_seconds is the first sample time), the
  !> interval wind maximum, and resets the sums.
  subroutine take_daily_means(this, means)
    class(radiation_daily_accumulator), intent(inout) :: this
    type(radiation_diagnostics), intent(out) :: means
    real(real64) :: inverse_count

    if (this%sample_count <= 0) error stop 'radiation daily accumulator is empty'
    means%ice_checks = this%ice_checks
    inverse_count = 1.0_real64/real(this%sample_count, real64)
    means%sea_ice_area = this%sea_ice_area_sum*inverse_count
    means%sea_ice_total_volume = this%sea_ice_total_volume_sum*inverse_count
    if (this%sea_ice_area_sum > 0.0_real64) &
      means%mean_sea_ice_thickness = this%sea_ice_total_volume_sum/this%sea_ice_area_sum
    means%time_seconds = this%first_time_seconds
    means%mean_atmospheric_temperature = this%atmospheric_temperature_sum*inverse_count
    means%mean_surface_temperature = this%surface_temperature_sum*inverse_count
    means%mean_deep_temperature = this%deep_temperature_sum*inverse_count
    means%mean_land_surface_temperature = this%land_surface_temperature_sum*inverse_count
    means%mean_ocean_surface_temperature = this%ocean_surface_temperature_sum*inverse_count
    means%mean_land_precipitation = this%land_precipitation_sum*inverse_count
    means%mean_ocean_precipitation = this%ocean_precipitation_sum*inverse_count
    means%mean_land_evaporation = this%land_evaporation_sum*inverse_count
    means%mean_ocean_evaporation = this%ocean_evaporation_sum*inverse_count
    means%mean_surface_water = this%surface_water_sum*inverse_count
    means%mean_surface_wetness = this%surface_wetness_sum*inverse_count
    means%mean_runoff = this%runoff_sum*inverse_count
    means%dry_land_fraction = this%dry_land_fraction_sum*inverse_count
    means%mean_water_budget_residual = this%water_budget_residual_sum*inverse_count
    means%maximum_water_budget_residual = this%maximum_water_budget_residual
    means%mean_kinetic_energy = this%kinetic_energy_sum*inverse_count
    means%mean_surface_pressure = this%surface_pressure_sum*inverse_count
    means%mean_incoming_shortwave = this%incoming_shortwave_sum*inverse_count
    means%mean_reflected_shortwave = this%reflected_shortwave_sum*inverse_count
    means%mean_outgoing_longwave = this%outgoing_longwave_sum*inverse_count
    means%mean_convective_precipitation = this%convective_precipitation_sum*inverse_count
    means%mean_large_scale_precipitation = this%large_scale_precipitation_sum*inverse_count
    means%mean_evaporation = this%evaporation_sum*inverse_count
    means%mean_latent_heat_flux = this%latent_heat_flux_sum*inverse_count
    means%mean_precipitable_water = this%precipitable_water_sum*inverse_count
    means%mean_signed_column_water = this%signed_column_water_sum*inverse_count
    means%mean_negative_column_water = this%negative_column_water_sum*inverse_count
    means%mean_cloud_cover = this%cloud_cover_sum*inverse_count
    means%maximum_wind_speed = max(this%maximum_wind_speed, 0.0_real64)
    means%maximum_wind_longitude_degrees = this%maximum_wind_longitude_degrees
    means%maximum_wind_latitude_degrees = this%maximum_wind_latitude_degrees
    means%maximum_wind_level = this%maximum_wind_level
    means%maximum_wind_eta = this%maximum_wind_eta
    call this%reset()
  end subroutine take_daily_means

  subroutine reset_daily_accumulator(this)
    class(radiation_daily_accumulator), intent(inout) :: this

    this%ice_checks = sea_ice_checks()
    this%sample_count = 0
    this%first_time_seconds = 0.0_real64
    this%sea_ice_area_sum = 0.0_real64
    this%sea_ice_total_volume_sum = 0.0_real64
    this%atmospheric_temperature_sum = 0.0_real64
    this%surface_temperature_sum = 0.0_real64
    this%deep_temperature_sum = 0.0_real64
    this%land_surface_temperature_sum = 0.0_real64
    this%ocean_surface_temperature_sum = 0.0_real64
    this%land_precipitation_sum = 0.0_real64
    this%ocean_precipitation_sum = 0.0_real64
    this%land_evaporation_sum = 0.0_real64
    this%ocean_evaporation_sum = 0.0_real64
    this%surface_water_sum = 0.0_real64
    this%surface_wetness_sum = 0.0_real64
    this%runoff_sum = 0.0_real64
    this%dry_land_fraction_sum = 0.0_real64
    this%water_budget_residual_sum = 0.0_real64
    this%maximum_water_budget_residual = 0.0_real64
    this%kinetic_energy_sum = 0.0_real64
    this%surface_pressure_sum = 0.0_real64
    this%incoming_shortwave_sum = 0.0_real64
    this%reflected_shortwave_sum = 0.0_real64
    this%outgoing_longwave_sum = 0.0_real64
    this%convective_precipitation_sum = 0.0_real64
    this%large_scale_precipitation_sum = 0.0_real64
    this%evaporation_sum = 0.0_real64
    this%latent_heat_flux_sum = 0.0_real64
    this%precipitable_water_sum = 0.0_real64
    this%signed_column_water_sum = 0.0_real64
    this%negative_column_water_sum = 0.0_real64
    this%cloud_cover_sum = 0.0_real64
    this%maximum_wind_speed = -1.0_real64
    this%maximum_wind_longitude_degrees = 0.0_real64
    this%maximum_wind_latitude_degrees = 0.0_real64
    this%maximum_wind_level = 0
    this%maximum_wind_eta = 0.0_real64
  end subroutine reset_daily_accumulator

  integer function daily_sample_count(this) result(count)
    class(radiation_daily_accumulator), intent(in) :: this

    count = this%sample_count
  end function daily_sample_count

end module radiation_diagnostics_collector
