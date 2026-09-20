!> Case-specific time aggregation of the radiation diagnostic samples.
!>
!> The physics module produces one instantaneous sample per step and knows
!> nothing about output calendars.  Deciding what a "day" or a "month" means,
!> and accumulating the samples over those intervals, belongs to the radiation
!> case and lives here.
module radiation_diagnostics_collector
  use iso_fortran_env, only: real64
  use dry_radiation, only: radiation_diagnostics
  implicit none
  private

  !> Online equal-weight mean of the zonal and surface fields over one output interval.
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

  subroutine take_monthly(this, surface_temperature, surface_pressure, zonal_temperature, &
                          zonal_u, zonal_v, eddy_uv, eddy_vt)
    class(radiation_case_diagnostics), intent(inout) :: this
    real(real64), allocatable, intent(out) :: surface_temperature(:, :), surface_pressure(:, :)
    real(real64), allocatable, intent(out) :: zonal_temperature(:, :), zonal_u(:, :), zonal_v(:, :)
    real(real64), allocatable, intent(out) :: eddy_uv(:, :), eddy_vt(:, :)
    call this%monthly%take(surface_temperature, surface_pressure, zonal_temperature, &
                           zonal_u, zonal_v, eddy_uv, eddy_vt)
  end subroutine take_monthly

  subroutine reset(this)
    class(radiation_case_diagnostics), intent(inout) :: this
    call this%daily%reset()
    call this%monthly%reset()
  end subroutine reset

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

end module radiation_diagnostics_collector
