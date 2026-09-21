!> Convective adjustment and large-scale condensation, evaluated column by column
!> on the RAW-filtered previous time level in the order of
!> docs/dynamics/moist.md ("physical processes"):
!>
!>   1. dry convective adjustment on (T, q^+),
!>   2. moist convective adjustment on the provisional field advanced by the dry
!>      tendency over the leapfrog width,
!>   3. large-scale condensation on the field advanced further by the moist
!>      convective tendency,
!>
!> so that the same instability or supersaturation is never removed twice.  In
!> a dry atmosphere only the first step is active.
module dry_convection_tendency
  use iso_fortran_env, only: real64
  use dry_physics_config, only: dry_model_physics_config
  use dry_convection, only: dry_convective_adjustment_tendency
  use moist_convection, only: moist_convective_adjustment_tendency
  use large_scale_condensation, only: large_scale_condensation_tendency
  use dry_tendency_workspace, only: dry_workspace_type
  implicit none
  private
  public :: add_convection_column_tendency, add_dry_convection_tendency

contains

  !> One column.  `interval` is the width the leapfrog advances by (2 dt, or dt
  !> and dt/2 in the two initialization steps).  `humidity` is the signed grid
  !> humidity; the processes see max(humidity, 0).
  subroutine add_convection_column_tendency(physics, interval, pressure_half, temperature, humidity, &
                                            temperature_rhs, humidity_rhs, convective_precipitation, &
                                            large_scale_precipitation)
    type(dry_model_physics_config), intent(in) :: physics
    real(real64), intent(in) :: interval
    real(real64), intent(in) :: pressure_half(0:), temperature(:), humidity(:)
    real(real64), intent(inout) :: temperature_rhs(:), humidity_rhs(:)
    real(real64), intent(out) :: convective_precipitation, large_scale_precipitation
    real(real64), dimension(size(temperature)) :: clipped_humidity, dry_temperature_tendency, dry_humidity_tendency
    real(real64), dimension(size(temperature)) :: provisional_temperature, provisional_humidity
    real(real64), dimension(size(temperature)) :: moist_temperature_tendency, moist_humidity_tendency
    real(real64), dimension(size(temperature)) :: condensation_temperature_tendency, condensation_humidity_tendency
    logical :: moist

    moist = physics%moisture%enabled
    convective_precipitation = 0.0_real64
    large_scale_precipitation = 0.0_real64
    clipped_humidity = max(humidity, 0.0_real64)
    dry_temperature_tendency = 0.0_real64
    dry_humidity_tendency = 0.0_real64
    if (physics%convection%enabled) then
      if (moist) then
        call dry_convective_adjustment_tendency(physics%convection, pressure_half, temperature, &
                                                dry_temperature_tendency, clipped_humidity, dry_humidity_tendency)
        humidity_rhs = humidity_rhs + dry_humidity_tendency
      else
        call dry_convective_adjustment_tendency(physics%convection, pressure_half, temperature, &
                                                dry_temperature_tendency)
      end if
      temperature_rhs = temperature_rhs + dry_temperature_tendency
    end if
    if (.not. moist) return

    provisional_temperature = temperature + interval*dry_temperature_tendency
    provisional_humidity = max(clipped_humidity + interval*dry_humidity_tendency, 0.0_real64)
    moist_temperature_tendency = 0.0_real64
    moist_humidity_tendency = 0.0_real64
    if (physics%moist_convection%enabled) then
      call moist_convective_adjustment_tendency(physics%moist_convection, pressure_half, provisional_temperature, &
                                                provisional_humidity, moist_temperature_tendency, &
                                                moist_humidity_tendency, convective_precipitation)
      temperature_rhs = temperature_rhs + moist_temperature_tendency
      humidity_rhs = humidity_rhs + moist_humidity_tendency
    end if
    if (physics%condensation%enabled) then
      provisional_temperature = provisional_temperature + interval*moist_temperature_tendency
      provisional_humidity = max(provisional_humidity + interval*moist_humidity_tendency, 0.0_real64)
      call large_scale_condensation_tendency(physics%condensation, pressure_half, provisional_temperature, &
                                             provisional_humidity, interval, condensation_temperature_tendency, &
                                             condensation_humidity_tendency, large_scale_precipitation)
      temperature_rhs = temperature_rhs + condensation_temperature_tendency
      humidity_rhs = humidity_rhs + condensation_humidity_tendency
    end if
  end subroutine add_convection_column_tendency

  !> Evaluated column by column on the RAW-filtered previous time level; the
  !> precipitation of each column is recorded for the diagnostics.
  subroutine add_dry_convection_tendency(physics, interval, workspace)
    type(dry_model_physics_config), intent(in) :: physics
    real(real64), intent(in) :: interval
    type(dry_workspace_type), intent(inout) :: workspace
    integer :: i, j

    do j = 1, workspace%ny
      do i = 1, workspace%ring_nlon(j)
        call add_convection_column_tendency(physics, interval, workspace%previous_pressure_half(i, j, :), &
          workspace%previous_temperature_grid(i, j, :), workspace%previous_humidity_grid(i, j, :), &
          workspace%forcing_temperature(i, j, :), workspace%forcing_humidity(i, j, :), &
          workspace%convective_precipitation(i, j), workspace%large_scale_precipitation(i, j))
      end do
    end do
  end subroutine add_dry_convection_tendency

end module dry_convection_tendency
