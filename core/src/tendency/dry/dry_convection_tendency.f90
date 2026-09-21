!> Convective adjustment and large-scale condensation, evaluated column by column
!> on the RAW-filtered previous time level in the order of
!> docs/dynamics/moist.md ("physical processes"):
!>
!>   1. dry convective adjustment on (T, q^+),
!>   2. moist convective adjustment on the provisional field advanced by the dry
!>      tendency over the leapfrog width,
!>   3. large-scale condensation on the field advanced further by the moist
!>      convective tendency,
!>   4. the effective column cloud cover from the field advanced by the
!>      condensation and from the convective precipitation (docs/tendency/cloud.md),
!>
!> so that the same instability or supersaturation is never removed twice.  In
!> a dry atmosphere only the first step is active.
!>
!> The full-level pressures come from the workspace and the Exner function is
!> evaluated once per column, shared by the three processes.
module dry_convection_tendency
  use iso_fortran_env, only: real64
  use dry_vertical_coordinate, only: dry_air_kappa, reference_surface_pressure
  use dry_physics_config, only: dry_model_physics_config
  use dry_convection, only: dry_convective_adjustment_from_levels
  use moist_convection, only: moist_convective_adjustment_from_levels
  use large_scale_condensation, only: large_scale_condensation_from_levels
  use cloud_diagnostics, only: diagnose_cloud_cover
  use moist_thermodynamics, only: full_level_pressures
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
                                            large_scale_precipitation, cloud_cover)
    type(dry_model_physics_config), intent(in) :: physics
    real(real64), intent(in) :: interval
    real(real64), intent(in) :: pressure_half(0:), temperature(:), humidity(:)
    real(real64), intent(inout) :: temperature_rhs(:), humidity_rhs(:)
    real(real64), intent(out) :: convective_precipitation, large_scale_precipitation, cloud_cover
    real(real64), dimension(size(temperature)) :: full_level_pressure, delta_pressure

    call full_level_pressures(pressure_half, full_level_pressure, delta_pressure)
    call add_convection_column_tendency_from_levels(physics, interval, full_level_pressure, delta_pressure, &
                                                    temperature, humidity, temperature_rhs, humidity_rhs, &
                                                    convective_precipitation, large_scale_precipitation, cloud_cover)
  end subroutine add_convection_column_tendency

  !> `cloud_cover` is the effective column cloud cover diagnosed from the field
  !> left by the three processes; zero unless the cloud diagnosis is enabled.
  subroutine add_convection_column_tendency_from_levels(physics, interval, full_level_pressure, delta_pressure, &
                                                        temperature, humidity, temperature_rhs, humidity_rhs, &
                                                        convective_precipitation, large_scale_precipitation, &
                                                        cloud_cover)
    type(dry_model_physics_config), intent(in) :: physics
    real(real64), intent(in) :: interval
    real(real64), intent(in) :: full_level_pressure(:), delta_pressure(:), temperature(:), humidity(:)
    real(real64), intent(inout) :: temperature_rhs(:), humidity_rhs(:)
    real(real64), intent(out) :: convective_precipitation, large_scale_precipitation, cloud_cover
    real(real64), dimension(size(temperature)) :: exner
    real(real64), dimension(size(temperature)) :: clipped_humidity, dry_temperature_tendency, dry_humidity_tendency
    real(real64), dimension(size(temperature)) :: provisional_temperature, provisional_humidity
    real(real64), dimension(size(temperature)) :: moist_temperature_tendency, moist_humidity_tendency
    real(real64), dimension(size(temperature)) :: condensation_temperature_tendency, condensation_humidity_tendency
    logical :: moist

    moist = physics%moisture%enabled
    convective_precipitation = 0.0_real64
    large_scale_precipitation = 0.0_real64
    cloud_cover = 0.0_real64
    exner = (full_level_pressure/reference_surface_pressure)**dry_air_kappa
    clipped_humidity = max(humidity, 0.0_real64)
    dry_temperature_tendency = 0.0_real64
    dry_humidity_tendency = 0.0_real64
    if (physics%convection%enabled) then
      if (moist) then
        call dry_convective_adjustment_from_levels(physics%convection, delta_pressure, exner, temperature, &
                                                   dry_temperature_tendency, clipped_humidity, dry_humidity_tendency)
        humidity_rhs = humidity_rhs + dry_humidity_tendency
      else
        call dry_convective_adjustment_from_levels(physics%convection, delta_pressure, exner, temperature, &
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
      call moist_convective_adjustment_from_levels(physics%moist_convection, full_level_pressure, delta_pressure, &
                                                   exner, provisional_temperature, provisional_humidity, &
                                                   moist_temperature_tendency, moist_humidity_tendency, &
                                                   convective_precipitation)
      temperature_rhs = temperature_rhs + moist_temperature_tendency
      humidity_rhs = humidity_rhs + moist_humidity_tendency
    end if
    provisional_temperature = provisional_temperature + interval*moist_temperature_tendency
    provisional_humidity = max(provisional_humidity + interval*moist_humidity_tendency, 0.0_real64)
    if (physics%condensation%enabled) then
      call large_scale_condensation_from_levels(physics%condensation, full_level_pressure, delta_pressure, &
                                                provisional_temperature, provisional_humidity, interval, &
                                                condensation_temperature_tendency, condensation_humidity_tendency, &
                                                large_scale_precipitation)
      temperature_rhs = temperature_rhs + condensation_temperature_tendency
      humidity_rhs = humidity_rhs + condensation_humidity_tendency
      provisional_temperature = provisional_temperature + interval*condensation_temperature_tendency
      provisional_humidity = max(provisional_humidity + interval*condensation_humidity_tendency, 0.0_real64)
    end if
    if (physics%cloud%enabled) then
      cloud_cover = diagnose_cloud_cover(physics%cloud, full_level_pressure, provisional_temperature, &
                                         provisional_humidity, convective_precipitation)
    end if
  end subroutine add_convection_column_tendency_from_levels

  !> Evaluated column by column on the RAW-filtered previous time level; the
  !> precipitation of each column is recorded for the diagnostics and the cloud
  !> cover for the radiation tendency that follows.
  subroutine add_dry_convection_tendency(physics, interval, workspace)
    type(dry_model_physics_config), intent(in) :: physics
    real(real64), intent(in) :: interval
    type(dry_workspace_type), intent(inout) :: workspace
    integer :: i, j

    !$omp parallel do default(shared) private(i, j) schedule(dynamic, 2)
    do j = 1, workspace%ny
      do i = 1, workspace%ring_nlon(j)
        call add_convection_column_tendency_from_levels(physics, interval, &
          workspace%previous_full_level_pressure(i, j, :), workspace%previous_delta_p(i, j, :), &
          workspace%previous_temperature_grid(i, j, :), workspace%previous_humidity_grid(i, j, :), &
          workspace%forcing_temperature(i, j, :), workspace%forcing_humidity(i, j, :), &
          workspace%convective_precipitation(i, j), workspace%large_scale_precipitation(i, j), &
          workspace%cloud_cover(i, j))
      end do
    end do
    !$omp end parallel do
  end subroutine add_dry_convection_tendency

end module dry_convection_tendency
