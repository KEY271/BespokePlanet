!> Precipitation recharge and overflow runoff of the land bucket.
module dry_bucket_tendency
  use iso_fortran_env, only: real64
  use dry_physics_config, only: bucket_config
  use dry_tendency_workspace, only: dry_workspace_type
  use land_bucket, only: advance_bucket
  implicit none
  private
  public :: add_dry_bucket_tendency

contains

  subroutine add_dry_bucket_tendency(config, interval, workspace)
    type(bucket_config), intent(in) :: config
    real(real64), intent(in) :: interval
    type(dry_workspace_type), intent(inout) :: workspace
    real(real64) :: precipitation
    integer :: i, j

    if (.not. config%enabled) return
    if (config%capacity <= 0.0_real64) error stop 'land-bucket capacity must be positive'
    if (interval <= 0.0_real64) error stop 'land-bucket interval must be positive'
    !$omp parallel do default(shared) private(i, j, precipitation) schedule(dynamic, 2)
    do j = 1, workspace%ny
      do i = 1, workspace%ring_nlon(j)
        if (workspace%land_fraction(i, j) > 0.0_real64) then
          precipitation = workspace%convective_precipitation(i, j) + &
                          workspace%large_scale_precipitation(i, j)
          call advance_bucket(workspace%previous_surface_water(i, j), precipitation, &
            workspace%land_evaporation(i, j), config%capacity, interval, &
            workspace%forcing_surface_water(i, j), workspace%runoff(i, j))
          workspace%water_budget_residual(i, j) = &
            workspace%forcing_surface_water(i, j) - &
            (precipitation - workspace%land_evaporation(i, j) - workspace%runoff(i, j))
        end if
      end do
    end do
    !$omp end parallel do
  end subroutine add_dry_bucket_tendency

end module dry_bucket_tendency
