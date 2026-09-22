!> Projects the accumulated grid-space right-hand side onto spectral space.
!>
!> Every tendency module adds into the shared grid accumulators, so the
!> spectral transforms happen exactly once here regardless of how many physical
!> processes are active.  This is what keeps the transform count of the split
!> tendency equal to that of the original single kernel.
module dry_tendency_projection
  use iso_fortran_env, only: real64
  use harmonics, only: harmonic_transform
  use planet_parameters, only: earth_radius
  use spectral_vector_operators, only: flux_curl_divergence
  use dry_state, only: dry_tendency_type, enforce_dry_spectral_field
  use dry_tendency_workspace, only: dry_workspace_type
  implicit none
  private
  public :: project_dry_tendency

contains

  !> Specific humidity is transformed only when the atmosphere carries moisture;
  !> a dry run keeps its (zero) humidity tendency without an extra transform.
  subroutine project_dry_tendency(transform, truncation, workspace, radiation_enabled, moisture_enabled, rhs)
    type(harmonic_transform), intent(inout) :: transform
    integer, intent(in) :: truncation
    type(dry_workspace_type), intent(inout) :: workspace
    logical, intent(in) :: radiation_enabled, moisture_enabled
    type(dry_tendency_type), intent(inout) :: rhs
    complex(real64), allocatable :: curl_spectral(:, :), divergence_spectral(:, :)
    complex(real64), allocatable :: temporary_spectral(:, :)
    integer :: k, n, m

    !$omp parallel do default(shared) private(k, n, m, curl_spectral, divergence_spectral, temporary_spectral) &
    !$omp   schedule(dynamic, 1)
    do k = 1, workspace%number_of_levels
      call flux_curl_divergence(transform, workspace%forcing_u(:, :, k), workspace%forcing_v(:, :, k), &
                                curl_spectral, divergence_spectral)
      rhs%zeta(:, :, k) = curl_spectral
      rhs%delta(:, :, k) = divergence_spectral

      call transform%grid_to_spectral(workspace%kinetic_geopotential(:, :, k), temporary_spectral)
      do m = 0, truncation
        do n = m, truncation
          rhs%delta(n, m, k) = rhs%delta(n, m, k) + &
            real(n*(n + 1), real64)*temporary_spectral(n, m)/earth_radius**2
        end do
      end do

      call transform%grid_to_spectral(workspace%forcing_temperature(:, :, k), temporary_spectral)
      rhs%temperature(:, :, k) = temporary_spectral
      if (moisture_enabled) then
        call transform%grid_to_spectral(workspace%forcing_humidity(:, :, k), temporary_spectral)
        rhs%specific_humidity(:, :, k) = temporary_spectral
      else
        rhs%specific_humidity(:, :, k) = 0.0_real64
      end if

      call enforce_dry_spectral_field(rhs%zeta(:, :, k), truncation, .true.)
      call enforce_dry_spectral_field(rhs%delta(:, :, k), truncation, .true.)
      call enforce_dry_spectral_field(rhs%temperature(:, :, k), truncation, .false.)
      call enforce_dry_spectral_field(rhs%specific_humidity(:, :, k), truncation, .false.)
    end do
    !$omp end parallel do

    call transform%grid_to_spectral(workspace%forcing_log_ps, temporary_spectral)
    rhs%log_surface_pressure = temporary_spectral
    call enforce_dry_spectral_field(rhs%log_surface_pressure, truncation, .false.)

    ! The surface temperatures stay on the grid: their tendencies are copied, not transformed,
    ! so a sharp land--sea contrast in the heat capacity cannot produce Gibbs ripples.
    if (any(shape(rhs%surface_temperature) /= shape(workspace%forcing_surface_temperature))) then
      error stop 'dry tendency surface temperature has an inconsistent grid shape'
    end if
    if (radiation_enabled) then
      rhs%surface_temperature = workspace%forcing_surface_temperature
      rhs%deep_temperature = workspace%forcing_deep_temperature
    else
      rhs%surface_temperature = 0.0_real64
      rhs%deep_temperature = 0.0_real64
    end if
  end subroutine project_dry_tendency

end module dry_tendency_projection
