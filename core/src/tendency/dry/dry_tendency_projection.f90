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

  subroutine project_dry_tendency(transform, truncation, workspace, radiation_enabled, rhs)
    type(harmonic_transform), intent(inout) :: transform
    integer, intent(in) :: truncation
    type(dry_workspace_type), intent(inout) :: workspace
    logical, intent(in) :: radiation_enabled
    type(dry_tendency_type), intent(inout) :: rhs
    complex(real64), allocatable :: curl_spectral(:, :), divergence_spectral(:, :)
    complex(real64), allocatable :: temporary_spectral(:, :)
    integer :: k, n, m

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

      call enforce_dry_spectral_field(rhs%zeta(:, :, k), truncation, .true.)
      call enforce_dry_spectral_field(rhs%delta(:, :, k), truncation, .true.)
      call enforce_dry_spectral_field(rhs%temperature(:, :, k), truncation, .false.)
    end do

    call transform%grid_to_spectral(workspace%forcing_log_ps, temporary_spectral)
    rhs%log_surface_pressure = temporary_spectral
    call enforce_dry_spectral_field(rhs%log_surface_pressure, truncation, .false.)

    if (radiation_enabled) then
      call transform%grid_to_spectral(workspace%forcing_surface_temperature, temporary_spectral)
      rhs%surface_temperature = temporary_spectral
      call transform%grid_to_spectral(workspace%forcing_deep_temperature, temporary_spectral)
      rhs%deep_temperature = temporary_spectral
      call enforce_dry_spectral_field(rhs%surface_temperature, truncation, .false.)
      call enforce_dry_spectral_field(rhs%deep_temperature, truncation, .false.)
    else
      rhs%surface_temperature = 0.0_real64
      rhs%deep_temperature = 0.0_real64
    end if
  end subroutine project_dry_tendency

end module dry_tendency_projection
