!> Radiative heating of the atmosphere together with the surface and deep
!> ground energy budget.
!>
!> The instantaneous diagnostics that the radiation case aggregates are a
!> by-product of this tendency, so they are produced here rather than being
!> recomputed later.  This module accumulates nothing over time and writes no
!> files; the case-side collector owns the daily and monthly means.
module dry_radiation_tendency
  use iso_fortran_env, only: real64
  use harmonics, only: harmonic_transform
  use dry_physics_config, only: radiation_config
  use dry_radiation, only: radiation_tendency, radiation_diagnostics
  use dry_tendency_workspace, only: dry_workspace_type
  implicit none
  private
  public :: add_dry_radiation_column_tendency, add_dry_radiation_tendency

contains

  subroutine add_dry_radiation_column_tendency(config, pressure_half, temperature, surface_temperature, &
                                                deep_temperature, lowest_u, lowest_v, sin_latitude, &
                                                longitude, time_seconds, temperature_rhs, &
                                                surface_temperature_rhs, deep_temperature_rhs, &
                                                incoming_shortwave, reflected_shortwave, outgoing_longwave)
    type(radiation_config), intent(in) :: config
    real(real64), intent(in) :: pressure_half(0:), temperature(:)
    real(real64), intent(in) :: surface_temperature, deep_temperature, lowest_u, lowest_v
    real(real64), intent(in) :: sin_latitude, longitude, time_seconds
    real(real64), intent(inout) :: temperature_rhs(:)
    real(real64), intent(inout) :: surface_temperature_rhs, deep_temperature_rhs
    real(real64), intent(out) :: incoming_shortwave, reflected_shortwave, outgoing_longwave
    real(real64) :: temperature_contribution(size(temperature))
    real(real64) :: surface_contribution, deep_contribution

    call radiation_tendency(config, pressure_half, temperature, surface_temperature, deep_temperature, &
      lowest_u, lowest_v, sin_latitude, longitude, time_seconds, temperature_contribution, &
      surface_contribution, deep_contribution, incoming_shortwave, reflected_shortwave, outgoing_longwave)
    temperature_rhs = temperature_rhs + temperature_contribution
    surface_temperature_rhs = surface_temperature_rhs + surface_contribution
    deep_temperature_rhs = deep_temperature_rhs + deep_contribution
  end subroutine add_dry_radiation_column_tendency

  !> Radiation and surface exchange use one consistent RAW-filtered
  !> previous-time column.  The diagnostics are sampled from the current state.
  subroutine add_dry_radiation_tendency(config, transform, workspace, diagnostics)
    type(radiation_config), intent(in) :: config
    type(harmonic_transform), intent(in) :: transform
    type(dry_workspace_type), intent(inout) :: workspace
    type(radiation_diagnostics), intent(out), optional :: diagnostics
    integer :: i, j, k, levels
    real(real64) :: longitude, incoming_shortwave, reflected_shortwave, outgoing_longwave
    real(real64) :: area_weight, atmospheric_mass, temperature_mass_sum, kinetic_energy_mass_sum

    levels = workspace%number_of_levels
    if (present(diagnostics)) then
      diagnostics%time_seconds = workspace%evaluation_time
      diagnostics%mean_atmospheric_temperature = 0.0_real64
      diagnostics%mean_surface_temperature = 0.0_real64
      diagnostics%mean_deep_temperature = 0.0_real64
      diagnostics%mean_kinetic_energy = 0.0_real64
      diagnostics%mean_surface_pressure = 0.0_real64
      diagnostics%mean_incoming_shortwave = 0.0_real64
      diagnostics%mean_reflected_shortwave = 0.0_real64
      diagnostics%mean_outgoing_longwave = 0.0_real64
      diagnostics%surface_temperature = workspace%surface_temperature_grid
      diagnostics%surface_pressure = workspace%ps
      allocate (diagnostics%zonal_temperature(workspace%ny, levels))
      allocate (diagnostics%zonal_u(workspace%ny, levels), diagnostics%zonal_v(workspace%ny, levels))
      allocate (diagnostics%zonal_uv(workspace%ny, levels), diagnostics%zonal_vt(workspace%ny, levels))
      diagnostics%zonal_temperature = 0.0_real64
      diagnostics%zonal_u = 0.0_real64
      diagnostics%zonal_v = 0.0_real64
      diagnostics%zonal_uv = 0.0_real64
      diagnostics%zonal_vt = 0.0_real64
      atmospheric_mass = 0.0_real64
      temperature_mass_sum = 0.0_real64
      kinetic_energy_mass_sum = 0.0_real64
    end if
    do j = 1, workspace%ny
      do i = 1, workspace%ring_nlon(j)
        longitude = 2.0_real64*acos(-1.0_real64)*real(i - 1, real64)/real(workspace%ring_nlon(j), real64)
        call add_dry_radiation_column_tendency(config, workspace%previous_pressure_half(i, j, :), &
          workspace%previous_temperature_grid(i, j, :), &
          workspace%previous_surface_temperature_grid(i, j), &
          workspace%previous_deep_temperature_grid(i, j), &
          workspace%previous_u(i, j, levels), workspace%previous_v(i, j, levels), transform%mu(j), &
          longitude, workspace%evaluation_time, &
          workspace%forcing_temperature(i, j, :), workspace%forcing_surface_temperature(i, j), &
          workspace%forcing_deep_temperature(i, j), incoming_shortwave, reflected_shortwave, &
          outgoing_longwave)
        if (present(diagnostics)) then
          area_weight = 0.5_real64*workspace%gaussian_weights(j)/real(workspace%ring_nlon(j), real64)
          diagnostics%mean_surface_temperature = diagnostics%mean_surface_temperature + &
            area_weight*workspace%surface_temperature_grid(i, j)
          diagnostics%mean_deep_temperature = diagnostics%mean_deep_temperature + &
            area_weight*workspace%deep_temperature_grid(i, j)
          diagnostics%mean_surface_pressure = diagnostics%mean_surface_pressure + &
            area_weight*workspace%ps(i, j)
          diagnostics%mean_incoming_shortwave = diagnostics%mean_incoming_shortwave + &
            area_weight*incoming_shortwave
          diagnostics%mean_reflected_shortwave = diagnostics%mean_reflected_shortwave + &
            area_weight*reflected_shortwave
          diagnostics%mean_outgoing_longwave = diagnostics%mean_outgoing_longwave + &
            area_weight*outgoing_longwave
          do k = 1, levels
            diagnostics%zonal_temperature(j, k) = diagnostics%zonal_temperature(j, k) + &
              workspace%temperature_grid(i, j, k)/real(workspace%ring_nlon(j), real64)
            diagnostics%zonal_u(j, k) = diagnostics%zonal_u(j, k) + &
              workspace%u(i, j, k)/real(workspace%ring_nlon(j), real64)
            diagnostics%zonal_v(j, k) = diagnostics%zonal_v(j, k) + &
              workspace%v(i, j, k)/real(workspace%ring_nlon(j), real64)
            diagnostics%zonal_uv(j, k) = diagnostics%zonal_uv(j, k) + &
              workspace%u(i, j, k)*workspace%v(i, j, k)/real(workspace%ring_nlon(j), real64)
            diagnostics%zonal_vt(j, k) = diagnostics%zonal_vt(j, k) + &
              workspace%v(i, j, k)*workspace%temperature_grid(i, j, k)/real(workspace%ring_nlon(j), real64)
            atmospheric_mass = atmospheric_mass + area_weight*workspace%delta_p(i, j, k)
            temperature_mass_sum = temperature_mass_sum + &
              area_weight*workspace%delta_p(i, j, k)*workspace%temperature_grid(i, j, k)
            kinetic_energy_mass_sum = kinetic_energy_mass_sum + area_weight*workspace%delta_p(i, j, k)* &
              0.5_real64*(workspace%u(i, j, k)**2 + workspace%v(i, j, k)**2)
          end do
        end if
      end do
    end do
    if (present(diagnostics)) then
      diagnostics%mean_atmospheric_temperature = temperature_mass_sum/atmospheric_mass
      diagnostics%mean_kinetic_energy = kinetic_energy_mass_sum/atmospheric_mass
    end if
  end subroutine add_dry_radiation_tendency

end module dry_radiation_tendency
