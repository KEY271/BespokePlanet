!> Adiabatic dry dynamics: the contribution of advection, the pressure-gradient
!> term and the thermodynamic equation to the grid-space right-hand side.
!>
!> This is the first tendency the evaluator adds, so it writes the dynamical
!> part into the accumulators that the physical tendencies then add to.
module dry_dynamics_tendency
  use iso_fortran_env, only: real64
  use harmonics, only: harmonic_transform
  use planet_parameters, only: earth_radius
  use dry_vertical_coordinate, only: dry_air_gas_constant, dry_air_kappa
  use dry_state, only: dry_state_type
  use dry_tendency_workspace, only: dry_workspace_type
  implicit none
  private
  public :: add_dry_dynamics_tendency

contains

  !> With the momentum forcing F = -(vertical advection) - R T (pressure-gradient term),
  !>   d(zeta)/dt = -div((zeta+f) u) + curl F,   d(delta)/dt = curl((zeta+f) u) + div F - lap(K+Phi).
  !> Since -div(A u, A v) = curl(A v, -A u) and curl(A u, A v) = div(A v, -A u), both
  !> tendencies are the curl and divergence of the single vector stored in
  !> forcing_u/forcing_v, which costs one pair of grid-to-spectral transforms.
  !>
  !> The grid loops are in the private procedures below, whose explicit-shape dummies
  !> let the compiler treat one level as the contiguous, non-aliasing arrays that the
  !> single kernel used before the split.
  subroutine add_dry_dynamics_tendency(transform, state, workspace, maximum_speed)
    type(harmonic_transform), intent(inout) :: transform
    type(dry_state_type), intent(in) :: state
    type(dry_workspace_type), intent(inout) :: workspace
    real(real64), intent(out) :: maximum_speed
    integer :: k, levels
    real(real64) :: maximum_speed_squared, level_maximum_speed_squared

    levels = workspace%number_of_levels
    maximum_speed_squared = 0.0_real64
    do k = 1, levels
      call add_momentum_forcing(workspace%nx, workspace%ny, workspace%ring_nlon, transform%mu, &
        workspace%active_rotation_rate, workspace%zeta_grid(:, :, k), workspace%u(:, :, k), &
        workspace%v(:, :, k), workspace%vertical_u(:, :, k), workspace%vertical_v(:, :, k), &
        workspace%temperature_grid(:, :, k), workspace%pressure_gradient_u(:, :, k), &
        workspace%pressure_gradient_v(:, :, k), workspace%geopotential(:, :, k), &
        workspace%forcing_u(:, :, k), workspace%forcing_v(:, :, k), &
        workspace%kinetic_geopotential(:, :, k), level_maximum_speed_squared)
      maximum_speed_squared = max(maximum_speed_squared, level_maximum_speed_squared)

      call transform%gradient_to_grid(state%temperature(:, :, k), workspace%dtdlambda, workspace%dtdphi)
      call add_thermodynamic_forcing(workspace%nx, workspace%ny, workspace%ring_nlon, transform%mu, &
        workspace%u(:, :, k), workspace%v(:, :, k), workspace%dtdlambda, workspace%dtdphi, &
        workspace%temperature_grid(:, :, k), workspace%pressure_gradient_u(:, :, k), &
        workspace%pressure_gradient_v(:, :, k), workspace%vertical_t(:, :, k), &
        workspace%layer_l(:, :, k), workspace%alpha(:, :, k), workspace%delta_p(:, :, k), &
        workspace%cumulative(:, :, k - 1), workspace%mass_divergence(:, :, k), &
        workspace%forcing_temperature(:, :, k))
    end do
    maximum_speed = sqrt(maximum_speed_squared)

    call add_surface_pressure_forcing(workspace%nx, workspace%ny, workspace%ring_nlon, &
      workspace%cumulative(:, :, levels), workspace%ps, workspace%forcing_log_ps)
  end subroutine add_dry_dynamics_tendency

  subroutine add_momentum_forcing(nx, ny, ring_nlon, mu, rotation_rate, zeta_grid, u, v, &
                                  vertical_u, vertical_v, temperature_grid, pressure_gradient_u, &
                                  pressure_gradient_v, geopotential, forcing_u, forcing_v, &
                                  kinetic_geopotential, maximum_speed_squared)
    integer, intent(in) :: nx, ny, ring_nlon(ny)
    real(real64), intent(in) :: mu(ny), rotation_rate
    real(real64), intent(in) :: zeta_grid(nx, ny), u(nx, ny), v(nx, ny)
    real(real64), intent(in) :: vertical_u(nx, ny), vertical_v(nx, ny), temperature_grid(nx, ny)
    real(real64), intent(in) :: pressure_gradient_u(nx, ny), pressure_gradient_v(nx, ny)
    real(real64), intent(in) :: geopotential(nx, ny)
    real(real64), intent(inout) :: forcing_u(nx, ny), forcing_v(nx, ny), kinetic_geopotential(nx, ny)
    real(real64), intent(out) :: maximum_speed_squared
    integer :: i, j
    real(real64) :: absolute_vorticity, kinetic

    maximum_speed_squared = 0.0_real64
    do j = 1, ny
      do i = 1, ring_nlon(j)
        absolute_vorticity = zeta_grid(i, j) + 2.0_real64*rotation_rate*mu(j)
        kinetic = 0.5_real64*(u(i, j)**2 + v(i, j)**2)
        maximum_speed_squared = max(maximum_speed_squared, 2.0_real64*kinetic)
        forcing_u(i, j) = forcing_u(i, j) + absolute_vorticity*v(i, j) - vertical_u(i, j) - &
          dry_air_gas_constant*temperature_grid(i, j)*pressure_gradient_u(i, j)
        forcing_v(i, j) = forcing_v(i, j) - absolute_vorticity*u(i, j) - vertical_v(i, j) - &
          dry_air_gas_constant*temperature_grid(i, j)*pressure_gradient_v(i, j)
        kinetic_geopotential(i, j) = kinetic_geopotential(i, j) + kinetic + geopotential(i, j)
      end do
    end do
  end subroutine add_momentum_forcing

  subroutine add_thermodynamic_forcing(nx, ny, ring_nlon, mu, u, v, dtdlambda, dtdphi, &
                                       temperature_grid, pressure_gradient_u, pressure_gradient_v, &
                                       vertical_t, layer_l, alpha, delta_p, cumulative_above, &
                                       mass_divergence, forcing_temperature)
    integer, intent(in) :: nx, ny, ring_nlon(ny)
    real(real64), intent(in) :: mu(ny)
    real(real64), intent(in) :: u(nx, ny), v(nx, ny), dtdlambda(nx, ny), dtdphi(nx, ny)
    real(real64), intent(in) :: temperature_grid(nx, ny)
    real(real64), intent(in) :: pressure_gradient_u(nx, ny), pressure_gradient_v(nx, ny)
    real(real64), intent(in) :: vertical_t(nx, ny), layer_l(nx, ny), alpha(nx, ny), delta_p(nx, ny)
    real(real64), intent(in) :: cumulative_above(nx, ny), mass_divergence(nx, ny)
    real(real64), intent(inout) :: forcing_temperature(nx, ny)
    integer :: i, j
    real(real64) :: cosphi, thermodynamic_q

    do j = 1, ny
      cosphi = sqrt(max(0.0_real64, 1.0_real64 - mu(j)**2))
      do i = 1, ring_nlon(j)
        thermodynamic_q = u(i, j)*pressure_gradient_u(i, j) + v(i, j)*pressure_gradient_v(i, j) - &
          (layer_l(i, j)*cumulative_above(i, j) + alpha(i, j)*mass_divergence(i, j))/delta_p(i, j)
        forcing_temperature(i, j) = forcing_temperature(i, j) + &
          (-u(i, j)*dtdlambda(i, j)/(earth_radius*cosphi) - &
           v(i, j)*dtdphi(i, j)/earth_radius - vertical_t(i, j) + &
           dry_air_kappa*temperature_grid(i, j)*thermodynamic_q)
      end do
    end do
  end subroutine add_thermodynamic_forcing

  subroutine add_surface_pressure_forcing(nx, ny, ring_nlon, cumulative_bottom, ps, forcing_log_ps)
    integer, intent(in) :: nx, ny, ring_nlon(ny)
    real(real64), intent(in) :: cumulative_bottom(nx, ny), ps(nx, ny)
    real(real64), intent(inout) :: forcing_log_ps(nx, ny)
    integer :: i, j

    do j = 1, ny
      do i = 1, ring_nlon(j)
        forcing_log_ps(i, j) = forcing_log_ps(i, j) - cumulative_bottom(i, j)/ps(i, j)
      end do
    end do
  end subroutine add_surface_pressure_forcing

end module dry_dynamics_tendency
