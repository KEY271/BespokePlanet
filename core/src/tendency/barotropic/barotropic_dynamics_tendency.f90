module barotropic_dynamics_tendency
  use iso_fortran_env, only: real64
  use harmonics, only: harmonic_transform
  use planet_parameters, only: earth_radius, earth_rotation_rate
  use barotropic_state, only: barotropic_state_type, barotropic_tendency_type, enforce_barotropic_constraints
  implicit none
  private

  public :: add_barotropic_dynamics_tendency

contains

  subroutine add_barotropic_dynamics_tendency(transform, truncation, state, rhs)
    type(harmonic_transform), intent(inout) :: transform
    integer, intent(in) :: truncation
    type(barotropic_state_type), intent(in) :: state
    type(barotropic_tendency_type), intent(inout) :: rhs
    complex(real64), allocatable :: psi(:, :), flux_u_spec(:, :), flux_v_spec(:, :), contribution(:, :)
    real(real64), allocatable :: zeta_grid(:, :), u(:, :), v(:, :), q(:, :)
    real(real64), allocatable :: dpsi_dlambda(:, :), dpsi_dphi(:, :)
    real(real64), allocatable :: flux_u(:, :), flux_v(:, :), tendency_grid(:, :)
    real(real64), allocatable :: dflux_u_dlambda(:, :), unused_u(:, :)
    real(real64), allocatable :: unused_v(:, :), dflux_v_dphi(:, :)
    integer, allocatable :: nlon(:)
    integer :: j, n, m
    real(real64) :: cosphi

    allocate (psi(0:truncation + 1, 0:truncation))
    psi = cmplx(0.0_real64, 0.0_real64, kind=real64)
    do m = 0, truncation
      do n = max(1, m), truncation
        psi(n, m) = -earth_radius**2*state%zeta(n, m)/real(n*(n + 1), real64)
      end do
    end do

    call transform%spectral_to_grid(state%zeta, zeta_grid)
    call transform%gradient_to_grid(psi, dpsi_dlambda, dpsi_dphi)
    call transform%allocate_field(u)
    call transform%allocate_field(v)
    call transform%allocate_field(q)
    call transform%allocate_field(flux_u)
    call transform%allocate_field(flux_v)
    u = 0.0_real64
    v = 0.0_real64
    q = 0.0_real64
    flux_u = 0.0_real64
    flux_v = 0.0_real64
    nlon = transform%get_nlon()

    do j = 1, size(nlon)
      cosphi = sqrt(max(0.0_real64, 1.0_real64 - transform%mu(j)**2))
      u(1:nlon(j), j) = -dpsi_dphi(1:nlon(j), j)/earth_radius
      v(1:nlon(j), j) = dpsi_dlambda(1:nlon(j), j)/(earth_radius*cosphi)
      q(1:nlon(j), j) = zeta_grid(1:nlon(j), j) + 2.0_real64*earth_rotation_rate*transform%mu(j)
      flux_u(1:nlon(j), j) = u(1:nlon(j), j)*q(1:nlon(j), j)
      flux_v(1:nlon(j), j) = v(1:nlon(j), j)*q(1:nlon(j), j)*cosphi
    end do

    call transform%grid_to_spectral(flux_u, flux_u_spec)
    call transform%grid_to_spectral(flux_v, flux_v_spec)
    call transform%gradient_to_grid(flux_u_spec, dflux_u_dlambda, unused_u)
    call transform%gradient_to_grid(flux_v_spec, unused_v, dflux_v_dphi)
    call transform%allocate_field(tendency_grid)
    tendency_grid = 0.0_real64
    do j = 1, size(nlon)
      cosphi = sqrt(max(0.0_real64, 1.0_real64 - transform%mu(j)**2))
      tendency_grid(1:nlon(j), j) = &
        -(dflux_u_dlambda(1:nlon(j), j) + dflux_v_dphi(1:nlon(j), j))/(earth_radius*cosphi)
    end do

    call transform%grid_to_spectral(tendency_grid, contribution)
    call enforce_barotropic_constraints(contribution, truncation)
    rhs%zeta = rhs%zeta + contribution
  end subroutine add_barotropic_dynamics_tendency

end module barotropic_dynamics_tendency
