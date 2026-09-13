module shallow_water_nonlinear
  use iso_fortran_env, only: real64
  use harmonics, only: harmonic_transform
  use barotropic_vorticity, only: earth_radius, rotation_rate
  implicit none
  private

  public :: compute_shallow_water_nonlinear_tendency
  public :: diagnose_shallow_water_velocity

contains

  subroutine compute_shallow_water_nonlinear_tendency(transform, truncation, zeta, delta, eta, &
                                                      rhs_zeta, rhs_delta, rhs_eta)
    type(harmonic_transform), intent(inout) :: transform
    integer, intent(in) :: truncation
    complex(real64), intent(in) :: zeta(0:, 0:), delta(0:, 0:), eta(0:, 0:)
    complex(real64), allocatable, intent(out) :: rhs_zeta(:, :), rhs_delta(:, :), rhs_eta(:, :)
    complex(real64), allocatable :: kinetic_energy_spec(:, :), flux_term(:, :)
    real(real64), allocatable :: zeta_grid(:, :), eta_grid(:, :), u(:, :), v(:, :)
    real(real64), allocatable :: q(:, :), kinetic_energy(:, :), flux_u(:, :), flux_v(:, :)
    integer, allocatable :: nlon(:)
    integer :: j, n, m

    call diagnose_shallow_water_velocity(transform, truncation, zeta, delta, u, v)
    call transform%spectral_to_grid(zeta, zeta_grid)
    call transform%spectral_to_grid(eta, eta_grid)
    call transform%allocate_field(q)
    call transform%allocate_field(kinetic_energy)
    call transform%allocate_field(flux_u)
    call transform%allocate_field(flux_v)
    q = 0.0_real64
    kinetic_energy = 0.0_real64
    flux_u = 0.0_real64
    flux_v = 0.0_real64
    nlon = transform%get_nlon()

    do j = 1, size(nlon)
      q(1:nlon(j), j) = zeta_grid(1:nlon(j), j) + 2.0_real64*rotation_rate*transform%mu(j)
      kinetic_energy(1:nlon(j), j) = 0.5_real64*(u(1:nlon(j), j)**2 + v(1:nlon(j), j)**2)
      flux_u(1:nlon(j), j) = q(1:nlon(j), j)*u(1:nlon(j), j)
      flux_v(1:nlon(j), j) = q(1:nlon(j), j)*v(1:nlon(j), j)
    end do
    call allocate_spectral_state(truncation, rhs_zeta, rhs_delta, rhs_eta)
    call flux_divergence(transform, flux_u, flux_v, flux_term)
    rhs_zeta = -flux_term

    do j = 1, size(nlon)
      flux_u(1:nlon(j), j) = -q(1:nlon(j), j)*v(1:nlon(j), j)
      flux_v(1:nlon(j), j) = q(1:nlon(j), j)*u(1:nlon(j), j)
    end do
    call flux_divergence(transform, flux_u, flux_v, flux_term)
    rhs_delta = -flux_term
    call transform%grid_to_spectral(kinetic_energy, kinetic_energy_spec)
    do m = 0, truncation
      do n = m, truncation
        rhs_delta(n, m) = rhs_delta(n, m) + &
                          real(n*(n + 1), real64)*kinetic_energy_spec(n, m)/earth_radius**2
      end do
    end do

    do j = 1, size(nlon)
      flux_u(1:nlon(j), j) = eta_grid(1:nlon(j), j)*u(1:nlon(j), j)
      flux_v(1:nlon(j), j) = eta_grid(1:nlon(j), j)*v(1:nlon(j), j)
    end do
    call flux_divergence(transform, flux_u, flux_v, flux_term)
    rhs_eta = -flux_term
    call enforce_tendency_constraints(truncation, rhs_zeta)
    call enforce_tendency_constraints(truncation, rhs_delta)
    call enforce_tendency_constraints(truncation, rhs_eta)
  end subroutine compute_shallow_water_nonlinear_tendency

  subroutine flux_divergence(transform, flux_u, flux_v, divergence)
    type(harmonic_transform), intent(inout) :: transform
    real(real64), intent(in) :: flux_u(:, :), flux_v(:, :)
    complex(real64), allocatable, intent(out) :: divergence(:, :)
    complex(real64), allocatable :: flux_u_spec(:, :), flux_v_cos_spec(:, :)
    real(real64), allocatable :: flux_v_cos(:, :), du_dlambda(:, :), dvcos_dphi(:, :)
    real(real64), allocatable :: unused_one(:, :), unused_two(:, :), divergence_grid(:, :)
    integer, allocatable :: nlon(:)
    integer :: j
    real(real64) :: cosphi

    call transform%allocate_field(flux_v_cos)
    call transform%allocate_field(divergence_grid)
    flux_v_cos = 0.0_real64
    divergence_grid = 0.0_real64
    nlon = transform%get_nlon()
    do j = 1, size(nlon)
      cosphi = sqrt(max(0.0_real64, 1.0_real64 - transform%mu(j)**2))
      flux_v_cos(1:nlon(j), j) = flux_v(1:nlon(j), j)*cosphi
    end do
    call transform%grid_to_spectral(flux_u, flux_u_spec)
    call transform%grid_to_spectral(flux_v_cos, flux_v_cos_spec)
    call transform%gradient_to_grid(flux_u_spec, du_dlambda, unused_one)
    call transform%gradient_to_grid(flux_v_cos_spec, unused_two, dvcos_dphi)
    do j = 1, size(nlon)
      cosphi = sqrt(max(0.0_real64, 1.0_real64 - transform%mu(j)**2))
      divergence_grid(1:nlon(j), j) = &
        (du_dlambda(1:nlon(j), j) + dvcos_dphi(1:nlon(j), j))/(earth_radius*cosphi)
    end do
    call transform%grid_to_spectral(divergence_grid, divergence)
  end subroutine flux_divergence

  subroutine diagnose_shallow_water_velocity(transform, truncation, zeta, delta, u, v)
    type(harmonic_transform), intent(inout) :: transform
    integer, intent(in) :: truncation
    complex(real64), intent(in) :: zeta(0:, 0:), delta(0:, 0:)
    real(real64), allocatable, intent(out) :: u(:, :), v(:, :)
    complex(real64), allocatable :: psi(:, :), chi(:, :)
    real(real64), allocatable :: dpsi_dlambda(:, :), dpsi_dphi(:, :)
    real(real64), allocatable :: dchi_dlambda(:, :), dchi_dphi(:, :)
    integer, allocatable :: nlon(:)
    integer :: j, n, m
    real(real64) :: cosphi

    allocate (psi(0:truncation + 1, 0:truncation), chi(0:truncation + 1, 0:truncation))
    psi = cmplx(0.0_real64, 0.0_real64, kind=real64)
    chi = cmplx(0.0_real64, 0.0_real64, kind=real64)
    do m = 0, truncation
      do n = max(1, m), truncation
        psi(n, m) = -earth_radius**2*zeta(n, m)/real(n*(n + 1), real64)
        chi(n, m) = -earth_radius**2*delta(n, m)/real(n*(n + 1), real64)
      end do
    end do
    call transform%gradient_to_grid(psi, dpsi_dlambda, dpsi_dphi)
    call transform%gradient_to_grid(chi, dchi_dlambda, dchi_dphi)
    call transform%allocate_field(u)
    call transform%allocate_field(v)
    u = 0.0_real64
    v = 0.0_real64
    nlon = transform%get_nlon()
    do j = 1, size(nlon)
      cosphi = sqrt(max(0.0_real64, 1.0_real64 - transform%mu(j)**2))
      u(1:nlon(j), j) = -dpsi_dphi(1:nlon(j), j)/earth_radius + &
                         dchi_dlambda(1:nlon(j), j)/(earth_radius*cosphi)
      v(1:nlon(j), j) = dpsi_dlambda(1:nlon(j), j)/(earth_radius*cosphi) + &
                         dchi_dphi(1:nlon(j), j)/earth_radius
    end do
  end subroutine diagnose_shallow_water_velocity

  subroutine allocate_spectral_state(truncation, zeta, delta, eta)
    integer, intent(in) :: truncation
    complex(real64), allocatable, intent(out) :: zeta(:, :), delta(:, :), eta(:, :)

    allocate (zeta(0:truncation + 1, 0:truncation), &
              delta(0:truncation + 1, 0:truncation), &
              eta(0:truncation + 1, 0:truncation))
  end subroutine allocate_spectral_state

  subroutine enforce_tendency_constraints(truncation, field)
    integer, intent(in) :: truncation
    complex(real64), intent(inout) :: field(0:, 0:)
    integer :: n, m

    field(0, 0) = cmplx(0.0_real64, 0.0_real64, kind=real64)
    field(truncation + 1, :) = cmplx(0.0_real64, 0.0_real64, kind=real64)
    do m = 1, truncation
      do n = 0, m - 1
        field(n, m) = cmplx(0.0_real64, 0.0_real64, kind=real64)
      end do
    end do
  end subroutine enforce_tendency_constraints

end module shallow_water_nonlinear
