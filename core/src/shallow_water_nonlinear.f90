module shallow_water_nonlinear
  use iso_fortran_env, only: real64
  use harmonics, only: harmonic_transform
  use planet_parameters, only: earth_radius, rotation_rate => earth_rotation_rate
  use spectral_vector_operators, only: diagnose_horizontal_velocity, &
                                       flux_divergence, flux_curl, flux_curl_divergence
  implicit none
  private

  public :: compute_shallow_water_nonlinear_tendency
  public :: flux_divergence
  public :: flux_curl
  public :: flux_curl_divergence

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

    call diagnose_horizontal_velocity(transform, truncation, zeta, delta, u, v)
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
    ! Both come from the same two analyses: d(zeta)/dt = -div(q u), d(delta)/dt = curl(q u) - ...
    call flux_curl_divergence(transform, flux_u, flux_v, rhs_delta, flux_term)
    rhs_zeta = -flux_term
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
