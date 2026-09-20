module spectral_vector_operators
  use iso_fortran_env, only: real64
  use harmonics, only: harmonic_transform
  use planet_parameters, only: earth_radius
  implicit none
  private
  public :: diagnose_horizontal_velocity, flux_curl_divergence, flux_divergence, flux_curl

contains

  subroutine flux_curl_divergence(transform, vector_u, vector_v, curl, divergence)
    type(harmonic_transform), intent(inout) :: transform
    real(real64), intent(in) :: vector_u(:, :), vector_v(:, :)
    complex(real64), allocatable, intent(out) :: curl(:, :), divergence(:, :)
    call transform%curl_divergence(vector_u, vector_v, curl, divergence)
    curl = curl/earth_radius
    divergence = divergence/earth_radius
  end subroutine flux_curl_divergence

  subroutine flux_divergence(transform, flux_u, flux_v, divergence)
    type(harmonic_transform), intent(inout) :: transform
    real(real64), intent(in) :: flux_u(:, :), flux_v(:, :)
    complex(real64), allocatable, intent(out) :: divergence(:, :)
    complex(real64), allocatable :: unused_curl(:, :)
    call flux_curl_divergence(transform, flux_u, flux_v, unused_curl, divergence)
  end subroutine flux_divergence

  subroutine flux_curl(transform, vector_u, vector_v, curl)
    type(harmonic_transform), intent(inout) :: transform
    real(real64), intent(in) :: vector_u(:, :), vector_v(:, :)
    complex(real64), allocatable, intent(out) :: curl(:, :)
    complex(real64), allocatable :: unused_divergence(:, :)
    call flux_curl_divergence(transform, vector_u, vector_v, curl, unused_divergence)
  end subroutine flux_curl

  subroutine diagnose_horizontal_velocity(transform, truncation, zeta, delta, u, v)
    type(harmonic_transform), intent(inout) :: transform
    integer, intent(in) :: truncation
    complex(real64), intent(in) :: zeta(0:, 0:), delta(0:, 0:)
    real(real64), allocatable, intent(out) :: u(:, :), v(:, :)
    complex(real64), allocatable :: psi(:, :), chi(:, :)
    integer :: n, m
    allocate (psi(0:truncation + 1, 0:truncation), chi(0:truncation + 1, 0:truncation))
    psi = cmplx(0.0_real64, 0.0_real64, kind=real64)
    chi = cmplx(0.0_real64, 0.0_real64, kind=real64)
    do m = 0, truncation
      do n = max(1, m), truncation
        psi(n, m) = -earth_radius*zeta(n, m)/real(n*(n + 1), real64)
        chi(n, m) = -earth_radius*delta(n, m)/real(n*(n + 1), real64)
      end do
    end do
    call transform%wind_to_grid(psi, chi, u, v)
  end subroutine diagnose_horizontal_velocity

end module spectral_vector_operators
