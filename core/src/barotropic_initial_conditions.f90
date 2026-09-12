module barotropic_initial_conditions
  use iso_fortran_env, only: real64, int64
  use harmonics, only: harmonic_transform
  use barotropic_vorticity, only: earth_radius
  implicit none
  private

  public :: single_harmonic_vorticity
  public :: rossby_haurwitz_vorticity
  public :: random_low_wavenumber_vorticity

  integer, parameter, public :: single_mode_degree = 4
  integer, parameter, public :: single_mode_order = 2
  real(real64), parameter, public :: prescribed_maximum_vorticity = 2.0e-5_real64
  integer, parameter, public :: rossby_haurwitz_wave_number = 4
  real(real64), parameter, public :: rossby_haurwitz_omega = 7.848e-6_real64
  real(real64), parameter, public :: rossby_haurwitz_wave_amplitude = 7.848e-6_real64
  integer, parameter, public :: random_minimum_degree = 8
  integer, parameter, public :: random_maximum_degree = 12

contains

  subroutine single_harmonic_vorticity(transform, truncation, zeta)
    type(harmonic_transform), intent(in) :: transform
    integer, intent(in) :: truncation
    complex(real64), allocatable, intent(out) :: zeta(:, :)

    if (truncation < single_mode_degree) then
      error stop 'single harmonic initial condition requires a larger T'
    end if
    allocate (zeta(0:truncation + 1, 0:truncation))
    zeta = cmplx(0.0_real64, 0.0_real64, kind=real64)
    zeta(single_mode_degree, single_mode_order) = &
      cmplx(1.0_real64, 0.0_real64, kind=real64)
    call scale_to_maximum(transform, zeta, prescribed_maximum_vorticity)
  end subroutine single_harmonic_vorticity

  subroutine rossby_haurwitz_vorticity(transform, truncation, zeta)
    type(harmonic_transform), intent(in) :: transform
    integer, intent(in) :: truncation
    complex(real64), allocatable, intent(out) :: zeta(:, :)
    real(real64), allocatable :: psi_grid(:, :)
    complex(real64), allocatable :: psi(:, :)
    integer, allocatable :: nlon(:)
    integer :: j, k, n, m
    real(real64) :: pi, lambda, mu, cosphi

    if (truncation < rossby_haurwitz_wave_number + 1) then
      error stop 'Rossby-Haurwitz R=4 initial condition requires T >= 5'
    end if

    call transform%allocate_field(psi_grid)
    psi_grid = 0.0_real64
    nlon = transform%get_nlon()
    pi = acos(-1.0_real64)
    do j = 1, size(nlon)
      mu = transform%mu(j)
      cosphi = sqrt(max(0.0_real64, 1.0_real64 - mu**2))
      do k = 1, nlon(j)
        lambda = 2.0_real64*pi*real(k - 1, real64)/real(nlon(j), real64)
        psi_grid(k, j) = -earth_radius**2*rossby_haurwitz_omega*mu + &
                         earth_radius**2*rossby_haurwitz_wave_amplitude* &
                         cosphi**rossby_haurwitz_wave_number*mu* &
                         cos(real(rossby_haurwitz_wave_number, real64)*lambda)
      end do
    end do

    call transform%grid_to_spectral(psi_grid, psi)
    allocate (zeta(0:truncation + 1, 0:truncation))
    zeta = cmplx(0.0_real64, 0.0_real64, kind=real64)
    do m = 0, truncation
      do n = max(1, m), truncation
        zeta(n, m) = -real(n*(n + 1), real64)*psi(n, m)/earth_radius**2
      end do
    end do
  end subroutine rossby_haurwitz_vorticity

  subroutine random_low_wavenumber_vorticity(transform, truncation, seed, zeta)
    type(harmonic_transform), intent(in) :: transform
    integer, intent(in) :: truncation
    integer(int64), intent(in) :: seed
    complex(real64), allocatable, intent(out) :: zeta(:, :)
    integer(int64) :: state
    integer :: n, m
    real(real64) :: real_part, imag_part

    if (truncation < random_maximum_degree) then
      error stop 'random low-wavenumber initial condition requires a larger T'
    end if
    if (seed <= 0_int64) error stop 'random low-wavenumber seed must be positive'

    allocate (zeta(0:truncation + 1, 0:truncation))
    zeta = cmplx(0.0_real64, 0.0_real64, kind=real64)
    state = modulo(seed, 2147483647_int64)
    if (state == 0_int64) state = 1_int64
    do n = random_minimum_degree, random_maximum_degree
      do m = 0, n
        real_part = 2.0_real64*next_uniform(state) - 1.0_real64
        if (m == 0) then
          zeta(n, m) = cmplx(real_part, 0.0_real64, kind=real64)
        else
          imag_part = 2.0_real64*next_uniform(state) - 1.0_real64
          zeta(n, m) = cmplx(real_part, imag_part, kind=real64)/sqrt(2.0_real64)
        end if
      end do
    end do
    call scale_to_maximum(transform, zeta, prescribed_maximum_vorticity)
  end subroutine random_low_wavenumber_vorticity

  subroutine scale_to_maximum(transform, coefficients, target_maximum)
    type(harmonic_transform), intent(in) :: transform
    complex(real64), intent(inout) :: coefficients(0:, 0:)
    real(real64), intent(in) :: target_maximum
    real(real64), allocatable :: field(:, :)
    integer, allocatable :: nlon(:)
    integer :: j
    real(real64) :: maximum

    call transform%spectral_to_grid(coefficients, field)
    nlon = transform%get_nlon()
    maximum = 0.0_real64
    do j = 1, size(nlon)
      maximum = max(maximum, maxval(abs(field(1:nlon(j), j))))
    end do
    if (maximum <= 0.0_real64) error stop 'cannot scale a zero initial condition'
    coefficients = coefficients*(target_maximum/maximum)
  end subroutine scale_to_maximum

  real(real64) function next_uniform(state) result(value)
    integer(int64), intent(inout) :: state

    state = modulo(16807_int64*state, 2147483647_int64)
    value = real(state, real64)/2147483647.0_real64
  end function next_uniform

end module barotropic_initial_conditions
