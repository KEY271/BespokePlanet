module shallow_water_initial_conditions
  use iso_fortran_env, only: real64
  use harmonics, only: harmonic_transform
  implicit none
  private

  public :: isolated_height_mountain
  public :: single_harmonic_height

  real(real64), parameter, public :: mountain_height_metres = 100.0_real64
  real(real64), parameter, public :: mountain_longitude_degrees = 0.0_real64
  real(real64), parameter, public :: mountain_latitude_degrees = 0.0_real64
  real(real64), parameter, public :: mountain_angular_radius_degrees = 10.0_real64
  integer, parameter, public :: height_mode_degree = 4
  integer, parameter, public :: height_mode_order = 2
  real(real64), parameter, public :: height_mode_amplitude_metres = 100.0_real64

contains

  subroutine isolated_height_mountain(transform, truncation, zeta, delta, eta)
    type(harmonic_transform), intent(in) :: transform
    integer, intent(in) :: truncation
    complex(real64), allocatable, intent(out) :: zeta(:, :), delta(:, :), eta(:, :)
    real(real64), allocatable :: eta_grid(:, :)
    integer, allocatable :: nlon(:)
    integer :: j, k
    real(real64) :: pi, longitude, centre_longitude, centre_latitude
    real(real64) :: cosphi, cosine_distance, angular_distance, angular_radius

    if (truncation < 1) error stop 'isolated mountain requires positive T'
    call zero_state(truncation, zeta, delta)
    call transform%allocate_field(eta_grid)
    eta_grid = 0.0_real64
    nlon = transform%get_nlon()
    pi = acos(-1.0_real64)
    centre_longitude = mountain_longitude_degrees*pi/180.0_real64
    centre_latitude = mountain_latitude_degrees*pi/180.0_real64
    angular_radius = mountain_angular_radius_degrees*pi/180.0_real64
    do j = 1, size(nlon)
      cosphi = sqrt(max(0.0_real64, 1.0_real64 - transform%mu(j)**2))
      do k = 1, nlon(j)
        longitude = 2.0_real64*pi*real(k - 1, real64)/real(nlon(j), real64)
        cosine_distance = sin(centre_latitude)*transform%mu(j) + &
                          cos(centre_latitude)*cosphi*cos(longitude - centre_longitude)
        angular_distance = acos(max(-1.0_real64, min(1.0_real64, cosine_distance)))
        eta_grid(k, j) = mountain_height_metres*exp(-(angular_distance/angular_radius)**2)
      end do
    end do
    call transform%grid_to_spectral(eta_grid, eta)
    eta(0, 0) = cmplx(0.0_real64, 0.0_real64, kind=real64)
    eta(truncation + 1, :) = cmplx(0.0_real64, 0.0_real64, kind=real64)
  end subroutine isolated_height_mountain

  subroutine single_harmonic_height(transform, truncation, zeta, delta, eta)
    type(harmonic_transform), intent(in) :: transform
    integer, intent(in) :: truncation
    complex(real64), allocatable, intent(out) :: zeta(:, :), delta(:, :), eta(:, :)
    real(real64), allocatable :: eta_grid(:, :)
    integer, allocatable :: nlon(:)
    integer :: j
    real(real64) :: maximum

    if (truncation < height_mode_degree) then
      error stop 'single harmonic height initial condition requires a larger T'
    end if
    call zero_state(truncation, zeta, delta)
    allocate (eta(0:truncation + 1, 0:truncation))
    eta = cmplx(0.0_real64, 0.0_real64, kind=real64)
    eta(height_mode_degree, height_mode_order) = cmplx(1.0_real64, 0.0_real64, kind=real64)
    call transform%spectral_to_grid(eta, eta_grid)
    nlon = transform%get_nlon()
    maximum = 0.0_real64
    do j = 1, size(nlon)
      maximum = max(maximum, maxval(abs(eta_grid(1:nlon(j), j))))
    end do
    if (maximum <= 0.0_real64) error stop 'single harmonic height initial condition is zero'
    eta = eta*(height_mode_amplitude_metres/maximum)
  end subroutine single_harmonic_height

  subroutine zero_state(truncation, zeta, delta)
    integer, intent(in) :: truncation
    complex(real64), allocatable, intent(out) :: zeta(:, :), delta(:, :)

    allocate (zeta(0:truncation + 1, 0:truncation), delta(0:truncation + 1, 0:truncation))
    zeta = cmplx(0.0_real64, 0.0_real64, kind=real64)
    delta = cmplx(0.0_real64, 0.0_real64, kind=real64)
  end subroutine zero_state

end module shallow_water_initial_conditions
