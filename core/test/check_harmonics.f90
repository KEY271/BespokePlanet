program check_harmonics
  use iso_fortran_env, only: real64
  use harmonics, only: harmonic_transform
  implicit none

  integer, parameter :: T = 63
  real(real64), parameter :: tolerance = 1.0e-10_real64
  type(harmonic_transform) :: transform
  integer, allocatable :: nlon(:)
  integer :: mmax_all

  ! Reinitialization must release and replace the FFTW plans cleanly.
  call transform%init(7)
  call transform%init(T)
  nlon = transform%get_nlon()
  mmax_all = minval(nlon)/2 - 1

  call check_constant()
  call check_mode(1, 0, cmplx(1.0_real64, 0.0_real64, kind=real64))
  call check_mode(2, 1, cmplx(0.75_real64, 0.0_real64, kind=real64))
  call check_mode(3, 2, cmplx(0.0_real64, -0.5_real64, kind=real64))
  call check_mode(8, 3, cmplx(0.25_real64, 0.5_real64, kind=real64))
  call check_mode(T, 0, cmplx(1.0_real64, 0.0_real64, kind=real64))
  call check_mode(T, mmax_all, cmplx(-0.25_real64, 0.5_real64, kind=real64))
  call check_mode(T, T, cmplx(0.125_real64, -0.25_real64, kind=real64))
  call check_spectral_extent()
  call check_extension_mode()
  call check_analytic_gradient(0, 'constant')
  call check_analytic_gradient(1, 'sin(latitude)')
  call check_analytic_gradient(2, 'cos(latitude) cos(longitude)')
  call check_analytic_gradient(3, 'cos(latitude) sin(longitude)')
  call check_analytic_gradient(4, 'cos(latitude)^2 cos(2 longitude)')
  call check_wind_curl_divergence(1, 0, 2, 1)
  call check_wind_curl_divergence(8, 3, 5, 0)
  call check_wind_curl_divergence(T, mmax_all, T - 1, 2)
  call check_wind_curl_divergence(T, T, T, 0)
  call check_solid_body_rotation()

contains

  !> wind_to_grid must agree with the explicit gradients of psi and chi, and
  !> curl_divergence of that wind must return the Laplacians -n(n+1) psi and -n(n+1) chi.
  subroutine check_wind_curl_divergence(n_psi, m_psi, n_chi, m_chi)
    integer, intent(in) :: n_psi, m_psi, n_chi, m_chi
    complex(real64), allocatable :: psi(:, :), chi(:, :), curl(:, :), divergence(:, :)
    complex(real64), allocatable :: expected_curl(:, :), expected_divergence(:, :)
    real(real64), allocatable :: u(:, :), v(:, :), dpsi_dlambda(:, :), dpsi_dphi(:, :)
    real(real64), allocatable :: dchi_dlambda(:, :), dchi_dphi(:, :)
    real(real64) :: cosphi, wind_scale, wind_err, curl_err, divergence_err
    integer :: j, nl

    allocate (psi(0:T + 1, 0:T), chi(0:T + 1, 0:T))
    psi = cmplx(0.0_real64, 0.0_real64, kind=real64)
    chi = cmplx(0.0_real64, 0.0_real64, kind=real64)
    ! A real grid field has real zonal-mean (m = 0) coefficients.
    psi(n_psi, m_psi) = cmplx(0.75_real64, merge(0.0_real64, -0.5_real64, m_psi == 0), kind=real64)
    chi(n_chi, m_chi) = cmplx(-0.25_real64, merge(0.0_real64, 0.5_real64, m_chi == 0), kind=real64)
    expected_curl = -real(n_psi*(n_psi + 1), real64)*psi
    expected_divergence = -real(n_chi*(n_chi + 1), real64)*chi

    call transform%wind_to_grid(psi, chi, u, v)
    call transform%gradient_to_grid(psi, dpsi_dlambda, dpsi_dphi)
    call transform%gradient_to_grid(chi, dchi_dlambda, dchi_dphi)
    wind_err = 0.0_real64
    wind_scale = 0.0_real64
    do j = 1, size(nlon)
      nl = nlon(j)
      cosphi = sqrt(1.0_real64 - transform%mu(j)**2)
      wind_err = max(wind_err, maxval(abs(u(1:nl, j) - &
        (dchi_dlambda(1:nl, j)/cosphi - dpsi_dphi(1:nl, j)))))
      wind_err = max(wind_err, maxval(abs(v(1:nl, j) - &
        (dpsi_dlambda(1:nl, j)/cosphi + dchi_dphi(1:nl, j)))))
      wind_scale = max(wind_scale, maxval(abs(u(1:nl, j))), maxval(abs(v(1:nl, j))))
    end do
    wind_err = wind_err/wind_scale

    call transform%curl_divergence(u, v, curl, divergence)
    curl_err = maxval(abs(curl - expected_curl))/maxval(abs(expected_curl))
    divergence_err = maxval(abs(divergence - expected_divergence))/maxval(abs(expected_divergence))

    write (*, '(4(a,i0),a,3(a,es12.4))') 'wind/curl/divergence psi(n=', n_psi, ', m=', m_psi, &
      ') chi(n=', n_chi, ', m=', m_chi, '):', ' relative wind error = ', wind_err, &
      ' curl error = ', curl_err, ' divergence error = ', divergence_err
    if (wind_err > tolerance .or. curl_err > tolerance .or. divergence_err > tolerance) then
      error stop 'wind_to_grid / curl_divergence test failed'
    end if
  end subroutine check_wind_curl_divergence

  !> Solid-body rotation u = cos(latitude), v = 0 on the unit sphere has
  !> curl = 2 sin(latitude) and zero divergence.
  subroutine check_solid_body_rotation()
    real(real64), allocatable :: u(:, :), v(:, :), sin_latitude(:, :)
    complex(real64), allocatable :: curl(:, :), divergence(:, :), expected_curl(:, :)
    real(real64) :: curl_err, divergence_err
    integer :: j

    call transform%allocate_field(u)
    call transform%allocate_field(v)
    call transform%allocate_field(sin_latitude)
    u = 0.0_real64
    v = 0.0_real64
    sin_latitude = 0.0_real64
    do j = 1, size(nlon)
      u(1:nlon(j), j) = sqrt(1.0_real64 - transform%mu(j)**2)
      sin_latitude(1:nlon(j), j) = transform%mu(j)
    end do
    call transform%grid_to_spectral(sin_latitude, expected_curl)
    expected_curl = 2.0_real64*expected_curl
    call transform%curl_divergence(u, v, curl, divergence)
    curl_err = maxval(abs(curl - expected_curl))
    divergence_err = maxval(abs(divergence))

    write (*, '(a,2(a,es12.4))') 'solid-body rotation:', ' curl error = ', curl_err, &
      ' divergence error = ', divergence_err
    if (curl_err > tolerance .or. divergence_err > tolerance) then
      error stop 'solid-body rotation curl/divergence test failed'
    end if
  end subroutine check_solid_body_rotation

  subroutine check_constant()
    real(real64), allocatable :: field(:, :), field2(:, :)
    complex(real64), allocatable :: a(:, :)
    real(real64) :: max_err, rel_rms_err

    call transform%allocate_field(field)
    field = 1.0_real64

    call transform%grid_to_spectral(field, a)
    call transform%spectral_to_grid(a, field2)
    call calc_field_error(field, field2, max_err, rel_rms_err)

    write (*, '(a,2(a,es12.4))') 'constant:', ' max error = ', max_err, &
      ' relative RMS error = ', rel_rms_err
    if (max_err > tolerance .or. rel_rms_err > tolerance) then
      error stop 'constant round-trip failed'
    end if
  end subroutine check_constant

  subroutine check_mode(n_mode, m_mode, amplitude)
    integer, intent(in) :: n_mode, m_mode
    complex(real64), intent(in) :: amplitude
    real(real64), allocatable :: field(:, :), field2(:, :)
    complex(real64), allocatable :: a(:, :), a2(:, :)
    real(real64) :: coefficient_err, field_err, rel_field_err

    allocate (a(0:T + 1, 0:T))
    a = cmplx(0.0_real64, 0.0_real64, kind=real64)
    a(n_mode, m_mode) = amplitude

    call transform%spectral_to_grid(a, field)
    call transform%grid_to_spectral(field, a2)
    call transform%spectral_to_grid(a2, field2)

    coefficient_err = maxval(abs(a2 - a))
    call calc_field_error(field, field2, field_err, rel_field_err)

    write (*, '(a,i0,a,i0,a,3(a,es12.4))') 'mode (n=', n_mode, ', m=', m_mode, '):', &
      ' coefficient error = ', coefficient_err, ' max field error = ', field_err, &
      ' relative RMS error = ', rel_field_err
    if (coefficient_err > tolerance .or. field_err > tolerance .or. rel_field_err > tolerance) then
      error stop 'spherical harmonic mode round-trip failed'
    end if
  end subroutine check_mode

  subroutine check_spectral_extent()
    real(real64), allocatable :: field(:, :)
    complex(real64), allocatable :: a(:, :)

    call transform%allocate_field(field)
    field = 1.0_real64
    call transform%grid_to_spectral(field, a)

    if (lbound(a, 1) /= 0 .or. ubound(a, 1) /= T + 1 .or. &
        lbound(a, 2) /= 0 .or. ubound(a, 2) /= T) then
      error stop 'spectral coefficient array has an inconsistent shape'
    end if
    if (maxval(abs(a(T + 1, :))) > 0.0_real64) then
      error stop 'spectral extension row must be initialized to zero'
    end if
  end subroutine check_spectral_extent

  subroutine check_extension_mode()
    real(real64), allocatable :: field(:, :)
    complex(real64), allocatable :: a(:, :)
    integer :: j
    real(real64) :: max_value

    allocate (a(0:T + 1, 0:T))
    a = cmplx(0.0_real64, 0.0_real64, kind=real64)
    a(T + 1, 0) = cmplx(1.0_real64, 0.0_real64, kind=real64)

    call transform%spectral_to_grid(a, field)
    max_value = 0.0_real64
    do j = 1, size(nlon)
      max_value = max(max_value, maxval(abs(field(1:nlon(j), j))))
    end do
    if (max_value <= tolerance) then
      error stop 'spectral extension mode was not transformed to the grid'
    end if
  end subroutine check_extension_mode

  subroutine check_analytic_gradient(case_id, case_name)
    integer, intent(in) :: case_id
    character(*), intent(in) :: case_name
    real(real64), allocatable :: field(:, :), dfdlambda(:, :), dfdphi(:, :)
    real(real64), allocatable :: expected_lambda(:, :), expected_phi(:, :)
    complex(real64), allocatable :: a(:, :)
    real(real64) :: pi, lambda, mu_j, cosphi
    real(real64) :: lambda_err, phi_err
    integer :: j, k

    call transform%allocate_field(field)
    call transform%allocate_field(expected_lambda)
    call transform%allocate_field(expected_phi)
    field = 0.0_real64
    expected_lambda = 0.0_real64
    expected_phi = 0.0_real64
    pi = acos(-1.0_real64)

    do j = 1, size(nlon)
      mu_j = transform%mu(j)
      cosphi = sqrt(max(0.0_real64, 1.0_real64 - mu_j**2))
      do k = 1, nlon(j)
        lambda = 2.0_real64*pi*real(k - 1, real64)/real(nlon(j), real64)
        select case (case_id)
        case (0)
          field(k, j) = 1.0_real64
        case (1)
          field(k, j) = mu_j
          expected_phi(k, j) = cosphi
        case (2)
          field(k, j) = cosphi*cos(lambda)
          expected_lambda(k, j) = -cosphi*sin(lambda)
          expected_phi(k, j) = -mu_j*cos(lambda)
        case (3)
          field(k, j) = cosphi*sin(lambda)
          expected_lambda(k, j) = cosphi*cos(lambda)
          expected_phi(k, j) = -mu_j*sin(lambda)
        case (4)
          field(k, j) = (1.0_real64 - mu_j**2)*cos(2.0_real64*lambda)
          expected_lambda(k, j) = -2.0_real64*(1.0_real64 - mu_j**2)*sin(2.0_real64*lambda)
          expected_phi(k, j) = -2.0_real64*mu_j*cosphi*cos(2.0_real64*lambda)
        case default
          error stop 'unknown analytic-gradient test case'
        end select
      end do
    end do

    call transform%grid_to_spectral(field, a)
    call transform%gradient_to_grid(a, dfdlambda, dfdphi)

    lambda_err = 0.0_real64
    phi_err = 0.0_real64
    do j = 1, size(nlon)
      lambda_err = max(lambda_err, &
                       maxval(abs(dfdlambda(1:nlon(j), j) - expected_lambda(1:nlon(j), j))))
      phi_err = max(phi_err, &
                    maxval(abs(dfdphi(1:nlon(j), j) - expected_phi(1:nlon(j), j))))
    end do

    write (*, '(3a,2(a,es12.4))') 'analytic gradient ', case_name, ':', &
      ' longitude max error = ', lambda_err, ' latitude max error = ', phi_err
    if (lambda_err > tolerance .or. phi_err > tolerance) then
      error stop 'analytic gradient test failed'
    end if
  end subroutine check_analytic_gradient

  subroutine calc_field_error(field, field2, max_err, rel_rms_err)
    real(real64), intent(in) :: field(:, :), field2(:, :)
    real(real64), intent(out) :: max_err, rel_rms_err
    integer :: j, k, npts
    real(real64) :: err, sum_err2, sum_ref2

    max_err = 0.0_real64
    sum_err2 = 0.0_real64
    sum_ref2 = 0.0_real64
    npts = 0

    do j = 1, 2*(T + 1)
      do k = 1, nlon(j)
        err = field2(k, j) - field(k, j)
        max_err = max(max_err, abs(err))
        sum_err2 = sum_err2 + err**2
        sum_ref2 = sum_ref2 + field(k, j)**2
        npts = npts + 1
      end do
    end do

    if (sum_ref2 > 0.0_real64) then
      rel_rms_err = sqrt(sum_err2/sum_ref2)
    else
      rel_rms_err = sqrt(sum_err2/real(npts, real64))
    end if
  end subroutine calc_field_error
end program check_harmonics
