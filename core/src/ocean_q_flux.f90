!> Prescribed ocean heat convergence (Q flux, docs/tendency/q-flux.md).
!>
!> The zonal profile Q_0 = -q_* (1 - 3 mu^2) is the convergence of the northward
!> transport (3 sqrt(3)/2) Phi_max sin(phi) cos(phi)^2.  On a grid with land the
!> ocean-area mean is removed so that the discrete ocean integral vanishes.  The
!> field is per unit ocean area and zero on pure land and in the grid padding.
module ocean_q_flux
  use iso_fortran_env, only: real64
  use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
  use dry_physics_config, only: q_flux_config
  implicit none
  private
  public :: validate_q_flux_config, q_flux_scale, q_flux_profile, build_ocean_q_flux

  real(real64), parameter :: pi = acos(-1.0_real64)

  !> Values recorded in the case metadata.
  type, public :: q_flux_diagnostics
    real(real64) :: scale = 0.0_real64
    real(real64) :: removed_ocean_mean = 0.0_real64
    !> sum_ij w_ij f_O,ij Q_ij in W m^-2 of planet area (zero up to rounding).
    real(real64) :: ocean_integral = 0.0_real64
    real(real64) :: ocean_area_fraction = 0.0_real64
    !> Extremes of the implied northward transport at the Gauss cell boundaries, W.
    real(real64) :: maximum_transport = 0.0_real64
    real(real64) :: maximum_transport_latitude = 0.0_real64
    real(real64) :: minimum_transport = 0.0_real64
    real(real64) :: minimum_transport_latitude = 0.0_real64
    !> Implied transport across the north pole, W; zero up to rounding.
    real(real64) :: polar_residual_transport = 0.0_real64
  end type q_flux_diagnostics

contains

  subroutine validate_q_flux_config(config)
    type(q_flux_config), intent(in) :: config
    if (.not. ieee_is_finite(config%maximum_transport)) error stop 'nonfinite Q-flux transport'
    if (config%maximum_transport < 0.0_real64) error stop 'Q-flux transport must be non-negative'
  end subroutine validate_q_flux_config

  !> q_* = 3 sqrt(3) Phi_max / (4 pi a^2), W m^-2.
  pure real(real64) function q_flux_scale(config, radius) result(scale)
    type(q_flux_config), intent(in) :: config
    real(real64), intent(in) :: radius
    scale = 3.0_real64*sqrt(3.0_real64)*config%maximum_transport/(4.0_real64*pi*radius**2)
  end function q_flux_scale

  !> Q_0(mu) = -q_* (1 - 3 mu^2) with mu = sin(latitude).
  pure real(real64) function q_flux_profile(config, radius, mu) result(q)
    type(q_flux_config), intent(in) :: config
    real(real64), intent(in) :: radius, mu
    q = -q_flux_scale(config, radius)*(1.0_real64 - 3.0_real64*mu**2)
  end function q_flux_profile

  !> Builds Q(lambda, phi) on the octahedral grid.  gaussian_weights sum to 2 and
  !> the area weight of a grid point is gaussian_weight/(2 nlon).
  subroutine build_ocean_q_flux(config, radius, mu, gaussian_weights, ring_nlon, land_fraction, q, diagnostics)
    type(q_flux_config), intent(in) :: config
    real(real64), intent(in) :: radius, mu(:), gaussian_weights(:), land_fraction(:, :)
    integer, intent(in) :: ring_nlon(:)
    real(real64), intent(out) :: q(:, :)
    type(q_flux_diagnostics), intent(out) :: diagnostics
    real(real64) :: weight, ocean, ring_heat(size(mu)), boundary_mu, transport, integral
    integer :: i, j, jj, ny
    integer :: order(size(mu))

    call validate_q_flux_config(config)
    ny = size(mu)
    if (size(gaussian_weights) /= ny .or. size(ring_nlon) /= ny .or. size(land_fraction, 2) /= ny .or. &
        any(shape(q) /= shape(land_fraction))) error stop 'Q-flux grid has an inconsistent shape'
    if (any(land_fraction < 0.0_real64) .or. any(land_fraction > 1.0_real64)) &
      error stop 'Q-flux land fraction is outside [0,1]'
    diagnostics = q_flux_diagnostics()
    diagnostics%scale = q_flux_scale(config, radius)
    q = 0.0_real64
    ocean = 0.0_real64
    integral = 0.0_real64
    do j = 1, ny
      do i = 1, ring_nlon(j)
        weight = 0.5_real64*gaussian_weights(j)/real(ring_nlon(j), real64)*(1.0_real64 - land_fraction(i, j))
        ocean = ocean + weight
        integral = integral + weight*q_flux_profile(config, radius, mu(j))
      end do
    end do
    diagnostics%ocean_area_fraction = ocean
    if (ocean <= 0.0_real64) error stop 'Q flux requires ocean on the grid'
    diagnostics%removed_ocean_mean = integral/ocean

    ring_heat = 0.0_real64
    integral = 0.0_real64
    do j = 1, ny
      do i = 1, ring_nlon(j)
        if (land_fraction(i, j) >= 1.0_real64) cycle
        q(i, j) = q_flux_profile(config, radius, mu(j)) - diagnostics%removed_ocean_mean
        weight = 0.5_real64*gaussian_weights(j)/real(ring_nlon(j), real64)*(1.0_real64 - land_fraction(i, j))
        ring_heat(j) = ring_heat(j) + weight*q(i, j)
      end do
      integral = integral + ring_heat(j)
    end do
    diagnostics%ocean_integral = integral

    ! Implied northward transport, accumulated from the south pole in increasing mu.
    order = sort_order(mu)
    transport = 0.0_real64
    boundary_mu = -1.0_real64
    diagnostics%maximum_transport = -huge(1.0_real64)
    diagnostics%minimum_transport = huge(1.0_real64)
    do jj = 1, ny - 1
      j = order(jj)
      transport = transport - 4.0_real64*pi*radius**2*ring_heat(j)
      boundary_mu = boundary_mu + gaussian_weights(j)
      if (transport > diagnostics%maximum_transport) then
        diagnostics%maximum_transport = transport
        diagnostics%maximum_transport_latitude = asin(max(-1.0_real64, min(1.0_real64, boundary_mu)))*180.0_real64/pi
      end if
      if (transport < diagnostics%minimum_transport) then
        diagnostics%minimum_transport = transport
        diagnostics%minimum_transport_latitude = asin(max(-1.0_real64, min(1.0_real64, boundary_mu)))*180.0_real64/pi
      end if
    end do
    diagnostics%polar_residual_transport = transport - 4.0_real64*pi*radius**2*ring_heat(order(ny))
  end subroutine build_ocean_q_flux

  pure function sort_order(values) result(order)
    real(real64), intent(in) :: values(:)
    integer :: order(size(values))
    integer :: i, j, key
    order = [(i, i=1, size(values))]
    do i = 2, size(values)
      key = order(i)
      j = i - 1
      do while (j >= 1)
        if (values(order(j)) <= values(key)) exit
        order(j + 1) = order(j)
        j = j - 1
      end do
      order(j + 1) = key
    end do
  end function sort_order
end module ocean_q_flux
