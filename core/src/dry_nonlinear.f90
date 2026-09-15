module dry_nonlinear
  use iso_fortran_env, only: real64
  use harmonics, only: harmonic_transform
  use barotropic_vorticity, only: earth_radius, rotation_rate
  use shallow_water_nonlinear, only: diagnose_shallow_water_velocity, flux_divergence, flux_curl
  use dry_vertical_coordinate, only: hybrid_sigma_coordinate, dry_air_gas_constant, dry_air_kappa
  implicit none
  private

  public :: compute_dry_nonlinear_tendency

contains

  subroutine compute_dry_nonlinear_tendency(transform, truncation, coordinate, &
                                            zeta, delta, temperature, log_surface_pressure, &
                                            surface_geopotential, &
                                            rhs_zeta, rhs_delta, rhs_temperature, &
                                            rhs_log_surface_pressure)
    type(harmonic_transform), intent(inout) :: transform
    integer, intent(in) :: truncation
    type(hybrid_sigma_coordinate), intent(in) :: coordinate
    complex(real64), intent(in) :: zeta(0:, 0:, :), delta(0:, 0:, :), temperature(0:, 0:, :)
    complex(real64), intent(in) :: log_surface_pressure(0:, 0:), surface_geopotential(0:, 0:)
    complex(real64), allocatable, intent(out) :: rhs_zeta(:, :, :), rhs_delta(:, :, :)
    complex(real64), allocatable, intent(out) :: rhs_temperature(:, :, :)
    complex(real64), allocatable, intent(out) :: rhs_log_surface_pressure(:, :)
    real(real64), allocatable :: zeta_grid(:, :, :), delta_grid(:, :, :), temperature_grid(:, :, :)
    real(real64), allocatable :: u(:, :, :), v(:, :, :), log_ps(:, :), ps(:, :)
    real(real64), allocatable :: surface_geopotential_grid(:, :)
    real(real64), allocatable :: dlogps_dlambda(:, :), dlogps_dphi(:, :)
    real(real64), allocatable :: pressure_half(:, :, :), delta_p(:, :, :), layer_l(:, :, :)
    real(real64), allocatable :: alpha(:, :, :), geopotential_half(:, :, :), geopotential(:, :, :)
    real(real64), allocatable :: mass_divergence(:, :, :), cumulative(:, :, :), mass_flux(:, :, :)
    real(real64), allocatable :: pressure_gradient_u(:, :, :), pressure_gradient_v(:, :, :)
    real(real64), allocatable :: vertical_u(:, :, :), vertical_v(:, :, :), vertical_t(:, :, :)
    real(real64), allocatable :: q(:, :), kinetic(:, :), vector_u(:, :), vector_v(:, :)
    real(real64), allocatable :: tendency_grid(:, :), dtdlambda(:, :), dtdphi(:, :)
    real(real64), allocatable :: temporary_grid(:, :), temporary_u(:, :), temporary_v(:, :)
    complex(real64), allocatable :: temporary_spectral(:, :), rotational_curl(:, :)
    integer, allocatable :: nlon(:)
    integer :: number_of_levels, nx, ny, i, j, k, n, m
    real(real64) :: cosphi, gradient_u, gradient_v, coefficient, thermodynamic_q

    number_of_levels = coordinate%number_of_levels
    if (size(zeta, 3) /= number_of_levels .or. size(delta, 3) /= number_of_levels .or. &
        size(temperature, 3) /= number_of_levels) then
      error stop 'dry nonlinear state has the wrong number of vertical levels'
    end if
    call transform%allocate_field(temporary_grid)
    nx = size(temporary_grid, 1)
    ny = size(temporary_grid, 2)
    nlon = transform%get_nlon()
    allocate (zeta_grid(nx, ny, number_of_levels), delta_grid(nx, ny, number_of_levels))
    allocate (temperature_grid(nx, ny, number_of_levels), u(nx, ny, number_of_levels))
    allocate (v(nx, ny, number_of_levels), log_ps(nx, ny), ps(nx, ny))
    allocate (pressure_half(nx, ny, 0:number_of_levels))
    allocate (delta_p(nx, ny, number_of_levels), layer_l(nx, ny, number_of_levels))
    allocate (alpha(nx, ny, number_of_levels), geopotential_half(nx, ny, 0:number_of_levels))
    allocate (geopotential(nx, ny, number_of_levels))
    allocate (mass_divergence(nx, ny, number_of_levels), cumulative(nx, ny, 0:number_of_levels))
    allocate (mass_flux(nx, ny, 0:number_of_levels))
    allocate (pressure_gradient_u(nx, ny, number_of_levels))
    allocate (pressure_gradient_v(nx, ny, number_of_levels))
    allocate (vertical_u(nx, ny, number_of_levels), vertical_v(nx, ny, number_of_levels))
    allocate (vertical_t(nx, ny, number_of_levels))
    allocate (rhs_zeta(0:truncation + 1, 0:truncation, number_of_levels))
    allocate (rhs_delta(0:truncation + 1, 0:truncation, number_of_levels))
    allocate (rhs_temperature(0:truncation + 1, 0:truncation, number_of_levels))
    rhs_zeta = 0.0_real64
    rhs_delta = 0.0_real64
    rhs_temperature = 0.0_real64

    zeta_grid = 0.0_real64
    delta_grid = 0.0_real64
    temperature_grid = 0.0_real64
    u = 0.0_real64
    v = 0.0_real64
    do k = 1, number_of_levels
      call transform%spectral_to_grid(zeta(:, :, k), temporary_grid)
      zeta_grid(:, :, k) = temporary_grid
      call transform%spectral_to_grid(delta(:, :, k), temporary_grid)
      delta_grid(:, :, k) = temporary_grid
      call transform%spectral_to_grid(temperature(:, :, k), temporary_grid)
      temperature_grid(:, :, k) = temporary_grid
      call diagnose_shallow_water_velocity(transform, truncation, zeta(:, :, k), delta(:, :, k), &
                                            temporary_u, temporary_v)
      u(:, :, k) = temporary_u
      v(:, :, k) = temporary_v
    end do
    call transform%spectral_to_grid(log_surface_pressure, log_ps)
    ps = exp(log_ps)
    call transform%gradient_to_grid(log_surface_pressure, dlogps_dlambda, dlogps_dphi)

    pressure_half = 0.0_real64
    delta_p = 0.0_real64
    layer_l = 0.0_real64
    alpha = 0.0_real64
    do k = 0, number_of_levels
      pressure_half(:, :, k) = coordinate%a_half(k) + coordinate%b_half(k)*ps
    end do
    do k = 1, number_of_levels
      delta_p(:, :, k) = pressure_half(:, :, k) - pressure_half(:, :, k - 1)
      layer_l(:, :, k) = log(pressure_half(:, :, k)/pressure_half(:, :, k - 1))
      alpha(:, :, k) = 1.0_real64 - pressure_half(:, :, k - 1)* &
                        layer_l(:, :, k)/delta_p(:, :, k)
    end do

    ! The surface geopotential is a fixed lower boundary condition, not a prognostic field.
    call transform%spectral_to_grid(surface_geopotential, surface_geopotential_grid)
    geopotential_half(:, :, number_of_levels) = surface_geopotential_grid
    do k = number_of_levels, 1, -1
      geopotential(:, :, k) = geopotential_half(:, :, k) + &
        alpha(:, :, k)*dry_air_gas_constant*temperature_grid(:, :, k)
      geopotential_half(:, :, k - 1) = geopotential_half(:, :, k) + &
        dry_air_gas_constant*temperature_grid(:, :, k)*layer_l(:, :, k)
    end do

    mass_divergence = 0.0_real64
    cumulative = 0.0_real64
    do k = 1, number_of_levels
      do j = 1, ny
        cosphi = sqrt(max(0.0_real64, 1.0_real64 - transform%mu(j)**2))
        do i = 1, nlon(j)
          gradient_u = dlogps_dlambda(i, j)/(earth_radius*cosphi)
          gradient_v = dlogps_dphi(i, j)/earth_radius
          mass_divergence(i, j, k) = delta_p(i, j, k)*delta_grid(i, j, k) + &
            ps(i, j)*coordinate%delta_b(k)*(u(i, j, k)*gradient_u + v(i, j, k)*gradient_v)
          cumulative(i, j, k) = cumulative(i, j, k - 1) + mass_divergence(i, j, k)
        end do
      end do
    end do
    do k = 0, number_of_levels
      mass_flux(:, :, k) = coordinate%b_half(k)*cumulative(:, :, number_of_levels) - &
                            cumulative(:, :, k)
    end do

    pressure_gradient_u = 0.0_real64
    pressure_gradient_v = 0.0_real64
    vertical_u = 0.0_real64
    vertical_v = 0.0_real64
    vertical_t = 0.0_real64
    do k = 1, number_of_levels
      do j = 1, ny
        cosphi = sqrt(max(0.0_real64, 1.0_real64 - transform%mu(j)**2))
        do i = 1, nlon(j)
          coefficient = ps(i, j)/delta_p(i, j, k)*(coordinate%b_half(k - 1)* &
            layer_l(i, j, k) + alpha(i, j, k)*coordinate%delta_b(k))
          pressure_gradient_u(i, j, k) = coefficient*dlogps_dlambda(i, j)/(earth_radius*cosphi)
          pressure_gradient_v(i, j, k) = coefficient*dlogps_dphi(i, j)/earth_radius
          if (k < number_of_levels) then
            vertical_u(i, j, k) = vertical_u(i, j, k) + mass_flux(i, j, k)* &
              (u(i, j, k + 1) - u(i, j, k))
            vertical_v(i, j, k) = vertical_v(i, j, k) + mass_flux(i, j, k)* &
              (v(i, j, k + 1) - v(i, j, k))
            vertical_t(i, j, k) = vertical_t(i, j, k) + mass_flux(i, j, k)* &
              (temperature_grid(i, j, k + 1) - temperature_grid(i, j, k))
          end if
          if (k > 1) then
            vertical_u(i, j, k) = vertical_u(i, j, k) + mass_flux(i, j, k - 1)* &
              (u(i, j, k) - u(i, j, k - 1))
            vertical_v(i, j, k) = vertical_v(i, j, k) + mass_flux(i, j, k - 1)* &
              (v(i, j, k) - v(i, j, k - 1))
            vertical_t(i, j, k) = vertical_t(i, j, k) + mass_flux(i, j, k - 1)* &
              (temperature_grid(i, j, k) - temperature_grid(i, j, k - 1))
          end if
          vertical_u(i, j, k) = vertical_u(i, j, k)/(2.0_real64*delta_p(i, j, k))
          vertical_v(i, j, k) = vertical_v(i, j, k)/(2.0_real64*delta_p(i, j, k))
          vertical_t(i, j, k) = vertical_t(i, j, k)/(2.0_real64*delta_p(i, j, k))
        end do
      end do
    end do

    call transform%allocate_field(q)
    call transform%allocate_field(kinetic)
    call transform%allocate_field(vector_u)
    call transform%allocate_field(vector_v)
    call transform%allocate_field(tendency_grid)
    do k = 1, number_of_levels
      q = 0.0_real64
      kinetic = 0.0_real64
      vector_u = 0.0_real64
      vector_v = 0.0_real64
      do j = 1, ny
        do i = 1, nlon(j)
          q(i, j) = zeta_grid(i, j, k) + 2.0_real64*rotation_rate*transform%mu(j)
          kinetic(i, j) = 0.5_real64*(u(i, j, k)**2 + v(i, j, k)**2)
          vector_u(i, j) = q(i, j)*u(i, j, k)
          vector_v(i, j) = q(i, j)*v(i, j, k)
        end do
      end do
      call flux_divergence(transform, vector_u, vector_v, temporary_spectral)
      rhs_zeta(:, :, k) = -temporary_spectral
      call flux_curl(transform, vector_u, vector_v, rotational_curl)
      rhs_delta(:, :, k) = rotational_curl

      tendency_grid = kinetic + geopotential(:, :, k)
      call transform%grid_to_spectral(tendency_grid, temporary_spectral)
      do m = 0, truncation
        do n = m, truncation
          rhs_delta(n, m, k) = rhs_delta(n, m, k) + &
            real(n*(n + 1), real64)*temporary_spectral(n, m)/earth_radius**2
        end do
      end do

      vector_u = -vertical_u(:, :, k) - dry_air_gas_constant*temperature_grid(:, :, k)* &
                                      pressure_gradient_u(:, :, k)
      vector_v = -vertical_v(:, :, k) - dry_air_gas_constant*temperature_grid(:, :, k)* &
                                      pressure_gradient_v(:, :, k)
      call flux_curl(transform, vector_u, vector_v, temporary_spectral)
      rhs_zeta(:, :, k) = rhs_zeta(:, :, k) + temporary_spectral
      call flux_divergence(transform, vector_u, vector_v, temporary_spectral)
      rhs_delta(:, :, k) = rhs_delta(:, :, k) + temporary_spectral

      call transform%gradient_to_grid(temperature(:, :, k), dtdlambda, dtdphi)
      tendency_grid = 0.0_real64
      do j = 1, ny
        cosphi = sqrt(max(0.0_real64, 1.0_real64 - transform%mu(j)**2))
        do i = 1, nlon(j)
          thermodynamic_q = u(i, j, k)*pressure_gradient_u(i, j, k) + &
                            v(i, j, k)*pressure_gradient_v(i, j, k) - &
            (layer_l(i, j, k)*cumulative(i, j, k - 1) + &
             alpha(i, j, k)*mass_divergence(i, j, k))/delta_p(i, j, k)
          tendency_grid(i, j) = -u(i, j, k)*dtdlambda(i, j)/(earth_radius*cosphi) - &
            v(i, j, k)*dtdphi(i, j)/earth_radius - vertical_t(i, j, k) + &
            dry_air_kappa*temperature_grid(i, j, k)*thermodynamic_q
        end do
      end do
      call transform%grid_to_spectral(tendency_grid, temporary_spectral)
      rhs_temperature(:, :, k) = temporary_spectral
      call enforce_tendency(truncation, rhs_zeta(:, :, k), .true.)
      call enforce_tendency(truncation, rhs_delta(:, :, k), .true.)
      call enforce_tendency(truncation, rhs_temperature(:, :, k), .false.)
    end do

    tendency_grid = 0.0_real64
    do j = 1, ny
      do i = 1, nlon(j)
        tendency_grid(i, j) = -cumulative(i, j, number_of_levels)/ps(i, j)
      end do
    end do
    call transform%grid_to_spectral(tendency_grid, rhs_log_surface_pressure)
    call enforce_tendency(truncation, rhs_log_surface_pressure, .false.)
  end subroutine compute_dry_nonlinear_tendency

  subroutine enforce_tendency(truncation, field, zero_mean)
    integer, intent(in) :: truncation
    complex(real64), intent(inout) :: field(0:, 0:)
    logical, intent(in) :: zero_mean
    integer :: n, m

    if (zero_mean) field(0, 0) = 0.0_real64
    field(truncation + 1, :) = 0.0_real64
    do m = 1, truncation
      do n = 0, m - 1
        field(n, m) = 0.0_real64
      end do
    end do
  end subroutine enforce_tendency

end module dry_nonlinear
