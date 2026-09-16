module dry_nonlinear
  use iso_fortran_env, only: real64
  use harmonics, only: harmonic_transform
  use barotropic_vorticity, only: earth_radius, rotation_rate
  use shallow_water_nonlinear, only: diagnose_shallow_water_velocity, flux_curl_divergence
  use dry_vertical_coordinate, only: hybrid_sigma_coordinate, dry_air_gas_constant, dry_air_kappa
  use dry_held_suarez, only: held_suarez_forcing
  use dry_radiation, only: radiation_tendency, radiation_diagnostics, planetary_rotation_rate, &
                           radiation_rayleigh_rate
  use dry_convection, only: dry_convective_adjustment_tendency
  implicit none
  private

  public :: compute_dry_nonlinear_tendency

contains

  subroutine compute_dry_nonlinear_tendency(transform, truncation, coordinate, &
                                            zeta, delta, temperature, log_surface_pressure, &
                                            previous_zeta, previous_delta, previous_temperature, &
                                            previous_log_surface_pressure, previous_surface_temperature, &
                                            previous_deep_temperature, &
                                            surface_geopotential, surface_temperature, deep_temperature, &
                                            rhs_zeta, rhs_delta, rhs_temperature, &
                                            rhs_log_surface_pressure, rhs_surface_temperature, &
                                            rhs_deep_temperature, maximum_speed, apply_held_suarez, &
                                            apply_radiation, evaluation_time, diagnostics)
    type(harmonic_transform), intent(inout) :: transform
    integer, intent(in) :: truncation
    type(hybrid_sigma_coordinate), intent(in) :: coordinate
    complex(real64), intent(in) :: zeta(0:, 0:, :), delta(0:, 0:, :), temperature(0:, 0:, :)
    complex(real64), intent(in) :: log_surface_pressure(0:, 0:), surface_geopotential(0:, 0:)
    !> RAW-filtered previous time level, used by all prescribed physical tendencies.
    complex(real64), intent(in) :: previous_zeta(0:, 0:, :), previous_delta(0:, 0:, :)
    complex(real64), intent(in) :: previous_temperature(0:, 0:, :)
    complex(real64), intent(in) :: previous_log_surface_pressure(0:, 0:)
    complex(real64), intent(in) :: previous_surface_temperature(0:, 0:), previous_deep_temperature(0:, 0:)
    complex(real64), intent(in) :: surface_temperature(0:, 0:), deep_temperature(0:, 0:)
    complex(real64), allocatable, intent(out) :: rhs_zeta(:, :, :), rhs_delta(:, :, :)
    complex(real64), allocatable, intent(out) :: rhs_temperature(:, :, :)
    complex(real64), allocatable, intent(out) :: rhs_log_surface_pressure(:, :)
    complex(real64), allocatable, intent(out) :: rhs_surface_temperature(:, :), rhs_deep_temperature(:, :)
    !> Largest horizontal wind speed (m/s) of the input state, a by-product of the grid winds.
    real(real64), intent(out), optional :: maximum_speed
    logical, intent(in), optional :: apply_held_suarez
    logical, intent(in), optional :: apply_radiation
    real(real64), intent(in), optional :: evaluation_time
    type(radiation_diagnostics), intent(out), optional :: diagnostics
    real(real64), allocatable :: zeta_grid(:, :, :), delta_grid(:, :, :), temperature_grid(:, :, :)
    real(real64), allocatable :: u(:, :, :), v(:, :, :), log_ps(:, :), ps(:, :)
    real(real64), allocatable :: surface_geopotential_grid(:, :)
    real(real64), allocatable :: dlogps_dlambda(:, :), dlogps_dphi(:, :)
    real(real64), allocatable :: pressure_half(:, :, :), delta_p(:, :, :), layer_l(:, :, :)
    real(real64), allocatable :: alpha(:, :, :), geopotential_half(:, :, :), geopotential(:, :, :)
    real(real64), allocatable :: mass_divergence(:, :, :), cumulative(:, :, :), mass_flux(:, :, :)
    real(real64), allocatable :: pressure_gradient_u(:, :, :), pressure_gradient_v(:, :, :)
    real(real64), allocatable :: vertical_u(:, :, :), vertical_v(:, :, :), vertical_t(:, :, :)
    real(real64), allocatable :: vector_u(:, :), vector_v(:, :)
    real(real64), allocatable :: held_temperature_tendency(:, :)
    real(real64), allocatable :: radiative_temperature_tendency(:, :, :)
    real(real64), allocatable :: convective_temperature_tendency(:, :, :)
    real(real64), allocatable :: previous_temperature_grid(:, :, :), previous_u(:, :, :), previous_v(:, :, :)
    real(real64), allocatable :: previous_log_ps(:, :), previous_ps(:, :), previous_pressure_half(:, :, :)
    real(real64), allocatable :: previous_alpha(:, :, :)
    real(real64), allocatable :: previous_surface_temperature_grid(:, :), previous_deep_temperature_grid(:, :)
    real(real64), allocatable :: surface_temperature_grid(:, :), deep_temperature_grid(:, :)
    real(real64), allocatable :: surface_temperature_tendency(:, :), deep_temperature_tendency(:, :)
    real(real64), allocatable :: tendency_grid(:, :), dtdlambda(:, :), dtdphi(:, :)
    real(real64), allocatable :: temporary_grid(:, :), temporary_u(:, :), temporary_v(:, :)
    complex(real64), allocatable :: temporary_spectral(:, :), curl_spectral(:, :), divergence_spectral(:, :)
    integer, allocatable :: nlon(:)
    real(real64), allocatable :: gaussian_weights(:)
    integer :: number_of_levels, nx, ny, i, j, k, n, m
    real(real64) :: cosphi, gradient_u, gradient_v, coefficient, thermodynamic_q, active_rotation_rate
    real(real64) :: absolute_vorticity, kinetic, maximum_speed_squared
    real(real64) :: full_level_pressure, forcing_u, forcing_v, forcing_temperature, top_rayleigh_rate
    real(real64) :: longitude, current_time, incoming_shortwave, reflected_shortwave, outgoing_longwave
    real(real64) :: area_weight, atmospheric_mass, temperature_mass_sum, kinetic_energy_mass_sum
    logical :: use_held_suarez, use_radiation

    number_of_levels = coordinate%number_of_levels
    use_held_suarez = .false.
    if (present(apply_held_suarez)) use_held_suarez = apply_held_suarez
    use_radiation = .false.
    if (present(apply_radiation)) use_radiation = apply_radiation
    active_rotation_rate = rotation_rate
    if (use_radiation) active_rotation_rate = planetary_rotation_rate
    current_time = 0.0_real64
    if (present(evaluation_time)) current_time = evaluation_time
    if (present(diagnostics) .and. .not. use_radiation) then
      error stop 'radiation diagnostics requested without radiation forcing'
    end if
    if (size(zeta, 3) /= number_of_levels .or. size(delta, 3) /= number_of_levels .or. &
        size(temperature, 3) /= number_of_levels .or. &
        size(previous_zeta, 3) /= number_of_levels .or. &
        size(previous_delta, 3) /= number_of_levels .or. &
        size(previous_temperature, 3) /= number_of_levels) then
      error stop 'dry nonlinear state has the wrong number of vertical levels'
    end if
    call transform%allocate_field(temporary_grid)
    nx = size(temporary_grid, 1)
    ny = size(temporary_grid, 2)
    nlon = transform%get_nlon()
    if (use_radiation) gaussian_weights = transform%get_gaussian_weights()
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
    if (use_held_suarez .or. use_radiation) then
      allocate (previous_temperature_grid(nx, ny, number_of_levels))
      allocate (previous_u(nx, ny, number_of_levels), previous_v(nx, ny, number_of_levels))
      allocate (previous_log_ps(nx, ny), previous_ps(nx, ny))
      allocate (previous_pressure_half(nx, ny, 0:number_of_levels))
      allocate (previous_alpha(nx, ny, number_of_levels))
    end if
    if (use_radiation) then
      allocate (previous_surface_temperature_grid(nx, ny), previous_deep_temperature_grid(nx, ny))
    end if
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
      if (use_held_suarez .or. use_radiation) then
        call transform%spectral_to_grid(previous_temperature(:, :, k), temporary_grid)
        previous_temperature_grid(:, :, k) = temporary_grid
        call diagnose_shallow_water_velocity(transform, truncation, previous_zeta(:, :, k), &
                                              previous_delta(:, :, k), temporary_u, temporary_v)
        previous_u(:, :, k) = temporary_u
        previous_v(:, :, k) = temporary_v
      end if
    end do
    call transform%spectral_to_grid(log_surface_pressure, log_ps)
    ps = exp(log_ps)
    if (use_held_suarez .or. use_radiation) then
      call transform%spectral_to_grid(previous_log_surface_pressure, previous_log_ps)
      previous_ps = exp(previous_log_ps)
      do k = 0, number_of_levels
        previous_pressure_half(:, :, k) = coordinate%a_half(k) + coordinate%b_half(k)*previous_ps
      end do
      do k = 1, number_of_levels
        previous_alpha(:, :, k) = 1.0_real64 - previous_pressure_half(:, :, k - 1)* &
          log(previous_pressure_half(:, :, k)/previous_pressure_half(:, :, k - 1))/ &
          (previous_pressure_half(:, :, k) - previous_pressure_half(:, :, k - 1))
      end do
    end if
    if (use_radiation) then
      call transform%spectral_to_grid(surface_temperature, surface_temperature_grid)
      call transform%spectral_to_grid(deep_temperature, deep_temperature_grid)
      call transform%spectral_to_grid(previous_surface_temperature, previous_surface_temperature_grid)
      call transform%spectral_to_grid(previous_deep_temperature, previous_deep_temperature_grid)
    end if
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

    call transform%allocate_field(vector_u)
    call transform%allocate_field(vector_v)
    call transform%allocate_field(tendency_grid)
    call transform%allocate_field(held_temperature_tendency)
    allocate (radiative_temperature_tendency(nx, ny, number_of_levels))
    allocate (convective_temperature_tendency(nx, ny, number_of_levels))
    radiative_temperature_tendency = 0.0_real64
    convective_temperature_tendency = 0.0_real64
    if (use_radiation) then
      do j = 1, ny
        do i = 1, nlon(j)
          call dry_convective_adjustment_tendency(previous_pressure_half(i, j, :), &
            previous_temperature_grid(i, j, :), convective_temperature_tendency(i, j, :))
        end do
      end do
      call transform%allocate_field(surface_temperature_tendency)
      call transform%allocate_field(deep_temperature_tendency)
      surface_temperature_tendency = 0.0_real64
      deep_temperature_tendency = 0.0_real64
      if (present(diagnostics)) then
        diagnostics%time_seconds = current_time
        diagnostics%mean_atmospheric_temperature = 0.0_real64
        diagnostics%mean_surface_temperature = 0.0_real64
        diagnostics%mean_deep_temperature = 0.0_real64
        diagnostics%mean_kinetic_energy = 0.0_real64
        diagnostics%mean_surface_pressure = 0.0_real64
        diagnostics%mean_incoming_shortwave = 0.0_real64
        diagnostics%mean_reflected_shortwave = 0.0_real64
        diagnostics%mean_outgoing_longwave = 0.0_real64
        diagnostics%surface_temperature = surface_temperature_grid
        diagnostics%surface_pressure = ps
        allocate (diagnostics%zonal_temperature(ny, number_of_levels))
        allocate (diagnostics%zonal_u(ny, number_of_levels), diagnostics%zonal_v(ny, number_of_levels))
        allocate (diagnostics%zonal_uv(ny, number_of_levels), diagnostics%zonal_vt(ny, number_of_levels))
        diagnostics%zonal_temperature = 0.0_real64
        diagnostics%zonal_u = 0.0_real64
        diagnostics%zonal_v = 0.0_real64
        diagnostics%zonal_uv = 0.0_real64
        diagnostics%zonal_vt = 0.0_real64
        atmospheric_mass = 0.0_real64
        temperature_mass_sum = 0.0_real64
        kinetic_energy_mass_sum = 0.0_real64
      end if
      do j = 1, ny
        do i = 1, nlon(j)
          longitude = 2.0_real64*acos(-1.0_real64)*real(i - 1, real64)/real(nlon(j), real64)
          ! Radiation and surface exchange use one consistent RAW-filtered previous-time column.
          call radiation_tendency(previous_pressure_half(i, j, :), previous_temperature_grid(i, j, :), &
            previous_surface_temperature_grid(i, j), previous_deep_temperature_grid(i, j), &
            previous_u(i, j, number_of_levels), previous_v(i, j, number_of_levels), transform%mu(j), &
            longitude, current_time, &
            radiative_temperature_tendency(i, j, :), surface_temperature_tendency(i, j), &
            deep_temperature_tendency(i, j), incoming_shortwave, reflected_shortwave, outgoing_longwave)
          if (present(diagnostics)) then
            area_weight = 0.5_real64*gaussian_weights(j)/real(nlon(j), real64)
            diagnostics%mean_surface_temperature = diagnostics%mean_surface_temperature + &
              area_weight*surface_temperature_grid(i, j)
            diagnostics%mean_deep_temperature = diagnostics%mean_deep_temperature + &
              area_weight*deep_temperature_grid(i, j)
            diagnostics%mean_surface_pressure = diagnostics%mean_surface_pressure + area_weight*ps(i, j)
            diagnostics%mean_incoming_shortwave = diagnostics%mean_incoming_shortwave + &
              area_weight*incoming_shortwave
            diagnostics%mean_reflected_shortwave = diagnostics%mean_reflected_shortwave + &
              area_weight*reflected_shortwave
            diagnostics%mean_outgoing_longwave = diagnostics%mean_outgoing_longwave + &
              area_weight*outgoing_longwave
            do k = 1, number_of_levels
              diagnostics%zonal_temperature(j, k) = diagnostics%zonal_temperature(j, k) + &
                temperature_grid(i, j, k)/real(nlon(j), real64)
              diagnostics%zonal_u(j, k) = diagnostics%zonal_u(j, k) + u(i, j, k)/real(nlon(j), real64)
              diagnostics%zonal_v(j, k) = diagnostics%zonal_v(j, k) + v(i, j, k)/real(nlon(j), real64)
              diagnostics%zonal_uv(j, k) = diagnostics%zonal_uv(j, k) + &
                u(i, j, k)*v(i, j, k)/real(nlon(j), real64)
              diagnostics%zonal_vt(j, k) = diagnostics%zonal_vt(j, k) + &
                v(i, j, k)*temperature_grid(i, j, k)/real(nlon(j), real64)
              atmospheric_mass = atmospheric_mass + area_weight*delta_p(i, j, k)
              temperature_mass_sum = temperature_mass_sum + &
                area_weight*delta_p(i, j, k)*temperature_grid(i, j, k)
              kinetic_energy_mass_sum = kinetic_energy_mass_sum + area_weight*delta_p(i, j, k)* &
                0.5_real64*(u(i, j, k)**2 + v(i, j, k)**2)
            end do
          end if
        end do
      end do
      if (present(diagnostics)) then
        diagnostics%mean_atmospheric_temperature = temperature_mass_sum/atmospheric_mass
        diagnostics%mean_kinetic_energy = kinetic_energy_mass_sum/atmospheric_mass
      end if
    end if
    maximum_speed_squared = 0.0_real64
    do k = 1, number_of_levels
      vector_u = 0.0_real64
      vector_v = 0.0_real64
      tendency_grid = 0.0_real64
      held_temperature_tendency = 0.0_real64
      top_rayleigh_rate = 0.0_real64
      if (use_radiation) top_rayleigh_rate = radiation_rayleigh_rate(k)
      ! With the momentum forcing F = -(vertical advection) - R T (pressure-gradient term),
      !   d(zeta)/dt = -div((zeta+f) u) + curl F,   d(delta)/dt = curl((zeta+f) u) + div F - lap(K+Phi).
      ! Since -div(A u, A v) = curl(A v, -A u) and curl(A u, A v) = div(A v, -A u), both
      ! tendencies are the curl and divergence of the single vector G = ((zeta+f) v + F_u,
      ! -(zeta+f) u + F_v), which costs two grid-to-spectral transforms.
      do j = 1, ny
        do i = 1, nlon(j)
          absolute_vorticity = zeta_grid(i, j, k) + 2.0_real64*active_rotation_rate*transform%mu(j)
          kinetic = 0.5_real64*(u(i, j, k)**2 + v(i, j, k)**2)
          maximum_speed_squared = max(maximum_speed_squared, 2.0_real64*kinetic)
          vector_u(i, j) = absolute_vorticity*v(i, j, k) - vertical_u(i, j, k) - &
            dry_air_gas_constant*temperature_grid(i, j, k)*pressure_gradient_u(i, j, k)
          vector_v(i, j) = -absolute_vorticity*u(i, j, k) - vertical_v(i, j, k) - &
            dry_air_gas_constant*temperature_grid(i, j, k)*pressure_gradient_v(i, j, k)
          if (use_held_suarez .or. use_radiation) then
            full_level_pressure = previous_pressure_half(i, j, k)*exp(-previous_alpha(i, j, k))
            call held_suarez_forcing(transform%mu(j), full_level_pressure, previous_ps(i, j), &
              previous_temperature_grid(i, j, k), previous_u(i, j, k), previous_v(i, j, k), &
              forcing_u, forcing_v, forcing_temperature)
            vector_u(i, j) = vector_u(i, j) + forcing_u
            vector_v(i, j) = vector_v(i, j) + forcing_v
            if (use_held_suarez) held_temperature_tendency(i, j) = forcing_temperature
          end if
          if (use_radiation) then
            vector_u(i, j) = vector_u(i, j) - top_rayleigh_rate*previous_u(i, j, k)
            vector_v(i, j) = vector_v(i, j) - top_rayleigh_rate*previous_v(i, j, k)
          end if
          tendency_grid(i, j) = kinetic + geopotential(i, j, k)
        end do
      end do
      call flux_curl_divergence(transform, vector_u, vector_v, curl_spectral, divergence_spectral)
      rhs_zeta(:, :, k) = curl_spectral
      rhs_delta(:, :, k) = divergence_spectral

      call transform%grid_to_spectral(tendency_grid, temporary_spectral)
      do m = 0, truncation
        do n = m, truncation
          rhs_delta(n, m, k) = rhs_delta(n, m, k) + &
            real(n*(n + 1), real64)*temporary_spectral(n, m)/earth_radius**2
        end do
      end do

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
            dry_air_kappa*temperature_grid(i, j, k)*thermodynamic_q + &
            held_temperature_tendency(i, j) + radiative_temperature_tendency(i, j, k) + &
            convective_temperature_tendency(i, j, k)
        end do
      end do
      call transform%grid_to_spectral(tendency_grid, temporary_spectral)
      rhs_temperature(:, :, k) = temporary_spectral
      call enforce_tendency(truncation, rhs_zeta(:, :, k), .true.)
      call enforce_tendency(truncation, rhs_delta(:, :, k), .true.)
      call enforce_tendency(truncation, rhs_temperature(:, :, k), .false.)
    end do
    if (present(maximum_speed)) maximum_speed = sqrt(maximum_speed_squared)

    tendency_grid = 0.0_real64
    do j = 1, ny
      do i = 1, nlon(j)
        tendency_grid(i, j) = -cumulative(i, j, number_of_levels)/ps(i, j)
      end do
    end do
    call transform%grid_to_spectral(tendency_grid, rhs_log_surface_pressure)
    call enforce_tendency(truncation, rhs_log_surface_pressure, .false.)
    if (use_radiation) then
      call transform%grid_to_spectral(surface_temperature_tendency, rhs_surface_temperature)
      call transform%grid_to_spectral(deep_temperature_tendency, rhs_deep_temperature)
      call enforce_tendency(truncation, rhs_surface_temperature, .false.)
      call enforce_tendency(truncation, rhs_deep_temperature, .false.)
    else
      allocate (rhs_surface_temperature(0:truncation + 1, 0:truncation))
      allocate (rhs_deep_temperature(0:truncation + 1, 0:truncation))
      rhs_surface_temperature = 0.0_real64
      rhs_deep_temperature = 0.0_real64
    end if
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
