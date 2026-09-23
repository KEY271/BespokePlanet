!> Acceptance checks of the prescribed ocean heat convergence (docs/tendency/q-flux.md).
program check_q_flux
  use iso_fortran_env, only: real64
  use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
  use planet_parameters, only: earth_radius
  use dry_physics_config, only: q_flux_config, sea_ice_config, radiation_config, dry_model_physics_config
  use ocean_q_flux, only: q_flux_diagnostics, build_ocean_q_flux, q_flux_profile, q_flux_scale
  use sea_ice, only: advance_sea_ice, sea_ice_budget
  use surface_tiles, only: tiled_surface_tendency
  use harmonics, only: harmonic_transform
  use dry_atmosphere, only: dry_atmosphere_solver
  use dry_case_initial_conditions, only: land_sea_case_physics, set_land_sea_case_state, radiation_case_planet
  implicit none
  real(real64), parameter :: co = 1.2558e8_real64, dt = 1200.0_real64
  type(q_flux_config) :: config
  type(sea_ice_config) :: ice
  real(real64) :: latent, tf

  config%enabled = .true.
  ice%enabled = .true.
  latent = ice%density*ice%latent_heat
  tf = ice%freezing_temperature
  call near(q_flux_scale(config, earth_radius), 15.280869504352504_real64, 1.0e-12_real64, 'q_* for 1.5 PW')
  call check_fields(31)
  call check_ocean_updates()
  call check_tile_path()
  call check_zero_transport_matches_disabled()
  write (*, '(a)') 'check_q_flux: all checks passed'

contains

  subroutine near(actual, target, tolerance, name)
    real(real64), intent(in) :: actual, target, tolerance
    character(*), intent(in) :: name
    if (.not. ieee_is_finite(actual) .or. abs(actual - target) > tolerance) then
      write (*, '(a,3es24.15)') name, actual, target, tolerance
      error stop 'Q-flux assertion failed'
    end if
  end subroutine near

  subroutine check_fields(truncation)
    integer, intent(in) :: truncation
    type(harmonic_transform) :: transform
    real(real64), allocatable :: land(:, :), q(:, :), weights(:)
    integer, allocatable :: nlon(:)
    type(q_flux_diagnostics) :: stats
    integer :: i, j
    real(real64) :: ocean_mean_q

    call transform%init(truncation)
    nlon = transform%get_nlon()
    weights = transform%get_gaussian_weights()
    call transform%allocate_field(land)
    allocate (q, mold=land)

    ! Aquaplanet: Gauss quadrature integrates the P_2 profile exactly.
    land = 0.0_real64
    call build_ocean_q_flux(config, earth_radius, transform%mu, weights, nlon, land, q, stats)
    call near(stats%removed_ocean_mean, 0.0_real64, 1.0e-12_real64, 'aquaplanet removed mean')
    call near(stats%ocean_integral, 0.0_real64, 1.0e-12_real64, 'aquaplanet ocean integral')
    call near(stats%ocean_area_fraction, 1.0_real64, 1.0e-13_real64, 'aquaplanet area')
    do j = 1, size(nlon)
      do i = 1, nlon(j)
        call near(q(i, j), q_flux_profile(config, earth_radius, transform%mu(j)), 1.0e-12_real64, 'aquaplanet profile')
      end do
    end do
    ! The implied transport approaches Phi_max near 35.3 degrees in both hemispheres.
    call near(stats%maximum_transport/config%maximum_transport, 1.0_real64, 0.02_real64, 'northern peak transport')
    call near(stats%minimum_transport/config%maximum_transport, -1.0_real64, 0.02_real64, 'southern peak transport')
    call near(stats%maximum_transport_latitude, 35.26_real64, 3.0_real64, 'northern peak latitude')
    call near(stats%minimum_transport_latitude, -35.26_real64, 3.0_real64, 'southern peak latitude')
    call near(stats%polar_residual_transport/config%maximum_transport, 0.0_real64, 1.0e-12_real64, 'polar closure')

    ! Land: pure land takes no flux and the ocean integral still vanishes.
    do j = 1, size(nlon)
      do i = 1, nlon(j)
        land(i, j) = 0.5_real64*(1.0_real64 + sin(3.0_real64*real(i, real64) + 2.0_real64*transform%mu(j)))
        if (transform%mu(j) < -0.8_real64) land(i, j) = 1.0_real64
        if (mod(i, 7) == 0) land(i, j) = 1.0_real64
        if (mod(i, 5) == 0) land(i, j) = 0.0_real64
      end do
    end do
    call build_ocean_q_flux(config, earth_radius, transform%mu, weights, nlon, land, q, stats)
    call near(stats%ocean_integral, 0.0_real64, 1.0e-12_real64, 'land ocean integral')
    if (abs(stats%removed_ocean_mean) < 1.0e-3_real64) error stop 'land should shift the ocean mean'
    call near(stats%polar_residual_transport/config%maximum_transport, 0.0_real64, 1.0e-12_real64, 'land polar closure')
    ocean_mean_q = 0.0_real64
    do j = 1, size(nlon)
      do i = 1, nlon(j)
        if (land(i, j) == 1.0_real64) then
          if (q(i, j) /= 0.0_real64) error stop 'Q flux on pure land'
        else
          call near(q(i, j), q_flux_profile(config, earth_radius, transform%mu(j)) - stats%removed_ocean_mean, &
            1.0e-12_real64, 'land-shifted profile')
        end if
      end do
    end do
  end subroutine check_fields

  subroutine check_ocean_updates()
    real(real64) :: t, a, v, q, f
    type(sea_ice_budget) :: budget

    ! Ice-free ocean: Co dT = dt (F + Q).
    f = 40.0_real64
    q = -15.0_real64
    call advance_sea_ice(ice, co, dt, tf + 5.0_real64, 0.0_real64, 0.0_real64, f, 0.0_real64, 0.0_real64, &
      t, a, v, budget, q)
    call near(co*(t - tf - 5.0_real64), dt*(f + q), 1.0e-6_real64, 'ice-free Q heating')
    if (a /= 0.0_real64 .or. v /= 0.0_real64) error stop 'ice-free Q made ice'
    call near(budget%energy_residual, 0.0_real64, 1.0e-6_real64, 'ice-free energy residual')

    ! Full cover, Q > 0: the open water adds nothing and Q melts dt Q / L from below.
    q = 25.0_real64
    call advance_sea_ice(ice, co, dt, tf, 1.0_real64, 1.0_real64, 1.0e6_real64, 0.0_real64, 0.0_real64, &
      t, a, v, budget, q)
    call near(v, 1.0_real64 - dt*q/latent, 1.0e-15_real64, 'full-cover basal Q melting')
    call near(t, tf, 0.0_real64, 'ice-covered water stays at freezing')
    call near(budget%energy_residual, 0.0_real64, 1.0e-6_real64, 'full-cover melting residual')

    ! Full cover, Q < 0: the area stays one and the volume grows by dt |Q| / L.
    q = -25.0_real64
    call advance_sea_ice(ice, co, dt, tf, 1.0_real64, 1.0_real64, -1.0e6_real64, 0.0_real64, 0.0_real64, &
      t, a, v, budget, q)
    call near(a, 1.0_real64, 0.0_real64, 'full-cover freezing area')
    call near(v, 1.0_real64 + dt*abs(q)/latent, 1.0e-15_real64, 'full-cover Q freezing')

    ! Partial cover: the energy change includes Q without an open-water factor.
    q = 30.0_real64
    call advance_sea_ice(ice, co, dt, tf, 0.4_real64, 0.2_real64, -50.0_real64, 10.0_real64, 0.0_real64, &
      t, a, v, budget, q)
    call near(co*(t - tf) - latent*(v - 0.2_real64), dt*(0.6_real64*(-50.0_real64) + 0.4_real64*(-10.0_real64) + q), &
      1.0e-6_real64, 'partial-cover energy change')
    call near(budget%energy_residual, 0.0_real64, 1.0e-6_real64, 'partial-cover residual')
  end subroutine check_ocean_updates

  subroutine check_tile_path()
    type(radiation_config) :: rad
    type(sea_ice_config) :: disabled
    real(real64) :: pressure(0:2), air(2), rhs(2), rhs_q(2), lr, dr, ocean_r, ocean_q, ar, vr, ti
    real(real64) :: incoming, reflected, outgoing
    type(sea_ice_budget) :: budget

    pressure = [1000.0_real64, 50000.0_real64, 100000.0_real64]
    air = [230.0_real64, 260.0_real64]
    call tiled_surface_tendency(rad, disabled, pressure, air, 2.0_real64, 1.0_real64, 0.0_real64, &
      0.0_real64, 0.0_real64, dt, 0.3_real64, 285.0_real64, 280.0_real64, 290.0_real64, 0.0_real64, 0.0_real64, &
      0.0_real64, 0.0_real64, 0.4_real64, rhs, lr, dr, ocean_r, ar, vr, ti, incoming, reflected, outgoing, budget)
    call tiled_surface_tendency(rad, disabled, pressure, air, 2.0_real64, 1.0_real64, 0.0_real64, &
      0.0_real64, 0.0_real64, dt, 0.3_real64, 285.0_real64, 280.0_real64, 290.0_real64, 0.0_real64, 0.0_real64, &
      0.0_real64, 0.0_real64, 0.4_real64, rhs_q, lr, dr, ocean_q, ar, vr, ti, incoming, reflected, outgoing, budget, &
      ocean_heat_convergence=20.0_real64)
    call near(ocean_q - ocean_r, 20.0_real64/co, 1.0e-18_real64, 'ice-disabled Q tendency')
    if (any(rhs_q /= rhs)) error stop 'Q flux changed the atmosphere'
  end subroutine check_tile_path

  !> Phi_max = 0 reproduces the run without Q flux bit for bit.
  subroutine check_zero_transport_matches_disabled()
    type(dry_model_physics_config) :: off, zero
    real(real64), allocatable :: t_off(:, :, :), t_zero(:, :, :), o_off(:, :), o_zero(:, :), a_off(:, :), a_zero(:, :)
    real(real64), allocatable :: q(:, :)
    type(q_flux_diagnostics) :: stats

    off = land_sea_case_physics()
    off%q_flux%enabled = .false.
    zero = land_sea_case_physics()
    zero%q_flux%maximum_transport = 0.0_real64
    call run(off, t_off, o_off, a_off, q, stats)
    if (any(q /= 0.0_real64)) error stop 'disabled Q flux is not zero'
    call run(zero, t_zero, o_zero, a_zero, q, stats)
    if (any(q /= 0.0_real64)) error stop 'zero-transport Q flux is not zero'
    if (any(t_off /= t_zero) .or. any(o_off /= o_zero) .or. any(a_off /= a_zero)) &
      error stop 'zero Q flux differs from disabled Q flux'
    zero%q_flux%maximum_transport = 1.5e15_real64
    call run(zero, t_zero, o_zero, a_zero, q, stats)
    if (all(o_off == o_zero)) error stop 'Q flux had no effect on the ocean'
    call near(stats%ocean_integral, 0.0_real64, 1.0e-12_real64, 'solver ocean integral')
  end subroutine check_zero_transport_matches_disabled

  subroutine run(physics, temperature, ocean, area, q, stats)
    type(dry_model_physics_config), intent(in) :: physics
    real(real64), allocatable, intent(out) :: temperature(:, :, :), ocean(:, :), area(:, :), q(:, :)
    type(q_flux_diagnostics), intent(out) :: stats
    type(harmonic_transform) :: transform
    type(dry_atmosphere_solver) :: solver
    complex(real64), allocatable :: phi(:, :)
    real(real64), allocatable :: land(:, :), z(:, :, :), d(:, :, :), ps(:, :), u(:, :, :), v(:, :, :)
    integer :: step

    call transform%init(7)
    call transform%allocate_field(land)
    land = 0.3_real64
    land(1, :) = 1.0_real64
    land(2, :) = 0.0_real64
    allocate (phi(0:8, 0:7))
    phi = 0.0_real64
    call solver%init(7, dt)
    call set_land_sea_case_state(solver, transform, physics, radiation_case_planet(physics), phi, land)
    call solver%get_ocean_q_flux(q, stats)
    do step = 1, 12
      call solver%advance()
    end do
    call solver%get_fields(z, d, temperature, ps, u, v, ocean_temperature=ocean, sea_ice_fraction=area)
  end subroutine run
end program check_q_flux
