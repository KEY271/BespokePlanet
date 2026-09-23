program check_sea_ice
  use iso_fortran_env, only: real64
  use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
  use dry_physics_config, only: sea_ice_config, radiation_config, dry_model_physics_config
  use sea_ice, only: advance_sea_ice, solve_ice_surface, equilibrate_sea_ice, sea_ice_budget, sea_ice_checks
  use radiation_case_output, only: radiation_output_options, write_radiation_monthly_output, &
    write_radiation_yearly_snapshot, write_radiation_metadata, initialize_radiation_daily_output, &
    append_radiation_daily_output
  use field_binary_writer, only: write_field
  use topography, only: topography_config, topography_diagnostics
  use surface_tiles, only: tiled_surface_tendency
  use harmonics, only: harmonic_transform
  use dry_atmosphere, only: dry_atmosphere_solver
  use dry_case_initial_conditions, only: land_sea_case_physics, set_land_sea_case_state, radiation_case_planet
  use dry_vertical_coordinate, only: hybrid_sigma_coordinate
  use dry_radiation, only: radiation_diagnostics, radiation_tendency, move_radiation_diagnostics
  use radiation_diagnostics_collector, only: radiation_monthly_accumulator, radiation_monthly_means, radiation_daily_accumulator
  use ocean_q_flux, only: q_flux_diagnostics
  implicit none
  type(sea_ice_config) :: ice
  real(real64), parameter :: co = 1.2558e8_real64, dt = 1200.0_real64, sigma = 5.670374419e-8_real64
  real(real64) :: latent, tf, t, a, v, ts, qc, melt, residual, expected, energy, year_a, year_v, half_a, half_v
  type(sea_ice_budget) :: budget
  type(sea_ice_checks) :: checks

  ice%enabled = .true.
  latent = ice%density*ice%latent_heat
  tf = ice%freezing_temperature
  call update(tf + 5.0_real64, 0.0_real64, 0.0_real64, 100.0_real64, 0.0_real64, 0.0_real64)
  call near(t, tf + 5.0_real64 + dt*100.0_real64/co, 1.0e-12_real64, 'warm open water')
  call near(v, 0.0_real64, 0.0_real64, 'no fictitious warm-water ice')
  call update(tf, 0.0_real64, 0.0_real64, -100.0_real64, 0.0_real64, 0.0_real64)
  call near(v, dt*100.0_real64/latent, 1.0e-15_real64, 'new ice volume')
  call near(a, v/ice%new_ice_thickness, 1.0e-15_real64, 'new ice area')
  call update(tf, 0.25_real64, 0.2_real64, -100.0_real64, 0.0_real64, 0.0_real64)
  call near(budget%open_water_volume, 90000.0_real64/latent, 1.0e-15_real64, 'open-water area factor')
  call update(tf, 0.99_real64, 0.5_real64, -1.0e8_real64, 0.0_real64, 0.0_real64)
  call near(a, 1.0_real64, 0.0_real64, 'area cap')
  if (v <= 0.505_real64) error stop 'area cap discarded volume'

  call solve_ice_surface(ice, 1.0_real64, 150.0_real64, sigma, 2.0_real64, 240.0_real64, ts, qc, melt, residual)
  if (ts >= tf .or. qc <= 0.0_real64 .or. melt /= 0.0_real64) error stop 'cold ice surface branch'
  call near(150.0_real64 - sigma*ts**4 - 2.0_real64*(ts - 240.0_real64) + &
    ice%conductivity*(tf - ts), 0.0_real64, 1.0e-5_real64, 'cold ice equilibrium')
  call update(tf, 1.0_real64, 1.0_real64, 1.0e6_real64, qc, melt)
  call near(budget%open_water_volume, 0.0_real64, 0.0_real64, 'full cover has no open-water flux')
  if (v <= 1.0_real64) error stop 'cold ice did not grow'
  expected = 272.0_real64
  call solve_ice_surface(ice, 1.0_real64, sigma*expected**4 + 2.0_real64*(expected - tf), &
    sigma, 0.0_real64, tf, ts, qc, melt)
  call near(ts, expected, 1.0e-7_real64, 'warm nonmelting surface temperature')
  if (qc >= 0.0_real64 .or. melt /= 0.0_real64) error stop 'basal melting sign'
  call solve_ice_surface(ice, 1.0_real64, 500.0_real64, sigma, 0.0_real64, tf, ts, qc, melt)
  call near(ts, ice%melting_temperature, 0.0_real64, 'surface temperature cap')
  if (melt <= 0.0_real64 .or. qc >= 0.0_real64) error stop 'surface melting branch'
  call near(melt - qc, 500.0_real64 - sigma*ts**4, 1.0e-12_real64, 'surface plus basal energy')
  call solve_ice_surface(ice, 1.0e-16_real64, 150.0_real64, sigma, 2.0_real64, 240.0_real64, ts, qc, melt)
  if (.not. all(ieee_is_finite([ts, qc, melt]))) error stop 'thin ice nonfinite'

  call update(tf, 0.4_real64, 0.2_real64, 0.0_real64, 0.0_real64, 0.1_real64*latent/(dt*0.4_real64))
  call near(a, 0.4_real64*sqrt(0.5_real64), 1.0e-14_real64, 'square-root area law')
  call near(v, 0.1_real64, 1.0e-14_real64, 'partial melt volume')
  call update(tf, 0.4_real64, 0.2_real64, 0.0_real64, 0.0_real64, 0.25_real64*latent/(dt*0.4_real64))
  call near(a, 0.0_real64, 0.0_real64, 'complete melt area')
  call near(t, tf + 0.05_real64*latent/co, 1.0e-12_real64, 'excess surface melt heat')
  call update(tf, 0.4_real64, 0.2_real64, 0.25_real64*latent/(dt*0.6_real64), 0.0_real64, 0.0_real64)
  call near(t, tf + 0.05_real64*latent/co, 1.0e-12_real64, 'excess open-water melt heat')
  call update(tf, 0.4_real64, 0.2_real64, 0.0_real64, -0.25_real64*latent/(dt*0.4_real64), 0.0_real64)
  call near(t, tf + 0.05_real64*latent/co, 1.0e-12_real64, 'excess basal melt heat')
  ! Simultaneous new ice (+0.01 m, +0.02 area) and surface melting (-0.11 m).
  call update(tf, 0.4_real64, 0.2_real64, -0.01_real64*latent/(dt*0.6_real64), &
    0.0_real64, 0.11_real64*latent/(dt*0.4_real64))
  call near(a, 0.42_real64*sqrt(0.5_real64), 1.0e-14_real64, 'simultaneous growth and melt')
  call update(tf, 0.4_real64, 0.2_real64, -0.01_real64*latent/(dt*0.6_real64), &
    0.0_real64, 0.01_real64*latent/(dt*0.4_real64))
  call near(a, 0.42_real64, 1.0e-14_real64, 'zero net volume with new ice')
  t = tf - 0.1_real64
  a = -0.01_real64
  v = -0.001_real64
  energy = co*(t - tf) - latent*v
  call equilibrate_sea_ice(ice, co, t, a, v, checks)
  call near(co*(t - tf) - latent*v, energy, 1.0e-7_real64, 'RAW cold projection energy')
  if (a <= 0.0_real64 .or. v <= 0.0_real64 .or. t /= tf) error stop 'RAW cold projection phase'
  t = tf + 1.0_real64
  v = 0.1_real64
  energy = co*(t - tf) - latent*v
  call equilibrate_sea_ice(ice, co, t, a, v, checks)
  call near(co*(t - tf) - latent*v, energy, 1.0e-5_real64, 'RAW warm projection energy')
  if (a /= 0.0_real64 .or. v /= 0.0_real64) error stop 'RAW warm projection phase'
  if (checks%projected_cells /= 2 .or. checks%checked_cells /= 2) error stop 'projection counts'
  call near(checks%projection_energy_residual, 0.0_real64, 1.0e-5_real64, 'projection residual diagnostic')
  call check_tile_energy()
  call check_monthly_weighting()
  call check_coupled_steps(dt, year_a, year_v)
  call check_coupled_steps(dt/2.0_real64, half_a, half_v)
  call near(year_a, half_a, 0.01_real64, 'coupled area timestep convergence')
  call near(year_v, half_v, 0.01_real64, 'coupled volume timestep convergence')
  call seasonal_column(3600.0_real64, year_a, year_v)
  call seasonal_column(1800.0_real64, half_a, half_v)
  call near(year_a, half_a, 0.01_real64, 'seasonal area timestep convergence')
  call near(year_v, half_v, 0.01_real64, 'seasonal volume timestep convergence')
  write (*, '(a,4f12.6)') 'seasonal mean A,V at 3600/1800 s: ', year_a, year_v, half_a, half_v
  write (*, '(a)') 'check_sea_ice: all checks passed'
contains
  subroutine near(actual, target, tolerance, name)
    real(real64), intent(in) :: actual, target, tolerance
    character(*), intent(in) :: name
    if (.not. ieee_is_finite(actual) .or. abs(actual - target) > tolerance) then
      write (*, '(a,3es24.15)') name, actual, target, tolerance
      error stop 'sea-ice assertion failed'
    end if
  end subroutine near

  subroutine update(t0, a0, v0, f, conduction, melting)
    real(real64), intent(in) :: t0, a0, v0, f, conduction, melting
    call advance_sea_ice(ice, co, dt, t0, a0, v0, f, conduction, melting, t, a, v, budget)
    call near(co*(t - t0) - latent*(v - v0), &
      dt*((1.0_real64 - a0)*f + a0*(melting - conduction)), 2.0e-5_real64, 'ocean energy closure')
    if (a < 0.0_real64 .or. a > 1.0_real64 .or. v < 0.0_real64) error stop 'ice update bounds'
    if ((a == 0.0_real64) .neqv. (v == 0.0_real64)) error stop 'ice update phase'
  end subroutine update

  subroutine check_tile_energy()
    type(radiation_config) :: rad
    real(real64) :: pressure(0:2), air(2), rhs(2), lr, dr, ocean_r, ar, vr, ti, incoming, reflected, outgoing
    real(real64) :: air_energy, surface_energy
    real(real64) :: old_rhs(2), old_surface, old_deep, old_in, old_ref, old_out
    type(sea_ice_config) :: disabled
    pressure = [1000.0_real64, 50000.0_real64, 100000.0_real64]
    air = [230.0_real64, 260.0_real64]
    call tiled_surface_tendency(rad, ice, pressure, air, 2.0_real64, 1.0_real64, 0.0_real64, &
      0.0_real64, 0.0_real64, dt, 0.3_real64, 285.0_real64, 280.0_real64, tf, 0.5_real64, 0.4_real64, &
      0.0_real64, 0.0_real64, 0.4_real64, rhs, lr, dr, ocean_r, ar, vr, ti, incoming, reflected, outgoing, budget)
    air_energy = rad%dry_air_specific_heat/rad%gravity_acceleration*sum((pressure(1:2) - pressure(0:1))*rhs)
    surface_energy = 0.3_real64*(rad%surface_heat_capacity*lr + rad%deep_ground_heat_capacity*dr) + &
      0.7_real64*(co*ocean_r - latent*vr)
    call near(air_energy + surface_energy, incoming - reflected - outgoing, 1.0e-7_real64, 'air/land/ocean/ice closure')
    ! Pure land must ignore dummy ocean/ice values, and reproduce the two-layer equations.
    rad%surface_shortwave_albedo = rad%land_shortwave_albedo
    call radiation_tendency(rad, pressure, air, 285.0_real64, 280.0_real64, 2.0_real64, 1.0_real64, &
      0.0_real64, 0.0_real64, 0.0_real64, old_rhs, old_surface, old_deep, old_in, old_ref, old_out)
    call tiled_surface_tendency(rad, disabled, pressure, air, 2.0_real64, 1.0_real64, 0.0_real64, &
      0.0_real64, 0.0_real64, dt, 1.0_real64, 285.0_real64, 280.0_real64, 0.0_real64, 0.0_real64, 0.0_real64, &
      0.0_real64, 0.0_real64, 0.0_real64, rhs, lr, dr, ocean_r, ar, vr, ti, incoming, reflected, outgoing, budget)
    call near(maxval(abs(rhs - old_rhs)), 0.0_real64, 1.0e-12_real64, 'pure land atmosphere')
    call near(lr, old_surface, 1.0e-12_real64, 'pure land surface')
    call near(dr, old_deep, 1.0e-12_real64, 'pure land deep layer')
    call near(ocean_r, 0.0_real64, 0.0_real64, 'absent ocean not updated')
    ! Without sea ice, pure ocean may supercool and must reproduce the slab equations.
    rad%slab_ocean_enabled = .true.
    rad%surface_shortwave_albedo = rad%ocean_shortwave_albedo
    call radiation_tendency(rad, pressure, air, 260.0_real64, 0.0_real64, 2.0_real64, 1.0_real64, &
      0.0_real64, 0.0_real64, 0.0_real64, old_rhs, old_surface, old_deep, old_in, old_ref, old_out)
    call tiled_surface_tendency(rad, disabled, pressure, air, 2.0_real64, 1.0_real64, 0.0_real64, &
      0.0_real64, 0.0_real64, dt, 0.0_real64, 0.0_real64, 0.0_real64, 260.0_real64, 0.0_real64, 0.0_real64, &
      0.0_real64, 0.0_real64, 0.0_real64, rhs, lr, dr, ocean_r, ar, vr, ti, incoming, reflected, outgoing, budget)
    call near(maxval(abs(rhs - old_rhs)), 0.0_real64, 1.0e-12_real64, 'pure ocean atmosphere')
    call near(ocean_r, old_surface, 1.0e-12_real64, 'ice-disabled slab')
    call near(abs(lr) + abs(dr) + abs(ar) + abs(vr), 0.0_real64, 0.0_real64, 'absent land and disabled ice')
  end subroutine check_tile_energy

  subroutine check_monthly_weighting()
    type(radiation_diagnostics) :: sample
    type(radiation_monthly_accumulator) :: accumulator
    type(radiation_monthly_means) :: means
    type(radiation_daily_accumulator) :: daily
    type(radiation_diagnostics) :: moved, day
    allocate (sample%surface_temperature(1,1), sample%deep_temperature(1,1), sample%surface_pressure(1,1))
    allocate (sample%zonal_temperature(1,1), sample%zonal_u(1,1), sample%zonal_v(1,1))
    allocate (sample%zonal_uv(1,1), sample%zonal_vt(1,1))
    allocate (sample%surface_water(1,1), sample%surface_wetness(1,1), sample%runoff(1,1))
    sample%surface_temperature = tf
    sample%deep_temperature = tf
    sample%surface_pressure = 1.0e5_real64
    sample%zonal_temperature = tf
    sample%zonal_u = 0.0_real64
    sample%zonal_v = 0.0_real64
    sample%zonal_uv = 0.0_real64
    sample%zonal_vt = 0.0_real64
    sample%surface_water = 0.0_real64
    sample%surface_wetness = 0.0_real64
    sample%runoff = 0.0_real64
    sample%land_temperature = sample%surface_temperature
    sample%ocean_temperature = sample%surface_temperature
    sample%sea_ice_fraction = 0.25_real64*sample%surface_wetness + 0.25_real64
    sample%sea_ice_volume = sample%sea_ice_fraction*2.0_real64
    sample%sea_ice_temperature = sample%surface_temperature - 20.0_real64
    call accumulator%add(sample)
    sample%sea_ice_fraction = 0.75_real64
    sample%sea_ice_volume = 0.75_real64
    sample%sea_ice_temperature = tf
    call accumulator%add(sample)
    call accumulator%take(means)
    call near(means%sea_ice_thickness(1,1), 1.25_real64, 1.0e-14_real64, 'monthly thickness weighting')
    call near(means%sea_ice_temperature(1,1), tf - 5.0_real64, 1.0e-12_real64, 'monthly ice temperature weighting')
    sample%ice_checks = checks
    call move_radiation_diagnostics(sample, moved)
    call daily%add(moved)
    call daily%add(moved)
    call daily%take(day)
    if (day%ice_checks%projected_cells /= 4) error stop 'daily projection counts lost'
    call near(day%ice_checks%projection_volume, checks%projection_volume, 0.0_real64, 'daily projection maximum')
    call daily%add(moved)
    call daily%take(day)
    if (day%ice_checks%projected_cells /= 2) error stop 'daily projection counts not reset'
  end subroutine check_monthly_weighting

  subroutine check_coupled_steps(interval, final_area, final_volume)
    real(real64), intent(in) :: interval
    real(real64), intent(out) :: final_area, final_volume
    type(harmonic_transform) :: transform
    type(hybrid_sigma_coordinate) :: coordinate
    type(dry_atmosphere_solver) :: solver
    type(dry_model_physics_config) :: physics
    complex(real64), allocatable :: phi(:,:)
    real(real64), allocatable :: land(:,:), surface(:,:), ocean(:,:), ai(:,:), vi(:,:), tl(:,:)
    real(real64), allocatable :: z(:,:,:), d(:,:,:), temp(:,:,:), ps(:,:), u(:,:,:), windv(:,:,:), skin(:,:)
    integer, allocatable :: nlon(:)
    integer :: step, j, output_status
    character(len=1024) :: output_directory
    type(radiation_monthly_means) :: output_tiles, monthly
    type(radiation_diagnostics) :: sample
    type(radiation_monthly_accumulator) :: monthly_accumulator
    type(radiation_output_options) :: options
    type(topography_config) :: terrain
    type(topography_diagnostics) :: terrain_stats
    complex(real64), allocatable :: zs(:,:,:), ds(:,:,:), temps(:,:,:), lps(:,:), qs(:,:,:)
    real(real64), allocatable :: humidity(:,:,:), deep(:,:), water_field(:,:), cloud(:,:)
    real(real64), allocatable :: ph(:), dp(:), layer(:), alpha(:), tref(:), ah(:), bh(:), q_flux(:,:)
    type(q_flux_diagnostics) :: q_flux_summary
    call transform%init(7)
    call coordinate%init_default()
    nlon = transform%get_nlon()
    call transform%allocate_field(land)
    land = 0.5_real64
    land(1,:) = 1.0_real64
    land(2,:) = 0.0_real64
    allocate (phi(0:8,0:7))
    phi = 0.0_real64
    physics = land_sea_case_physics()
    call solver%init(7, interval)
    call set_land_sea_case_state(solver, transform, physics, radiation_case_planet(physics), phi, land)
    ! The ocean starts at T_p + (T_e - T_p) cos^2(phi) >= T_f, so no initial ice forms.
    call solver%get_fields(z, d, temp, ps, u, windv, ocean_temperature=ocean, sea_ice_fraction=ai, sea_ice_volume=vi)
    do j = 1, size(nlon)
      call near(maxval(abs(ocean(1:nlon(j),j) - (physics%radiation%initial_ocean_pole_temperature + &
        (physics%radiation%initial_ocean_equator_temperature - physics%radiation%initial_ocean_pole_temperature)* &
        (1.0_real64 - transform%mu(j)**2)))), 0.0_real64, 1.0e-12_real64, 'initial ocean temperature profile')
      if (any(ai(1:nlon(j),j) /= 0.0_real64) .or. any(vi(1:nlon(j),j) /= 0.0_real64)) error stop 'initial ice formed'
    end do
    surface = land*0.0_real64 + 280.0_real64
    ocean = land*0.0_real64 + tf
    ai = land*0.0_real64 + 0.4_real64
    vi = land*0.0_real64 + 0.2_real64
    where (land == 1.0_real64)
      ai = 0.0_real64
      vi = 0.0_real64
    end where
    call solver%set_surface_state(surface, surface, ocean_temperature=ocean, sea_ice_fraction=ai, sea_ice_volume=vi)
    do step = 1, nint(120.0_real64*dt/interval)
      call solver%advance()
      call solver%take_latest_diagnostics(sample)
      if (sample%ice_checks%checked_cells /= 2*sum(nlon - 1)) error stop 'both RAW states not checked'
      call near(sample%ice_checks%energy_residual, 0.0_real64, 2.0e-5_real64, 'coupled physical energy residual')
      call near(sample%ice_checks%projection_energy_residual, 0.0_real64, 2.0e-5_real64, 'coupled projection energy')
    end do
    call solver%get_fields(z, d, temp, ps, u, windv, land_temperature=tl, ocean_temperature=ocean, &
      sea_ice_fraction=ai, sea_ice_volume=vi, surface_temperature=surface, sea_ice_temperature=skin)
    final_area = 0.0_real64
    final_volume = 0.0_real64
    do j = 1, size(nlon)
      final_area = final_area + sum(ai(1:nlon(j),j))/real(sum(nlon), real64)
      final_volume = final_volume + sum(vi(1:nlon(j),j))/real(sum(nlon), real64)
      if (.not. all(ieee_is_finite(surface(1:nlon(j),j)))) error stop 'coupled surface is nonfinite'
      if (any(ai(1:nlon(j),j) < 0.0_real64) .or. any(ai(1:nlon(j),j) > 1.0_real64) .or. &
          any(vi(1:nlon(j),j) < 0.0_real64) .or. any(ocean(1:nlon(j),j) < tf)) error stop 'coupled ice bounds'
      if (any(skin(1:nlon(j),j) > ice%melting_temperature)) error stop 'coupled skin melting cap'
      if (all(abs(tl(1:nlon(j),j) - ocean(1:nlon(j),j)) < 0.1_real64)) error stop 'land/ocean not separated'
    end do
    ! Optional integration fixture: use an externally created temporary directory.
    call get_environment_variable('BESPOKE_SEA_ICE_OUTPUT', output_directory, status=output_status)
    if (output_status == 0 .and. len_trim(output_directory) > 0 .and. interval == dt) then
      call monthly_accumulator%add(sample)
      call monthly_accumulator%take(monthly)
      options%include_moisture = .true.
      options%include_land_sea = .true.
      options%include_surface_tiles = .true.
      options%include_sea_ice = .true.
      call write_radiation_monthly_output(trim(output_directory), 1, nlon, monthly, options)
      call initialize_radiation_daily_output(trim(output_directory), options)
      call append_radiation_daily_output(trim(output_directory), sample, 1.0_real64, 1, 4, 2, options)
      call solver%get_spectral_state(zs, ds, temps, lps, specific_humidity=qs)
      call solver%get_fields(z, d, temp, ps, u, windv, surface_temperature=surface, deep_temperature=deep, &
        specific_humidity=humidity, surface_water=water_field, land_temperature=output_tiles%land_temperature, &
        ocean_temperature=output_tiles%ocean_temperature, sea_ice_fraction=output_tiles%sea_ice_fraction, &
        sea_ice_volume=output_tiles%sea_ice_volume, sea_ice_temperature=output_tiles%sea_ice_temperature, &
        sea_ice_thickness=output_tiles%sea_ice_thickness)
      call solver%get_cloud_cover(cloud)
      call write_radiation_yearly_snapshot(trim(output_directory), 1, nlon, zs, ds, temps, lps, z, d, temp, u, windv, &
        log(ps), surface, deep, options, humidity_spectral=qs, humidity=humidity, time_seconds=solver%get_time(), &
        step=solver%get_step(), cloud_cover=cloud, surface_water=water_field, tiles=output_tiles)
      call solver%get_reference_atmosphere(ph, dp, layer, alpha, tref, ah, bh)
      call solver%get_ocean_q_flux(q_flux, q_flux_summary)
      call write_radiation_metadata(trim(output_directory), 'sea_ice_test', physics, radiation_case_planet(physics), &
        7, dt, 120.0_real64*dt, 120, 0.0_real64, 0.0_real64, nlon, transform%mu, ph, dp, layer, alpha, tref, ah, bh, &
        options, terrain=terrain, terrain_diagnostics=terrain_stats, q_flux=q_flux_summary)
      call write_field(trim(output_directory)//'/ocean_q_flux.bin', nlon, q_flux)
      call write_field(trim(output_directory)//'/land_fraction.bin', nlon, land)
    end if
  end subroutine check_coupled_steps

  subroutine seasonal_column(interval, final_area, final_volume)
    real(real64), intent(in) :: interval
    real(real64), intent(out) :: final_area, final_volume
    real(real64) :: water, fraction, volume, time, forcing, skin, conduction, melting, nt, na, nv
    integer :: step
    type(sea_ice_budget) :: fluxes
    water = tf
    fraction = 0.4_real64
    volume = 0.2_real64
    final_area = 0.0_real64
    final_volume = 0.0_real64
    do step = 0, nint(360.0_real64*86400.0_real64/interval) - 1
      time = real(step,real64)*interval
      forcing = 70.0_real64*cos(2.0_real64*acos(-1.0_real64)*time/(360.0_real64*86400.0_real64))
      conduction = 0.0_real64
      melting = 0.0_real64
      if (fraction > 0.0_real64) call solve_ice_surface(ice, volume/fraction, sigma*tf**4 + forcing, &
        sigma, 2.0_real64, tf, skin, conduction, melting)
      call advance_sea_ice(ice, co, interval, water, fraction, volume, forcing, conduction, melting, nt, na, nv, fluxes)
      call near(fluxes%energy_residual, 0.0_real64, 2.0e-5_real64, 'seasonal energy closure')
      water = nt
      fraction = na
      volume = nv
      final_area = final_area + fraction*interval/(360.0_real64*86400.0_real64)
      final_volume = final_volume + volume*interval/(360.0_real64*86400.0_real64)
    end do
    if (final_area < 0.01_real64 .or. final_volume < 0.01_real64) error stop 'seasonal ice cycle is trivial'
  end subroutine seasonal_column
end program check_sea_ice
