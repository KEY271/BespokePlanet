!> Checks of docs/tendency/snow.md: the snowfall of the large-scale condensation,
!> its melting in warm layers, the land snowpack and its coupling to the tiled
!> surface, and a short land-sea integration with snow.
program check_snow
  use iso_fortran_env, only: real64
  use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
  use planet_parameters, only: earth_gravity, planet_config
  use dry_physics_config, only: condensation_config, snow_config, radiation_config, sea_ice_config, &
                                dry_model_physics_config
  use moist_thermodynamics, only: dry_air_specific_heat, latent_heat_of_condensation, saturation_specific_humidity
  use large_scale_condensation, only: large_scale_condensation_from_levels, saturation_condensate
  use land_snow, only: snow_cover_fraction, advance_land_snow
  use surface_tiles, only: tiled_surface_tendency
  use surface_exchange, only: surface_sensible_heat_flux
  use sea_ice, only: sea_ice_budget
  use harmonics, only: harmonic_transform
  use topography, only: topography_config, topography_diagnostics, generate_topography
  use dry_atmosphere, only: dry_atmosphere_solver
  use dry_case_initial_conditions, only: land_sea_case_physics, moist_case_physics, radiation_case_planet, &
                                         set_land_sea_case_state
  use numerics_config, only: model_numerics_config
  use dry_radiation, only: radiation_diagnostics
  implicit none

  real(real64), parameter :: interval = 2400.0_real64

  call check_column_snowfall()
  call check_disabled_snow_is_rain()
  call check_warm_column_melts_everything()
  call check_land_snowpack()
  call check_tile_coupling()
  call check_case_configuration()
  call check_land_integration()
  write (*, '(a)') 'check_snow: all checks passed'

contains

  subroutine near(actual, expected, tolerance, label)
    real(real64), intent(in) :: actual, expected, tolerance
    character(*), intent(in) :: label
    if (.not. ieee_is_finite(actual) .or. abs(actual - expected) > tolerance) then
      write (*, '(a,3es24.15)') 'snow check failed: '//label//' actual, expected, tolerance = ', &
        actual, expected, tolerance
      error stop 'snow assertion failed'
    end if
  end subroutine near

  !> A column of 5 layers (top first): two supersaturated layers colder than
  !> -10 C make snow, a supersaturated layer at -5 C makes rain, and two dry
  !> layers warmer than 5 C melt part of the falling snow.
  subroutine column(pressure, thickness, temperature, humidity)
    real(real64), intent(out) :: pressure(5), thickness(5), temperature(5), humidity(5)
    integer :: k
    pressure = [3.0e4_real64, 4.5e4_real64, 6.0e4_real64, 7.5e4_real64, 9.0e4_real64]
    thickness = 1.5e4_real64
    temperature = [240.0_real64, 255.0_real64, 268.15_real64, 278.17_real64, 278.19_real64]
    do k = 1, 3
      humidity(k) = 1.2_real64*saturation_specific_humidity(temperature(k), pressure(k))
    end do
    do k = 4, 5
      humidity(k) = 0.5_real64*saturation_specific_humidity(temperature(k), pressure(k))
    end do
  end subroutine column

  subroutine check_column_snowfall()
    type(condensation_config) :: condensation
    type(snow_config) :: snow
    real(real64) :: pressure(5), thickness(5), temperature(5), humidity(5)
    real(real64) :: temperature_tendency(5), humidity_tendency(5), precipitation, snowfall, melt
    real(real64) :: produced_snow, layer_mass(5), column_enthalpy, final_temperature(5), condensate, fusion
    real(real64) :: sublimation_heat
    integer :: k

    snow%enabled = .true.
    fusion = snow%latent_heat_of_fusion
    sublimation_heat = latent_heat_of_condensation + fusion
    call column(pressure, thickness, temperature, humidity)
    call large_scale_condensation_from_levels(condensation, pressure, thickness, temperature, humidity, interval, &
      temperature_tendency, humidity_tendency, precipitation, snow, snowfall, melt)
    layer_mass = thickness/earth_gravity
    produced_snow = -sum(layer_mass(1:2)*humidity_tendency(1:2))
    call near(precipitation, -sum(layer_mass*humidity_tendency), 1.0e-15_real64, 'precipitation is the condensate')
    call near(snowfall + melt, produced_snow, 1.0e-15_real64, 'snow mass is conserved while falling')
    if (snowfall <= 0.0_real64 .or. melt <= 0.0_real64) error stop 'column should both melt and deliver snow'
    if (any(humidity_tendency(4:5) /= 0.0_real64)) error stop 'unsaturated layers condensed'
    ! The snow layers release L_v + L_f and end saturated with that heating.
    do k = 1, 2
      condensate = -interval*humidity_tendency(k)
      call near(temperature_tendency(k), sublimation_heat/dry_air_specific_heat*condensate/interval, &
                1.0e-15_real64, 'snow layer heating')
      call near(humidity(k) - condensate, saturation_specific_humidity(temperature(k) + &
        sublimation_heat/dry_air_specific_heat*condensate, pressure(k)), &
        condensation%saturation_tolerance*saturation_specific_humidity(temperature(k), pressure(k)), &
        'snow layer saturation')
      call near(condensate, saturation_condensate(condensation, temperature(k), humidity(k), pressure(k), &
                sublimation_heat), 0.0_real64, 'snow condensate uses L_v + L_f')
    end do
    ! The rain layer at -5 C keeps L_v.
    condensate = -interval*humidity_tendency(3)
    call near(temperature_tendency(3), latent_heat_of_condensation/dry_air_specific_heat*condensate/interval, &
              1.0e-15_real64, 'rain layer heating')
    ! Melting layers cool, but never below 5 C; the upper one melts to exactly 5 C.
    final_temperature = temperature + interval*temperature_tendency
    call near(final_temperature(4), snow%atmospheric_melting_temperature, 1.0e-10_real64, 'melting layer at 5 C')
    call near(final_temperature(5), snow%atmospheric_melting_temperature, 1.0e-10_real64, 'second melting layer')
    ! The column gains L_f of the snow reaching the ground and nothing else.
    column_enthalpy = sum(layer_mass*(dry_air_specific_heat*temperature_tendency + &
      latent_heat_of_condensation*humidity_tendency))
    call near(column_enthalpy, fusion*snowfall, 1.0e-9_real64*fusion*snowfall, 'column enthalpy gains L_f snowfall')
  end subroutine check_column_snowfall

  !> Without an enabled snow configuration the old rain-only path is reproduced bit for bit.
  subroutine check_disabled_snow_is_rain()
    type(condensation_config) :: condensation
    type(snow_config) :: snow
    real(real64) :: pressure(5), thickness(5), temperature(5), humidity(5)
    real(real64) :: t_old(5), q_old(5), p_old, t_new(5), q_new(5), p_new, snowfall, melt

    call column(pressure, thickness, temperature, humidity)
    call large_scale_condensation_from_levels(condensation, pressure, thickness, temperature, humidity, interval, &
      t_old, q_old, p_old)
    call large_scale_condensation_from_levels(condensation, pressure, thickness, temperature, humidity, interval, &
      t_new, q_new, p_new, snow, snowfall, melt)
    if (any(t_new /= t_old) .or. any(q_new /= q_old) .or. p_new /= p_old .or. snowfall /= 0.0_real64 .or. &
        melt /= 0.0_real64) error stop 'disabled snow changed the rain-only condensation'
    call near(sum(thickness/earth_gravity*(dry_air_specific_heat*t_old + latent_heat_of_condensation*q_old)), &
              0.0_real64, 1.0e-9_real64, 'rain-only column enthalpy')
  end subroutine check_disabled_snow_is_rain

  !> A deep warm layer melts all the snow; the rest of the column is untouched.
  subroutine check_warm_column_melts_everything()
    type(condensation_config) :: condensation
    type(snow_config) :: snow
    real(real64) :: pressure(2), thickness(2), temperature(2), humidity(2)
    real(real64) :: temperature_tendency(2), humidity_tendency(2), precipitation, snowfall, melt, fusion

    snow%enabled = .true.
    fusion = snow%latent_heat_of_fusion
    pressure = [4.0e4_real64, 8.5e4_real64]
    thickness = [2.0e4_real64, 3.0e4_real64]
    temperature = [245.0_real64, 290.0_real64]
    humidity = [1.5_real64*saturation_specific_humidity(temperature(1), pressure(1)), 0.0_real64]
    call large_scale_condensation_from_levels(condensation, pressure, thickness, temperature, humidity, interval, &
      temperature_tendency, humidity_tendency, precipitation, snow, snowfall, melt)
    call near(snowfall, 0.0_real64, 0.0_real64, 'no snow reaches the warm ground')
    call near(melt, precipitation, 1.0e-15_real64, 'all snow melted in the warm layer')
    call near(thickness(2)/earth_gravity*dry_air_specific_heat*temperature_tendency(2), -fusion*melt, &
              1.0e-9_real64, 'melting cools by L_f')
    if (temperature(2) + interval*temperature_tendency(2) <= snow%atmospheric_melting_temperature) &
      error stop 'warm layer cooled below the melting temperature'
  end subroutine check_warm_column_melts_everything

  subroutine check_land_snowpack()
    type(snow_config) :: snow
    real(real64), parameter :: capacity = 2.0e6_real64
    real(real64) :: rhs, snow_rhs, melt, new_temperature, fusion

    snow%enabled = .true.
    fusion = snow%latent_heat_of_fusion
    call near(snow_cover_fraction(snow, 0.0_real64), 0.0_real64, 0.0_real64, 'no snow, no cover')
    call near(snow_cover_fraction(snow, 50.0_real64), 0.5_real64, 1.0e-15_real64, 'S = S_0 gives half cover')
    call near(snow_cover_fraction(snow, 150.0_real64), 0.75_real64, 1.0e-15_real64, 'S = 3 S_0 gives 3/4')
    ! A cold surface keeps its temperature tendency and accumulates the snowfall.
    rhs = -1.0e-4_real64
    call advance_land_snow(snow, capacity, interval, 270.0_real64, 10.0_real64, 1.0e-3_real64, rhs, snow_rhs, melt)
    call near(rhs, -1.0e-4_real64, 0.0_real64, 'cold land temperature unchanged')
    call near(snow_rhs, 1.0e-3_real64, 1.0e-16_real64, 'cold snowpack grows by the snowfall')
    call near(melt, 0.0_real64, 0.0_real64, 'no melt below 275 K')
    ! Heat above 275 K melts snow and brings the land back to 275 K.
    rhs = 1.0e-3_real64
    call advance_land_snow(snow, capacity, interval, 274.0_real64, 100.0_real64, 0.0_real64, rhs, snow_rhs, melt)
    new_temperature = 274.0_real64 + interval*rhs
    call near(new_temperature, snow%surface_melting_temperature, 1.0e-10_real64, 'melting holds land at 275 K')
    call near(melt*interval, capacity*(274.0_real64 + interval*1.0e-3_real64 - 275.0_real64)/fusion, 1.0e-10_real64, &
              'melt equals excess heat over L_f')
    call near(snow_rhs, -melt, 1.0e-16_real64, 'snowpack loses the melt')
    ! A thin snowpack melts away and the remaining heat warms the land.
    rhs = 1.0e-2_real64
    call advance_land_snow(snow, capacity, interval, 280.0_real64, 1.0_real64, 0.0_real64, rhs, snow_rhs, melt)
    call near(1.0_real64 + interval*snow_rhs, 0.0_real64, 1.0e-12_real64, 'thin snowpack melts completely')
    call near(capacity*(280.0_real64 + interval*rhs) - capacity*(280.0_real64 + interval*1.0e-2_real64), &
              -fusion*1.0_real64, 1.0e-6_real64, 'only the melted snow takes heat')
  end subroutine check_land_snowpack

  !> The tiled surface: albedo, (1 - f) sensible heat, land melt and ocean snowfall all conserve energy.
  subroutine check_tile_coupling()
    type(radiation_config) :: rad
    type(sea_ice_config) :: ice
    type(snow_config) :: snow
    type(sea_ice_budget) :: budget
    real(real64) :: pressure(0:2), air(2), rhs(2), lr, dr, ocean_r, ar, vr, ti, incoming, reflected, outgoing
    real(real64) :: bare_rhs(2), bare_lr, bare_reflected, snow_rhs, melt, air_energy, surface_energy
    real(real64) :: cover, fusion, snowfall, bare_ocean_r, sensible_bare, sensible_snow
    real(real64), parameter :: latent = 30.0_real64

    snow%enabled = .true.
    fusion = snow%latent_heat_of_fusion
    pressure = [1000.0_real64, 50000.0_real64, 100000.0_real64]
    air = [230.0_real64, 260.0_real64]
    snowfall = 2.0e-4_real64
    cover = snow_cover_fraction(snow, 50.0_real64)
    ! Pure land below the melting temperature: brighter, less sensible heat, and snowfall accumulates.
    call tiled_surface_tendency(rad, ice, pressure, air, 2.0_real64, 1.0_real64, 0.0_real64, 0.0_real64, 0.0_real64, &
      interval, 1.0_real64, 268.0_real64, 270.0_real64, 0.0_real64, 0.0_real64, 0.0_real64, latent, 0.0_real64, &
      0.0_real64, bare_rhs, bare_lr, dr, ocean_r, ar, vr, ti, incoming, bare_reflected, outgoing, budget)
    call tiled_surface_tendency(rad, ice, pressure, air, 2.0_real64, 1.0_real64, 0.0_real64, 0.0_real64, 0.0_real64, &
      interval, 1.0_real64, 268.0_real64, 270.0_real64, 0.0_real64, 0.0_real64, 0.0_real64, latent, 0.0_real64, &
      0.0_real64, rhs, lr, dr, ocean_r, ar, vr, ti, incoming, reflected, outgoing, budget, snow=snow, &
      snow_water=50.0_real64, snowfall=snowfall, snow_rhs=snow_rhs, snow_melt=melt)
    call near(reflected/bare_reflected, (rad%land_shortwave_albedo + snow%albedo_increase*cover)/ &
              rad%land_shortwave_albedo, 1.0e-12_real64, 'snow brightens the land by 0.4 f')
    ! Only the sensible heat into the lowest layer differs, by -f H of the bare land.
    sensible_bare = surface_sensible_heat_flux(rad, pressure, air(2), 268.0_real64, 2.0_real64, 1.0_real64)
    sensible_snow = (pressure(2) - pressure(1))*rad%dry_air_specific_heat/rad%gravity_acceleration*(rhs(2) - bare_rhs(2))
    call near(sensible_snow, -cover*sensible_bare, 1.0e-9_real64*abs(sensible_bare), 'sensible heat times 1 - f')
    call near(rhs(1), bare_rhs(1), 1.0e-15_real64, 'upper layer unaffected by the snowpack')
    call near(melt, 0.0_real64, 0.0_real64, 'cold land does not melt')
    call near(snow_rhs, snowfall, 1.0e-16_real64, 'snowfall accumulates on cold land')
    ! Energy of air, land and snowpack: the snowfall brings -L_f of snow to the surface.
    air_energy = rad%dry_air_specific_heat/rad%gravity_acceleration*sum((pressure(1:2) - pressure(0:1))*rhs)
    surface_energy = rad%surface_heat_capacity*lr + rad%deep_ground_heat_capacity*dr - fusion*snow_rhs
    call near(air_energy + surface_energy, incoming - reflected - outgoing - latent - fusion*snowfall, 1.0e-7_real64, &
              'cold land energy closure')
    ! A warm land surface melts snow at 275 K with the same closure.
    call tiled_surface_tendency(rad, ice, pressure, air, 2.0_real64, 1.0_real64, 0.0_real64, 0.0_real64, 0.0_real64, &
      interval, 1.0_real64, 276.0_real64, 270.0_real64, 0.0_real64, 0.0_real64, 0.0_real64, latent, 0.0_real64, &
      0.0_real64, rhs, lr, dr, ocean_r, ar, vr, ti, incoming, reflected, outgoing, budget, snow=snow, &
      snow_water=50.0_real64, snowfall=snowfall, snow_rhs=snow_rhs, snow_melt=melt)
    call near(276.0_real64 + interval*lr, snow%surface_melting_temperature, 1.0e-10_real64, 'warm land melts to 275 K')
    call near(snow_rhs, snowfall - melt, 1.0e-16_real64, 'snowpack budget')
    air_energy = rad%dry_air_specific_heat/rad%gravity_acceleration*sum((pressure(1:2) - pressure(0:1))*rhs)
    surface_energy = rad%surface_heat_capacity*lr + rad%deep_ground_heat_capacity*dr - fusion*snow_rhs
    call near(air_energy + surface_energy, incoming - reflected - outgoing - latent - fusion*snowfall, 1.0e-7_real64, &
              'melting land energy closure')
    ! Snow falling on an ice-free ocean melts with the heat of the mixed layer.
    call tiled_surface_tendency(rad, ice, pressure, air, 2.0_real64, 1.0_real64, 0.0_real64, 0.0_real64, 0.0_real64, &
      interval, 0.0_real64, 0.0_real64, 0.0_real64, 280.0_real64, 0.0_real64, 0.0_real64, 0.0_real64, latent, &
      0.0_real64, rhs, lr, dr, bare_ocean_r, ar, vr, ti, incoming, reflected, outgoing, budget)
    call tiled_surface_tendency(rad, ice, pressure, air, 2.0_real64, 1.0_real64, 0.0_real64, 0.0_real64, 0.0_real64, &
      interval, 0.0_real64, 0.0_real64, 0.0_real64, 280.0_real64, 0.0_real64, 0.0_real64, 0.0_real64, latent, &
      0.0_real64, rhs, lr, dr, ocean_r, ar, vr, ti, incoming, reflected, outgoing, budget, snow=snow, &
      snow_water=0.0_real64, snowfall=snowfall, snow_rhs=snow_rhs, snow_melt=melt)
    call near(rad%seawater_density*rad%seawater_specific_heat*rad%slab_ocean_depth*(ocean_r - bare_ocean_r), &
              -fusion*snowfall, 1.0e-9_real64, 'ocean pays L_f of the snowfall')
    call near(snow_rhs + melt, 0.0_real64, 0.0_real64, 'no snowpack on the ocean')
  end subroutine check_tile_coupling

  subroutine check_case_configuration()
    type(dry_model_physics_config) :: land, moist
    land = land_sea_case_physics()
    moist = moist_case_physics()
    if (.not. land%snow%enabled .or. moist%snow%enabled) error stop 'snow must be enabled only in the land cases'
    call near(land%snow%formation_temperature, 263.15_real64, 0.0_real64, 'snow formation at -10 C')
    call near(land%snow%atmospheric_melting_temperature, 278.15_real64, 0.0_real64, 'atmospheric melting at 5 C')
    call near(land%snow%surface_melting_temperature, 275.0_real64, 0.0_real64, 'surface melting at 275 K')
    call near(land%snow%masking_water_equivalent, 50.0_real64, 0.0_real64, 'S_0 = 0.05 m water equivalent')
    call near(land%snow%albedo_increase, 0.4_real64, 0.0_real64, 'albedo increase')
  end subroutine check_case_configuration

  !> A short land-sea run starting from a thick snowpack on every land point.
  subroutine check_land_integration()
    type(harmonic_transform) :: transform
    type(topography_config) :: terrain
    type(topography_diagnostics) :: terrain_diagnostics
    type(dry_model_physics_config) :: physics
    type(planet_config) :: planet
    type(model_numerics_config) :: numerics
    type(dry_atmosphere_solver) :: solver
    type(radiation_diagnostics) :: diagnostics
    real(real64), allocatable :: land_fraction(:, :), analytic_height(:, :), height(:, :)
    real(real64), allocatable :: zeta(:, :, :), delta(:, :, :), temperature(:, :, :)
    real(real64), allocatable :: surface_pressure(:, :), u(:, :, :), v(:, :, :), snow_water(:, :), cover(:, :)
    real(real64), allocatable :: surface_water(:, :)
    complex(real64), allocatable :: surface_geopotential(:, :)
    integer :: step

    call transform%init(31)
    call generate_topography(transform, terrain, land_fraction, analytic_height, surface_geopotential, &
                             height, terrain_diagnostics)
    physics = land_sea_case_physics()
    physics%snow%initial_water_equivalent = 100.0_real64
    planet = radiation_case_planet(physics)
    numerics = model_numerics_config()
    numerics%truncation = 31
    numerics%time_step = 1200.0_real64
    call solver%init_with_config(numerics)
    call set_land_sea_case_state(solver, transform, physics, planet, surface_geopotential, land_fraction)
    do step = 1, 6
      call solver%advance()
      call solver%take_latest_diagnostics(diagnostics)
    end do
    call solver%get_fields(zeta, delta, temperature, surface_pressure, u, v, snow_water=snow_water, &
                           snow_fraction=cover, surface_water=surface_water)
    if (.not. all(ieee_is_finite(snow_water)) .or. .not. all(ieee_is_finite(cover))) &
      error stop 'snow integration produced a non-finite snowpack'
    if (minval(snow_water) < 0.0_real64 .or. any(cover < 0.0_real64) .or. any(cover >= 1.0_real64)) &
      error stop 'snow integration snowpack is out of range'
    if (any(land_fraction == 0.0_real64 .and. snow_water /= 0.0_real64)) error stop 'snowpack on pure ocean'
    if (maxval(snow_water, mask=land_fraction > 0.0_real64) <= 0.0_real64) error stop 'the initial snowpack vanished'
    if (minval(surface_water) < 0.0_real64 .or. maxval(surface_water) > physics%bucket%capacity) &
      error stop 'snow integration bucket is out of range'
    if (.not. allocated(diagnostics%snow_water) .or. .not. allocated(diagnostics%snowfall) .or. &
        .not. allocated(diagnostics%snow_melt) .or. .not. allocated(diagnostics%snow_fraction)) &
      error stop 'snow integration did not produce snow diagnostics'
    if (diagnostics%mean_snow_water <= 0.0_real64 .or. diagnostics%mean_snow_water > 100.0_real64 + 1.0_real64 .or. &
        diagnostics%mean_snow_fraction <= 0.0_real64 .or. diagnostics%mean_snow_fraction >= 1.0_real64 .or. &
        diagnostics%mean_snowfall < 0.0_real64 .or. diagnostics%mean_snow_melt < 0.0_real64) &
      error stop 'snow integration diagnostics are out of range'
    call near(diagnostics%maximum_snow_budget_residual, 0.0_real64, 1.0e-15_real64, 'snowpack budget residual')
    write (*, '(a,4es13.5)') 'snow run: land S, f, snowfall, melt = ', diagnostics%mean_snow_water, &
      diagnostics%mean_snow_fraction, diagnostics%mean_snowfall, diagnostics%mean_snow_melt
  end subroutine check_land_integration

end program check_snow
