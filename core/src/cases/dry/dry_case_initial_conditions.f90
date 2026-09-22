!> Case-side composition of the dry-atmosphere initial conditions.
!>
!> The solver offers generic operations (install a spectral state, set the
!> ground temperatures, select the planet and the active physics).  Choosing
!> Jablonowski-Williamson, Held-Suarez or the radiative-equilibrium start, and
!> which physical processes and planet each case runs with, is a case decision,
!> so it lives here instead of in the model.
module dry_case_initial_conditions
  use iso_fortran_env, only: real64
  use harmonics, only: harmonic_transform
  use dry_vertical_coordinate, only: hybrid_sigma_coordinate
  use dry_atmosphere, only: dry_atmosphere_solver
  use dry_initial_conditions, only: jablonowski_williamson_initial_state, &
                                    jablonowski_williamson_topographic_initial_state
  use dry_held_suarez, only: held_suarez_initial_state
  use dry_physics_config, only: dry_model_physics_config, radiation_config, moisture_config, &
                                radiation_planet_rotation_rate
  use moist_thermodynamics, only: saturation_specific_humidity, full_level_pressures
  use planet_parameters, only: planet_config
  implicit none
  private

  public :: set_jablonowski_williamson_case_state
  public :: held_suarez_case_physics, set_held_suarez_case_state
  public :: radiation_case_physics, slab_ocean_case_physics
  public :: radiation_case_planet, set_radiation_case_state
  public :: moist_case_physics
  public :: land_sea_case_physics, set_land_sea_case_state

contains

  !> Balanced baroclinic base state over the Jablonowski-Williamson terrain.
  !> No physical forcing is enabled.
  subroutine set_jablonowski_williamson_case_state(solver, transform, include_perturbation)
    type(dry_atmosphere_solver), intent(inout) :: solver
    type(harmonic_transform), intent(inout) :: transform
    logical, intent(in), optional :: include_perturbation
    complex(real64), allocatable :: zeta(:, :, :), delta(:, :, :), temperature(:, :, :)
    complex(real64), allocatable :: log_ps(:, :), surface_geopotential(:, :)
    type(hybrid_sigma_coordinate) :: coordinate
    logical :: perturb

    perturb = .true.
    if (present(include_perturbation)) perturb = include_perturbation
    coordinate = solver%get_coordinate()
    call jablonowski_williamson_initial_state(transform, solver%get_truncation(), coordinate, &
                                              perturb, zeta, delta, temperature, log_ps, &
                                              surface_geopotential)
    call solver%set_initial_state(zeta, delta, temperature, log_ps, surface_geopotential)
  end subroutine set_jablonowski_williamson_case_state

  !> Held-Suarez forcing: thermal relaxation together with the boundary-layer drag.
  function held_suarez_case_physics() result(physics)
    type(dry_model_physics_config) :: physics

    physics = dry_model_physics_config()
    physics%held_suarez%enabled = .true.
    physics%surface_friction%enabled = .true.
  end function held_suarez_case_physics

  !> Resting isothermal atmosphere integrated with the given (Held-Suarez) physics.
  subroutine set_held_suarez_case_state(solver, transform, physics)
    type(dry_atmosphere_solver), intent(inout) :: solver
    type(harmonic_transform), intent(inout) :: transform
    type(dry_model_physics_config), intent(in) :: physics
    complex(real64), allocatable :: zeta(:, :, :), delta(:, :, :), temperature(:, :, :)
    complex(real64), allocatable :: log_ps(:, :), surface_geopotential(:, :)
    type(hybrid_sigma_coordinate) :: coordinate

    coordinate = solver%get_coordinate()
    call held_suarez_initial_state(transform, solver%get_truncation(), coordinate, &
                                   zeta, delta, temperature, log_ps, surface_geopotential)
    call solver%set_initial_state(zeta, delta, temperature, log_ps, surface_geopotential)
    call solver%set_physics(physics)
  end subroutine set_held_suarez_case_state

  !> Radiation case physics: radiation with the ground budget, convective
  !> adjustment, the Held-Suarez boundary drag and a top-level sponge.
  function radiation_case_physics() result(physics)
    type(dry_model_physics_config) :: physics

    physics = dry_model_physics_config()
    physics%radiation%enabled = .true.
    physics%convection%enabled = .true.
    physics%surface_friction%enabled = .true.
    physics%rayleigh_friction%enabled = .true.
  end function radiation_case_physics

  !> The radiation case with its two-layer ground replaced by a 30 m slab
  !> ocean and with zero obliquity.  Every other physical process is unchanged.
  function slab_ocean_case_physics() result(physics)
    type(dry_model_physics_config) :: physics

    physics = radiation_case_physics()
    physics%radiation%slab_ocean_enabled = .true.
    physics%radiation%axial_tilt = 0.0_real64
  end function slab_ocean_case_physics

  !> The slab-ocean case with water vapour (docs/cases/moist.md): specific
  !> humidity, virtual temperature, evaporation from the saturated ocean, moist
  !> convective adjustment and large-scale condensation.  The grey longwave
  !> optical depth follows the prognostic water vapour, and the seasonal cycle of
  !> the radiation case (Earth's obliquity) replaces the zero obliquity of the
  !> slab-ocean case.  The upper Rayleigh friction is inherited.  Clouds are
  !> diagnosed and reflect shortwave, so the ocean surface takes its own albedo in
  !> place of the planetary value with clouds folded in (docs/tendency/cloud.md).
  function moist_case_physics() result(physics)
    type(dry_model_physics_config) :: physics
    type(radiation_config) :: seasonal_radiation

    physics = slab_ocean_case_physics()
    physics%radiation%axial_tilt = seasonal_radiation%axial_tilt
    physics%radiation%surface_shortwave_albedo = seasonal_radiation%ocean_shortwave_albedo
    physics%moisture%enabled = .true.
    physics%evaporation%enabled = .true.
    physics%evaporation%surface_wetness = 1.0_real64
    physics%moist_convection%enabled = .true.
    physics%condensation%enabled = .true.
    physics%cloud%enabled = .true.
  end function moist_case_physics

  !> Moist physics with one land--ocean surface budget mixed at every grid point.
  function land_sea_case_physics() result(physics)
    type(dry_model_physics_config) :: physics

    physics = moist_case_physics()
    physics%radiation%land_sea_mixing_enabled = .true.
    physics%evaporation%ocean_surface_wetness = 1.0_real64
    physics%bucket%enabled = .true.
  end function land_sea_case_physics

  !> The planet of the radiation case rotates with the calendar of its radiation
  !> configuration, so that a solar day is exactly solar_day seconds.
  function radiation_case_planet(physics) result(planet)
    type(dry_model_physics_config), intent(in) :: physics
    type(planet_config) :: planet

    planet = planet_config()
    planet%rotation_rate = radiation_planet_rotation_rate(physics%radiation)
  end function radiation_case_planet

  !> Jablonowski-Williamson temperature over flat terrain, with the unperturbed
  !> zonal wind rebalanced for Phi_s = 0 on the radiation case planet.  The
  !> active ground or ocean surface starts at the lowest model-level temperature.
  !> With prognostic water vapour the troposphere starts at the configured
  !> relative humidity (docs/cases/moist.md).
  subroutine set_radiation_case_state(solver, transform, physics, planet)
    type(dry_atmosphere_solver), intent(inout) :: solver
    type(harmonic_transform), intent(inout) :: transform
    type(dry_model_physics_config), intent(in) :: physics
    type(planet_config), intent(in) :: planet
    complex(real64), allocatable :: zeta(:, :, :), delta(:, :, :), temperature(:, :, :)
    complex(real64), allocatable :: log_ps(:, :), unused_surface_geopotential(:, :)
    complex(real64), allocatable :: state_zeta(:, :, :), state_delta(:, :, :)
    complex(real64), allocatable :: state_temperature(:, :, :), state_log_ps(:, :)
    complex(real64), allocatable :: humidity(:, :, :)
    real(real64), allocatable :: lowest_temperature(:, :)
    type(hybrid_sigma_coordinate) :: coordinate
    integer :: number_of_levels

    coordinate = solver%get_coordinate()
    number_of_levels = coordinate%number_of_levels
    call jablonowski_williamson_initial_state(transform, solver%get_truncation(), coordinate, &
                                              .false., zeta, delta, temperature, log_ps, &
                                              unused_surface_geopotential, planet%rotation_rate, &
                                              flat_terrain=.true.)
    call solver%set_planet(planet)
    call solver%set_initial_state(zeta, delta, temperature, log_ps)
    ! Read the state back so that the ground starts from the constrained
    ! spectral temperature the solver actually integrates.
    call solver%get_spectral_state(state_zeta, state_delta, state_temperature, state_log_ps)
    if (physics%moisture%enabled) then
      call initial_humidity_state(transform, coordinate, physics%moisture, state_temperature, state_log_ps, humidity)
      call solver%set_initial_state(state_zeta, state_delta, state_temperature, state_log_ps, &
                                    specific_humidity=humidity)
      call solver%get_spectral_state(state_zeta, state_delta, state_temperature, state_log_ps)
    end if
    ! T_s = T_d = T_N on the grid (docs/cases/radiation.md).
    call transform%spectral_to_grid(state_temperature(:, :, number_of_levels), lowest_temperature)
    call solver%set_surface_state(lowest_temperature, lowest_temperature)
    call solver%set_physics(physics)
  end subroutine set_radiation_case_state

  !> Pressure-coordinate radiation basic state placed over the supplied terrain.
  subroutine set_land_sea_case_state(solver, transform, physics, planet, surface_geopotential, land_fraction)
    type(dry_atmosphere_solver), intent(inout) :: solver
    type(harmonic_transform), intent(inout) :: transform
    type(dry_model_physics_config), intent(in) :: physics
    type(planet_config), intent(in) :: planet
    complex(real64), intent(in) :: surface_geopotential(0:, 0:)
    real(real64), intent(in) :: land_fraction(:, :)
    complex(real64), allocatable :: zeta(:, :, :), delta(:, :, :), temperature(:, :, :), log_ps(:, :)
    complex(real64), allocatable :: state_zeta(:, :, :), state_delta(:, :, :), state_temperature(:, :, :)
    complex(real64), allocatable :: state_log_ps(:, :), humidity(:, :, :)
    real(real64), allocatable :: surface_water(:, :), lowest_temperature(:, :)
    type(hybrid_sigma_coordinate) :: coordinate
    integer :: number_of_levels

    coordinate = solver%get_coordinate()
    number_of_levels = coordinate%number_of_levels
    call jablonowski_williamson_topographic_initial_state(transform, solver%get_truncation(), coordinate, &
      surface_geopotential, planet%rotation_rate, zeta, delta, temperature, log_ps)
    call solver%set_planet(planet)
    call solver%set_initial_state(zeta, delta, temperature, log_ps, surface_geopotential, &
                                  land_fraction=land_fraction)
    call solver%get_spectral_state(state_zeta, state_delta, state_temperature, state_log_ps)
    call initial_humidity_state(transform, coordinate, physics%moisture, state_temperature, state_log_ps, humidity)
    call solver%set_initial_state(state_zeta, state_delta, state_temperature, state_log_ps, surface_geopotential, &
                                  specific_humidity=humidity, land_fraction=land_fraction)
    call solver%get_spectral_state(state_zeta, state_delta, state_temperature, state_log_ps)
    allocate (surface_water, mold=land_fraction)
    surface_water = 0.0_real64
    where (land_fraction > 0.0_real64) surface_water = physics%bucket%initial_water
    call transform%spectral_to_grid(state_temperature(:, :, number_of_levels), lowest_temperature)
    call solver%set_surface_state(lowest_temperature, lowest_temperature, surface_water)
    call solver%set_physics(physics)
  end subroutine set_land_sea_case_state

  !> Spectral specific humidity q = RH q_s(T_k, p_k) on the levels whose full-level
  !> pressure is at least the configured top pressure, and zero above.  The
  !> Jablonowski-Williamson stratosphere is warm and thin, where the saturation
  !> humidity is not meaningful (q_s -> 1 as p -> e_s), so it starts dry.
  subroutine initial_humidity_state(transform, coordinate, moisture, temperature, log_surface_pressure, humidity)
    type(harmonic_transform), intent(inout) :: transform
    type(hybrid_sigma_coordinate), intent(in) :: coordinate
    type(moisture_config), intent(in) :: moisture
    complex(real64), intent(in) :: temperature(0:, 0:, :), log_surface_pressure(0:, 0:)
    complex(real64), allocatable, intent(out) :: humidity(:, :, :)
    real(real64), allocatable :: grid(:, :), log_ps_grid(:, :)
    real(real64), allocatable :: temperature_grid(:, :, :), humidity_grid(:, :, :)
    complex(real64), allocatable :: spectral(:, :)
    real(real64) :: pressure_half(0:coordinate%number_of_levels)
    real(real64) :: full_level_pressure(coordinate%number_of_levels), delta_pressure(coordinate%number_of_levels)
    real(real64) :: surface_pressure
    integer, allocatable :: nlon(:)
    integer :: i, j, k, levels

    levels = coordinate%number_of_levels
    nlon = transform%get_nlon()
    call transform%allocate_field(grid)
    call transform%allocate_field(log_ps_grid)
    allocate (temperature_grid(size(grid, 1), size(grid, 2), levels))
    allocate (humidity_grid(size(grid, 1), size(grid, 2), levels))
    do k = 1, levels
      call transform%spectral_to_grid(temperature(:, :, k), grid)
      temperature_grid(:, :, k) = grid
    end do
    call transform%spectral_to_grid(log_surface_pressure, log_ps_grid)
    humidity_grid = 0.0_real64
    do j = 1, size(grid, 2)
      do i = 1, nlon(j)
        surface_pressure = exp(log_ps_grid(i, j))
        do k = 0, levels
          pressure_half(k) = coordinate%a_half(k) + coordinate%b_half(k)*surface_pressure
        end do
        call full_level_pressures(pressure_half, full_level_pressure, delta_pressure)
        do k = 1, levels
          if (full_level_pressure(k) >= moisture%initial_humidity_top_pressure) then
            humidity_grid(i, j, k) = moisture%initial_relative_humidity* &
              saturation_specific_humidity(temperature_grid(i, j, k), full_level_pressure(k))
          end if
        end do
      end do
    end do
    allocate (humidity(0:ubound(temperature, 1), 0:ubound(temperature, 2), levels))
    do k = 1, levels
      call transform%grid_to_spectral(humidity_grid(:, :, k), spectral)
      humidity(:, :, k) = spectral(0:ubound(temperature, 1), 0:ubound(temperature, 2))
    end do
  end subroutine initial_humidity_state

end module dry_case_initial_conditions
