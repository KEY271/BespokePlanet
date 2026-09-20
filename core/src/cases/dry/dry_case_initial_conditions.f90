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
  use dry_initial_conditions, only: jablonowski_williamson_initial_state
  use dry_held_suarez, only: held_suarez_initial_state
  use dry_physics_config, only: dry_model_physics_config, radiation_planet_rotation_rate
  use planet_parameters, only: planet_config
  implicit none
  private

  public :: set_jablonowski_williamson_case_state
  public :: held_suarez_case_physics, set_held_suarez_case_state
  public :: radiation_case_physics, slab_ocean_case_physics
  public :: radiation_case_planet, set_radiation_case_state

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
  subroutine set_radiation_case_state(solver, transform, physics, planet)
    type(dry_atmosphere_solver), intent(inout) :: solver
    type(harmonic_transform), intent(inout) :: transform
    type(dry_model_physics_config), intent(in) :: physics
    type(planet_config), intent(in) :: planet
    complex(real64), allocatable :: zeta(:, :, :), delta(:, :, :), temperature(:, :, :)
    complex(real64), allocatable :: log_ps(:, :), unused_surface_geopotential(:, :)
    complex(real64), allocatable :: state_zeta(:, :, :), state_delta(:, :, :)
    complex(real64), allocatable :: state_temperature(:, :, :), state_log_ps(:, :)
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
    call solver%set_surface_state(state_temperature(:, :, number_of_levels), &
                                  state_temperature(:, :, number_of_levels))
    call solver%set_physics(physics)
  end subroutine set_radiation_case_state

end module dry_case_initial_conditions
