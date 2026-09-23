!> Composes the explicit right-hand side of the dry or moist atmosphere from the
!> individual dynamical and physical tendencies.
!>
!> Each tendency adds into the shared grid-space accumulators of the workspace;
!> the projection then performs the spectral transforms once.  Adding a new
!> physical process means adding one call here, not editing a kernel.
!>
!> The order of the calls is the order in which the contributions are summed.
!> Radiation, surface fluxes, friction and the convective processes are all
!> evaluated independently on the RAW-filtered previous time level; only the
!> convective adjustments, the condensation and the cloud diagnosis are chained
!> through provisional fields (dry_convection_tendency).  The evaporation is
!> evaluated before the radiation so that the surface budget can take the same
!> latent heat flux, and the convective processes before the radiation so that
!> the shortwave reflection can take the cloud cover of the same evaluation
!> (docs/dynamics/moist.md, "physical processes").  The land bucket follows the
!> radiation, which advances the land snowpack, so that the bucket takes the
!> rain and the snowmelt of the same evaluation (docs/tendency/snow.md).
module dry_tendency_evaluator
  use iso_fortran_env, only: real64
  use harmonics, only: harmonic_transform
  use planet_parameters, only: planet_config
  use dry_vertical_coordinate, only: hybrid_sigma_coordinate
  use dry_state, only: dry_state_type, dry_tendency_type, zero_dry_tendency
  use dry_physics_config, only: dry_model_physics_config
  use dry_radiation, only: radiation_diagnostics
  use dry_tendency_workspace, only: dry_workspace_type
  use dry_dynamics_tendency, only: add_dry_dynamics_tendency
  use dry_surface_friction_tendency, only: add_dry_surface_friction_tendency
  use dry_rayleigh_friction_tendency, only: add_dry_rayleigh_friction_tendency
  use dry_held_suarez_tendency, only: add_dry_held_suarez_tendency
  use dry_evaporation_tendency, only: add_dry_evaporation_tendency
  use dry_bucket_tendency, only: add_dry_bucket_tendency
  use sea_ice, only: validate_sea_ice_config
  use land_snow, only: validate_snow_config
  use dry_radiation_tendency, only: add_dry_radiation_tendency, diagnose_surface_tiles
  use dry_convection_tendency, only: add_dry_convection_tendency
  use dry_tendency_diagnostics, only: collect_dry_diagnostics
  use dry_tendency_projection, only: project_dry_tendency
  implicit none
  private
  public :: evaluate_dry_tendency

contains

  !> `interval` is the width the leapfrog advances the state by with this
  !> tendency (2 dt in a regular step); the chained convective processes build
  !> their provisional fields with it.
  subroutine evaluate_dry_tendency(transform, truncation, coordinate, planet, state, physics_state, &
                                   surface_geopotential, physics, workspace, evaluation_time, interval, rhs, &
                                   maximum_speed, diagnostics, land_fraction)
    type(harmonic_transform), intent(inout) :: transform
    integer, intent(in) :: truncation
    type(hybrid_sigma_coordinate), intent(in) :: coordinate
    type(planet_config), intent(in) :: planet
    !> state advances the dynamics; physics_state is the RAW-filtered previous
    !> time level that every prescribed physical tendency is evaluated on.
    type(dry_state_type), intent(in) :: state, physics_state
    complex(real64), intent(in) :: surface_geopotential(0:, 0:)
    real(real64), intent(in), optional :: land_fraction(:, :)
    type(dry_model_physics_config), intent(in) :: physics
    type(dry_workspace_type), intent(inout) :: workspace
    real(real64), intent(in) :: evaluation_time, interval
    type(dry_tendency_type), intent(inout) :: rhs
    real(real64), intent(out) :: maximum_speed
    type(radiation_diagnostics), intent(out), optional :: diagnostics

    if (present(diagnostics) .and. .not. physics%radiation%enabled) then
      error stop 'radiation diagnostics requested without radiation forcing'
    end if
    if ((physics%evaporation%enabled .or. physics%moist_convection%enabled .or. physics%condensation%enabled) &
        .and. .not. physics%moisture%enabled) then
      error stop 'moist processes require the moisture (specific humidity) prognostic variable'
    end if
    if (physics%evaporation%enabled .and. .not. physics%radiation%enabled) then
      error stop 'evaporation requires the radiation surface energy budget'
    end if
    if (physics%bucket%enabled .and. (.not. physics%evaporation%enabled .or. &
        .not. physics%moisture%enabled .or. .not. physics%radiation%land_sea_mixing_enabled)) then
      error stop 'land bucket requires moisture, evaporation and land-sea mixing'
    end if
    if (physics%bucket%enabled .and. physics%bucket%capacity <= 0.0_real64) then
      error stop 'land bucket capacity must be positive'
    end if
    if (physics%bucket%enabled .and. (physics%bucket%initial_water < 0.0_real64 .or. &
        physics%bucket%initial_water > physics%bucket%capacity)) then
      error stop 'land bucket initial water must lie between zero and capacity'
    end if
    if (physics%bucket%enabled .and. (physics%bucket%dry_threshold_fraction < 0.0_real64 .or. &
        physics%bucket%dry_threshold_fraction > 1.0_real64)) then
      error stop 'land bucket dry threshold fraction must lie between zero and one'
    end if
    if (physics%cloud%enabled .and. .not. physics%moisture%enabled) then
      error stop 'the cloud diagnosis requires the moisture (specific humidity) prognostic variable'
    end if

    if (physics%sea_ice%enabled) then
      call validate_sea_ice_config(physics%sea_ice)
      if (.not. physics%radiation%enabled .or. &
          .not. (physics%radiation%slab_ocean_enabled .or. physics%radiation%land_sea_mixing_enabled)) &
        error stop 'sea ice requires a radiative ocean surface'
    end if
    if (physics%snow%enabled) call validate_snow_config(physics)
    if (physics%q_flux%enabled .and. (.not. physics%radiation%enabled .or. &
        .not. (physics%sea_ice%enabled .or. physics%radiation%land_sea_mixing_enabled))) &
      error stop 'Q flux requires the tiled ocean surface (sea ice or land-sea mixing)'
    call zero_dry_tendency(rhs)
    call workspace%zero_forcing()
    if (present(land_fraction)) then
      call workspace%prepare(transform, coordinate, planet%rotation_rate, state, physics_state, &
                             surface_geopotential, physics, evaluation_time, land_fraction)
    else
      call workspace%prepare(transform, coordinate, planet%rotation_rate, state, physics_state, &
                             surface_geopotential, physics, evaluation_time)
    end if

    call add_dry_dynamics_tendency(transform, state, workspace, maximum_speed)

    if (physics%surface_friction%enabled) then
      call add_dry_surface_friction_tendency(physics%surface_friction, workspace)
    end if
    if (physics%rayleigh_friction%enabled) then
      call add_dry_rayleigh_friction_tendency(physics%rayleigh_friction, workspace)
    end if
    if (physics%held_suarez%enabled) then
      call add_dry_held_suarez_tendency(physics%held_suarez, transform, workspace)
    end if
    if (physics%evaporation%enabled) then
      call add_dry_evaporation_tendency(physics%evaporation, physics%radiation, physics%bucket, physics%snow, &
        interval, workspace)
    end if
    if (physics%convection%enabled .or. physics%moist_convection%enabled .or. physics%condensation%enabled .or. &
        physics%cloud%enabled) then
      call add_dry_convection_tendency(physics, interval, workspace)
    end if
    if (physics%radiation%enabled) then
      call add_dry_radiation_tendency(physics%radiation, physics%sea_ice, physics%snow, interval, &
        physics%moisture%enabled, transform%mu, workspace)
    end if
    if (physics%bucket%enabled) call add_dry_bucket_tendency(physics%bucket, interval, workspace)

    if (present(diagnostics)) then
      call diagnose_surface_tiles(physics%radiation, physics%sea_ice, physics%moisture%enabled, transform%mu, workspace)
      call collect_dry_diagnostics(transform%mu, workspace, diagnostics)
    end if
    call project_dry_tendency(transform, truncation, workspace, physics%radiation%enabled, &
                              physics%moisture%enabled, rhs)
  end subroutine evaluate_dry_tendency

end module dry_tendency_evaluator
