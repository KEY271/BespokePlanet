!> Composes the explicit right-hand side of the dry atmosphere from the
!> individual dynamical and physical tendencies.
!>
!> Each tendency adds into the shared grid-space accumulators of the workspace;
!> the projection then performs the spectral transforms once.  Adding a new
!> physical process means adding one call here, not editing a kernel.
!>
!> The order of the calls is the order in which the contributions are summed.
!> It reproduces the accumulation order of the original single kernel, so the
!> split does not change the rounding of the right-hand side.
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
  use dry_radiation_tendency, only: add_dry_radiation_tendency
  use dry_convection_tendency, only: add_dry_convection_tendency
  use dry_tendency_projection, only: project_dry_tendency
  implicit none
  private
  public :: evaluate_dry_tendency

contains

  subroutine evaluate_dry_tendency(transform, truncation, coordinate, planet, state, physics_state, &
                                   surface_geopotential, physics, workspace, evaluation_time, rhs, &
                                   maximum_speed, diagnostics)
    type(harmonic_transform), intent(inout) :: transform
    integer, intent(in) :: truncation
    type(hybrid_sigma_coordinate), intent(in) :: coordinate
    type(planet_config), intent(in) :: planet
    !> state advances the dynamics; physics_state is the RAW-filtered previous
    !> time level that every prescribed physical tendency is evaluated on.
    type(dry_state_type), intent(in) :: state, physics_state
    complex(real64), intent(in) :: surface_geopotential(0:, 0:)
    type(dry_model_physics_config), intent(in) :: physics
    type(dry_workspace_type), intent(inout) :: workspace
    real(real64), intent(in) :: evaluation_time
    type(dry_tendency_type), intent(inout) :: rhs
    real(real64), intent(out) :: maximum_speed
    type(radiation_diagnostics), intent(out), optional :: diagnostics

    if (present(diagnostics) .and. .not. physics%radiation%enabled) then
      error stop 'radiation diagnostics requested without radiation forcing'
    end if

    call zero_dry_tendency(rhs)
    call workspace%zero_forcing()
    call workspace%prepare(transform, coordinate, planet%rotation_rate, state, physics_state, &
                           surface_geopotential, physics, evaluation_time)

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
    if (physics%radiation%enabled) then
      if (present(diagnostics)) then
        call add_dry_radiation_tendency(physics%radiation, transform, workspace, diagnostics)
      else
        call add_dry_radiation_tendency(physics%radiation, transform, workspace)
      end if
    end if
    if (physics%convection%enabled) then
      call add_dry_convection_tendency(physics%convection, workspace)
    end if

    call project_dry_tendency(transform, truncation, workspace, physics%radiation%enabled, rhs)
  end subroutine evaluate_dry_tendency

end module dry_tendency_evaluator
