module dry_stepper
  use iso_fortran_env, only: real64
  use harmonics, only: harmonic_transform
  use planet_parameters, only: planet_config
  use dry_vertical_coordinate, only: hybrid_sigma_coordinate
  use dry_gravity_wave, only: dry_gravity_wave_solver
  use dry_gravity_wave_operator, only: solve_dry_gravity_wave
  use dry_state, only: dry_state_type, dry_tendency_type, allocate_dry_state, allocate_dry_surface_water, &
                       allocate_dry_tendency, copy_dry_state, swap_dry_states, enforce_dry_state_constraints
  use dry_physics_config, only: dry_model_physics_config
  use dry_tendency_evaluator, only: evaluate_dry_tendency
  use dry_radiation, only: radiation_diagnostics
  use spectral_hyperdiffusion, only: apply_spectral_hyperdiffusion
  use raw_filter, only: apply_raw_filter
  use numerics_config, only: model_numerics_config, dry_hyperdiffusion_config
  use dry_tendency_workspace, only: dry_workspace_type
  implicit none
  private

  !> Semi-implicit RAW-filtered leapfrog integration of the dry atmosphere.  The
  !> work storage is allocated once and exchanged with the model state by swapping.
  type, public :: dry_stepper_type
    private
    integer :: truncation = -1
    integer :: number_of_levels = 0
    type(dry_tendency_type) :: rhs
    type(dry_state_type) :: half, candidate, next, filtered
  contains
    procedure, public :: initialize => initialize_dry_stepper
    procedure, public :: advance => advance_dry
  end type dry_stepper_type

contains

  subroutine initialize_dry_stepper(this, transform, truncation, number_of_levels)
    class(dry_stepper_type), intent(inout) :: this
    type(harmonic_transform), intent(inout) :: transform
    integer, intent(in) :: truncation, number_of_levels
    real(real64), allocatable :: grid(:, :)

    this%truncation = truncation
    this%number_of_levels = number_of_levels
    call allocate_dry_tendency(this%rhs, truncation, number_of_levels)
    call allocate_dry_state(this%half, truncation, number_of_levels)
    call allocate_dry_state(this%candidate, truncation, number_of_levels)
    call allocate_dry_state(this%next, truncation, number_of_levels)
    call allocate_dry_state(this%filtered, truncation, number_of_levels)
    call transform%allocate_field(grid)
    call allocate_dry_surface_water(this%half, size(grid, 1), size(grid, 2))
    call allocate_dry_surface_water(this%candidate, size(grid, 1), size(grid, 2))
    call allocate_dry_surface_water(this%next, size(grid, 1), size(grid, 2))
    call allocate_dry_surface_water(this%filtered, size(grid, 1), size(grid, 2))
  end subroutine initialize_dry_stepper

  !> Advances one step.  When radiation is enabled, diagnostics receives the
  !> sample evaluated at the state this step started from.
  subroutine advance_dry(this, transform, coordinate, gravity_wave, numerics, hyperdiffusion, planet, &
                         physics, workspace, surface_geopotential, step_number, previous, current, &
                         land_fraction, maximum_speed, diagnostics)
    class(dry_stepper_type), intent(inout) :: this
    type(harmonic_transform), intent(inout) :: transform
    type(hybrid_sigma_coordinate), intent(in) :: coordinate
    type(dry_gravity_wave_solver), intent(in) :: gravity_wave
    type(model_numerics_config), intent(in) :: numerics
    type(dry_hyperdiffusion_config), intent(in) :: hyperdiffusion
    type(planet_config), intent(in) :: planet
    type(dry_model_physics_config), intent(in) :: physics
    type(dry_workspace_type), intent(inout) :: workspace
    complex(real64), intent(in) :: surface_geopotential(0:, 0:)
    real(real64), intent(in) :: land_fraction(:, :)
    integer, intent(inout) :: step_number
    type(dry_state_type), intent(inout) :: previous, current
    real(real64), intent(out) :: maximum_speed
    type(radiation_diagnostics), intent(out) :: diagnostics
    real(real64) :: half_step_maximum_speed, time_step
    logical :: collect_diagnostics

    if (this%truncation /= numerics%truncation .or. this%number_of_levels /= coordinate%number_of_levels) then
      error stop 'dry stepper is not initialized for this model shape'
    end if
    time_step = numerics%time_step
    collect_diagnostics = physics%radiation%enabled

    if (step_number == 0) then
      call integration_step(transform, coordinate, gravity_wave, numerics, hyperdiffusion, planet, &
        0.25_real64*time_step, current, current, 0.0_real64, physics, workspace, surface_geopotential, &
        land_fraction, collect_diagnostics, .false., this%rhs, this%candidate, this%half, this%filtered, maximum_speed, &
        diagnostics)
      call integration_step(transform, coordinate, gravity_wave, numerics, hyperdiffusion, planet, &
        0.5_real64*time_step, current, this%half, 0.5_real64*time_step, physics, workspace, &
        surface_geopotential, land_fraction, .false., .true., this%rhs, this%candidate, this%next, this%filtered, &
        half_step_maximum_speed)
      call copy_dry_state(current, previous)
    else
      call integration_step(transform, coordinate, gravity_wave, numerics, hyperdiffusion, planet, &
        time_step, previous, current, real(step_number, real64)*time_step, physics, workspace, &
        surface_geopotential, land_fraction, collect_diagnostics, .true., this%rhs, this%candidate, this%next, &
        this%filtered, maximum_speed, diagnostics)
      call swap_dry_states(previous, this%filtered)
    end if
    call swap_dry_states(current, this%next)
    step_number = step_number + 1
  end subroutine advance_dry

  subroutine integration_step(transform, coordinate, gravity_wave, numerics, hyperdiffusion, planet, interval, &
                              previous, current, evaluation_time, physics, workspace, surface_geopotential, &
                              land_fraction, collect_diagnostics, apply_raw, rhs, candidate, next, filtered, &
                              maximum_speed, diagnostics)
    type(harmonic_transform), intent(inout) :: transform
    type(hybrid_sigma_coordinate), intent(in) :: coordinate
    type(dry_gravity_wave_solver), intent(in) :: gravity_wave
    type(model_numerics_config), intent(in) :: numerics
    type(dry_hyperdiffusion_config), intent(in) :: hyperdiffusion
    type(planet_config), intent(in) :: planet
    real(real64), intent(in) :: interval, evaluation_time
    type(dry_state_type), intent(in) :: previous, current
    type(dry_model_physics_config), intent(in) :: physics
    type(dry_workspace_type), intent(inout) :: workspace
    complex(real64), intent(in) :: surface_geopotential(0:, 0:)
    real(real64), intent(in) :: land_fraction(:, :)
    logical, intent(in) :: collect_diagnostics, apply_raw
    type(dry_tendency_type), intent(inout) :: rhs
    type(dry_state_type), intent(inout) :: candidate, next, filtered
    real(real64), intent(out) :: maximum_speed
    type(radiation_diagnostics), intent(out), optional :: diagnostics
    real(real64) :: centered_interval, order
    integer :: k, truncation

    truncation = numerics%truncation
    centered_interval = 2.0_real64*interval
    order = real(hyperdiffusion%order, real64)
    if (collect_diagnostics) then
      if (.not. present(diagnostics)) error stop 'dry stepper diagnostics output is required'
      call evaluate_dry_tendency(transform, truncation, coordinate, planet, current, previous, &
        surface_geopotential, physics, workspace, evaluation_time, centered_interval, rhs, maximum_speed, &
        diagnostics, land_fraction)
    else
      call evaluate_dry_tendency(transform, truncation, coordinate, planet, current, previous, &
        surface_geopotential, physics, workspace, evaluation_time, centered_interval, rhs, maximum_speed, &
        land_fraction=land_fraction)
    end if

    ! Vorticity, specific humidity and the surface temperatures have no gravity-wave part
    ! and are advanced explicitly; the semi-implicit solve below handles the rest.
    candidate%zeta(:, :, :) = previous%zeta + centered_interval*rhs%zeta
    candidate%specific_humidity(:, :, :) = previous%specific_humidity + centered_interval*rhs%specific_humidity
    candidate%surface_temperature(:, :) = previous%surface_temperature + centered_interval*rhs%surface_temperature
    candidate%deep_temperature(:, :) = previous%deep_temperature + centered_interval*rhs%deep_temperature
    candidate%surface_water = previous%surface_water + centered_interval*workspace%forcing_surface_water
    if (physics%bucket%enabled) then
      candidate%surface_water = min(physics%bucket%capacity, max(0.0_real64, candidate%surface_water))
    end if
    call solve_dry_gravity_wave(gravity_wave, centered_interval, previous, current, rhs, candidate)
    do k = 1, coordinate%number_of_levels
      call apply_spectral_hyperdiffusion(truncation, centered_interval, candidate%zeta(:, :, k), &
                                         hyperdiffusion%vorticity_timescale_seconds, order)
      call apply_spectral_hyperdiffusion(truncation, centered_interval, candidate%delta(:, :, k), &
                                         hyperdiffusion%divergence_timescale_seconds, order)
      call apply_spectral_hyperdiffusion(truncation, centered_interval, candidate%temperature(:, :, k), &
                                         hyperdiffusion%temperature_timescale_seconds, order)
      call apply_spectral_hyperdiffusion(truncation, centered_interval, candidate%specific_humidity(:, :, k), &
                                         hyperdiffusion%humidity_timescale_seconds, order)
    end do

    if (apply_raw) then
      do k = 1, coordinate%number_of_levels
        call apply_raw_filter(previous%zeta(:, :, k), current%zeta(:, :, k), candidate%zeta(:, :, k), &
                              filtered%zeta(:, :, k), next%zeta(:, :, k), numerics%raw_filter)
        call apply_raw_filter(previous%delta(:, :, k), current%delta(:, :, k), candidate%delta(:, :, k), &
                              filtered%delta(:, :, k), next%delta(:, :, k), numerics%raw_filter)
        call apply_raw_filter(previous%temperature(:, :, k), current%temperature(:, :, k), &
                              candidate%temperature(:, :, k), filtered%temperature(:, :, k), &
                              next%temperature(:, :, k), numerics%raw_filter)
        call apply_raw_filter(previous%specific_humidity(:, :, k), current%specific_humidity(:, :, k), &
                              candidate%specific_humidity(:, :, k), filtered%specific_humidity(:, :, k), &
                              next%specific_humidity(:, :, k), numerics%raw_filter)
      end do
      call apply_raw_filter(previous%log_surface_pressure, current%log_surface_pressure, &
        candidate%log_surface_pressure, filtered%log_surface_pressure, next%log_surface_pressure, &
        numerics%raw_filter)
      call apply_raw_filter(previous%surface_temperature, current%surface_temperature, &
        candidate%surface_temperature, filtered%surface_temperature, next%surface_temperature, &
        numerics%raw_filter)
      call apply_raw_filter(previous%deep_temperature, current%deep_temperature, &
        candidate%deep_temperature, filtered%deep_temperature, next%deep_temperature, &
        numerics%raw_filter)
      call apply_raw_filter(previous%surface_water, current%surface_water, candidate%surface_water, &
                            filtered%surface_water, next%surface_water, numerics%raw_filter)
      if (physics%bucket%enabled) then
        filtered%surface_water = min(physics%bucket%capacity, max(0.0_real64, filtered%surface_water))
        next%surface_water = min(physics%bucket%capacity, max(0.0_real64, next%surface_water))
      end if
    else
      call copy_dry_state(current, filtered)
      call copy_dry_state(candidate, next)
    end if
    call enforce_dry_state_constraints(filtered, truncation)
    call enforce_dry_state_constraints(next, truncation)
  end subroutine integration_step

end module dry_stepper
