program check_refactoring_types
  use iso_fortran_env, only: real64
  use numerics_config, only: model_numerics_config, dry_hyperdiffusion_config, &
                             default_raw_filter_epsilon, default_raw_filter_alpha, &
                             default_hyperdiffusion_order, default_hyperdiffusion_timescale_seconds
  use shallow_water_config, only: shallow_water_equation_config
  use dry_physics_config, only: dry_model_physics_config
  use barotropic_state, only: barotropic_state_type, barotropic_tendency_type, &
                              allocate_barotropic_state, allocate_barotropic_tendency, &
                              copy_barotropic_state, swap_barotropic_states, enforce_barotropic_constraints
  use shallow_water_state, only: shallow_water_state_type, allocate_shallow_water_state, &
                                 enforce_shallow_water_constraints
  use dry_state, only: dry_state_type, dry_tendency_type, allocate_dry_state, &
                       allocate_dry_tendency, copy_dry_state, swap_dry_states, enforce_dry_state_constraints
  implicit none

  type(model_numerics_config) :: numerics
  type(dry_hyperdiffusion_config) :: dry_hyperdiffusion
  type(shallow_water_equation_config) :: equation
  type(dry_model_physics_config) :: physics
  type(barotropic_state_type) :: barotropic, barotropic_other
  type(barotropic_tendency_type) :: barotropic_rhs
  type(shallow_water_state_type) :: shallow
  type(dry_state_type) :: dry, dry_other
  type(dry_tendency_type) :: dry_rhs
  integer, parameter :: truncation = 4, levels = 3

  numerics%truncation = truncation
  numerics%time_step = 300.0_real64
  if (numerics%raw_filter%epsilon /= default_raw_filter_epsilon .or. &
      numerics%raw_filter%alpha /= default_raw_filter_alpha .or. &
      numerics%hyperdiffusion%order /= default_hyperdiffusion_order .or. &
      numerics%hyperdiffusion%timescale_seconds /= default_hyperdiffusion_timescale_seconds) then
    error stop 'numerics defaults are not taken from numerics_config'
  end if
  if (dry_hyperdiffusion%vorticity_timescale_seconds /= 4.0_real64*3600.0_real64 .or. &
      dry_hyperdiffusion%divergence_timescale_seconds /= 1.0_real64*3600.0_real64 .or. &
      dry_hyperdiffusion%temperature_timescale_seconds /= 4.0_real64*3600.0_real64 .or. &
      dry_hyperdiffusion%order /= default_hyperdiffusion_order) then
    error stop 'dry hyperdiffusion defaults changed'
  end if
  if (equation%gravity /= 9.8_real64 .or. equation%mean_depth /= 1.0e4_real64 .or. &
      equation%gravity_wave_implicitness /= 0.5_real64) then
    error stop 'shallow-water equation defaults changed'
  end if
  if (physics%held_suarez%enabled .or. physics%surface_friction%enabled .or. physics%rayleigh_friction%enabled .or. &
      physics%radiation%enabled .or. physics%convection%enabled) then
    error stop 'dry physics must default to disabled'
  end if

  call allocate_barotropic_state(barotropic, truncation)
  call allocate_barotropic_state(barotropic_other, truncation)
  call allocate_barotropic_tendency(barotropic_rhs, truncation)
  barotropic%zeta = cmplx(1.0_real64, 0.0_real64, kind=real64)
  call enforce_barotropic_constraints(barotropic%zeta, truncation)
  if (barotropic%zeta(0, 0) /= 0.0_real64 .or. any(barotropic%zeta(truncation + 1, :) /= 0.0_real64)) then
    error stop 'barotropic constraints failed'
  end if
  if (any(barotropic_rhs%zeta /= 0.0_real64)) error stop 'barotropic tendency did not initialize to zero'
  call copy_barotropic_state(barotropic, barotropic_other)
  if (any(barotropic_other%zeta /= barotropic%zeta)) error stop 'barotropic state copy failed'
  barotropic_other%zeta = cmplx(2.0_real64, 0.0_real64, kind=real64)
  call swap_barotropic_states(barotropic, barotropic_other)
  if (any(barotropic%zeta /= cmplx(2.0_real64, 0.0_real64, kind=real64)) .or. &
      barotropic_other%zeta(1, 1) /= cmplx(1.0_real64, 0.0_real64, kind=real64) .or. &
      lbound(barotropic%zeta, 1) /= 0 .or. lbound(barotropic_other%zeta, 2) /= 0) then
    error stop 'barotropic state swap failed'
  end if

  call allocate_shallow_water_state(shallow, truncation)
  shallow%zeta = cmplx(1.0_real64, 0.0_real64, kind=real64)
  shallow%delta = shallow%zeta
  shallow%eta = shallow%zeta
  call enforce_shallow_water_constraints(shallow, truncation)
  if (shallow%zeta(0, 0) /= 0.0_real64 .or. shallow%delta(0, 0) /= 0.0_real64 .or. &
      shallow%eta(0, 0) /= 0.0_real64) error stop 'shallow-water constraints failed'

  call allocate_dry_state(dry, truncation, levels)
  call allocate_dry_state(dry_other, truncation, levels)
  call allocate_dry_tendency(dry_rhs, truncation, levels)
  dry%zeta = cmplx(1.0_real64, 0.0_real64, kind=real64)
  dry%delta = dry%zeta
  dry%temperature = dry%zeta
  dry%log_surface_pressure = dry%zeta(:, :, 1)
  dry%deep_temperature = cmplx(3.0_real64, 0.0_real64, kind=real64)
  call enforce_dry_state_constraints(dry, truncation)
  if (any(dry%zeta(0, 0, :) /= 0.0_real64) .or. any(dry%delta(0, 0, :) /= 0.0_real64)) then
    error stop 'dry zero-mean constraints failed'
  end if
  if (any(dry_rhs%temperature /= 0.0_real64)) error stop 'dry tendency did not initialize to zero'
  call copy_dry_state(dry, dry_other)
  if (any(dry_other%temperature /= dry%temperature) .or. any(dry_other%deep_temperature /= dry%deep_temperature)) then
    error stop 'dry state copy failed'
  end if
  dry_other%deep_temperature = cmplx(5.0_real64, 0.0_real64, kind=real64)
  call swap_dry_states(dry, dry_other)
  if (dry%deep_temperature(0, 0) /= cmplx(5.0_real64, 0.0_real64, kind=real64) .or. &
      dry_other%deep_temperature(0, 0) /= cmplx(3.0_real64, 0.0_real64, kind=real64)) then
    error stop 'dry state swap failed'
  end if

end program check_refactoring_types
