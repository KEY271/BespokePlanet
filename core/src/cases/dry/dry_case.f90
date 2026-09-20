!> The dry-atmosphere cases without radiation: Jablonowski-Williamson (steady and
!> perturbed) and Held-Suarez.  Chooses the initial condition, the physics, the
!> run length and the snapshot schedule.
module dry_case
  use iso_fortran_env, only: real64, int64
  use numerics_config, only: model_numerics_config, dry_hyperdiffusion_config
  use dry_atmosphere, only: dry_atmosphere_solver
  use dry_case_initial_conditions, only: set_jablonowski_williamson_case_state, &
                                         held_suarez_case_physics, set_held_suarez_case_state
  use filesystem, only: make_directory
  use dry_case_output, only: write_dry_snapshot, write_dry_metadata
  use case_runtime, only: case_context, ensure_context, dt, duration, output_interval_steps, &
                          write_case_header, is_log_step, write_progress, elapsed_seconds
  implicit none
  private
  public :: run_dry_case

  real(real64), parameter :: held_suarez_duration = 200.0_real64*24.0_real64*3600.0_real64
  integer, parameter :: held_suarez_output_interval_steps = nint(5.0_real64*24.0_real64*3600.0_real64/dt)

contains

  subroutine run_dry_case(context, case_name, include_perturbation, held_suarez)
    type(case_context), intent(inout) :: context
    character(*), intent(in) :: case_name
    !> .false. keeps the balanced zonal base state (should stay steady);
    !> .true. adds the localized wind perturbation that triggers the baroclinic wave.
    logical, intent(in) :: include_perturbation
    logical, intent(in), optional :: held_suarez
    type(dry_atmosphere_solver) :: solver
    type(model_numerics_config) :: numerics
    type(dry_hyperdiffusion_config) :: hyperdiffusion
    real(real64), allocatable :: zeta(:, :, :), delta(:, :, :), temperature(:, :, :)
    real(real64), allocatable :: surface_pressure(:, :), u(:, :, :), v(:, :, :)
    real(real64), allocatable :: pressure_half(:), delta_pressure(:), layer_l(:), alpha(:), reference_temperature(:)
    real(real64), allocatable :: a_half(:), b_half(:)
    character(len=:), allocatable :: case_directory, initial_condition
    integer :: step, number_of_steps, snapshot_interval_steps
    integer(int64) :: start_count
    real(real64) :: cfl, maximum_cfl, elapsed_wall_seconds, case_duration
    logical :: use_held_suarez

    call ensure_context(context)
    numerics = context%numerics
    hyperdiffusion = dry_hyperdiffusion_config()
    call system_clock(start_count)
    use_held_suarez = .false.
    if (present(held_suarez)) use_held_suarez = held_suarez
    if (use_held_suarez) then
      initial_condition = 'Held-Suarez resting isothermal atmosphere with thermal and Rayleigh forcing'
      case_duration = held_suarez_duration
      snapshot_interval_steps = held_suarez_output_interval_steps
    else if (include_perturbation) then
      initial_condition = 'Jablonowski-Williamson with localized wind perturbation'
      case_duration = duration
      snapshot_interval_steps = output_interval_steps
    else
      initial_condition = 'Jablonowski-Williamson balanced base state without perturbation'
      case_duration = duration
      snapshot_interval_steps = output_interval_steps
    end if
    call write_case_header(case_name, case_duration, numerics%time_step)
    case_directory = context%output_root//'/'//case_name
    call make_directory(case_directory)
    call solver%init_with_config(numerics, hyperdiffusion=hyperdiffusion)
    if (use_held_suarez) then
      call set_held_suarez_case_state(solver, context%transform, held_suarez_case_physics())
    else
      call set_jablonowski_williamson_case_state(solver, context%transform, include_perturbation)
    end if
    number_of_steps = nint(case_duration/numerics%time_step)
    maximum_cfl = 0.0_real64
    do step = 0, number_of_steps
      ! Grid fields are synthesized only for snapshots.  On every other step the CFL of
      ! the current state comes from the grid winds that advance already builds.
      if (mod(step, snapshot_interval_steps) == 0 .or. step == number_of_steps) then
        call solver%get_fields(zeta, delta, temperature, surface_pressure, u, v, cfl)
        maximum_cfl = max(maximum_cfl, cfl)
        call write_dry_snapshot(case_directory, step, context%nlon, zeta, delta, temperature, &
                                surface_pressure, u, v)
      end if
      if (step < number_of_steps) then
        call solver%advance()
        cfl = solver%get_last_advance_cfl()
        maximum_cfl = max(maximum_cfl, cfl)
      end if
      ! Logged after advance, so the elapsed time already covers step + 1 completed steps.
      if (is_log_step(step, number_of_steps)) then
        call write_progress(step, number_of_steps, min(step + 1, number_of_steps), &
                            'advective CFL', cfl, start_count, numerics%time_step)
      end if
    end do
    elapsed_wall_seconds = elapsed_seconds(start_count)
    call solver%get_reference_atmosphere(pressure_half, delta_pressure, layer_l, alpha, reference_temperature, &
                                         a_half, b_half)
    call write_dry_metadata(case_directory, initial_condition, numerics%truncation, numerics%time_step, &
                            case_duration, number_of_steps, snapshot_interval_steps, &
                            maximum_cfl, elapsed_wall_seconds, context%nlon, context%transform%mu, hyperdiffusion, &
                            pressure_half, delta_pressure, layer_l, alpha, reference_temperature, a_half, b_half)
  end subroutine run_dry_case

end module dry_case
