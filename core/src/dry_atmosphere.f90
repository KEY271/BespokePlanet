module dry_atmosphere
  use iso_fortran_env, only: real64
  use harmonics, only: harmonic_transform
  use planet_parameters, only: earth_radius, planet_config
  use spectral_vector_operators, only: diagnose_horizontal_velocity
  use dry_vertical_coordinate, only: hybrid_sigma_coordinate
  use dry_gravity_wave, only: dry_gravity_wave_solver
  use dry_radiation, only: radiation_diagnostics, move_radiation_diagnostics
  use dry_tendency_evaluator, only: evaluate_dry_tendency
  use dry_state, only: dry_tendency_type, allocate_dry_tendency, dry_state_type, allocate_dry_state, allocate_dry_surface_fields, copy_dry_state, &
                       enforce_dry_state_constraints, enforce_dry_spectral_field
  use dry_physics_config, only: dry_model_physics_config, dry_reference_temperature
  use sea_ice, only: equilibrate_sea_ice, validate_sea_ice_config
  use ocean_q_flux, only: q_flux_diagnostics, build_ocean_q_flux, validate_q_flux_config
  use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
  use dry_stepper, only: dry_stepper_type
  use numerics_config, only: model_numerics_config, dry_hyperdiffusion_config
  use dry_tendency_workspace, only: dry_workspace_type
  implicit none
  private

  real(real64), parameter, public :: dry_gravity_wave_reference_temperature = dry_reference_temperature

  type, public :: dry_atmosphere_solver
    private
    type(harmonic_transform) :: transform
    type(hybrid_sigma_coordinate) :: coordinate
    type(dry_gravity_wave_solver) :: gravity_wave
    integer :: truncation = -1
    integer :: number_of_levels = 0
    integer :: step_number = -1
    type(model_numerics_config) :: numerics
    type(dry_hyperdiffusion_config) :: hyperdiffusion
    !> Planet the dynamics rotate with; a case may override it (set_planet).
    type(planet_config) :: planet
    !> Active physical processes and their coefficients; the only record of which forcing is enabled.
    type(dry_model_physics_config) :: physics
    type(dry_stepper_type) :: stepper
    type(dry_workspace_type) :: workspace
    type(radiation_diagnostics) :: latest_radiation_diagnostics
    type(q_flux_diagnostics) :: q_flux_diagnostics
    logical :: latest_diagnostics_available = .false.
    real(real64), allocatable :: gravity_wave_reference_temperature(:)
    !> Advective CFL of the state that the most recent advance started from.
    real(real64) :: last_advance_cfl = 0.0_real64
    type(dry_state_type) :: previous
    type(dry_state_type) :: current
    ! Time-independent lower boundary condition; it is never advanced or filtered.
    complex(real64), allocatable :: surface_geopotential(:, :)
    !> Time-independent grid land fraction used only by surface physics.
    real(real64), allocatable :: land_fraction(:, :)
  contains
    procedure, public :: init => initialize_dry_solver
    procedure, public :: init_with_config => initialize_dry_solver_with_config
    procedure, public :: set_initial_state
    !> Sets the grid surface and deep ground temperature on both time levels.
    !> Which initial condition they come from is a case decision, not a model one.
    procedure, public :: set_surface_state
    !> Selects the active physical processes; set_initial_state disables them all.
    procedure, public :: set_physics
    !> Selects the planet parameters the dynamics use.
    procedure, public :: set_planet
    procedure, public :: advance
    procedure, public :: get_fields
    procedure, public :: get_spectral_state
    procedure, public :: get_reference_atmosphere
    procedure, public :: get_truncation
    procedure, public :: get_coordinate
    procedure, public :: get_step
    procedure, public :: get_time
    procedure, public :: get_last_advance_cfl
    procedure, public :: take_latest_diagnostics
    !> Effective column cloud cover of the most recent tendency evaluation
    !> (docs/tendency/cloud.md); zero before the first advance or without clouds.
    procedure, public :: get_cloud_cover
    !> Prescribed Q flux per ocean area and its diagnostics (docs/tendency/q-flux.md).
    procedure, public :: get_ocean_q_flux
  end type dry_atmosphere_solver

contains

  subroutine initialize_dry_solver(this, truncation, dt, a_half, b_half)
    class(dry_atmosphere_solver), intent(inout) :: this
    integer, intent(in) :: truncation
    real(real64), intent(in) :: dt
    real(real64), intent(in), optional :: a_half(0:), b_half(0:)

    if (truncation < 1) error stop 'dry atmosphere truncation must be positive'
    if (dt <= 0.0_real64) error stop 'dry atmosphere time step must be positive'
    if (present(a_half) .neqv. present(b_half)) then
      error stop 'dry atmosphere custom A and B must be supplied together'
    end if
    this%truncation = truncation
    this%numerics = model_numerics_config()
    this%numerics%truncation = truncation
    this%numerics%time_step = dt
    this%hyperdiffusion = dry_hyperdiffusion_config()
    this%planet = planet_config()
    this%step_number = -1
    this%last_advance_cfl = 0.0_real64
    this%physics = dry_model_physics_config()
    this%latest_diagnostics_available = .false.
    call this%transform%init(truncation)
    if (present(a_half)) then
      call this%coordinate%init(a_half, b_half)
    else
      call this%coordinate%init_default()
    end if
    this%number_of_levels = this%coordinate%number_of_levels
    call this%workspace%initialize(this%transform, truncation, this%number_of_levels)
    allocate (this%gravity_wave_reference_temperature(this%number_of_levels))
    this%gravity_wave_reference_temperature = dry_gravity_wave_reference_temperature
    call this%gravity_wave%init(truncation, dt, this%coordinate, this%gravity_wave_reference_temperature)
    call allocate_dry_state(this%previous, truncation, this%number_of_levels)
    call allocate_dry_state(this%current, truncation, this%number_of_levels)
    call allocate_dry_surface_fields(this%previous, this%workspace%nx, this%workspace%ny)
    call allocate_dry_surface_fields(this%current, this%workspace%nx, this%workspace%ny)
    call this%stepper%initialize(this%transform, truncation, this%number_of_levels)
    if (allocated(this%surface_geopotential)) deallocate (this%surface_geopotential)
    allocate (this%surface_geopotential(0:truncation + 1, 0:truncation))
    this%surface_geopotential = 0.0_real64
    if (allocated(this%land_fraction)) deallocate (this%land_fraction)
    call this%transform%allocate_field(this%land_fraction)
    this%land_fraction = 0.0_real64
  end subroutine initialize_dry_solver

  subroutine initialize_dry_solver_with_config(this, config, a_half, b_half, hyperdiffusion, planet)
    class(dry_atmosphere_solver), intent(inout) :: this
    type(model_numerics_config), intent(in) :: config
    real(real64), intent(in), optional :: a_half(0:), b_half(0:)
    type(dry_hyperdiffusion_config), intent(in), optional :: hyperdiffusion
    type(planet_config), intent(in), optional :: planet
    if (present(a_half) .neqv. present(b_half)) then
      error stop 'dry atmosphere custom A and B must be supplied together'
    end if
    if (present(a_half)) then
      call this%init(config%truncation, config%time_step, a_half, b_half)
    else
      call this%init(config%truncation, config%time_step)
    end if
    this%numerics = config
    if (present(hyperdiffusion)) this%hyperdiffusion = hyperdiffusion
    if (present(planet)) this%planet = planet
  end subroutine initialize_dry_solver_with_config

  !> Installs the spectral state on both time levels.  Specific humidity starts
  !> at zero (a dry atmosphere) unless supplied.
  subroutine set_initial_state(this, zeta, delta, temperature, log_surface_pressure, surface_geopotential, &
                               specific_humidity, land_fraction)
    class(dry_atmosphere_solver), intent(inout) :: this
    complex(real64), intent(in) :: zeta(0:, 0:, :), delta(0:, 0:, :), temperature(0:, 0:, :)
    complex(real64), intent(in) :: log_surface_pressure(0:, 0:)
    !> Fixed surface geopotential; a flat surface (zero) is used when absent.
    complex(real64), intent(in), optional :: surface_geopotential(0:, 0:)
    complex(real64), intent(in), optional :: specific_humidity(0:, 0:, :)
    real(real64), intent(in), optional :: land_fraction(:, :)

    call check_initialized(this)
    call check_state_shape(this, zeta, 'zeta')
    call check_state_shape(this, delta, 'delta')
    call check_state_shape(this, temperature, 'temperature')
    if (ubound(log_surface_pressure, 1) < this%truncation + 1 .or. &
        ubound(log_surface_pressure, 2) < this%truncation) then
      error stop 'dry atmosphere initial log surface pressure has an inconsistent shape'
    end if
    this%current%zeta = zeta(0:this%truncation + 1, 0:this%truncation, 1:this%number_of_levels)
    this%current%delta = delta(0:this%truncation + 1, 0:this%truncation, 1:this%number_of_levels)
    this%current%temperature = temperature(0:this%truncation + 1, 0:this%truncation, 1:this%number_of_levels)
    this%current%log_surface_pressure = log_surface_pressure(0:this%truncation + 1, 0:this%truncation)
    this%current%specific_humidity = 0.0_real64
    if (present(specific_humidity)) then
      call check_state_shape(this, specific_humidity, 'specific humidity')
      this%current%specific_humidity = &
        specific_humidity(0:this%truncation + 1, 0:this%truncation, 1:this%number_of_levels)
    end if
    this%surface_geopotential = 0.0_real64
    if (present(surface_geopotential)) then
      if (ubound(surface_geopotential, 1) < this%truncation + 1 .or. &
          ubound(surface_geopotential, 2) < this%truncation) then
        error stop 'dry atmosphere surface geopotential has an inconsistent shape'
      end if
      this%surface_geopotential = surface_geopotential(0:this%truncation + 1, 0:this%truncation)
    end if
    call enforce_dry_spectral_field(this%surface_geopotential, this%truncation, .false.)
    this%land_fraction = 0.0_real64
    if (present(land_fraction)) then
      if (any(shape(land_fraction) /= shape(this%land_fraction))) then
        error stop 'dry atmosphere land fraction has an inconsistent shape'
      end if
      if (any(land_fraction < 0.0_real64) .or. any(land_fraction > 1.0_real64)) then
        error stop 'dry atmosphere land fraction is outside [0,1]'
      end if
      this%land_fraction = land_fraction
    end if
    this%current%surface_temperature = 0.0_real64
    this%current%deep_temperature = 0.0_real64
    this%current%surface_water = 0.0_real64
    this%current%land_temperature = 0.0_real64
    this%current%ocean_temperature = 0.0_real64
    this%current%sea_ice_fraction = 0.0_real64
    this%current%sea_ice_volume = 0.0_real64
    call enforce_dry_state_constraints(this%current, this%truncation)
    call copy_dry_state(this%current, this%previous)
    ! A plain initial state runs without forcing; the case enables the processes it needs.
    this%physics = dry_model_physics_config()
    this%workspace%ocean_q_flux = 0.0_real64
    this%q_flux_diagnostics = q_flux_diagnostics()
    this%latest_diagnostics_available = .false.
    this%step_number = 0
  end subroutine set_initial_state

  !> Overrides the grid surface and deep ground temperature on both time levels.
  !> set_initial_state leaves them at zero; cases that carry a ground energy
  !> budget supply the initial values here.  The fields have the padded grid
  !> shape of the transform (transform%allocate_field).
  subroutine set_surface_state(this, surface_temperature, deep_temperature, surface_water, &
                               land_temperature, ocean_temperature, sea_ice_fraction, sea_ice_volume)
    class(dry_atmosphere_solver), intent(inout) :: this
    real(real64), intent(in) :: surface_temperature(:, :), deep_temperature(:, :)
    real(real64), intent(in), optional :: surface_water(:, :)
    real(real64), intent(in), optional :: land_temperature(:, :)
    real(real64), intent(in), optional :: ocean_temperature(:, :)
    real(real64), intent(in), optional :: sea_ice_fraction(:, :)
    real(real64), intent(in), optional :: sea_ice_volume(:, :)

    call check_initialized(this)
    if (any(shape(surface_temperature) /= shape(this%current%surface_temperature)) .or. &
        any(shape(deep_temperature) /= shape(this%current%deep_temperature))) then
      error stop 'dry atmosphere surface state has an inconsistent shape'
    end if
    this%current%surface_temperature = surface_temperature
    this%current%deep_temperature = deep_temperature
    this%current%land_temperature = surface_temperature
    this%current%ocean_temperature = surface_temperature
    this%current%sea_ice_fraction = 0.0_real64
    this%current%sea_ice_volume = 0.0_real64
    if (present(land_temperature)) then
      if (any(shape(land_temperature) /= shape(surface_temperature))) error stop 'land_temperature shape mismatch'
      this%current%land_temperature = land_temperature
    end if
    if (present(ocean_temperature)) then
      if (any(shape(ocean_temperature) /= shape(surface_temperature))) error stop 'ocean_temperature shape mismatch'
      this%current%ocean_temperature = ocean_temperature
    end if
    if (present(sea_ice_fraction)) then
      if (any(shape(sea_ice_fraction) /= shape(surface_temperature))) error stop 'sea_ice_fraction shape mismatch'
      this%current%sea_ice_fraction = sea_ice_fraction
    end if
    if (present(sea_ice_volume)) then
      if (any(shape(sea_ice_volume) /= shape(surface_temperature))) error stop 'sea_ice_volume shape mismatch'
      this%current%sea_ice_volume = sea_ice_volume
    end if
    if (present(sea_ice_fraction) .neqv. present(sea_ice_volume)) error stop 'supply ice area and volume together'
    if (present(sea_ice_fraction)) call validate_initial_ice(this)
    if (this%physics%sea_ice%enabled) call normalize_initial_ice(this)
    this%previous%land_temperature = this%current%land_temperature
    this%previous%ocean_temperature = this%current%ocean_temperature
    this%previous%sea_ice_fraction = this%current%sea_ice_fraction
    this%previous%sea_ice_volume = this%current%sea_ice_volume
    this%previous%surface_temperature = this%current%surface_temperature
    this%previous%deep_temperature = this%current%deep_temperature
    if (present(surface_water)) then
      if (any(shape(surface_water) /= shape(this%current%surface_water))) then
        error stop 'dry atmosphere surface water has an inconsistent shape'
      end if
      this%current%surface_water = surface_water
      this%previous%surface_water = surface_water
    end if
  end subroutine set_surface_state

  !> Selects which physical processes are active.  set_initial_state disables
  !> them all, so this is called after the initial state is installed.
  subroutine set_physics(this, physics)
    class(dry_atmosphere_solver), intent(inout) :: this
    type(dry_model_physics_config), intent(in) :: physics

    call check_initialized(this)
    this%physics = physics
    if (physics%sea_ice%enabled) then
      call validate_sea_ice_config(physics%sea_ice)
      if (.not. physics%radiation%enabled .or. &
          .not. (physics%radiation%slab_ocean_enabled .or. physics%radiation%land_sea_mixing_enabled)) &
        error stop 'sea ice requires a radiative ocean surface'
      call normalize_initial_ice(this)
    end if
    this%workspace%ocean_q_flux = 0.0_real64
    this%q_flux_diagnostics = q_flux_diagnostics()
    if (physics%q_flux%enabled) then
      call validate_q_flux_config(physics%q_flux)
      if (.not. physics%radiation%enabled .or. &
          .not. (physics%sea_ice%enabled .or. physics%radiation%land_sea_mixing_enabled)) &
        error stop 'Q flux requires the tiled ocean surface (sea ice or land-sea mixing)'
      ! The land fraction is fixed by set_initial_state, which precedes set_physics.
      call build_ocean_q_flux(physics%q_flux, earth_radius, this%transform%mu, this%workspace%gaussian_weights, &
        this%workspace%ring_nlon, this%land_fraction, this%workspace%ocean_q_flux, this%q_flux_diagnostics)
    end if
  end subroutine set_physics

  subroutine get_ocean_q_flux(this, q_flux, diagnostics)
    class(dry_atmosphere_solver), intent(in) :: this
    real(real64), allocatable, intent(out) :: q_flux(:, :)
    type(q_flux_diagnostics), intent(out), optional :: diagnostics
    call check_initialized(this)
    q_flux = this%workspace%ocean_q_flux
    if (present(diagnostics)) diagnostics = this%q_flux_diagnostics
  end subroutine get_ocean_q_flux

  subroutine validate_initial_ice(this)
    class(dry_atmosphere_solver), intent(in) :: this
    integer :: i, j
    real(real64) :: a, v, t
    do j = 1, this%workspace%ny
      do i = 1, this%workspace%ring_nlon(j)
        a = this%current%sea_ice_fraction(i, j)
        v = this%current%sea_ice_volume(i, j)
        t = this%current%ocean_temperature(i, j)
        if (.not. all(ieee_is_finite([a, v, t]))) error stop 'nonfinite initial sea ice'
        if (a < 0.0_real64 .or. a > 1.0_real64 .or. v < 0.0_real64) error stop 'invalid initial sea ice'
        if (.not. this%physics%sea_ice%enabled .and. a > 0.0_real64) error stop 'initial ice requires sea ice enabled'
        if ((a == 0.0_real64) .neqv. (v == 0.0_real64)) error stop 'initial sea-ice area/volume mismatch'
        if (a > 0.0_real64 .and. this%land_fraction(i, j) == 1.0_real64) error stop 'sea ice on land'
        if (a > 0.0_real64 .and. t /= this%physics%sea_ice%freezing_temperature) &
          error stop 'initial ice must be in equilibrium with the ocean'
        if (this%physics%sea_ice%enabled .and. this%land_fraction(i, j) < 1.0_real64 .and. &
            t < this%physics%sea_ice%freezing_temperature) error stop 'explicit initial ocean must not be supercooled'
      end do
    end do
  end subroutine validate_initial_ice

  subroutine normalize_initial_ice(this)
    class(dry_atmosphere_solver), intent(inout) :: this
    integer :: i, j
    real(real64) :: capacity
    capacity = this%physics%radiation%seawater_density*this%physics%radiation%seawater_specific_heat* &
      this%physics%radiation%slab_ocean_depth
    do j = 1, this%workspace%ny
      do i = 1, this%workspace%ring_nlon(j)
        if (this%land_fraction(i, j) == 1.0_real64) cycle
        ! set_physics may be called before the surface initial condition.
        if (this%current%ocean_temperature(i, j) <= 0.0_real64) cycle
        call equilibrate_sea_ice(this%physics%sea_ice, capacity, this%current%ocean_temperature(i, j), &
          this%current%sea_ice_fraction(i, j), this%current%sea_ice_volume(i, j))
      end do
    end do
    this%previous%ocean_temperature = this%current%ocean_temperature
    this%previous%sea_ice_fraction = this%current%sea_ice_fraction
    this%previous%sea_ice_volume = this%current%sea_ice_volume
  end subroutine normalize_initial_ice

  subroutine set_planet(this, planet)
    class(dry_atmosphere_solver), intent(inout) :: this
    type(planet_config), intent(in) :: planet

    call check_initialized(this)
    this%planet = planet
  end subroutine set_planet

  subroutine advance(this)
    class(dry_atmosphere_solver), intent(inout) :: this
    real(real64) :: maximum_speed

    call check_ready(this)
    call this%stepper%advance(this%transform, this%coordinate, this%gravity_wave, this%numerics, &
      this%hyperdiffusion, this%planet, this%physics, this%workspace, this%surface_geopotential, this%step_number, &
      this%previous, this%current, this%land_fraction, maximum_speed, this%latest_radiation_diagnostics)
    this%last_advance_cfl = advective_cfl(this, maximum_speed)
    this%latest_diagnostics_available = this%physics%radiation%enabled
  end subroutine advance

  subroutine get_fields(this, zeta, delta, temperature, surface_pressure, u, v, cfl, &
                        surface_temperature, deep_temperature, specific_humidity, surface_water, &
                        land_temperature, ocean_temperature, sea_ice_fraction, sea_ice_volume, &
                        sea_ice_temperature, sea_ice_thickness)
    class(dry_atmosphere_solver), intent(inout) :: this
    real(real64), allocatable, intent(out) :: zeta(:, :, :), delta(:, :, :), temperature(:, :, :)
    real(real64), allocatable, intent(out) :: surface_pressure(:, :), u(:, :, :), v(:, :, :)
    real(real64), intent(out), optional :: cfl
    real(real64), allocatable, intent(out), optional :: surface_temperature(:, :), deep_temperature(:, :)
    !> Signed grid specific humidity (it may be slightly negative where the spectral truncation overshoots).
    real(real64), allocatable, intent(out), optional :: specific_humidity(:, :, :)
    real(real64), allocatable, intent(out), optional :: surface_water(:, :)
    real(real64), allocatable, intent(out), optional :: land_temperature(:, :)
    real(real64), allocatable, intent(out), optional :: ocean_temperature(:, :)
    real(real64), allocatable, intent(out), optional :: sea_ice_fraction(:, :)
    real(real64), allocatable, intent(out), optional :: sea_ice_volume(:, :)
    real(real64), allocatable, intent(out), optional :: sea_ice_temperature(:, :)
    real(real64), allocatable, intent(out), optional :: sea_ice_thickness(:, :)
    type(radiation_diagnostics) :: surface_sample
    type(dry_tendency_type) :: snapshot_rhs
    real(real64), allocatable :: temporary(:, :), temporary_u(:, :), temporary_v(:, :), log_ps(:, :)
    integer, allocatable :: nlon(:)
    integer :: nx, ny, j, k
    real(real64) :: maximum_speed

    call check_ready(this)
    call this%transform%allocate_field(temporary)
    nx = size(temporary, 1)
    ny = size(temporary, 2)
    allocate (zeta(nx, ny, this%number_of_levels), delta(nx, ny, this%number_of_levels))
    allocate (temperature(nx, ny, this%number_of_levels))
    allocate (u(nx, ny, this%number_of_levels), v(nx, ny, this%number_of_levels))
    do k = 1, this%number_of_levels
      call this%transform%spectral_to_grid(this%current%zeta(:, :, k), temporary)
      zeta(:, :, k) = temporary
      call this%transform%spectral_to_grid(this%current%delta(:, :, k), temporary)
      delta(:, :, k) = temporary
      call this%transform%spectral_to_grid(this%current%temperature(:, :, k), temporary)
      temperature(:, :, k) = temporary
      call diagnose_horizontal_velocity(this%transform, this%truncation, &
        this%current%zeta(:, :, k), this%current%delta(:, :, k), temporary_u, temporary_v)
      u(:, :, k) = temporary_u
      v(:, :, k) = temporary_v
      if (present(specific_humidity)) then
        if (k == 1) allocate (specific_humidity(nx, ny, this%number_of_levels))
        call this%transform%spectral_to_grid(this%current%specific_humidity(:, :, k), temporary)
        specific_humidity(:, :, k) = temporary
      end if
    end do
    call this%transform%spectral_to_grid(this%current%log_surface_pressure, log_ps)
    surface_pressure = exp(log_ps)
    ! The ground temperatures are grid prognostic fields of every dry run; they stay
    ! at their initial values unless a process drives them.
    if (present(surface_temperature)) surface_temperature = this%current%surface_temperature
    if (this%physics%radiation%land_sea_mixing_enabled .or. this%physics%sea_ice%enabled) then
      if (present(surface_temperature) .or. present(land_temperature) .or. present(ocean_temperature) .or. &
          present(sea_ice_fraction) .or. present(sea_ice_volume) .or. present(sea_ice_temperature) .or. &
          present(sea_ice_thickness)) then
        call allocate_dry_tendency(snapshot_rhs, this%truncation, this%number_of_levels, nx, ny)
        call evaluate_dry_tendency(this%transform, this%truncation, this%coordinate, this%planet, &
          this%current, this%current, this%surface_geopotential, this%physics, this%workspace, &
          this%get_time(), this%numerics%time_step, snapshot_rhs, maximum_speed, surface_sample, this%land_fraction)
        if (present(surface_temperature)) surface_temperature = surface_sample%surface_temperature
        if (present(land_temperature)) land_temperature = surface_sample%land_temperature
        if (present(ocean_temperature)) ocean_temperature = surface_sample%ocean_temperature
        if (present(sea_ice_fraction)) sea_ice_fraction = surface_sample%sea_ice_fraction
        if (present(sea_ice_volume)) sea_ice_volume = surface_sample%sea_ice_volume
        if (present(sea_ice_temperature)) sea_ice_temperature = surface_sample%sea_ice_temperature
        if (present(sea_ice_thickness)) sea_ice_thickness = surface_sample%sea_ice_thickness
      end if
    else
      if (present(land_temperature)) then
        allocate (land_temperature(nx, ny))
        land_temperature = 0.0_real64
      end if
      if (present(ocean_temperature)) then
        allocate (ocean_temperature(nx, ny))
        ocean_temperature = 0.0_real64
      end if
      if (present(sea_ice_fraction)) then
        allocate (sea_ice_fraction(nx, ny))
        sea_ice_fraction = 0.0_real64
      end if
      if (present(sea_ice_volume)) then
        allocate (sea_ice_volume(nx, ny))
        sea_ice_volume = 0.0_real64
      end if
      if (present(sea_ice_temperature)) then
        allocate (sea_ice_temperature(nx, ny))
        sea_ice_temperature = 0.0_real64
      end if
      if (present(sea_ice_thickness)) then
        allocate (sea_ice_thickness(nx, ny))
        sea_ice_thickness = 0.0_real64
      end if
      if (this%physics%radiation%slab_ocean_enabled) then
        if (present(ocean_temperature)) ocean_temperature = this%current%surface_temperature
      else
        if (present(land_temperature)) land_temperature = this%current%surface_temperature
      end if
    end if
    if (present(deep_temperature)) deep_temperature = this%current%deep_temperature
    if (present(surface_water)) surface_water = this%current%surface_water
    if (present(cfl)) then
      nlon = this%transform%get_nlon()
      maximum_speed = 0.0_real64
      do k = 1, this%number_of_levels
        do j = 1, size(nlon)
          maximum_speed = max(maximum_speed, maxval(sqrt( &
            u(1:nlon(j), j, k)**2 + v(1:nlon(j), j, k)**2)))
        end do
      end do
      cfl = advective_cfl(this, maximum_speed)
    end if
  end subroutine get_fields

  !> The spectral prognostic fields.  The surface temperatures are grid fields and
  !> are read with get_fields.
  subroutine get_spectral_state(this, zeta, delta, temperature, log_surface_pressure, specific_humidity)
    class(dry_atmosphere_solver), intent(in) :: this
    complex(real64), allocatable, intent(out) :: zeta(:, :, :), delta(:, :, :), temperature(:, :, :)
    complex(real64), allocatable, intent(out) :: log_surface_pressure(:, :)
    complex(real64), allocatable, intent(out), optional :: specific_humidity(:, :, :)

    call check_ready(this)
    zeta = this%current%zeta
    delta = this%current%delta
    temperature = this%current%temperature
    log_surface_pressure = this%current%log_surface_pressure
    if (present(specific_humidity)) specific_humidity = this%current%specific_humidity
  end subroutine get_spectral_state

  !> Model time of the current state, in seconds since the initial state.
  real(real64) function get_time(this) result(time_seconds)
    class(dry_atmosphere_solver), intent(in) :: this
    time_seconds = real(max(this%step_number, 0), real64)*this%numerics%time_step
  end function get_time

  subroutine get_reference_atmosphere(this, pressure_half, delta_pressure, layer_l, alpha, temperature, &
                                      a_half, b_half)
    class(dry_atmosphere_solver), intent(in) :: this
    real(real64), allocatable, intent(out) :: pressure_half(:), delta_pressure(:), layer_l(:), alpha(:), temperature(:)
    real(real64), allocatable, intent(out), optional :: a_half(:), b_half(:)

    call check_initialized(this)
    pressure_half = this%coordinate%reference_p_half
    delta_pressure = this%coordinate%reference_delta_p
    layer_l = this%coordinate%reference_l
    alpha = this%coordinate%reference_alpha
    temperature = this%gravity_wave_reference_temperature
    if (present(a_half)) a_half = this%coordinate%a_half
    if (present(b_half)) b_half = this%coordinate%b_half
  end subroutine get_reference_atmosphere

  integer function get_truncation(this) result(truncation)
    class(dry_atmosphere_solver), intent(in) :: this

    call check_initialized(this)
    truncation = this%truncation
  end function get_truncation

  !> The vertical coordinate the solver was initialized with, so that cases can
  !> build an initial condition on exactly the levels the model integrates.
  function get_coordinate(this) result(coordinate)
    class(dry_atmosphere_solver), intent(in) :: this
    type(hybrid_sigma_coordinate) :: coordinate

    call check_initialized(this)
    coordinate = this%coordinate
  end function get_coordinate

  integer function get_step(this) result(step)
    class(dry_atmosphere_solver), intent(in) :: this
    step = this%step_number
  end function get_step

  !> Advective CFL of the state the most recent advance started from.  It is a
  !> by-product of the tendency evaluation, so no extra transforms are needed.
  real(real64) function get_last_advance_cfl(this) result(cfl)
    class(dry_atmosphere_solver), intent(in) :: this
    if (this%step_number < 1) error stop 'dry atmosphere solver: advance before requesting its CFL'
    cfl = this%last_advance_cfl
  end function get_last_advance_cfl

  !> Hands over the radiation sample of the most recent advance without copying its
  !> grid arrays.  Each sample can be taken once.
  subroutine take_latest_diagnostics(this, diagnostics)
    class(dry_atmosphere_solver), intent(inout) :: this
    type(radiation_diagnostics), intent(out) :: diagnostics
    if (.not. this%latest_diagnostics_available) then
      error stop 'dry atmosphere solver: no untaken radiation diagnostic sample is available'
    end if
    call move_radiation_diagnostics(this%latest_radiation_diagnostics, diagnostics)
    this%latest_diagnostics_available = .false.
  end subroutine take_latest_diagnostics

  subroutine get_cloud_cover(this, cloud_cover)
    class(dry_atmosphere_solver), intent(in) :: this
    real(real64), allocatable, intent(out) :: cloud_cover(:, :)
    call check_ready(this)
    cloud_cover = this%workspace%cloud_cover
  end subroutine get_cloud_cover

  pure real(real64) function advective_cfl(this, maximum_speed) result(cfl)
    class(dry_atmosphere_solver), intent(in) :: this
    real(real64), intent(in) :: maximum_speed
    cfl = maximum_speed*this%numerics%time_step/earth_radius* &
          sqrt(real(this%truncation*(this%truncation + 1), real64))
  end function advective_cfl

  subroutine check_initialized(this)
    class(dry_atmosphere_solver), intent(in) :: this
    if (this%truncation < 1 .or. .not. allocated(this%current%zeta)) then
      error stop 'dry atmosphere solver: call init before use'
    end if
  end subroutine check_initialized

  subroutine check_ready(this)
    class(dry_atmosphere_solver), intent(in) :: this
    call check_initialized(this)
    if (this%step_number < 0) error stop 'dry atmosphere solver: set initial state before integration'
  end subroutine check_ready

  subroutine check_state_shape(this, field, name)
    class(dry_atmosphere_solver), intent(in) :: this
    complex(real64), intent(in) :: field(0:, 0:, :)
    character(*), intent(in) :: name

    if (ubound(field, 1) < this%truncation + 1 .or. ubound(field, 2) < this%truncation .or. &
        size(field, 3) /= this%number_of_levels) then
      error stop 'dry atmosphere initial '//name//' has an inconsistent shape'
    end if
  end subroutine check_state_shape

end module dry_atmosphere
