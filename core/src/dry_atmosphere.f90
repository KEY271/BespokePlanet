module dry_atmosphere
  use iso_fortran_env, only: real64
  use harmonics, only: harmonic_transform
  use barotropic_vorticity, only: earth_radius
  use shallow_water_nonlinear, only: diagnose_shallow_water_velocity
  use spectral_hyperdiffusion, only: apply_spectral_hyperdiffusion
  use raw_filter, only: apply_raw_filter
  use dry_vertical_coordinate, only: hybrid_sigma_coordinate
  use dry_gravity_wave, only: dry_gravity_wave_solver
  use dry_nonlinear, only: compute_dry_nonlinear_tendency
  use dry_initial_conditions, only: jablonowski_williamson_initial_state
  use dry_held_suarez, only: held_suarez_initial_state
  implicit none
  private

  real(real64), parameter, public :: dry_vorticity_diffusion_time = 4.0_real64*3600.0_real64
  real(real64), parameter, public :: dry_divergence_diffusion_time = 1.0_real64*3600.0_real64
  real(real64), parameter, public :: dry_temperature_diffusion_time = 4.0_real64*3600.0_real64

  type, public :: dry_atmosphere_solver
    private
    type(harmonic_transform) :: transform
    type(hybrid_sigma_coordinate) :: coordinate
    type(dry_gravity_wave_solver) :: gravity_wave
    integer :: truncation = -1
    integer :: number_of_levels = 0
    integer :: step_number = -1
    real(real64) :: dt = 0.0_real64
    logical :: held_suarez_forcing_enabled = .false.
    !> Advective CFL of the state that the most recent advance started from.
    real(real64) :: last_advance_cfl = 0.0_real64
    complex(real64), allocatable :: previous_zeta(:, :, :), previous_delta(:, :, :)
    complex(real64), allocatable :: previous_temperature(:, :, :), previous_log_ps(:, :)
    complex(real64), allocatable :: current_zeta(:, :, :), current_delta(:, :, :)
    complex(real64), allocatable :: current_temperature(:, :, :), current_log_ps(:, :)
    ! Time-independent lower boundary condition; it is never advanced or filtered.
    complex(real64), allocatable :: surface_geopotential(:, :)
  contains
    procedure, public :: init => initialize_dry_solver
    procedure, public :: set_initial_state
    procedure, public :: set_jablonowski_williamson_state
    procedure, public :: set_held_suarez_state
    procedure, public :: advance
    procedure, public :: get_fields
    procedure, public :: get_spectral_state
    procedure, public :: get_reference_atmosphere
    procedure, public :: get_step
    procedure, public :: get_last_advance_cfl
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
    this%dt = dt
    this%step_number = -1
    this%last_advance_cfl = 0.0_real64
    this%held_suarez_forcing_enabled = .false.
    call this%transform%init(truncation)
    if (present(a_half)) then
      call this%coordinate%init(a_half, b_half)
    else
      call this%coordinate%init_default()
    end if
    this%number_of_levels = this%coordinate%number_of_levels
    call this%gravity_wave%init(truncation, dt, this%coordinate)
    call allocate_state(this, this%previous_zeta, this%previous_delta, &
                        this%previous_temperature, this%previous_log_ps)
    call allocate_state(this, this%current_zeta, this%current_delta, &
                        this%current_temperature, this%current_log_ps)
    this%previous_zeta = 0.0_real64
    this%previous_delta = 0.0_real64
    this%previous_temperature = 0.0_real64
    this%previous_log_ps = 0.0_real64
    this%current_zeta = 0.0_real64
    this%current_delta = 0.0_real64
    this%current_temperature = 0.0_real64
    this%current_log_ps = 0.0_real64
    if (allocated(this%surface_geopotential)) deallocate (this%surface_geopotential)
    allocate (this%surface_geopotential(0:truncation + 1, 0:truncation))
    this%surface_geopotential = 0.0_real64
  end subroutine initialize_dry_solver

  subroutine set_initial_state(this, zeta, delta, temperature, log_surface_pressure, surface_geopotential)
    class(dry_atmosphere_solver), intent(inout) :: this
    complex(real64), intent(in) :: zeta(0:, 0:, :), delta(0:, 0:, :), temperature(0:, 0:, :)
    complex(real64), intent(in) :: log_surface_pressure(0:, 0:)
    !> Fixed surface geopotential; a flat surface (zero) is used when absent.
    complex(real64), intent(in), optional :: surface_geopotential(0:, 0:)

    call check_initialized(this)
    call check_state_shape(this, zeta, 'zeta')
    call check_state_shape(this, delta, 'delta')
    call check_state_shape(this, temperature, 'temperature')
    if (ubound(log_surface_pressure, 1) < this%truncation + 1 .or. &
        ubound(log_surface_pressure, 2) < this%truncation) then
      error stop 'dry atmosphere initial log surface pressure has an inconsistent shape'
    end if
    this%current_zeta = zeta(0:this%truncation + 1, 0:this%truncation, 1:this%number_of_levels)
    this%current_delta = delta(0:this%truncation + 1, 0:this%truncation, 1:this%number_of_levels)
    this%current_temperature = temperature(0:this%truncation + 1, 0:this%truncation, 1:this%number_of_levels)
    this%current_log_ps = log_surface_pressure(0:this%truncation + 1, 0:this%truncation)
    this%surface_geopotential = 0.0_real64
    if (present(surface_geopotential)) then
      if (ubound(surface_geopotential, 1) < this%truncation + 1 .or. &
          ubound(surface_geopotential, 2) < this%truncation) then
        error stop 'dry atmosphere surface geopotential has an inconsistent shape'
      end if
      this%surface_geopotential = surface_geopotential(0:this%truncation + 1, 0:this%truncation)
    end if
    call enforce_spectral_field(this%truncation, this%surface_geopotential, .false.)
    call enforce_state_constraints(this, this%current_zeta, this%current_delta, &
                                   this%current_temperature, this%current_log_ps)
    this%previous_zeta = this%current_zeta
    this%previous_delta = this%current_delta
    this%previous_temperature = this%current_temperature
    this%previous_log_ps = this%current_log_ps
    this%held_suarez_forcing_enabled = .false.
    this%step_number = 0
  end subroutine set_initial_state

  subroutine set_jablonowski_williamson_state(this, include_perturbation)
    class(dry_atmosphere_solver), intent(inout) :: this
    logical, intent(in), optional :: include_perturbation
    complex(real64), allocatable :: zeta(:, :, :), delta(:, :, :), temperature(:, :, :), log_ps(:, :)
    complex(real64), allocatable :: surface_geopotential(:, :)
    logical :: perturb

    call check_initialized(this)
    perturb = .true.
    if (present(include_perturbation)) perturb = include_perturbation
    call jablonowski_williamson_initial_state(this%transform, this%truncation, this%coordinate, &
                                              perturb, zeta, delta, temperature, log_ps, &
                                              surface_geopotential)
    call this%set_initial_state(zeta, delta, temperature, log_ps, surface_geopotential)
  end subroutine set_jablonowski_williamson_state

  subroutine set_held_suarez_state(this)
    class(dry_atmosphere_solver), intent(inout) :: this
    complex(real64), allocatable :: zeta(:, :, :), delta(:, :, :), temperature(:, :, :), log_ps(:, :)
    complex(real64), allocatable :: surface_geopotential(:, :)

    call check_initialized(this)
    call held_suarez_initial_state(this%transform, this%truncation, this%coordinate, &
                                   zeta, delta, temperature, log_ps, surface_geopotential)
    call this%set_initial_state(zeta, delta, temperature, log_ps, surface_geopotential)
    this%held_suarez_forcing_enabled = .true.
  end subroutine set_held_suarez_state

  subroutine advance(this)
    class(dry_atmosphere_solver), intent(inout) :: this
    complex(real64), allocatable :: next_zeta(:, :, :), next_delta(:, :, :), next_temperature(:, :, :)
    complex(real64), allocatable :: next_log_ps(:, :)
    complex(real64), allocatable :: filtered_zeta(:, :, :), filtered_delta(:, :, :)
    complex(real64), allocatable :: filtered_temperature(:, :, :), filtered_log_ps(:, :)
    complex(real64), allocatable :: half_zeta(:, :, :), half_delta(:, :, :), half_temperature(:, :, :)
    complex(real64), allocatable :: half_log_ps(:, :)
    real(real64) :: maximum_speed, half_step_maximum_speed

    call check_ready(this)
    if (this%step_number == 0) then
      call integration_step(this, 0.25_real64*this%dt, &
                            this%current_zeta, this%current_delta, this%current_temperature, this%current_log_ps, &
                            this%current_zeta, this%current_delta, this%current_temperature, this%current_log_ps, &
                            .false., half_zeta, half_delta, half_temperature, half_log_ps, &
                            filtered_zeta, filtered_delta, filtered_temperature, filtered_log_ps, &
                            maximum_speed)
      call integration_step(this, 0.5_real64*this%dt, &
                            this%current_zeta, this%current_delta, this%current_temperature, this%current_log_ps, &
                            half_zeta, half_delta, half_temperature, half_log_ps, .true., &
                            next_zeta, next_delta, next_temperature, next_log_ps, &
                            filtered_zeta, filtered_delta, filtered_temperature, filtered_log_ps, &
                            half_step_maximum_speed)
      this%previous_zeta = this%current_zeta
      this%previous_delta = this%current_delta
      this%previous_temperature = this%current_temperature
      this%previous_log_ps = this%current_log_ps
    else
      call integration_step(this, this%dt, &
                            this%previous_zeta, this%previous_delta, this%previous_temperature, this%previous_log_ps, &
                            this%current_zeta, this%current_delta, this%current_temperature, this%current_log_ps, &
                            .true., next_zeta, next_delta, next_temperature, next_log_ps, &
                            filtered_zeta, filtered_delta, filtered_temperature, filtered_log_ps, &
                            maximum_speed)
      this%previous_zeta = filtered_zeta
      this%previous_delta = filtered_delta
      this%previous_temperature = filtered_temperature
      this%previous_log_ps = filtered_log_ps
    end if
    ! Both branches evaluated the first tendency at the state this step started from.
    this%last_advance_cfl = advective_cfl(this, maximum_speed)
    this%current_zeta = next_zeta
    this%current_delta = next_delta
    this%current_temperature = next_temperature
    this%current_log_ps = next_log_ps
    this%step_number = this%step_number + 1
  end subroutine advance

  subroutine integration_step(this, interval, previous_zeta, previous_delta, previous_temperature, previous_log_ps, &
                              current_zeta, current_delta, current_temperature, current_log_ps, apply_raw, &
                              next_zeta, next_delta, next_temperature, next_log_ps, &
                              filtered_zeta, filtered_delta, filtered_temperature, filtered_log_ps, &
                              maximum_speed)
    class(dry_atmosphere_solver), intent(inout) :: this
    real(real64), intent(in) :: interval
    complex(real64), intent(in) :: previous_zeta(0:, 0:, :), previous_delta(0:, 0:, :)
    complex(real64), intent(in) :: previous_temperature(0:, 0:, :), previous_log_ps(0:, 0:)
    complex(real64), intent(in) :: current_zeta(0:, 0:, :), current_delta(0:, 0:, :)
    complex(real64), intent(in) :: current_temperature(0:, 0:, :), current_log_ps(0:, 0:)
    logical, intent(in) :: apply_raw
    complex(real64), allocatable, intent(out) :: next_zeta(:, :, :), next_delta(:, :, :)
    complex(real64), allocatable, intent(out) :: next_temperature(:, :, :), next_log_ps(:, :)
    complex(real64), allocatable, intent(out) :: filtered_zeta(:, :, :), filtered_delta(:, :, :)
    complex(real64), allocatable, intent(out) :: filtered_temperature(:, :, :), filtered_log_ps(:, :)
    !> Largest wind speed (m/s) of the current state used for the tendency.
    real(real64), intent(out) :: maximum_speed
    complex(real64), allocatable :: rhs_zeta(:, :, :), rhs_delta(:, :, :), rhs_temperature(:, :, :), rhs_log_ps(:, :)
    complex(real64), allocatable :: candidate_zeta(:, :, :), candidate_delta(:, :, :)
    complex(real64), allocatable :: candidate_temperature(:, :, :), candidate_log_ps(:, :)
    complex(real64), allocatable :: filtered_level(:, :), next_level(:, :)
    real(real64) :: centered_interval
    integer :: k

    centered_interval = 2.0_real64*interval
    call compute_dry_nonlinear_tendency(this%transform, this%truncation, this%coordinate, &
                                        current_zeta, current_delta, current_temperature, current_log_ps, &
                                        this%surface_geopotential, &
                                        rhs_zeta, rhs_delta, rhs_temperature, rhs_log_ps, maximum_speed, &
                                        this%held_suarez_forcing_enabled)
    call allocate_state(this, candidate_zeta, candidate_delta, candidate_temperature, candidate_log_ps)
    candidate_zeta = previous_zeta + centered_interval*rhs_zeta
    call this%gravity_wave%solve(centered_interval, &
      previous_log_ps, previous_delta, previous_temperature, &
      current_log_ps, current_delta, current_temperature, &
      rhs_log_ps, rhs_delta, rhs_temperature, candidate_log_ps, candidate_delta, candidate_temperature)
    do k = 1, this%number_of_levels
      call apply_spectral_hyperdiffusion(this%truncation, centered_interval, candidate_zeta(:, :, k), &
                                         dry_vorticity_diffusion_time)
      call apply_spectral_hyperdiffusion(this%truncation, centered_interval, candidate_delta(:, :, k), &
                                         dry_divergence_diffusion_time)
      call apply_spectral_hyperdiffusion(this%truncation, centered_interval, candidate_temperature(:, :, k), &
                                         dry_temperature_diffusion_time)
    end do
    call allocate_state(this, next_zeta, next_delta, next_temperature, next_log_ps)
    call allocate_state(this, filtered_zeta, filtered_delta, filtered_temperature, filtered_log_ps)
    if (apply_raw) then
      do k = 1, this%number_of_levels
        call apply_raw_filter(previous_zeta(:, :, k), current_zeta(:, :, k), candidate_zeta(:, :, k), &
                              filtered_level, next_level)
        filtered_zeta(:, :, k) = filtered_level
        next_zeta(:, :, k) = next_level
        call apply_raw_filter(previous_delta(:, :, k), current_delta(:, :, k), candidate_delta(:, :, k), &
                              filtered_level, next_level)
        filtered_delta(:, :, k) = filtered_level
        next_delta(:, :, k) = next_level
        call apply_raw_filter(previous_temperature(:, :, k), current_temperature(:, :, k), &
                              candidate_temperature(:, :, k), filtered_level, next_level)
        filtered_temperature(:, :, k) = filtered_level
        next_temperature(:, :, k) = next_level
      end do
      call apply_raw_filter(previous_log_ps, current_log_ps, candidate_log_ps, filtered_log_ps, next_log_ps)
    else
      filtered_zeta = current_zeta
      filtered_delta = current_delta
      filtered_temperature = current_temperature
      filtered_log_ps = current_log_ps
      next_zeta = candidate_zeta
      next_delta = candidate_delta
      next_temperature = candidate_temperature
      next_log_ps = candidate_log_ps
    end if
    call enforce_state_constraints(this, filtered_zeta, filtered_delta, filtered_temperature, filtered_log_ps)
    call enforce_state_constraints(this, next_zeta, next_delta, next_temperature, next_log_ps)
  end subroutine integration_step

  subroutine get_fields(this, zeta, delta, temperature, surface_pressure, u, v, cfl)
    class(dry_atmosphere_solver), intent(inout) :: this
    real(real64), allocatable, intent(out) :: zeta(:, :, :), delta(:, :, :), temperature(:, :, :)
    real(real64), allocatable, intent(out) :: surface_pressure(:, :), u(:, :, :), v(:, :, :)
    real(real64), intent(out), optional :: cfl
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
      call this%transform%spectral_to_grid(this%current_zeta(:, :, k), temporary)
      zeta(:, :, k) = temporary
      call this%transform%spectral_to_grid(this%current_delta(:, :, k), temporary)
      delta(:, :, k) = temporary
      call this%transform%spectral_to_grid(this%current_temperature(:, :, k), temporary)
      temperature(:, :, k) = temporary
      call diagnose_shallow_water_velocity(this%transform, this%truncation, &
        this%current_zeta(:, :, k), this%current_delta(:, :, k), temporary_u, temporary_v)
      u(:, :, k) = temporary_u
      v(:, :, k) = temporary_v
    end do
    call this%transform%spectral_to_grid(this%current_log_ps, log_ps)
    surface_pressure = exp(log_ps)
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

  subroutine get_spectral_state(this, zeta, delta, temperature, log_surface_pressure)
    class(dry_atmosphere_solver), intent(in) :: this
    complex(real64), allocatable, intent(out) :: zeta(:, :, :), delta(:, :, :), temperature(:, :, :)
    complex(real64), allocatable, intent(out) :: log_surface_pressure(:, :)

    call check_ready(this)
    zeta = this%current_zeta
    delta = this%current_delta
    temperature = this%current_temperature
    log_surface_pressure = this%current_log_ps
  end subroutine get_spectral_state

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
    temperature = this%coordinate%reference_temperature
    if (present(a_half)) a_half = this%coordinate%a_half
    if (present(b_half)) b_half = this%coordinate%b_half
  end subroutine get_reference_atmosphere

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

  pure real(real64) function advective_cfl(this, maximum_speed) result(cfl)
    class(dry_atmosphere_solver), intent(in) :: this
    real(real64), intent(in) :: maximum_speed
    cfl = maximum_speed*this%dt/earth_radius* &
          sqrt(real(this%truncation*(this%truncation + 1), real64))
  end function advective_cfl

  subroutine allocate_state(this, zeta, delta, temperature, log_ps)
    class(dry_atmosphere_solver), intent(in) :: this
    complex(real64), allocatable, intent(out) :: zeta(:, :, :), delta(:, :, :), temperature(:, :, :)
    complex(real64), allocatable, intent(out) :: log_ps(:, :)

    allocate (zeta(0:this%truncation + 1, 0:this%truncation, this%number_of_levels))
    allocate (delta(0:this%truncation + 1, 0:this%truncation, this%number_of_levels))
    allocate (temperature(0:this%truncation + 1, 0:this%truncation, this%number_of_levels))
    allocate (log_ps(0:this%truncation + 1, 0:this%truncation))
  end subroutine allocate_state

  subroutine enforce_state_constraints(this, zeta, delta, temperature, log_ps)
    class(dry_atmosphere_solver), intent(in) :: this
    complex(real64), intent(inout) :: zeta(0:, 0:, :), delta(0:, 0:, :), temperature(0:, 0:, :)
    complex(real64), intent(inout) :: log_ps(0:, 0:)
    integer :: k

    do k = 1, this%number_of_levels
      call enforce_spectral_field(this%truncation, zeta(:, :, k), .true.)
      call enforce_spectral_field(this%truncation, delta(:, :, k), .true.)
      call enforce_spectral_field(this%truncation, temperature(:, :, k), .false.)
    end do
    call enforce_spectral_field(this%truncation, log_ps, .false.)
  end subroutine enforce_state_constraints

  subroutine enforce_spectral_field(truncation, field, zero_mean)
    integer, intent(in) :: truncation
    complex(real64), intent(inout) :: field(0:, 0:)
    logical, intent(in) :: zero_mean
    integer :: n, m

    if (zero_mean) field(0, 0) = 0.0_real64
    field(truncation + 1, :) = 0.0_real64
    do m = 1, truncation
      do n = 0, m - 1
        field(n, m) = 0.0_real64
      end do
    end do
  end subroutine enforce_spectral_field

  subroutine check_initialized(this)
    class(dry_atmosphere_solver), intent(in) :: this
    if (this%truncation < 1 .or. .not. allocated(this%current_zeta)) then
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
