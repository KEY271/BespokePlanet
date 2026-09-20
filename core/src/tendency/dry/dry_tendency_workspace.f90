!> Shared diagnostic fields and right-hand-side accumulators of one dry tendency
!> evaluation.
!>
!> Every array is allocated once, when the model is initialized, and reused on
!> each step.  Splitting the dry tendency into separate modules must not add
!> spectral transforms or large allocations, so the grid diagnostics that more
!> than one tendency needs are computed here exactly once per evaluation.
!>
!> The forcing_* components are the grid-space right-hand side.  Tendency
!> modules add into them; dry_tendency_projection turns them into spectral
!> tendencies with the same number of transforms the single kernel used.
module dry_tendency_workspace
  use iso_fortran_env, only: real64
  use harmonics, only: harmonic_transform
  use planet_parameters, only: earth_radius
  use spectral_vector_operators, only: diagnose_horizontal_velocity
  use dry_vertical_coordinate, only: hybrid_sigma_coordinate, dry_air_gas_constant
  use dry_state, only: dry_state_type
  use dry_physics_config, only: dry_model_physics_config
  implicit none
  private

  type, public :: dry_workspace_type
    integer :: truncation = -1
    integer :: number_of_levels = 0
    integer :: nx = 0
    integer :: ny = 0
    integer, allocatable :: ring_nlon(:)
    real(real64), allocatable :: gaussian_weights(:)

    !> Current state on the grid.
    real(real64), allocatable :: zeta_grid(:, :, :), delta_grid(:, :, :), temperature_grid(:, :, :)
    real(real64), allocatable :: u(:, :, :), v(:, :, :)
    real(real64), allocatable :: log_ps(:, :), ps(:, :)
    real(real64), allocatable :: dlogps_dlambda(:, :), dlogps_dphi(:, :)
    real(real64), allocatable :: surface_geopotential_grid(:, :)
    real(real64), allocatable :: pressure_half(:, :, :), delta_p(:, :, :)
    real(real64), allocatable :: layer_l(:, :, :), alpha(:, :, :)
    real(real64), allocatable :: geopotential(:, :, :), geopotential_half(:, :, :)
    real(real64), allocatable :: mass_divergence(:, :, :), cumulative(:, :, :), mass_flux(:, :, :)
    real(real64), allocatable :: pressure_gradient_u(:, :, :), pressure_gradient_v(:, :, :)
    real(real64), allocatable :: vertical_u(:, :, :), vertical_v(:, :, :), vertical_t(:, :, :)

    !> RAW-filtered previous time level; every prescribed physical tendency
    !> evaluates on this state, which is why it is diagnosed once here.
    real(real64), allocatable :: previous_temperature_grid(:, :, :)
    real(real64), allocatable :: previous_u(:, :, :), previous_v(:, :, :)
    real(real64), allocatable :: previous_log_ps(:, :), previous_ps(:, :)
    real(real64), allocatable :: previous_pressure_half(:, :, :), previous_alpha(:, :, :)
    !> Full-level pressure of the previous time level, shared by the Held-Suarez
    !> drag and relaxation so that neither recomputes it.
    real(real64), allocatable :: previous_full_level_pressure(:, :, :)
    real(real64), allocatable :: previous_surface_temperature_grid(:, :)
    real(real64), allocatable :: previous_deep_temperature_grid(:, :)
    real(real64), allocatable :: surface_temperature_grid(:, :), deep_temperature_grid(:, :)

    !> Grid-space right-hand side that the tendency modules add into.
    !> forcing_u/forcing_v hold the combined vector whose curl and divergence
    !> give the vorticity and divergence tendencies in one pair of transforms.
    real(real64), allocatable :: forcing_u(:, :, :), forcing_v(:, :, :)
    real(real64), allocatable :: forcing_temperature(:, :, :)
    real(real64), allocatable :: kinetic_geopotential(:, :, :)
    real(real64), allocatable :: forcing_log_ps(:, :)
    real(real64), allocatable :: forcing_surface_temperature(:, :), forcing_deep_temperature(:, :)

    !> Scratch reused within a level.
    real(real64), allocatable :: temporary_grid(:, :), temporary_u(:, :), temporary_v(:, :)
    real(real64), allocatable :: dtdlambda(:, :), dtdphi(:, :)

    real(real64) :: evaluation_time = 0.0_real64
    !> Planet rotation rate the dynamics use, supplied by the model configuration.
    real(real64) :: active_rotation_rate = 0.0_real64
    logical :: physics_grids_ready = .false.
    logical :: held_suarez_grids_ready = .false.
    logical :: radiation_grids_ready = .false.
  contains
    procedure, public :: initialize => initialize_workspace
    procedure, public :: prepare => prepare_workspace
    procedure, public :: zero_forcing
  end type dry_workspace_type

contains

  subroutine initialize_workspace(this, transform, truncation, number_of_levels)
    class(dry_workspace_type), intent(inout) :: this
    type(harmonic_transform), intent(inout) :: transform
    integer, intent(in) :: truncation, number_of_levels
    real(real64), allocatable :: field(:, :)
    integer :: nx, ny, levels

    call transform%allocate_field(field)
    nx = size(field, 1)
    ny = size(field, 2)
    levels = number_of_levels
    this%truncation = truncation
    this%number_of_levels = levels
    this%nx = nx
    this%ny = ny
    this%ring_nlon = transform%get_nlon()
    this%gaussian_weights = transform%get_gaussian_weights()

    allocate (this%zeta_grid(nx, ny, levels), this%delta_grid(nx, ny, levels))
    allocate (this%temperature_grid(nx, ny, levels))
    allocate (this%u(nx, ny, levels), this%v(nx, ny, levels))
    allocate (this%log_ps(nx, ny), this%ps(nx, ny))
    allocate (this%dlogps_dlambda(nx, ny), this%dlogps_dphi(nx, ny))
    allocate (this%surface_geopotential_grid(nx, ny))
    allocate (this%pressure_half(nx, ny, 0:levels))
    allocate (this%delta_p(nx, ny, levels), this%layer_l(nx, ny, levels))
    allocate (this%alpha(nx, ny, levels))
    allocate (this%geopotential(nx, ny, levels), this%geopotential_half(nx, ny, 0:levels))
    allocate (this%mass_divergence(nx, ny, levels), this%cumulative(nx, ny, 0:levels))
    allocate (this%mass_flux(nx, ny, 0:levels))
    allocate (this%pressure_gradient_u(nx, ny, levels), this%pressure_gradient_v(nx, ny, levels))
    allocate (this%vertical_u(nx, ny, levels), this%vertical_v(nx, ny, levels))
    allocate (this%vertical_t(nx, ny, levels))

    allocate (this%previous_temperature_grid(nx, ny, levels))
    allocate (this%previous_u(nx, ny, levels), this%previous_v(nx, ny, levels))
    allocate (this%previous_log_ps(nx, ny), this%previous_ps(nx, ny))
    allocate (this%previous_pressure_half(nx, ny, 0:levels))
    allocate (this%previous_alpha(nx, ny, levels))
    allocate (this%previous_full_level_pressure(nx, ny, levels))
    allocate (this%previous_surface_temperature_grid(nx, ny))
    allocate (this%previous_deep_temperature_grid(nx, ny))
    allocate (this%surface_temperature_grid(nx, ny), this%deep_temperature_grid(nx, ny))

    allocate (this%forcing_u(nx, ny, levels), this%forcing_v(nx, ny, levels))
    allocate (this%forcing_temperature(nx, ny, levels))
    allocate (this%kinetic_geopotential(nx, ny, levels))
    allocate (this%forcing_log_ps(nx, ny))
    allocate (this%forcing_surface_temperature(nx, ny), this%forcing_deep_temperature(nx, ny))

    allocate (this%temporary_grid(nx, ny), this%temporary_u(nx, ny), this%temporary_v(nx, ny))
    allocate (this%dtdlambda(nx, ny), this%dtdphi(nx, ny))
  end subroutine initialize_workspace

  !> Zeroes the grid-space right-hand side.  Points outside a reduced-grid ring
  !> are never written by a tendency, so they must start at zero.
  subroutine zero_forcing(this)
    class(dry_workspace_type), intent(inout) :: this

    this%forcing_u = 0.0_real64
    this%forcing_v = 0.0_real64
    this%forcing_temperature = 0.0_real64
    this%kinetic_geopotential = 0.0_real64
    this%forcing_log_ps = 0.0_real64
    this%forcing_surface_temperature = 0.0_real64
    this%forcing_deep_temperature = 0.0_real64
  end subroutine zero_forcing

  !> Diagnoses every grid field that more than one tendency needs.
  subroutine prepare_workspace(this, transform, coordinate, rotation_rate, state, physics_state, &
                               surface_geopotential, physics, evaluation_time)
    class(dry_workspace_type), intent(inout) :: this
    type(harmonic_transform), intent(inout) :: transform
    type(hybrid_sigma_coordinate), intent(in) :: coordinate
    real(real64), intent(in) :: rotation_rate
    type(dry_state_type), intent(in) :: state, physics_state
    complex(real64), intent(in) :: surface_geopotential(0:, 0:)
    type(dry_model_physics_config), intent(in) :: physics
    real(real64), intent(in) :: evaluation_time
    integer :: k, levels

    levels = this%number_of_levels
    if (coordinate%number_of_levels /= levels) then
      error stop 'dry tendency workspace is not initialized for this vertical coordinate'
    end if
    if (size(state%zeta, 3) /= levels .or. size(physics_state%zeta, 3) /= levels) then
      error stop 'dry nonlinear state has the wrong number of vertical levels'
    end if
    this%evaluation_time = evaluation_time
    this%physics_grids_ready = physics%held_suarez%enabled .or. physics%surface_friction%enabled .or. &
                               physics%rayleigh_friction%enabled .or. physics%radiation%enabled .or. &
                               physics%convection%enabled
    this%held_suarez_grids_ready = physics%held_suarez%enabled .or. physics%surface_friction%enabled
    this%radiation_grids_ready = physics%radiation%enabled
    this%active_rotation_rate = rotation_rate

    this%zeta_grid = 0.0_real64
    this%delta_grid = 0.0_real64
    this%temperature_grid = 0.0_real64
    this%u = 0.0_real64
    this%v = 0.0_real64
    do k = 1, levels
      call transform%spectral_to_grid(state%zeta(:, :, k), this%temporary_grid)
      this%zeta_grid(:, :, k) = this%temporary_grid
      call transform%spectral_to_grid(state%delta(:, :, k), this%temporary_grid)
      this%delta_grid(:, :, k) = this%temporary_grid
      call transform%spectral_to_grid(state%temperature(:, :, k), this%temporary_grid)
      this%temperature_grid(:, :, k) = this%temporary_grid
      call diagnose_horizontal_velocity(transform, this%truncation, state%zeta(:, :, k), &
                                        state%delta(:, :, k), this%temporary_u, this%temporary_v)
      this%u(:, :, k) = this%temporary_u
      this%v(:, :, k) = this%temporary_v
      if (this%physics_grids_ready) then
        call transform%spectral_to_grid(physics_state%temperature(:, :, k), this%temporary_grid)
        this%previous_temperature_grid(:, :, k) = this%temporary_grid
        call diagnose_horizontal_velocity(transform, this%truncation, physics_state%zeta(:, :, k), &
                                          physics_state%delta(:, :, k), this%temporary_u, this%temporary_v)
        this%previous_u(:, :, k) = this%temporary_u
        this%previous_v(:, :, k) = this%temporary_v
      end if
    end do
    call transform%spectral_to_grid(state%log_surface_pressure, this%log_ps)
    this%ps = exp(this%log_ps)
    if (this%physics_grids_ready) then
      call transform%spectral_to_grid(physics_state%log_surface_pressure, this%previous_log_ps)
      this%previous_ps = exp(this%previous_log_ps)
      do k = 0, levels
        this%previous_pressure_half(:, :, k) = coordinate%a_half(k) + coordinate%b_half(k)*this%previous_ps
      end do
      do k = 1, levels
        this%previous_alpha(:, :, k) = 1.0_real64 - this%previous_pressure_half(:, :, k - 1)* &
          log(this%previous_pressure_half(:, :, k)/this%previous_pressure_half(:, :, k - 1))/ &
          (this%previous_pressure_half(:, :, k) - this%previous_pressure_half(:, :, k - 1))
      end do
      if (this%held_suarez_grids_ready) then
        do k = 1, levels
          this%previous_full_level_pressure(:, :, k) = this%previous_pressure_half(:, :, k)* &
                                                       exp(-this%previous_alpha(:, :, k))
        end do
      end if
    end if
    if (this%radiation_grids_ready) then
      call transform%spectral_to_grid(state%surface_temperature, this%surface_temperature_grid)
      call transform%spectral_to_grid(state%deep_temperature, this%deep_temperature_grid)
      call transform%spectral_to_grid(physics_state%surface_temperature, &
                                      this%previous_surface_temperature_grid)
      call transform%spectral_to_grid(physics_state%deep_temperature, &
                                      this%previous_deep_temperature_grid)
    end if
    call transform%gradient_to_grid(state%log_surface_pressure, this%dlogps_dlambda, this%dlogps_dphi)

    this%pressure_half = 0.0_real64
    this%delta_p = 0.0_real64
    this%layer_l = 0.0_real64
    this%alpha = 0.0_real64
    do k = 0, levels
      this%pressure_half(:, :, k) = coordinate%a_half(k) + coordinate%b_half(k)*this%ps
    end do
    do k = 1, levels
      this%delta_p(:, :, k) = this%pressure_half(:, :, k) - this%pressure_half(:, :, k - 1)
      this%layer_l(:, :, k) = log(this%pressure_half(:, :, k)/this%pressure_half(:, :, k - 1))
      this%alpha(:, :, k) = 1.0_real64 - this%pressure_half(:, :, k - 1)* &
                            this%layer_l(:, :, k)/this%delta_p(:, :, k)
    end do

    ! The surface geopotential is a fixed lower boundary condition, not a prognostic field.
    call transform%spectral_to_grid(surface_geopotential, this%surface_geopotential_grid)
    this%geopotential_half(:, :, levels) = this%surface_geopotential_grid
    do k = levels, 1, -1
      this%geopotential(:, :, k) = this%geopotential_half(:, :, k) + &
        this%alpha(:, :, k)*dry_air_gas_constant*this%temperature_grid(:, :, k)
      this%geopotential_half(:, :, k - 1) = this%geopotential_half(:, :, k) + &
        dry_air_gas_constant*this%temperature_grid(:, :, k)*this%layer_l(:, :, k)
    end do

    call accumulate_mass_divergence(this%nx, this%ny, levels, this%ring_nlon, transform%mu, &
                                    coordinate%delta_b, coordinate%b_half, this%delta_p, this%delta_grid, &
                                    this%u, this%v, this%ps, this%dlogps_dlambda, this%dlogps_dphi, &
                                    this%mass_divergence, this%cumulative, this%mass_flux)

    call compute_pressure_gradient_and_vertical_advection(this%nx, this%ny, levels, this%ring_nlon, &
                                    transform%mu, coordinate%delta_b, coordinate%b_half, this%delta_p, &
                                    this%layer_l, this%alpha, this%ps, this%dlogps_dlambda, this%dlogps_dphi, &
                                    this%u, this%v, this%temperature_grid, this%mass_flux, &
                                    this%pressure_gradient_u, this%pressure_gradient_v, &
                                    this%vertical_u, this%vertical_v, this%vertical_t)
  end subroutine prepare_workspace

  !> The grid loops below are separate procedures with explicit-shape dummies so that the
  !> compiler sees contiguous arrays that cannot alias, as it did when these were local
  !> arrays of one large kernel.  Accessing them as components of the workspace through
  !> the type-bound procedure costs about a factor of two here.
  subroutine accumulate_mass_divergence(nx, ny, levels, ring_nlon, mu, delta_b, b_half, delta_p, &
                                        delta_grid, u, v, ps, dlogps_dlambda, dlogps_dphi, &
                                        mass_divergence, cumulative, mass_flux)
    integer, intent(in) :: nx, ny, levels, ring_nlon(ny)
    real(real64), intent(in) :: mu(ny), delta_b(levels), b_half(0:levels)
    real(real64), intent(in) :: delta_p(nx, ny, levels), delta_grid(nx, ny, levels)
    real(real64), intent(in) :: u(nx, ny, levels), v(nx, ny, levels)
    real(real64), intent(in) :: ps(nx, ny), dlogps_dlambda(nx, ny), dlogps_dphi(nx, ny)
    real(real64), intent(out) :: mass_divergence(nx, ny, levels), cumulative(nx, ny, 0:levels)
    real(real64), intent(out) :: mass_flux(nx, ny, 0:levels)
    integer :: i, j, k
    real(real64) :: cosphi, gradient_u, gradient_v

    mass_divergence = 0.0_real64
    cumulative = 0.0_real64
    do k = 1, levels
      do j = 1, ny
        cosphi = sqrt(max(0.0_real64, 1.0_real64 - mu(j)**2))
        do i = 1, ring_nlon(j)
          gradient_u = dlogps_dlambda(i, j)/(earth_radius*cosphi)
          gradient_v = dlogps_dphi(i, j)/earth_radius
          mass_divergence(i, j, k) = delta_p(i, j, k)*delta_grid(i, j, k) + &
            ps(i, j)*delta_b(k)*(u(i, j, k)*gradient_u + v(i, j, k)*gradient_v)
          cumulative(i, j, k) = cumulative(i, j, k - 1) + mass_divergence(i, j, k)
        end do
      end do
    end do
    do k = 0, levels
      mass_flux(:, :, k) = b_half(k)*cumulative(:, :, levels) - cumulative(:, :, k)
    end do
  end subroutine accumulate_mass_divergence

  subroutine compute_pressure_gradient_and_vertical_advection(nx, ny, levels, ring_nlon, mu, delta_b, &
                                        b_half, delta_p, layer_l, alpha, ps, dlogps_dlambda, dlogps_dphi, &
                                        u, v, temperature_grid, mass_flux, pressure_gradient_u, &
                                        pressure_gradient_v, vertical_u, vertical_v, vertical_t)
    integer, intent(in) :: nx, ny, levels, ring_nlon(ny)
    real(real64), intent(in) :: mu(ny), delta_b(levels), b_half(0:levels)
    real(real64), intent(in) :: delta_p(nx, ny, levels), layer_l(nx, ny, levels), alpha(nx, ny, levels)
    real(real64), intent(in) :: ps(nx, ny), dlogps_dlambda(nx, ny), dlogps_dphi(nx, ny)
    real(real64), intent(in) :: u(nx, ny, levels), v(nx, ny, levels), temperature_grid(nx, ny, levels)
    real(real64), intent(in) :: mass_flux(nx, ny, 0:levels)
    real(real64), intent(out) :: pressure_gradient_u(nx, ny, levels), pressure_gradient_v(nx, ny, levels)
    real(real64), intent(out) :: vertical_u(nx, ny, levels), vertical_v(nx, ny, levels)
    real(real64), intent(out) :: vertical_t(nx, ny, levels)
    integer :: i, j, k
    real(real64) :: cosphi, coefficient

    pressure_gradient_u = 0.0_real64
    pressure_gradient_v = 0.0_real64
    vertical_u = 0.0_real64
    vertical_v = 0.0_real64
    vertical_t = 0.0_real64
    do k = 1, levels
      do j = 1, ny
        cosphi = sqrt(max(0.0_real64, 1.0_real64 - mu(j)**2))
        do i = 1, ring_nlon(j)
          coefficient = ps(i, j)/delta_p(i, j, k)*(b_half(k - 1)* &
            layer_l(i, j, k) + alpha(i, j, k)*delta_b(k))
          pressure_gradient_u(i, j, k) = coefficient*dlogps_dlambda(i, j)/(earth_radius*cosphi)
          pressure_gradient_v(i, j, k) = coefficient*dlogps_dphi(i, j)/earth_radius
          if (k < levels) then
            vertical_u(i, j, k) = vertical_u(i, j, k) + mass_flux(i, j, k)* &
              (u(i, j, k + 1) - u(i, j, k))
            vertical_v(i, j, k) = vertical_v(i, j, k) + mass_flux(i, j, k)* &
              (v(i, j, k + 1) - v(i, j, k))
            vertical_t(i, j, k) = vertical_t(i, j, k) + mass_flux(i, j, k)* &
              (temperature_grid(i, j, k + 1) - temperature_grid(i, j, k))
          end if
          if (k > 1) then
            vertical_u(i, j, k) = vertical_u(i, j, k) + mass_flux(i, j, k - 1)* &
              (u(i, j, k) - u(i, j, k - 1))
            vertical_v(i, j, k) = vertical_v(i, j, k) + mass_flux(i, j, k - 1)* &
              (v(i, j, k) - v(i, j, k - 1))
            vertical_t(i, j, k) = vertical_t(i, j, k) + mass_flux(i, j, k - 1)* &
              (temperature_grid(i, j, k) - temperature_grid(i, j, k - 1))
          end if
          vertical_u(i, j, k) = vertical_u(i, j, k)/(2.0_real64*delta_p(i, j, k))
          vertical_v(i, j, k) = vertical_v(i, j, k)/(2.0_real64*delta_p(i, j, k))
          vertical_t(i, j, k) = vertical_t(i, j, k)/(2.0_real64*delta_p(i, j, k))
        end do
      end do
    end do
  end subroutine compute_pressure_gradient_and_vertical_advection

end module dry_tendency_workspace
