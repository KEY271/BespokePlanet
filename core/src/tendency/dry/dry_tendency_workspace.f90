!> Shared diagnostic fields and right-hand-side accumulators of one dry (or
!> moist) tendency evaluation.
!>
!> Every array is allocated once, when the model is initialized, and reused on
!> each step.  Splitting the dry tendency into separate modules must not add
!> spectral transforms or large allocations, so the grid diagnostics that more
!> than one tendency needs are computed here exactly once per evaluation.
!>
!> The forcing_* components are the grid-space right-hand side.  Tendency
!> modules add into them; dry_tendency_projection turns them into spectral
!> tendencies with the same number of transforms the single kernel used.
!>
!> With moisture enabled the virtual temperature T_v = (1 + delta_v q^+) T
!> replaces T in the hydrostatic relation, the pressure-gradient term and the
!> adiabatic heating; without it virtual_temperature_grid is a copy of
!> temperature_grid and the dry equations are reproduced exactly.
module dry_tendency_workspace
  use iso_fortran_env, only: real64
  use harmonics, only: harmonic_transform
  use planet_parameters, only: earth_radius
  use spectral_vector_operators, only: diagnose_horizontal_velocity
  use dry_vertical_coordinate, only: hybrid_sigma_coordinate, dry_air_gas_constant
  use dry_state, only: dry_state_type
  use dry_physics_config, only: dry_model_physics_config
  use moist_thermodynamics, only: virtual_temperature_coefficient
  implicit none
  private

  type, public :: dry_workspace_type
    integer :: truncation = -1
    integer :: number_of_levels = 0
    integer :: nx = 0
    integer :: ny = 0
    integer, allocatable :: ring_nlon(:)
    real(real64), allocatable :: gaussian_weights(:)
    !> Full-level eta of the vertical coordinate, for the wind-maximum diagnostic.
    real(real64), allocatable :: full_level_eta(:)

    !> Current state on the grid.
    real(real64), allocatable :: zeta_grid(:, :, :), delta_grid(:, :, :), temperature_grid(:, :, :)
    real(real64), allocatable :: humidity_grid(:, :, :), virtual_temperature_grid(:, :, :)
    real(real64), allocatable :: u(:, :, :), v(:, :, :)
    real(real64), allocatable :: log_ps(:, :), ps(:, :)
    real(real64), allocatable :: dlogps_dlambda(:, :), dlogps_dphi(:, :)
    real(real64), allocatable :: surface_geopotential_grid(:, :)
    real(real64), allocatable :: land_fraction(:, :)
    real(real64), allocatable :: pressure_half(:, :, :), delta_p(:, :, :)
    real(real64), allocatable :: layer_l(:, :, :), alpha(:, :, :)
    real(real64), allocatable :: geopotential(:, :, :), geopotential_half(:, :, :)
    real(real64), allocatable :: mass_divergence(:, :, :), cumulative(:, :, :), mass_flux(:, :, :)
    real(real64), allocatable :: pressure_gradient_u(:, :, :), pressure_gradient_v(:, :, :)
    real(real64), allocatable :: vertical_u(:, :, :), vertical_v(:, :, :), vertical_t(:, :, :)
    real(real64), allocatable :: vertical_q(:, :, :)

    !> RAW-filtered previous time level; every prescribed physical tendency
    !> evaluates on this state, which is why it is diagnosed once here.
    real(real64), allocatable :: previous_temperature_grid(:, :, :)
    real(real64), allocatable :: previous_humidity_grid(:, :, :)
    real(real64), allocatable :: previous_u(:, :, :), previous_v(:, :, :)
    real(real64), allocatable :: previous_log_ps(:, :), previous_ps(:, :)
    real(real64), allocatable :: previous_pressure_half(:, :, :), previous_alpha(:, :, :)
    real(real64), allocatable :: previous_delta_p(:, :, :), previous_layer_l(:, :, :)
    !> Full-level pressure of the previous time level, shared by the Held-Suarez
    !> drag and relaxation so that neither recomputes it.
    real(real64), allocatable :: previous_full_level_pressure(:, :, :)
    real(real64), allocatable :: previous_surface_temperature_grid(:, :)
    real(real64), allocatable :: previous_deep_temperature_grid(:, :)
    real(real64), allocatable :: surface_temperature_grid(:, :), deep_temperature_grid(:, :)
    real(real64), allocatable :: previous_surface_water(:, :), surface_water(:, :)

    !> Grid-space right-hand side that the tendency modules add into.
    !> forcing_u/forcing_v hold the combined vector whose curl and divergence
    !> give the vorticity and divergence tendencies in one pair of transforms.
    real(real64), allocatable :: forcing_u(:, :, :), forcing_v(:, :, :)
    real(real64), allocatable :: forcing_temperature(:, :, :)
    real(real64), allocatable :: forcing_humidity(:, :, :)
    real(real64), allocatable :: kinetic_geopotential(:, :, :)
    real(real64), allocatable :: forcing_log_ps(:, :)
    real(real64), allocatable :: forcing_surface_temperature(:, :), forcing_deep_temperature(:, :)
    real(real64), allocatable :: forcing_surface_water(:, :)

    !> Column fluxes of the physical processes, kept for the diagnostics sample.
    !> The latent heat flux L E is what the radiation tendency takes from the surface.
    real(real64), allocatable :: incoming_shortwave(:, :), reflected_shortwave(:, :), outgoing_longwave(:, :)
    real(real64), allocatable :: evaporation(:, :), land_evaporation(:, :), ocean_evaporation(:, :)
    real(real64), allocatable :: latent_heat_flux(:, :), surface_wetness(:, :)
    real(real64), allocatable :: runoff(:, :), water_budget_residual(:, :)
    real(real64), allocatable :: convective_precipitation(:, :), large_scale_precipitation(:, :)
    !> Effective column cloud cover diagnosed after the convective processes and
    !> used by the shortwave reflection of the same evaluation (docs/tendency/cloud.md).
    real(real64), allocatable :: cloud_cover(:, :)

    !> Scratch reused within a level.
    real(real64), allocatable :: temporary_grid(:, :), temporary_u(:, :), temporary_v(:, :)
    real(real64), allocatable :: dtdlambda(:, :), dtdphi(:, :)

    real(real64) :: evaluation_time = 0.0_real64
    !> Planet rotation rate the dynamics use, supplied by the model configuration.
    real(real64) :: active_rotation_rate = 0.0_real64
    !> Wetness below this configured fraction is counted as dry land.
    real(real64) :: bucket_dry_threshold_fraction = 0.1_real64
    logical :: surface_tiles_enabled = .false.
    logical :: sea_ice_enabled = .false.
    real(real64), allocatable :: land_temperature(:, :), previous_land_temperature(:, :), forcing_land_temperature(:, :)
    real(real64), allocatable :: ocean_temperature(:, :), previous_ocean_temperature(:, :), forcing_ocean_temperature(:, :)
    real(real64), allocatable :: sea_ice_fraction(:, :), previous_sea_ice_fraction(:, :), forcing_sea_ice_fraction(:, :)
    real(real64), allocatable :: sea_ice_volume(:, :), previous_sea_ice_volume(:, :), forcing_sea_ice_volume(:, :)
    real(real64), allocatable :: sea_ice_temperature(:, :), sea_ice_thickness(:, :)
    real(real64), allocatable :: ice_energy_residual(:, :), ice_surface_residual(:, :)
    logical :: physics_grids_ready = .false.
    logical :: held_suarez_grids_ready = .false.
    logical :: radiation_grids_ready = .false.
    logical :: moisture_grids_ready = .false.
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
    allocate (this%full_level_eta(levels))
    this%full_level_eta = 0.0_real64

    allocate (this%zeta_grid(nx, ny, levels), this%delta_grid(nx, ny, levels))
    allocate (this%temperature_grid(nx, ny, levels))
    allocate (this%humidity_grid(nx, ny, levels), this%virtual_temperature_grid(nx, ny, levels))
    allocate (this%u(nx, ny, levels), this%v(nx, ny, levels))
    allocate (this%log_ps(nx, ny), this%ps(nx, ny))
    allocate (this%dlogps_dlambda(nx, ny), this%dlogps_dphi(nx, ny))
    allocate (this%surface_geopotential_grid(nx, ny))
    allocate (this%land_fraction(nx, ny))
    allocate (this%pressure_half(nx, ny, 0:levels))
    allocate (this%delta_p(nx, ny, levels), this%layer_l(nx, ny, levels))
    allocate (this%alpha(nx, ny, levels))
    allocate (this%geopotential(nx, ny, levels), this%geopotential_half(nx, ny, 0:levels))
    allocate (this%mass_divergence(nx, ny, levels), this%cumulative(nx, ny, 0:levels))
    allocate (this%mass_flux(nx, ny, 0:levels))
    allocate (this%pressure_gradient_u(nx, ny, levels), this%pressure_gradient_v(nx, ny, levels))
    allocate (this%vertical_u(nx, ny, levels), this%vertical_v(nx, ny, levels))
    allocate (this%vertical_t(nx, ny, levels), this%vertical_q(nx, ny, levels))

    allocate (this%previous_temperature_grid(nx, ny, levels))
    allocate (this%previous_humidity_grid(nx, ny, levels))
    allocate (this%previous_u(nx, ny, levels), this%previous_v(nx, ny, levels))
    allocate (this%previous_log_ps(nx, ny), this%previous_ps(nx, ny))
    allocate (this%previous_pressure_half(nx, ny, 0:levels))
    allocate (this%previous_alpha(nx, ny, levels))
    allocate (this%previous_delta_p(nx, ny, levels), this%previous_layer_l(nx, ny, levels))
    allocate (this%previous_full_level_pressure(nx, ny, levels))
    allocate (this%previous_surface_temperature_grid(nx, ny))
    allocate (this%previous_deep_temperature_grid(nx, ny))
    allocate (this%surface_temperature_grid(nx, ny), this%deep_temperature_grid(nx, ny))
    allocate (this%previous_surface_water(nx, ny), this%surface_water(nx, ny))
    this%previous_surface_temperature_grid = 0.0_real64
    this%previous_deep_temperature_grid = 0.0_real64
    this%surface_temperature_grid = 0.0_real64
    this%deep_temperature_grid = 0.0_real64

    allocate (this%forcing_u(nx, ny, levels), this%forcing_v(nx, ny, levels))
    allocate (this%forcing_temperature(nx, ny, levels), this%forcing_humidity(nx, ny, levels))
    allocate (this%kinetic_geopotential(nx, ny, levels))
    allocate (this%forcing_log_ps(nx, ny))
    allocate (this%forcing_surface_temperature(nx, ny), this%forcing_deep_temperature(nx, ny))
    allocate (this%forcing_surface_water(nx, ny))

    allocate (this%incoming_shortwave(nx, ny), this%reflected_shortwave(nx, ny), this%outgoing_longwave(nx, ny))
    allocate (this%evaporation(nx, ny), this%land_evaporation(nx, ny), this%ocean_evaporation(nx, ny))
    allocate (this%latent_heat_flux(nx, ny), this%surface_wetness(nx, ny))
    allocate (this%runoff(nx, ny), this%water_budget_residual(nx, ny))
    allocate (this%convective_precipitation(nx, ny), this%large_scale_precipitation(nx, ny))
    allocate (this%cloud_cover(nx, ny))
    this%cloud_cover = 0.0_real64

    allocate (this%temporary_grid(nx, ny), this%temporary_u(nx, ny), this%temporary_v(nx, ny))
    allocate (this%dtdlambda(nx, ny), this%dtdphi(nx, ny))
    this%humidity_grid = 0.0_real64
    this%previous_humidity_grid = 0.0_real64
    this%vertical_q = 0.0_real64
    ! The ring-limited geometry loops leave the padding untouched, so define it once.
    this%ps = 1.0_real64
    this%previous_ps = 1.0_real64
    this%pressure_half = 0.0_real64
    this%delta_p = 1.0_real64
    this%layer_l = 0.0_real64
    this%alpha = 0.0_real64
    this%previous_pressure_half = 0.0_real64
    this%previous_alpha = 0.0_real64
    this%previous_delta_p = 1.0_real64
    this%previous_layer_l = 0.0_real64
    this%previous_full_level_pressure = 1.0_real64
    allocate (this%land_temperature(nx, ny), this%previous_land_temperature(nx, ny), this%forcing_land_temperature(nx, ny))
    this%land_temperature = 0.0_real64
    this%previous_land_temperature = 0.0_real64
    allocate (this%ocean_temperature(nx, ny), this%previous_ocean_temperature(nx, ny), this%forcing_ocean_temperature(nx, ny))
    this%ocean_temperature = 0.0_real64
    this%previous_ocean_temperature = 0.0_real64
    allocate (this%sea_ice_fraction(nx, ny), this%previous_sea_ice_fraction(nx, ny), this%forcing_sea_ice_fraction(nx, ny))
    this%sea_ice_fraction = 0.0_real64
    this%previous_sea_ice_fraction = 0.0_real64
    allocate (this%sea_ice_volume(nx, ny), this%previous_sea_ice_volume(nx, ny), this%forcing_sea_ice_volume(nx, ny))
    this%sea_ice_volume = 0.0_real64
    this%previous_sea_ice_volume = 0.0_real64
    allocate (this%sea_ice_temperature(nx, ny), this%sea_ice_thickness(nx, ny))
    allocate (this%ice_energy_residual(nx, ny), this%ice_surface_residual(nx, ny))
    this%sea_ice_temperature = 0.0_real64
    this%sea_ice_thickness = 0.0_real64
    this%previous_surface_water = 0.0_real64
    this%surface_water = 0.0_real64
  end subroutine initialize_workspace

  !> Zeroes the grid-space right-hand side and the flux records.  Points outside a
  !> reduced-grid ring are never written by a tendency, so they must start at zero.
  subroutine zero_forcing(this)
    class(dry_workspace_type), intent(inout) :: this

    this%forcing_u = 0.0_real64
    this%forcing_v = 0.0_real64
    this%forcing_temperature = 0.0_real64
    this%forcing_humidity = 0.0_real64
    this%kinetic_geopotential = 0.0_real64
    this%forcing_log_ps = 0.0_real64
    this%forcing_surface_temperature = 0.0_real64
    this%forcing_deep_temperature = 0.0_real64
    this%forcing_surface_water = 0.0_real64
    this%forcing_land_temperature = 0.0_real64
    this%forcing_ocean_temperature = 0.0_real64
    this%forcing_sea_ice_fraction = 0.0_real64
    this%forcing_sea_ice_volume = 0.0_real64
    this%ice_energy_residual = 0.0_real64
    this%ice_surface_residual = 0.0_real64
    this%incoming_shortwave = 0.0_real64
    this%reflected_shortwave = 0.0_real64
    this%outgoing_longwave = 0.0_real64
    this%evaporation = 0.0_real64
    this%land_evaporation = 0.0_real64
    this%ocean_evaporation = 0.0_real64
    this%latent_heat_flux = 0.0_real64
    this%surface_wetness = 0.0_real64
    this%runoff = 0.0_real64
    this%water_budget_residual = 0.0_real64
    this%convective_precipitation = 0.0_real64
    this%large_scale_precipitation = 0.0_real64
    this%cloud_cover = 0.0_real64
  end subroutine zero_forcing

  !> Diagnoses every grid field that more than one tendency needs.
  subroutine prepare_workspace(this, transform, coordinate, rotation_rate, state, physics_state, &
                               surface_geopotential, physics, evaluation_time, land_fraction)
    class(dry_workspace_type), intent(inout) :: this
    type(harmonic_transform), intent(inout) :: transform
    type(hybrid_sigma_coordinate), intent(in) :: coordinate
    real(real64), intent(in) :: rotation_rate
    type(dry_state_type), intent(in) :: state, physics_state
    complex(real64), intent(in) :: surface_geopotential(0:, 0:)
    real(real64), intent(in), optional :: land_fraction(:, :)
    type(dry_model_physics_config), intent(in) :: physics
    real(real64), intent(in) :: evaluation_time
    real(real64), allocatable :: grid(:, :), grid_u(:, :), grid_v(:, :)
    integer :: k, levels

    levels = this%number_of_levels
    if (coordinate%number_of_levels /= levels) then
      error stop 'dry tendency workspace is not initialized for this vertical coordinate'
    end if
    if (size(state%zeta, 3) /= levels .or. size(physics_state%zeta, 3) /= levels) then
      error stop 'dry nonlinear state has the wrong number of vertical levels'
    end if
    this%evaluation_time = evaluation_time
    this%surface_tiles_enabled = physics%radiation%land_sea_mixing_enabled .or. physics%sea_ice%enabled
    this%sea_ice_enabled = physics%sea_ice%enabled
    call copy_surface_field(state%land_temperature, this%land_temperature, 'land_temperature')
    call copy_surface_field(physics_state%land_temperature, this%previous_land_temperature, 'land_temperature')
    call copy_surface_field(state%ocean_temperature, this%ocean_temperature, 'ocean_temperature')
    call copy_surface_field(physics_state%ocean_temperature, this%previous_ocean_temperature, 'ocean_temperature')
    call copy_surface_field(state%sea_ice_fraction, this%sea_ice_fraction, 'sea_ice_fraction')
    call copy_surface_field(physics_state%sea_ice_fraction, this%previous_sea_ice_fraction, 'sea_ice_fraction')
    call copy_surface_field(state%sea_ice_volume, this%sea_ice_volume, 'sea_ice_volume')
    call copy_surface_field(physics_state%sea_ice_volume, this%previous_sea_ice_volume, 'sea_ice_volume')
    this%full_level_eta = coordinate%full_level_eta
    this%physics_grids_ready = physics%held_suarez%enabled .or. physics%surface_friction%enabled .or. &
                               physics%rayleigh_friction%enabled .or. physics%radiation%enabled .or. &
                               physics%convection%enabled .or. physics%moist_convection%enabled .or. &
                               physics%condensation%enabled .or. physics%evaporation%enabled
    this%held_suarez_grids_ready = physics%held_suarez%enabled .or. physics%surface_friction%enabled .or. &
                                   physics%convection%enabled .or. physics%moist_convection%enabled .or. &
                                   physics%condensation%enabled
    this%radiation_grids_ready = physics%radiation%enabled
    this%moisture_grids_ready = physics%moisture%enabled
    this%active_rotation_rate = rotation_rate
    this%bucket_dry_threshold_fraction = physics%bucket%dry_threshold_fraction
    this%land_fraction = 0.0_real64
    if (present(land_fraction)) then
      if (any(shape(land_fraction) /= shape(this%land_fraction))) then
        error stop 'dry tendency land fraction has an inconsistent shape'
      end if
      this%land_fraction = land_fraction
    end if
    this%surface_water = 0.0_real64
    this%previous_surface_water = 0.0_real64
    if (allocated(state%surface_water)) then
      if (any(shape(state%surface_water) /= shape(this%surface_water))) then
        error stop 'dry tendency surface water has an inconsistent shape'
      end if
      this%surface_water = state%surface_water
    end if
    if (allocated(physics_state%surface_water)) then
      if (any(shape(physics_state%surface_water) /= shape(this%previous_surface_water))) then
        error stop 'dry tendency previous surface water has an inconsistent shape'
      end if
      this%previous_surface_water = physics_state%surface_water
    end if

    ! Every level is transformed independently; the scratch arrays are private to the thread.
    !$omp parallel do default(shared) private(k, grid, grid_u, grid_v) schedule(dynamic, 1)
    do k = 1, levels
      call transform%spectral_to_grid(state%zeta(:, :, k), grid)
      this%zeta_grid(:, :, k) = grid
      call transform%spectral_to_grid(state%delta(:, :, k), grid)
      this%delta_grid(:, :, k) = grid
      call transform%spectral_to_grid(state%temperature(:, :, k), grid)
      this%temperature_grid(:, :, k) = grid
      call diagnose_horizontal_velocity(transform, this%truncation, state%zeta(:, :, k), &
                                        state%delta(:, :, k), grid_u, grid_v)
      this%u(:, :, k) = grid_u
      this%v(:, :, k) = grid_v
      if (this%moisture_grids_ready) then
        call transform%spectral_to_grid(state%specific_humidity(:, :, k), grid)
        this%humidity_grid(:, :, k) = grid
        this%virtual_temperature_grid(:, :, k) = (1.0_real64 + virtual_temperature_coefficient* &
          max(this%humidity_grid(:, :, k), 0.0_real64))*this%temperature_grid(:, :, k)
      else
        this%virtual_temperature_grid(:, :, k) = this%temperature_grid(:, :, k)
      end if
      if (this%physics_grids_ready) then
        call transform%spectral_to_grid(physics_state%temperature(:, :, k), grid)
        this%previous_temperature_grid(:, :, k) = grid
        call diagnose_horizontal_velocity(transform, this%truncation, physics_state%zeta(:, :, k), &
                                          physics_state%delta(:, :, k), grid_u, grid_v)
        this%previous_u(:, :, k) = grid_u
        this%previous_v(:, :, k) = grid_v
        if (this%moisture_grids_ready) then
          call transform%spectral_to_grid(physics_state%specific_humidity(:, :, k), grid)
          this%previous_humidity_grid(:, :, k) = grid
        end if
      end if
    end do
    !$omp end parallel do
    call transform%spectral_to_grid(state%log_surface_pressure, this%log_ps)
    call compute_layer_geometry(this%nx, this%ny, levels, this%ring_nlon, coordinate%a_half, coordinate%b_half, &
                                this%log_ps, this%ps, this%pressure_half, this%delta_p, this%layer_l, this%alpha)
    if (this%physics_grids_ready) then
      call transform%spectral_to_grid(physics_state%log_surface_pressure, this%previous_log_ps)
      call compute_layer_geometry(this%nx, this%ny, levels, this%ring_nlon, coordinate%a_half, coordinate%b_half, &
                                  this%previous_log_ps, this%previous_ps, this%previous_pressure_half, &
                                  this%previous_delta_p, this%previous_layer_l, this%previous_alpha, &
                                  full_level_pressure=this%previous_full_level_pressure)
    end if
    ! The surface temperatures are grid prognostic fields; no transform is involved.
    call copy_surface_field(state%surface_temperature, this%surface_temperature_grid, 'surface temperature')
    call copy_surface_field(state%deep_temperature, this%deep_temperature_grid, 'deep temperature')
    call copy_surface_field(physics_state%surface_temperature, this%previous_surface_temperature_grid, &
                            'previous surface temperature')
    call copy_surface_field(physics_state%deep_temperature, this%previous_deep_temperature_grid, &
                            'previous deep temperature')
    call transform%gradient_to_grid(state%log_surface_pressure, this%dlogps_dlambda, this%dlogps_dphi)

    ! The surface geopotential is a fixed lower boundary condition, not a prognostic field.
    ! The hydrostatic integration uses the virtual temperature.
    call transform%spectral_to_grid(surface_geopotential, this%surface_geopotential_grid)
    this%geopotential_half(:, :, levels) = this%surface_geopotential_grid
    do k = levels, 1, -1
      this%geopotential(:, :, k) = this%geopotential_half(:, :, k) + &
        this%alpha(:, :, k)*dry_air_gas_constant*this%virtual_temperature_grid(:, :, k)
      this%geopotential_half(:, :, k - 1) = this%geopotential_half(:, :, k) + &
        dry_air_gas_constant*this%virtual_temperature_grid(:, :, k)*this%layer_l(:, :, k)
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
    if (this%moisture_grids_ready) then
      call compute_vertical_scalar_advection(this%nx, this%ny, levels, this%ring_nlon, this%delta_p, &
                                             this%humidity_grid, this%mass_flux, this%vertical_q)
    end if
  end subroutine prepare_workspace

  !> Copies a grid surface field of the state; an unallocated field reads as zero.
  subroutine copy_surface_field(source, destination, name)
    real(real64), allocatable, intent(in) :: source(:, :)
    real(real64), intent(inout) :: destination(:, :)
    character(*), intent(in) :: name
    if (.not. allocated(source)) then
      destination = 0.0_real64
      return
    end if
    if (any(shape(source) /= shape(destination))) then
      error stop 'dry tendency '//name//' has an inconsistent shape'
    end if
    destination = source
  end subroutine copy_surface_field

  !> Surface pressure, half-level pressures and the layer geometry (dp, ln(p_k/p_{k-1}),
  !> alpha and optionally the full-level pressure) at the grid points of every ring.
  !> The padding beyond ring_nlon(j) is left untouched: it is never transformed or read.
  !> Only these loops evaluate exp/log per grid point, so they stay ring-limited.
  subroutine compute_layer_geometry(nx, ny, levels, ring_nlon, a_half, b_half, log_ps, ps, pressure_half, &
                                    delta_p, layer_l, alpha, full_level_pressure)
    integer, intent(in) :: nx, ny, levels, ring_nlon(ny)
    real(real64), intent(in) :: a_half(0:levels), b_half(0:levels), log_ps(nx, ny)
    real(real64), intent(inout) :: ps(nx, ny), pressure_half(nx, ny, 0:levels)
    real(real64), intent(inout) :: delta_p(nx, ny, levels), layer_l(nx, ny, levels), alpha(nx, ny, levels)
    real(real64), intent(inout), optional :: full_level_pressure(nx, ny, levels)
    integer :: i, j, k
    real(real64) :: thickness, log_ratio
    logical :: with_full_level

    with_full_level = present(full_level_pressure)

    do j = 1, ny
      do i = 1, ring_nlon(j)
        ps(i, j) = exp(log_ps(i, j))
      end do
    end do
    do k = 0, levels
      do j = 1, ny
        do i = 1, ring_nlon(j)
          pressure_half(i, j, k) = a_half(k) + b_half(k)*ps(i, j)
        end do
      end do
    end do
    !$omp parallel do default(shared) private(i, j, k, thickness, log_ratio) schedule(static)
    do k = 1, levels
      do j = 1, ny
        do i = 1, ring_nlon(j)
          thickness = pressure_half(i, j, k) - pressure_half(i, j, k - 1)
          log_ratio = log(pressure_half(i, j, k)/pressure_half(i, j, k - 1))
          delta_p(i, j, k) = thickness
          layer_l(i, j, k) = log_ratio
          alpha(i, j, k) = 1.0_real64 - pressure_half(i, j, k - 1)*log_ratio/thickness
          if (with_full_level) full_level_pressure(i, j, k) = pressure_half(i, j, k)*exp(-alpha(i, j, k))
        end do
      end do
    end do
    !$omp end parallel do
  end subroutine compute_layer_geometry

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

  !> W_k(X) for one more advected scalar (specific humidity), with the same
  !> discretization as the temperature above.
  subroutine compute_vertical_scalar_advection(nx, ny, levels, ring_nlon, delta_p, scalar, mass_flux, vertical)
    integer, intent(in) :: nx, ny, levels, ring_nlon(ny)
    real(real64), intent(in) :: delta_p(nx, ny, levels), scalar(nx, ny, levels), mass_flux(nx, ny, 0:levels)
    real(real64), intent(out) :: vertical(nx, ny, levels)
    integer :: i, j, k

    vertical = 0.0_real64
    do k = 1, levels
      do j = 1, ny
        do i = 1, ring_nlon(j)
          if (k < levels) then
            vertical(i, j, k) = vertical(i, j, k) + mass_flux(i, j, k)*(scalar(i, j, k + 1) - scalar(i, j, k))
          end if
          if (k > 1) then
            vertical(i, j, k) = vertical(i, j, k) + mass_flux(i, j, k - 1)*(scalar(i, j, k) - scalar(i, j, k - 1))
          end if
          vertical(i, j, k) = vertical(i, j, k)/(2.0_real64*delta_p(i, j, k))
        end do
      end do
    end do
  end subroutine compute_vertical_scalar_advection

end module dry_tendency_workspace
