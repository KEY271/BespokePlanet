!> Prognostic state of the hydrostatic atmosphere.  Vorticity, divergence,
!> temperature, specific humidity and log surface pressure are spectral fields
!> on every level.  The surface (ground or ocean) temperature, the deep-ground
!> temperature and the land-bucket water have no horizontal dynamics and are
!> carried on the grid, so that the land--sea contrast in their tendencies never
!> passes through a spectral truncation (docs/tendency/land-sea-surface.md).
!> Specific humidity is always carried; it stays zero when the atmosphere is dry.
module dry_state
  use iso_fortran_env, only: real64
  implicit none
  private

  type, public :: dry_state_type
    complex(real64), allocatable :: zeta(:, :, :)
    complex(real64), allocatable :: delta(:, :, :)
    complex(real64), allocatable :: temperature(:, :, :)
    complex(real64), allocatable :: specific_humidity(:, :, :)
    complex(real64), allocatable :: log_surface_pressure(:, :)
    !> Grid-space surface fields (allocate_dry_surface_fields): surface and
    !> deep-ground temperature in K, land-bucket water in kg m^-2 of land area.
    real(real64), allocatable :: surface_temperature(:, :)
    real(real64), allocatable :: deep_temperature(:, :)
    real(real64), allocatable :: surface_water(:, :)
    real(real64), allocatable :: land_temperature(:, :)
    real(real64), allocatable :: ocean_temperature(:, :)
    real(real64), allocatable :: sea_ice_fraction(:, :)
    real(real64), allocatable :: sea_ice_volume(:, :)
  end type dry_state_type

  type, public :: dry_tendency_type
    complex(real64), allocatable :: zeta(:, :, :)
    complex(real64), allocatable :: delta(:, :, :)
    complex(real64), allocatable :: temperature(:, :, :)
    complex(real64), allocatable :: specific_humidity(:, :, :)
    complex(real64), allocatable :: log_surface_pressure(:, :)
    !> Grid-space tendencies of the surface temperatures (K s^-1).
    real(real64), allocatable :: surface_temperature(:, :)
    real(real64), allocatable :: deep_temperature(:, :)
  end type dry_tendency_type

  public :: allocate_dry_state, allocate_dry_surface_fields, allocate_dry_tendency, zero_dry_tendency
  public :: copy_dry_state, swap_dry_states
  public :: enforce_dry_state_constraints, enforce_dry_spectral_field

contains

  !> Allocates the spectral fields.  The grid surface fields are released here and
  !> allocated separately with allocate_dry_surface_fields, which needs the grid shape.
  subroutine allocate_dry_state(state, truncation, number_of_levels)
    type(dry_state_type), intent(inout) :: state
    integer, intent(in) :: truncation, number_of_levels
    if (allocated(state%zeta)) then
      deallocate (state%zeta, state%delta, state%temperature, state%specific_humidity, state%log_surface_pressure)
    end if
    if (allocated(state%surface_temperature)) deallocate (state%surface_temperature)
    if (allocated(state%deep_temperature)) deallocate (state%deep_temperature)
    if (allocated(state%surface_water)) deallocate (state%surface_water)
    if (allocated(state%land_temperature)) deallocate (state%land_temperature)
    if (allocated(state%ocean_temperature)) deallocate (state%ocean_temperature)
    if (allocated(state%sea_ice_fraction)) deallocate (state%sea_ice_fraction)
    if (allocated(state%sea_ice_volume)) deallocate (state%sea_ice_volume)
    allocate (state%zeta(0:truncation + 1, 0:truncation, number_of_levels))
    allocate (state%delta(0:truncation + 1, 0:truncation, number_of_levels))
    allocate (state%temperature(0:truncation + 1, 0:truncation, number_of_levels))
    allocate (state%specific_humidity(0:truncation + 1, 0:truncation, number_of_levels))
    allocate (state%log_surface_pressure(0:truncation + 1, 0:truncation))
    state%zeta = cmplx(0.0_real64, 0.0_real64, kind=real64)
    state%delta = cmplx(0.0_real64, 0.0_real64, kind=real64)
    state%temperature = cmplx(0.0_real64, 0.0_real64, kind=real64)
    state%specific_humidity = cmplx(0.0_real64, 0.0_real64, kind=real64)
    state%log_surface_pressure = cmplx(0.0_real64, 0.0_real64, kind=real64)
  end subroutine allocate_dry_state

  !> Allocates the grid surface fields (surface temperature, deep temperature and
  !> bucket water) with the padded grid shape of the transform and zeroes them.
  subroutine allocate_dry_surface_fields(state, nx, ny)
    type(dry_state_type), intent(inout) :: state
    integer, intent(in) :: nx, ny
    if (nx < 1 .or. ny < 1) error stop 'dry surface field shape must be positive'
    if (allocated(state%surface_temperature)) deallocate (state%surface_temperature)
    if (allocated(state%deep_temperature)) deallocate (state%deep_temperature)
    if (allocated(state%surface_water)) deallocate (state%surface_water)
    if (allocated(state%land_temperature)) deallocate (state%land_temperature)
    if (allocated(state%ocean_temperature)) deallocate (state%ocean_temperature)
    if (allocated(state%sea_ice_fraction)) deallocate (state%sea_ice_fraction)
    if (allocated(state%sea_ice_volume)) deallocate (state%sea_ice_volume)
    allocate (state%surface_temperature(nx, ny), state%deep_temperature(nx, ny), state%surface_water(nx, ny))
    state%surface_temperature = 0.0_real64
    state%deep_temperature = 0.0_real64
    state%surface_water = 0.0_real64
    allocate (state%land_temperature(nx, ny))
    state%land_temperature = 0.0_real64
    allocate (state%ocean_temperature(nx, ny))
    state%ocean_temperature = 0.0_real64
    allocate (state%sea_ice_fraction(nx, ny))
    state%sea_ice_fraction = 0.0_real64
    allocate (state%sea_ice_volume(nx, ny))
    state%sea_ice_volume = 0.0_real64
  end subroutine allocate_dry_surface_fields

  subroutine allocate_dry_tendency(tendency, truncation, number_of_levels, nx, ny)
    type(dry_tendency_type), intent(inout) :: tendency
    integer, intent(in) :: truncation, number_of_levels, nx, ny
    if (nx < 1 .or. ny < 1) error stop 'dry tendency grid shape must be positive'
    if (allocated(tendency%zeta)) then
      deallocate (tendency%zeta, tendency%delta, tendency%temperature, tendency%specific_humidity)
      deallocate (tendency%log_surface_pressure, tendency%surface_temperature, tendency%deep_temperature)
    end if
    allocate (tendency%zeta(0:truncation + 1, 0:truncation, number_of_levels))
    allocate (tendency%delta(0:truncation + 1, 0:truncation, number_of_levels))
    allocate (tendency%temperature(0:truncation + 1, 0:truncation, number_of_levels))
    allocate (tendency%specific_humidity(0:truncation + 1, 0:truncation, number_of_levels))
    allocate (tendency%log_surface_pressure(0:truncation + 1, 0:truncation))
    allocate (tendency%surface_temperature(nx, ny))
    allocate (tendency%deep_temperature(nx, ny))
    call zero_dry_tendency(tendency)
  end subroutine allocate_dry_tendency

  subroutine zero_dry_tendency(tendency)
    type(dry_tendency_type), intent(inout) :: tendency
    if (.not. allocated(tendency%zeta)) error stop 'dry tendency is not allocated'
    tendency%zeta = cmplx(0.0_real64, 0.0_real64, kind=real64)
    tendency%delta = cmplx(0.0_real64, 0.0_real64, kind=real64)
    tendency%temperature = cmplx(0.0_real64, 0.0_real64, kind=real64)
    tendency%specific_humidity = cmplx(0.0_real64, 0.0_real64, kind=real64)
    tendency%log_surface_pressure = cmplx(0.0_real64, 0.0_real64, kind=real64)
    tendency%surface_temperature = 0.0_real64
    tendency%deep_temperature = 0.0_real64
  end subroutine zero_dry_tendency

  !> Copies values into already allocated storage of the same shape.  Grid surface
  !> fields that the source carries are allocated in the destination when absent.
  subroutine copy_dry_state(source, destination)
    type(dry_state_type), intent(in) :: source
    type(dry_state_type), intent(inout) :: destination
    if (any(shape(source%zeta) /= shape(destination%zeta))) error stop 'dry state copy shape mismatch'
    destination%zeta(:, :, :) = source%zeta
    destination%delta(:, :, :) = source%delta
    destination%temperature(:, :, :) = source%temperature
    destination%specific_humidity(:, :, :) = source%specific_humidity
    destination%log_surface_pressure(:, :) = source%log_surface_pressure
    call copy_grid_field(source%surface_temperature, destination%surface_temperature, 'surface temperature')
    call copy_grid_field(source%deep_temperature, destination%deep_temperature, 'deep temperature')
    call copy_grid_field(source%surface_water, destination%surface_water, 'surface water')
    call copy_grid_field(source%land_temperature, destination%land_temperature, 'land_temperature')
    call copy_grid_field(source%ocean_temperature, destination%ocean_temperature, 'ocean_temperature')
    call copy_grid_field(source%sea_ice_fraction, destination%sea_ice_fraction, 'sea_ice_fraction')
    call copy_grid_field(source%sea_ice_volume, destination%sea_ice_volume, 'sea_ice_volume')
  end subroutine copy_dry_state

  subroutine copy_grid_field(source, destination, name)
    real(real64), allocatable, intent(in) :: source(:, :)
    real(real64), allocatable, intent(inout) :: destination(:, :)
    character(*), intent(in) :: name
    if (.not. allocated(source)) return
    if (.not. allocated(destination)) then
      allocate (destination, mold=source)
    else if (any(shape(source) /= shape(destination))) then
      error stop 'dry '//name//' copy shape mismatch'
    end if
    destination = source
  end subroutine copy_grid_field

  !> Exchanges storage without copying array data.
  subroutine swap_dry_states(first, second)
    type(dry_state_type), intent(inout) :: first, second
    call swap_level_field(first%zeta, second%zeta)
    call swap_level_field(first%delta, second%delta)
    call swap_level_field(first%temperature, second%temperature)
    call swap_level_field(first%specific_humidity, second%specific_humidity)
    call swap_surface_field(first%log_surface_pressure, second%log_surface_pressure)
    call swap_grid_field(first%surface_temperature, second%surface_temperature)
    call swap_grid_field(first%deep_temperature, second%deep_temperature)
    call swap_grid_field(first%surface_water, second%surface_water)
    call swap_grid_field(first%land_temperature, second%land_temperature)
    call swap_grid_field(first%ocean_temperature, second%ocean_temperature)
    call swap_grid_field(first%sea_ice_fraction, second%sea_ice_fraction)
    call swap_grid_field(first%sea_ice_volume, second%sea_ice_volume)
  end subroutine swap_dry_states

  subroutine swap_level_field(first, second)
    complex(real64), allocatable, intent(inout) :: first(:, :, :), second(:, :, :)
    complex(real64), allocatable :: temporary(:, :, :)
    call move_alloc(first, temporary)
    call move_alloc(second, first)
    call move_alloc(temporary, second)
  end subroutine swap_level_field

  subroutine swap_surface_field(first, second)
    complex(real64), allocatable, intent(inout) :: first(:, :), second(:, :)
    complex(real64), allocatable :: temporary(:, :)
    call move_alloc(first, temporary)
    call move_alloc(second, first)
    call move_alloc(temporary, second)
  end subroutine swap_surface_field

  subroutine swap_grid_field(first, second)
    real(real64), allocatable, intent(inout) :: first(:, :), second(:, :)
    real(real64), allocatable :: temporary(:, :)
    call move_alloc(first, temporary)
    call move_alloc(second, first)
    call move_alloc(temporary, second)
  end subroutine swap_grid_field

  !> Zeroes the global mean of vorticity and divergence and the unused spectral
  !> entries of every spectral field.  The global mean of specific humidity is not
  !> constrained, and the grid surface fields have no spectral entries to constrain.
  subroutine enforce_dry_state_constraints(state, truncation)
    type(dry_state_type), intent(inout) :: state
    integer, intent(in) :: truncation
    integer :: k
    do k = 1, size(state%zeta, 3)
      call enforce_dry_spectral_field(state%zeta(:, :, k), truncation, .true.)
      call enforce_dry_spectral_field(state%delta(:, :, k), truncation, .true.)
      call enforce_dry_spectral_field(state%temperature(:, :, k), truncation, .false.)
      call enforce_dry_spectral_field(state%specific_humidity(:, :, k), truncation, .false.)
    end do
    call enforce_dry_spectral_field(state%log_surface_pressure, truncation, .false.)
  end subroutine enforce_dry_state_constraints

  subroutine enforce_dry_spectral_field(field, truncation, zero_mean)
    complex(real64), intent(inout) :: field(0:, 0:)
    integer, intent(in) :: truncation
    logical, intent(in) :: zero_mean
    integer :: m, n
    if (zero_mean) field(0, 0) = cmplx(0.0_real64, 0.0_real64, kind=real64)
    field(truncation + 1, :) = cmplx(0.0_real64, 0.0_real64, kind=real64)
    do m = 1, truncation
      do n = 0, m - 1
        field(n, m) = cmplx(0.0_real64, 0.0_real64, kind=real64)
      end do
    end do
  end subroutine enforce_dry_spectral_field

end module dry_state
