!> Spectral prognostic state of the hydrostatic atmosphere: vorticity, divergence,
!> temperature, specific humidity and log surface pressure on every level, plus
!> the surface (ground or ocean) and deep-ground temperatures.  Specific
!> humidity is always carried; it stays zero when the atmosphere is dry.
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
    complex(real64), allocatable :: surface_temperature(:, :)
    complex(real64), allocatable :: deep_temperature(:, :)
    !> Land-bucket water in grid space (kg m^-2 of land area).  It has no
    !> horizontal dynamics and is therefore not represented spectrally.
    real(real64), allocatable :: surface_water(:, :)
  end type dry_state_type

  type, public :: dry_tendency_type
    complex(real64), allocatable :: zeta(:, :, :)
    complex(real64), allocatable :: delta(:, :, :)
    complex(real64), allocatable :: temperature(:, :, :)
    complex(real64), allocatable :: specific_humidity(:, :, :)
    complex(real64), allocatable :: log_surface_pressure(:, :)
    complex(real64), allocatable :: surface_temperature(:, :)
    complex(real64), allocatable :: deep_temperature(:, :)
  end type dry_tendency_type

  public :: allocate_dry_state, allocate_dry_surface_water, allocate_dry_tendency, zero_dry_tendency
  public :: copy_dry_state, swap_dry_states
  public :: enforce_dry_state_constraints, enforce_dry_spectral_field

contains

  subroutine allocate_dry_state(state, truncation, number_of_levels)
    type(dry_state_type), intent(inout) :: state
    integer, intent(in) :: truncation, number_of_levels
    if (allocated(state%zeta)) then
      deallocate (state%zeta, state%delta, state%temperature, state%specific_humidity, state%log_surface_pressure)
      deallocate (state%surface_temperature, state%deep_temperature)
    end if
    if (allocated(state%surface_water)) deallocate (state%surface_water)
    allocate (state%zeta(0:truncation + 1, 0:truncation, number_of_levels))
    allocate (state%delta(0:truncation + 1, 0:truncation, number_of_levels))
    allocate (state%temperature(0:truncation + 1, 0:truncation, number_of_levels))
    allocate (state%specific_humidity(0:truncation + 1, 0:truncation, number_of_levels))
    allocate (state%log_surface_pressure(0:truncation + 1, 0:truncation))
    allocate (state%surface_temperature(0:truncation + 1, 0:truncation))
    allocate (state%deep_temperature(0:truncation + 1, 0:truncation))
    state%zeta = cmplx(0.0_real64, 0.0_real64, kind=real64)
    state%delta = cmplx(0.0_real64, 0.0_real64, kind=real64)
    state%temperature = cmplx(0.0_real64, 0.0_real64, kind=real64)
    state%specific_humidity = cmplx(0.0_real64, 0.0_real64, kind=real64)
    state%log_surface_pressure = cmplx(0.0_real64, 0.0_real64, kind=real64)
    state%surface_temperature = cmplx(0.0_real64, 0.0_real64, kind=real64)
    state%deep_temperature = cmplx(0.0_real64, 0.0_real64, kind=real64)
  end subroutine allocate_dry_state

  subroutine allocate_dry_surface_water(state, nx, ny)
    type(dry_state_type), intent(inout) :: state
    integer, intent(in) :: nx, ny
    if (nx < 1 .or. ny < 1) error stop 'dry surface-water shape must be positive'
    if (allocated(state%surface_water)) deallocate (state%surface_water)
    allocate (state%surface_water(nx, ny))
    state%surface_water = 0.0_real64
  end subroutine allocate_dry_surface_water

  subroutine allocate_dry_tendency(tendency, truncation, number_of_levels)
    type(dry_tendency_type), intent(inout) :: tendency
    integer, intent(in) :: truncation, number_of_levels
    if (allocated(tendency%zeta)) then
      deallocate (tendency%zeta, tendency%delta, tendency%temperature, tendency%specific_humidity)
      deallocate (tendency%log_surface_pressure, tendency%surface_temperature, tendency%deep_temperature)
    end if
    allocate (tendency%zeta(0:truncation + 1, 0:truncation, number_of_levels))
    allocate (tendency%delta(0:truncation + 1, 0:truncation, number_of_levels))
    allocate (tendency%temperature(0:truncation + 1, 0:truncation, number_of_levels))
    allocate (tendency%specific_humidity(0:truncation + 1, 0:truncation, number_of_levels))
    allocate (tendency%log_surface_pressure(0:truncation + 1, 0:truncation))
    allocate (tendency%surface_temperature(0:truncation + 1, 0:truncation))
    allocate (tendency%deep_temperature(0:truncation + 1, 0:truncation))
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
    tendency%surface_temperature = cmplx(0.0_real64, 0.0_real64, kind=real64)
    tendency%deep_temperature = cmplx(0.0_real64, 0.0_real64, kind=real64)
  end subroutine zero_dry_tendency

  !> Copies values into already allocated storage of the same shape.
  subroutine copy_dry_state(source, destination)
    type(dry_state_type), intent(in) :: source
    type(dry_state_type), intent(inout) :: destination
    if (any(shape(source%zeta) /= shape(destination%zeta))) error stop 'dry state copy shape mismatch'
    destination%zeta(:, :, :) = source%zeta
    destination%delta(:, :, :) = source%delta
    destination%temperature(:, :, :) = source%temperature
    destination%specific_humidity(:, :, :) = source%specific_humidity
    destination%log_surface_pressure(:, :) = source%log_surface_pressure
    destination%surface_temperature(:, :) = source%surface_temperature
    destination%deep_temperature(:, :) = source%deep_temperature
    if (allocated(source%surface_water)) then
      if (.not. allocated(destination%surface_water)) then
        allocate (destination%surface_water, mold=source%surface_water)
      else if (any(shape(source%surface_water) /= shape(destination%surface_water))) then
        error stop 'dry surface-water copy shape mismatch'
      end if
      destination%surface_water = source%surface_water
    end if
  end subroutine copy_dry_state

  !> Exchanges storage without copying array data.
  subroutine swap_dry_states(first, second)
    type(dry_state_type), intent(inout) :: first, second
    call swap_level_field(first%zeta, second%zeta)
    call swap_level_field(first%delta, second%delta)
    call swap_level_field(first%temperature, second%temperature)
    call swap_level_field(first%specific_humidity, second%specific_humidity)
    call swap_surface_field(first%log_surface_pressure, second%log_surface_pressure)
    call swap_surface_field(first%surface_temperature, second%surface_temperature)
    call swap_surface_field(first%deep_temperature, second%deep_temperature)
    call swap_grid_field(first%surface_water, second%surface_water)
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
  !> entries of every field.  The global mean of specific humidity is not constrained.
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
    call enforce_dry_spectral_field(state%surface_temperature, truncation, .false.)
    call enforce_dry_spectral_field(state%deep_temperature, truncation, .false.)
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
