module shallow_water_state
  use iso_fortran_env, only: real64
  implicit none
  private

  type, public :: shallow_water_state_type
    complex(real64), allocatable :: zeta(:, :)
    complex(real64), allocatable :: delta(:, :)
    complex(real64), allocatable :: eta(:, :)
  end type shallow_water_state_type

  type, public :: shallow_water_tendency_type
    complex(real64), allocatable :: zeta(:, :)
    complex(real64), allocatable :: delta(:, :)
    complex(real64), allocatable :: eta(:, :)
  end type shallow_water_tendency_type

  public :: allocate_shallow_water_state, allocate_shallow_water_tendency
  public :: zero_shallow_water_tendency, enforce_shallow_water_constraints
  public :: copy_shallow_water_state, swap_shallow_water_states

contains

  subroutine allocate_shallow_water_state(state, truncation)
    type(shallow_water_state_type), intent(inout) :: state
    integer, intent(in) :: truncation
    if (allocated(state%zeta)) deallocate (state%zeta, state%delta, state%eta)
    allocate (state%zeta(0:truncation + 1, 0:truncation))
    allocate (state%delta(0:truncation + 1, 0:truncation))
    allocate (state%eta(0:truncation + 1, 0:truncation))
    state%zeta = cmplx(0.0_real64, 0.0_real64, kind=real64)
    state%delta = cmplx(0.0_real64, 0.0_real64, kind=real64)
    state%eta = cmplx(0.0_real64, 0.0_real64, kind=real64)
  end subroutine allocate_shallow_water_state

  subroutine allocate_shallow_water_tendency(tendency, truncation)
    type(shallow_water_tendency_type), intent(inout) :: tendency
    integer, intent(in) :: truncation
    if (allocated(tendency%zeta)) deallocate (tendency%zeta, tendency%delta, tendency%eta)
    allocate (tendency%zeta(0:truncation + 1, 0:truncation))
    allocate (tendency%delta(0:truncation + 1, 0:truncation))
    allocate (tendency%eta(0:truncation + 1, 0:truncation))
    call zero_shallow_water_tendency(tendency)
  end subroutine allocate_shallow_water_tendency

  subroutine zero_shallow_water_tendency(tendency)
    type(shallow_water_tendency_type), intent(inout) :: tendency
    if (.not. allocated(tendency%zeta)) error stop 'shallow-water tendency is not allocated'
    tendency%zeta = cmplx(0.0_real64, 0.0_real64, kind=real64)
    tendency%delta = cmplx(0.0_real64, 0.0_real64, kind=real64)
    tendency%eta = cmplx(0.0_real64, 0.0_real64, kind=real64)
  end subroutine zero_shallow_water_tendency

  !> Copies values into already allocated storage of the same shape.
  subroutine copy_shallow_water_state(source, destination)
    type(shallow_water_state_type), intent(in) :: source
    type(shallow_water_state_type), intent(inout) :: destination
    if (any(shape(source%zeta) /= shape(destination%zeta))) error stop 'shallow-water state copy shape mismatch'
    destination%zeta(:, :) = source%zeta
    destination%delta(:, :) = source%delta
    destination%eta(:, :) = source%eta
  end subroutine copy_shallow_water_state

  !> Exchanges storage without copying array data.
  subroutine swap_shallow_water_states(first, second)
    type(shallow_water_state_type), intent(inout) :: first, second
    call swap_field(first%zeta, second%zeta)
    call swap_field(first%delta, second%delta)
    call swap_field(first%eta, second%eta)
  end subroutine swap_shallow_water_states

  subroutine swap_field(first, second)
    complex(real64), allocatable, intent(inout) :: first(:, :), second(:, :)
    complex(real64), allocatable :: temporary(:, :)
    call move_alloc(first, temporary)
    call move_alloc(second, first)
    call move_alloc(temporary, second)
  end subroutine swap_field

  subroutine enforce_shallow_water_constraints(state, truncation)
    type(shallow_water_state_type), intent(inout) :: state
    integer, intent(in) :: truncation
    call enforce_field(state%zeta, truncation)
    call enforce_field(state%delta, truncation)
    call enforce_field(state%eta, truncation)
  end subroutine enforce_shallow_water_constraints

  subroutine enforce_field(field, truncation)
    complex(real64), intent(inout) :: field(0:, 0:)
    integer, intent(in) :: truncation
    integer :: m, n
    field(0, 0) = cmplx(0.0_real64, 0.0_real64, kind=real64)
    field(truncation + 1, :) = cmplx(0.0_real64, 0.0_real64, kind=real64)
    do m = 1, truncation
      do n = 0, m - 1
        field(n, m) = cmplx(0.0_real64, 0.0_real64, kind=real64)
      end do
    end do
  end subroutine enforce_field

end module shallow_water_state
