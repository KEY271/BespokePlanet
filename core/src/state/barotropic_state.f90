module barotropic_state
  use iso_fortran_env, only: real64
  implicit none
  private

  type, public :: barotropic_state_type
    complex(real64), allocatable :: zeta(:, :)
  end type barotropic_state_type

  type, public :: barotropic_tendency_type
    complex(real64), allocatable :: zeta(:, :)
  end type barotropic_tendency_type

  public :: allocate_barotropic_state, allocate_barotropic_tendency
  public :: zero_barotropic_tendency, enforce_barotropic_constraints
  public :: copy_barotropic_state, swap_barotropic_states

contains

  subroutine allocate_barotropic_state(state, truncation)
    type(barotropic_state_type), intent(inout) :: state
    integer, intent(in) :: truncation
    if (allocated(state%zeta)) deallocate (state%zeta)
    allocate (state%zeta(0:truncation + 1, 0:truncation))
    state%zeta = cmplx(0.0_real64, 0.0_real64, kind=real64)
  end subroutine allocate_barotropic_state

  subroutine allocate_barotropic_tendency(tendency, truncation)
    type(barotropic_tendency_type), intent(inout) :: tendency
    integer, intent(in) :: truncation
    if (allocated(tendency%zeta)) deallocate (tendency%zeta)
    allocate (tendency%zeta(0:truncation + 1, 0:truncation))
    call zero_barotropic_tendency(tendency)
  end subroutine allocate_barotropic_tendency

  subroutine zero_barotropic_tendency(tendency)
    type(barotropic_tendency_type), intent(inout) :: tendency
    if (.not. allocated(tendency%zeta)) error stop 'barotropic tendency is not allocated'
    tendency%zeta = cmplx(0.0_real64, 0.0_real64, kind=real64)
  end subroutine zero_barotropic_tendency

  !> Copies values into already allocated storage of the same shape.
  subroutine copy_barotropic_state(source, destination)
    type(barotropic_state_type), intent(in) :: source
    type(barotropic_state_type), intent(inout) :: destination
    if (any(shape(source%zeta) /= shape(destination%zeta))) error stop 'barotropic state copy shape mismatch'
    destination%zeta(:, :) = source%zeta
  end subroutine copy_barotropic_state

  !> Exchanges storage without copying array data.
  subroutine swap_barotropic_states(first, second)
    type(barotropic_state_type), intent(inout) :: first, second
    complex(real64), allocatable :: temporary(:, :)
    call move_alloc(first%zeta, temporary)
    call move_alloc(second%zeta, first%zeta)
    call move_alloc(temporary, second%zeta)
  end subroutine swap_barotropic_states

  subroutine enforce_barotropic_constraints(field, truncation)
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
  end subroutine enforce_barotropic_constraints

end module barotropic_state
