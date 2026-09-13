module raw_filter
  use iso_fortran_env, only: real64
  implicit none
  private

  real(real64), parameter, public :: raw_filter_epsilon = 0.1_real64
  real(real64), parameter, public :: raw_filter_alpha = 0.53_real64

  public :: apply_raw_filter

contains

  subroutine apply_raw_filter(previous_filtered, current, candidate, filtered_current, next)
    complex(real64), intent(in) :: previous_filtered(0:, 0:), current(0:, 0:), candidate(0:, 0:)
    complex(real64), allocatable, intent(out) :: filtered_current(:, :), next(:, :)
    complex(real64), allocatable :: change(:, :)
    integer :: maximum_degree, maximum_order

    if (any(shape(previous_filtered) /= shape(current)) .or. &
        any(shape(previous_filtered) /= shape(candidate))) then
      error stop 'RAW filter fields must have identical shapes'
    end if
    maximum_degree = ubound(previous_filtered, 1)
    maximum_order = ubound(previous_filtered, 2)
    allocate (change(0:maximum_degree, 0:maximum_order))
    allocate (filtered_current(0:maximum_degree, 0:maximum_order))
    allocate (next(0:maximum_degree, 0:maximum_order))
    change = 0.5_real64*raw_filter_epsilon* &
             (previous_filtered - 2.0_real64*current + candidate)
    filtered_current = current + raw_filter_alpha*change
    next = candidate - (1.0_real64 - raw_filter_alpha)*change
  end subroutine apply_raw_filter

end module raw_filter
