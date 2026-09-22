module raw_filter
  use iso_fortran_env, only: real64
  use numerics_config, only: raw_filter_config
  implicit none
  private

  public :: apply_raw_filter

  interface apply_raw_filter
    module procedure apply_raw_filter_spectral
    module procedure apply_raw_filter_grid
  end interface apply_raw_filter

contains

  !> Robert-Asselin-Williams filter of one spectral field.
  !>
  !> The outputs are written into caller-owned arrays so that steppers can reuse
  !> their work storage instead of allocating new arrays every time step.
  subroutine apply_raw_filter_spectral(previous_filtered, current, candidate, filtered_current, next, config)
    complex(real64), intent(in) :: previous_filtered(0:, 0:), current(0:, 0:), candidate(0:, 0:)
    complex(real64), intent(inout) :: filtered_current(0:, 0:), next(0:, 0:)
    type(raw_filter_config), intent(in) :: config
    complex(real64) :: change
    integer :: n, m

    if (any(shape(previous_filtered) /= shape(current)) .or. &
        any(shape(previous_filtered) /= shape(candidate)) .or. &
        any(shape(previous_filtered) /= shape(filtered_current)) .or. &
        any(shape(previous_filtered) /= shape(next))) then
      error stop 'RAW filter fields must have identical shapes'
    end if
    if (config%epsilon < 0.0_real64 .or. config%alpha < 0.0_real64 .or. config%alpha > 1.0_real64) then
      error stop 'RAW filter configuration is invalid'
    end if
    do m = 0, ubound(previous_filtered, 2)
      do n = 0, ubound(previous_filtered, 1)
        change = 0.5_real64*config%epsilon* &
                 (previous_filtered(n, m) - 2.0_real64*current(n, m) + candidate(n, m))
        filtered_current(n, m) = current(n, m) + config%alpha*change
        next(n, m) = candidate(n, m) - (1.0_real64 - config%alpha)*change
      end do
    end do
  end subroutine apply_raw_filter_spectral

  !> RAW filter for a grid-space surface field.
  subroutine apply_raw_filter_grid(previous_filtered, current, candidate, filtered_current, next, config)
    real(real64), intent(in) :: previous_filtered(:, :), current(:, :), candidate(:, :)
    real(real64), intent(inout) :: filtered_current(:, :), next(:, :)
    type(raw_filter_config), intent(in) :: config
    real(real64) :: change
    integer :: i, j

    if (any(shape(previous_filtered) /= shape(current)) .or. &
        any(shape(previous_filtered) /= shape(candidate)) .or. &
        any(shape(previous_filtered) /= shape(filtered_current)) .or. &
        any(shape(previous_filtered) /= shape(next))) then
      error stop 'RAW filter grid fields must have identical shapes'
    end if
    if (config%epsilon < 0.0_real64 .or. config%alpha < 0.0_real64 .or. config%alpha > 1.0_real64) then
      error stop 'RAW filter configuration is invalid'
    end if
    do j = 1, size(previous_filtered, 2)
      do i = 1, size(previous_filtered, 1)
        change = 0.5_real64*config%epsilon* &
          (previous_filtered(i, j) - 2.0_real64*current(i, j) + candidate(i, j))
        filtered_current(i, j) = current(i, j) + config%alpha*change
        next(i, j) = candidate(i, j) - (1.0_real64 - config%alpha)*change
      end do
    end do
  end subroutine apply_raw_filter_grid

end module raw_filter
