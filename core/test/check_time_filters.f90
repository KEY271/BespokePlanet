program check_time_filters
  use iso_fortran_env, only: real64
  use raw_filter, only: raw_filter_epsilon, raw_filter_alpha, apply_raw_filter
  use spectral_hyperdiffusion, only: hyperdiffusion_order, &
                                       hyperdiffusion_timescale_seconds, &
                                       apply_spectral_hyperdiffusion
  implicit none

  integer, parameter :: T = 4

  call check_raw_formula()
  call check_hyperdiffusion_formula()

contains

  subroutine check_raw_formula()
    complex(real64), allocatable :: previous(:, :), current(:, :), candidate(:, :)
    complex(real64), allocatable :: filtered(:, :), next(:, :)
    real(real64) :: change, expected_filtered, expected_next

    call allocate_test_field(previous)
    call allocate_test_field(current)
    call allocate_test_field(candidate)
    previous = cmplx(0.0_real64, 0.0_real64, kind=real64)
    current = cmplx(0.0_real64, 0.0_real64, kind=real64)
    candidate = cmplx(0.0_real64, 0.0_real64, kind=real64)
    previous(2, 1) = cmplx(1.0_real64, 0.0_real64, kind=real64)
    current(2, 1) = cmplx(2.0_real64, 0.0_real64, kind=real64)
    candidate(2, 1) = cmplx(4.0_real64, 0.0_real64, kind=real64)

    call apply_raw_filter(previous, current, candidate, filtered, next)
    change = 0.5_real64*raw_filter_epsilon*(1.0_real64 - 4.0_real64 + 4.0_real64)
    expected_filtered = 2.0_real64 + raw_filter_alpha*change
    expected_next = 4.0_real64 - (1.0_real64 - raw_filter_alpha)*change
    if (abs(real(filtered(2, 1), real64) - expected_filtered) > 1.0e-15_real64) then
      error stop 'shared RAW filter returned the wrong filtered value'
    end if
    if (abs(real(next(2, 1), real64) - expected_next) > 1.0e-15_real64) then
      error stop 'shared RAW filter returned the wrong next value'
    end if
    if (lbound(filtered, 1) /= 0 .or. lbound(filtered, 2) /= 0) then
      error stop 'shared RAW filter did not preserve spectral index bounds'
    end if
  end subroutine check_raw_formula

  subroutine check_hyperdiffusion_formula()
    complex(real64), allocatable :: field(:, :)
    real(real64) :: ratio, expected

    call allocate_test_field(field)
    field = cmplx(0.0_real64, 0.0_real64, kind=real64)
    field(T, 0) = cmplx(1.0_real64, 0.0_real64, kind=real64)
    field(2, 1) = cmplx(1.0_real64, 0.0_real64, kind=real64)
    call apply_spectral_hyperdiffusion(T, hyperdiffusion_timescale_seconds, field)
    if (abs(real(field(T, 0), real64) - 0.5_real64) > 1.0e-15_real64) then
      error stop 'shared hyperdiffusion has the wrong truncation-scale damping'
    end if
    ratio = real(2*3, real64)/real(T*(T + 1), real64)
    expected = 1.0_real64/(1.0_real64 + ratio**hyperdiffusion_order)
    if (abs(real(field(2, 1), real64) - expected) > 1.0e-15_real64) then
      error stop 'shared hyperdiffusion has the wrong modal damping'
    end if
  end subroutine check_hyperdiffusion_formula

  subroutine allocate_test_field(field)
    complex(real64), allocatable, intent(out) :: field(:, :)
    allocate (field(0:T + 1, 0:T))
  end subroutine allocate_test_field

end program check_time_filters
