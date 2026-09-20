program check_time_filters
  use iso_fortran_env, only: real64
  use numerics_config, only: raw_filter_config, &
                             hyperdiffusion_order => default_hyperdiffusion_order, &
                             hyperdiffusion_timescale_seconds => default_hyperdiffusion_timescale_seconds
  use raw_filter, only: apply_raw_filter
  use spectral_hyperdiffusion, only: apply_spectral_hyperdiffusion
  implicit none

  integer, parameter :: T = 4

  call check_raw_formula()
  call check_custom_raw_config()
  call check_hyperdiffusion_formula()
  call check_custom_hyperdiffusion_timescale()
  call check_custom_hyperdiffusion_order()

contains

  subroutine check_raw_formula()
    complex(real64), allocatable :: previous(:, :), current(:, :), candidate(:, :)
    complex(real64), allocatable :: filtered(:, :), next(:, :)
    real(real64) :: change, expected_filtered, expected_next
    type(raw_filter_config) :: config

    call allocate_test_field(previous)
    call allocate_test_field(current)
    call allocate_test_field(candidate)
    call allocate_test_field(filtered)
    call allocate_test_field(next)
    previous = cmplx(0.0_real64, 0.0_real64, kind=real64)
    current = cmplx(0.0_real64, 0.0_real64, kind=real64)
    candidate = cmplx(0.0_real64, 0.0_real64, kind=real64)
    previous(2, 1) = cmplx(1.0_real64, 0.0_real64, kind=real64)
    current(2, 1) = cmplx(2.0_real64, 0.0_real64, kind=real64)
    candidate(2, 1) = cmplx(4.0_real64, 0.0_real64, kind=real64)

    call apply_raw_filter(previous, current, candidate, filtered, next, config)
    change = 0.5_real64*config%epsilon*(1.0_real64 - 4.0_real64 + 4.0_real64)
    expected_filtered = 2.0_real64 + config%alpha*change
    expected_next = 4.0_real64 - (1.0_real64 - config%alpha)*change
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

  subroutine check_custom_raw_config()
    complex(real64), allocatable :: previous(:, :), current(:, :), candidate(:, :)
    complex(real64), allocatable :: filtered(:, :), next(:, :)
    type(raw_filter_config) :: config

    call allocate_test_field(previous)
    call allocate_test_field(current)
    call allocate_test_field(candidate)
    call allocate_test_field(filtered)
    call allocate_test_field(next)
    previous = cmplx(1.0_real64, 0.0_real64, kind=real64)
    current = cmplx(0.0_real64, 0.0_real64, kind=real64)
    candidate = cmplx(1.0_real64, 0.0_real64, kind=real64)
    config%epsilon = 0.2_real64
    config%alpha = 1.0_real64
    ! change = 0.5*0.2*(1 - 0 + 1) = 0.2; alpha = 1 moves the whole change to the filtered level.
    call apply_raw_filter(previous, current, candidate, filtered, next, config)
    if (any(abs(filtered - cmplx(0.2_real64, 0.0_real64, kind=real64)) > 1.0e-15_real64) .or. &
        any(abs(next - cmplx(1.0_real64, 0.0_real64, kind=real64)) > 1.0e-15_real64)) then
      error stop 'shared RAW filter ignored a custom configuration'
    end if
  end subroutine check_custom_raw_config

  subroutine check_hyperdiffusion_formula()
    complex(real64), allocatable :: field(:, :)
    real(real64) :: ratio, expected

    call allocate_test_field(field)
    field = cmplx(0.0_real64, 0.0_real64, kind=real64)
    field(T, 0) = cmplx(1.0_real64, 0.0_real64, kind=real64)
    field(2, 1) = cmplx(1.0_real64, 0.0_real64, kind=real64)
    call apply_spectral_hyperdiffusion(T, hyperdiffusion_timescale_seconds, field, &
                                       hyperdiffusion_timescale_seconds, real(hyperdiffusion_order, real64))
    if (abs(real(field(T, 0), real64) - 0.5_real64) > 1.0e-15_real64) then
      error stop 'shared hyperdiffusion has the wrong truncation-scale damping'
    end if
    ratio = real(2*3, real64)/real(T*(T + 1), real64)
    expected = 1.0_real64/(1.0_real64 + ratio**hyperdiffusion_order)
    if (abs(real(field(2, 1), real64) - expected) > 1.0e-15_real64) then
      error stop 'shared hyperdiffusion has the wrong modal damping'
    end if
  end subroutine check_hyperdiffusion_formula

  subroutine check_custom_hyperdiffusion_timescale()
    complex(real64), allocatable :: field(:, :)
    real(real64), parameter :: custom_timescale = 3600.0_real64

    call allocate_test_field(field)
    field = cmplx(0.0_real64, 0.0_real64, kind=real64)
    field(T, 0) = cmplx(1.0_real64, 0.0_real64, kind=real64)
    call apply_spectral_hyperdiffusion(T, custom_timescale, field, custom_timescale, &
                                       real(hyperdiffusion_order, real64))
    if (abs(real(field(T, 0), real64) - 0.5_real64) > 1.0e-15_real64) then
      error stop 'shared hyperdiffusion ignored a custom e-folding timescale'
    end if
  end subroutine check_custom_hyperdiffusion_timescale

  subroutine check_custom_hyperdiffusion_order()
    complex(real64), allocatable :: field(:, :)
    real(real64), parameter :: custom_order = 2.5_real64
    real(real64) :: ratio, expected

    call allocate_test_field(field)
    field = cmplx(0.0_real64, 0.0_real64, kind=real64)
    field(2, 1) = cmplx(1.0_real64, 0.0_real64, kind=real64)
    call apply_spectral_hyperdiffusion(T, hyperdiffusion_timescale_seconds, field, &
                                       hyperdiffusion_timescale_seconds, custom_order)
    ratio = real(2*3, real64)/real(T*(T + 1), real64)
    expected = 1.0_real64/(1.0_real64 + ratio**custom_order)
    if (abs(real(field(2, 1), real64) - expected) > 1.0e-15_real64) then
      error stop 'shared hyperdiffusion ignored a custom order'
    end if
  end subroutine check_custom_hyperdiffusion_order

  subroutine allocate_test_field(field)
    complex(real64), allocatable, intent(out) :: field(:, :)
    allocate (field(0:T + 1, 0:T))
  end subroutine allocate_test_field

end program check_time_filters
