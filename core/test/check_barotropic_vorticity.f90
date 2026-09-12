program check_barotropic_vorticity
  use iso_fortran_env, only: real64, int64
  use harmonics, only: harmonic_transform
  use barotropic_vorticity, only: barotropic_solver
  use barotropic_initial_conditions, only: single_harmonic_vorticity, &
                                            rossby_haurwitz_vorticity, &
                                            random_low_wavenumber_vorticity
  implicit none

  integer, parameter :: T = 12
  type(harmonic_transform) :: transform

  call transform%init(T)
  call check_zero_solution()
  call check_initial_conditions()

contains

  subroutine check_zero_solution()
    type(barotropic_solver) :: solver
    complex(real64), allocatable :: initial(:, :)
    real(real64), allocatable :: zeta(:, :), u(:, :), v(:, :)

    allocate (initial(0:T + 1, 0:T))
    initial = cmplx(0.0_real64, 0.0_real64, kind=real64)
    call solver%init(T, 900.0_real64)
    call solver%set_initial_vorticity(initial)
    call solver%advance()
    call solver%advance()
    call solver%get_fields(zeta, u, v)
    if (maxval(abs(zeta)) > 1.0e-14_real64) error stop 'zero-vorticity state changed'
    if (maxval(abs(u)) > 1.0e-14_real64) error stop 'zero-vorticity zonal wind is nonzero'
    if (maxval(abs(v)) > 1.0e-14_real64) error stop 'zero-vorticity meridional wind is nonzero'
  end subroutine check_zero_solution

  subroutine check_initial_conditions()
    complex(real64), allocatable :: first(:, :), second(:, :)

    call single_harmonic_vorticity(transform, T, first)
    if (count(abs(first) > 0.0_real64) /= 1) then
      error stop 'single-harmonic initial condition contains extra modes'
    end if

    call rossby_haurwitz_vorticity(transform, T, first)
    if (maxval(abs(first)) <= 0.0_real64) error stop 'Rossby-Haurwitz initial condition is zero'

    call random_low_wavenumber_vorticity(transform, T, 20260913_int64, first)
    call random_low_wavenumber_vorticity(transform, T, 20260913_int64, second)
    if (maxval(abs(first - second)) /= 0.0_real64) then
      error stop 'fixed-seed random initial condition is not reproducible'
    end if
    if (maxval(abs(first(0:7, :))) > 0.0_real64 .or. &
        maxval(abs(first(13:T + 1, :))) > 0.0_real64) then
      error stop 'random initial condition escaped n=8:12'
    end if
  end subroutine check_initial_conditions

end program check_barotropic_vorticity
