program check_shallow_water
  use iso_fortran_env, only: real64
  use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
  use harmonics, only: harmonic_transform
  use shallow_water, only: shallow_water_solver
  use shallow_water_initial_conditions, only: isolated_height_mountain, &
                                                single_harmonic_height, &
                                                height_mode_degree, height_mode_order, &
                                                height_mode_amplitude_metres
  implicit none

  integer, parameter :: T = 8
  real(real64), parameter :: dt = 300.0_real64
  type(harmonic_transform) :: transform

  call transform%init(T)
  call check_zero_solution()
  call check_single_harmonic()
  call check_mountain_evolves()

contains

  subroutine check_zero_solution()
    type(shallow_water_solver) :: solver
    complex(real64), allocatable :: initial_zeta(:, :), initial_delta(:, :), initial_eta(:, :)
    complex(real64), allocatable :: zeta_spec(:, :), delta_spec(:, :), eta_spec(:, :)
    real(real64), allocatable :: zeta(:, :), delta(:, :), eta(:, :), u(:, :), v(:, :)

    call allocate_zero_state(initial_zeta, initial_delta, initial_eta)
    call solver%init(T, dt)
    call solver%set_initial_state(initial_zeta, initial_delta, initial_eta)
    call solver%advance()
    call solver%advance()
    call solver%get_spectral_state(zeta_spec, delta_spec, eta_spec)
    call solver%get_fields(zeta, delta, eta, u, v)
    if (maxval(abs(zeta_spec)) > 1.0e-14_real64) error stop 'zero-state zeta changed'
    if (maxval(abs(delta_spec)) > 1.0e-14_real64) error stop 'zero-state delta changed'
    if (maxval(abs(eta_spec)) > 1.0e-12_real64) error stop 'zero-state eta changed'
    if (maxval(abs(u)) > 1.0e-12_real64 .or. maxval(abs(v)) > 1.0e-12_real64) then
      error stop 'zero-state velocity changed'
    end if
  end subroutine check_zero_solution

  subroutine check_single_harmonic()
    complex(real64), allocatable :: zeta(:, :), delta(:, :), eta(:, :)
    real(real64), allocatable :: eta_grid(:, :)
    integer, allocatable :: nlon(:)
    integer :: j
    real(real64) :: maximum

    call single_harmonic_height(transform, T, zeta, delta, eta)
    if (maxval(abs(zeta)) /= 0.0_real64 .or. maxval(abs(delta)) /= 0.0_real64) then
      error stop 'single height mode initialized velocity fields'
    end if
    if (count(abs(eta) > 0.0_real64) /= 1) error stop 'single height mode contains extra modes'
    if (abs(eta(height_mode_degree, height_mode_order)) <= 0.0_real64) then
      error stop 'requested single height mode is absent'
    end if
    call transform%spectral_to_grid(eta, eta_grid)
    nlon = transform%get_nlon()
    maximum = 0.0_real64
    do j = 1, size(nlon)
      maximum = max(maximum, maxval(abs(eta_grid(1:nlon(j), j))))
    end do
    if (abs(maximum - height_mode_amplitude_metres) > 1.0e-10_real64) then
      error stop 'single height mode has the wrong amplitude'
    end if
  end subroutine check_single_harmonic

  subroutine check_mountain_evolves()
    type(shallow_water_solver) :: solver
    complex(real64), allocatable :: zeta(:, :), delta(:, :), eta(:, :), initial_eta(:, :)
    real(real64), allocatable :: zeta_grid(:, :), delta_grid(:, :), eta_grid(:, :), u(:, :), v(:, :)
    real(real64) :: cfl

    call isolated_height_mountain(transform, T, zeta, delta, eta)
    initial_eta = eta
    if (abs(eta(0, 0)) > 0.0_real64) error stop 'mountain height has a nonzero global mean mode'
    call solver%init(T, dt)
    call solver%set_initial_state(zeta, delta, eta)
    call solver%advance()
    call solver%get_spectral_state(zeta, delta, eta)
    call solver%get_fields(zeta_grid, delta_grid, eta_grid, u, v, cfl)
    if (maxval(abs(delta)) <= 1.0e-14_real64) error stop 'height mountain did not launch a gravity wave'
    if (maxval(abs(eta - initial_eta)) <= 1.0e-10_real64) error stop 'height mountain did not evolve'
    if (abs(zeta(0, 0)) > 0.0_real64 .or. abs(delta(0, 0)) > 0.0_real64 .or. &
        abs(eta(0, 0)) > 0.0_real64) then
      error stop 'mean spectral constraints were not enforced'
    end if
    if (.not. all(ieee_is_finite(zeta_grid)) .or. .not. all(ieee_is_finite(delta_grid)) .or. &
        .not. all(ieee_is_finite(eta_grid)) .or. .not. all(ieee_is_finite(u)) .or. &
        .not. all(ieee_is_finite(v)) .or. .not. ieee_is_finite(cfl)) then
      error stop 'shallow-water evolution produced a non-finite value'
    end if
  end subroutine check_mountain_evolves

  subroutine allocate_zero_state(zeta, delta, eta)
    complex(real64), allocatable, intent(out) :: zeta(:, :), delta(:, :), eta(:, :)

    allocate (zeta(0:T + 1, 0:T), delta(0:T + 1, 0:T), eta(0:T + 1, 0:T))
    zeta = cmplx(0.0_real64, 0.0_real64, kind=real64)
    delta = cmplx(0.0_real64, 0.0_real64, kind=real64)
    eta = cmplx(0.0_real64, 0.0_real64, kind=real64)
  end subroutine allocate_zero_state

end program check_shallow_water
