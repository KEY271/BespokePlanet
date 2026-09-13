module shallow_water_gravity_wave
  use iso_fortran_env, only: real64
  use barotropic_vorticity, only: earth_radius
  implicit none
  private

  real(real64), parameter, public :: gravity_acceleration = 9.8_real64
  real(real64), parameter, public :: mean_depth = 1.0e4_real64
  real(real64), parameter, public :: gravity_wave_implicitness = 0.5_real64

  public :: solve_implicit_gravity_wave

contains

  subroutine solve_implicit_gravity_wave(truncation, centered_interval, &
                                         previous_delta, previous_eta, &
                                         nonlinear_delta, nonlinear_eta, &
                                         next_delta, next_eta)
    integer, intent(in) :: truncation
    real(real64), intent(in) :: centered_interval
    complex(real64), intent(in) :: previous_delta(0:, 0:), previous_eta(0:, 0:)
    complex(real64), intent(in) :: nonlinear_delta(0:, 0:), nonlinear_eta(0:, 0:)
    complex(real64), intent(out) :: next_delta(0:, 0:), next_eta(0:, 0:)
    complex(real64) :: delta_rhs, eta_rhs
    integer :: n, m
    real(real64) :: gravity_laplacian, coupling, determinant

    next_delta = cmplx(0.0_real64, 0.0_real64, kind=real64)
    next_eta = cmplx(0.0_real64, 0.0_real64, kind=real64)
    coupling = centered_interval*gravity_wave_implicitness
    do m = 0, truncation
      do n = m, truncation
        gravity_laplacian = gravity_acceleration*real(n*(n + 1), real64)/earth_radius**2
        determinant = 1.0_real64 + coupling**2*mean_depth*gravity_laplacian
        delta_rhs = previous_delta(n, m) + centered_interval*(nonlinear_delta(n, m) + &
                    (1.0_real64 - gravity_wave_implicitness)*gravity_laplacian*previous_eta(n, m))
        eta_rhs = previous_eta(n, m) + centered_interval*(nonlinear_eta(n, m) - &
                  (1.0_real64 - gravity_wave_implicitness)*mean_depth*previous_delta(n, m))
        next_delta(n, m) = (delta_rhs + coupling*gravity_laplacian*eta_rhs)/determinant
        next_eta(n, m) = (eta_rhs - coupling*mean_depth*delta_rhs)/determinant
      end do
    end do
  end subroutine solve_implicit_gravity_wave

end module shallow_water_gravity_wave
