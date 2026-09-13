module fftw_interfaces
  use iso_c_binding, only: c_double, c_double_complex, c_int, c_ptr
  implicit none
  private

  integer(c_int), parameter, public :: FFTW_UNALIGNED = 2_c_int
  integer(c_int), parameter, public :: FFTW_ESTIMATE = 64_c_int

  public :: fftw_plan_dft_r2c_1d
  public :: fftw_plan_dft_c2r_1d
  public :: fftw_execute_dft_r2c
  public :: fftw_execute_dft_c2r
  public :: fftw_destroy_plan

  interface
    type(c_ptr) function fftw_plan_dft_r2c_1d(n, input, output, flags) &
      bind(C, name='fftw_plan_dft_r2c_1d')
      import
      integer(c_int), value :: n
      real(c_double), dimension(*), intent(out) :: input
      complex(c_double_complex), dimension(*), intent(out) :: output
      integer(c_int), value :: flags
    end function fftw_plan_dft_r2c_1d

    type(c_ptr) function fftw_plan_dft_c2r_1d(n, input, output, flags) &
      bind(C, name='fftw_plan_dft_c2r_1d')
      import
      integer(c_int), value :: n
      complex(c_double_complex), dimension(*), intent(out) :: input
      real(c_double), dimension(*), intent(out) :: output
      integer(c_int), value :: flags
    end function fftw_plan_dft_c2r_1d

    subroutine fftw_execute_dft_r2c(plan, input, output) bind(C, name='fftw_execute_dft_r2c')
      import
      type(c_ptr), value :: plan
      real(c_double), dimension(*), intent(inout) :: input
      complex(c_double_complex), dimension(*), intent(out) :: output
    end subroutine fftw_execute_dft_r2c

    subroutine fftw_execute_dft_c2r(plan, input, output) bind(C, name='fftw_execute_dft_c2r')
      import
      type(c_ptr), value :: plan
      complex(c_double_complex), dimension(*), intent(inout) :: input
      real(c_double), dimension(*), intent(out) :: output
    end subroutine fftw_execute_dft_c2r

    subroutine fftw_destroy_plan(plan) bind(C, name='fftw_destroy_plan')
      import
      type(c_ptr), value :: plan
    end subroutine fftw_destroy_plan
  end interface
end module fftw_interfaces
