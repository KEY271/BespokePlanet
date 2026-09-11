module lapack_interfaces
  use iso_fortran_env, only: real64
  implicit none
  private

  public :: dstev

  interface
    subroutine dstev(jobz, n, d, e, z, ldz, work, info)
      import :: real64
      implicit none

      character(len=1), intent(in) :: jobz
      integer, intent(in) :: n
      integer, intent(in) :: ldz

      real(real64), intent(inout) :: d(*)
      real(real64), intent(inout) :: e(*)
      real(real64), intent(out)   :: z(ldz, *)
      real(real64), intent(out)   :: work(*)

      integer, intent(out) :: info
    end subroutine dstev
  end interface
end module lapack_interfaces
