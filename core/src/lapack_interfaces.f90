module lapack_interfaces
  use iso_fortran_env, only: real64
  implicit none
  private

  public :: dstev, dgetrf, dgetri

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

    subroutine dgetrf(m, n, a, lda, ipiv, info)
      import :: real64
      implicit none
      integer, intent(in) :: m, n, lda
      real(real64), intent(inout) :: a(lda, *)
      integer, intent(out) :: ipiv(*)
      integer, intent(out) :: info
    end subroutine dgetrf

    subroutine dgetri(n, a, lda, ipiv, work, lwork, info)
      import :: real64
      implicit none
      integer, intent(in) :: n, lda, lwork
      real(real64), intent(inout) :: a(lda, *)
      integer, intent(in) :: ipiv(*)
      real(real64), intent(out) :: work(*)
      integer, intent(out) :: info
    end subroutine dgetri
  end interface
end module lapack_interfaces
