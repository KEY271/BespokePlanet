program main
  use iso_fortran_env, only: real64
  use harmonics, only: init_harmonics
  implicit none

  integer :: T = 63

  call init_harmonics(T)
end program main
