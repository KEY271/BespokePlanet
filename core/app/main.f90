program main
  use harmonics, only: harmonic_transform
  implicit none

  integer :: T = 63
  type(harmonic_transform) :: transform

  call transform%init(T)
end program main
