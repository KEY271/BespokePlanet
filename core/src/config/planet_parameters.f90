!> Physical constants of the planet shared by every equation system.
module planet_parameters
  use iso_fortran_env, only: real64
  implicit none
  private

  real(real64), parameter, public :: earth_radius = 6.371e6_real64
  real(real64), parameter, public :: earth_rotation_rate = 7.2921159e-5_real64
  !> Gravity used by the hydrostatic dry atmosphere.  The shallow-water system
  !> deliberately keeps its own idealized value in shallow_water_config.
  real(real64), parameter, public :: earth_gravity = 9.80616_real64

  !> Planet values a case may select per run.  The radius is a fixed constant
  !> because the spectral operators scale with it; only the rotation rate varies
  !> between the cases at present.
  type, public :: planet_config
    real(real64) :: rotation_rate = earth_rotation_rate
  end type planet_config

end module planet_parameters
