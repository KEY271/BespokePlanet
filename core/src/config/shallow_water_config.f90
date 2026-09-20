!> Parameters of the shallow-water equation system.
module shallow_water_config
  use iso_fortran_env, only: real64
  implicit none
  private

  real(real64), parameter, public :: default_shallow_water_mean_depth = 1.0e4_real64
  !> Idealized test-case gravity; intentionally distinct from planet_parameters.
  real(real64), parameter, public :: default_shallow_water_gravity = 9.8_real64
  real(real64), parameter, public :: default_shallow_water_implicitness = 0.5_real64

  type, public :: shallow_water_equation_config
    real(real64) :: mean_depth = default_shallow_water_mean_depth
    real(real64) :: gravity = default_shallow_water_gravity
    real(real64) :: gravity_wave_implicitness = default_shallow_water_implicitness
  end type shallow_water_equation_config

end module shallow_water_config
