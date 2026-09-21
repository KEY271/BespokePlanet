!> Numerical-method settings and their defaults.
!>
!> The named constants below are the only place where the default numerical
!> coefficients are written.  The configuration types use them as default
!> initializers, and the integration utilities receive configuration values
!> explicitly instead of carrying their own copies.
module numerics_config
  use iso_fortran_env, only: real64
  implicit none
  private

  real(real64), parameter, public :: default_raw_filter_epsilon = 0.1_real64
  real(real64), parameter, public :: default_raw_filter_alpha = 0.53_real64
  integer, parameter, public :: default_hyperdiffusion_order = 4
  real(real64), parameter, public :: default_hyperdiffusion_timescale_seconds = 6.0_real64*3600.0_real64
  real(real64), parameter, public :: default_dry_vorticity_diffusion_seconds = 4.0_real64*3600.0_real64
  real(real64), parameter, public :: default_dry_divergence_diffusion_seconds = 1.0_real64*3600.0_real64
  real(real64), parameter, public :: default_dry_temperature_diffusion_seconds = 4.0_real64*3600.0_real64
  real(real64), parameter, public :: default_dry_humidity_diffusion_seconds = 4.0_real64*3600.0_real64

  type, public :: raw_filter_config
    real(real64) :: epsilon = default_raw_filter_epsilon
    real(real64) :: alpha = default_raw_filter_alpha
  end type raw_filter_config

  !> Hyperdiffusion with one e-folding time at the truncation wavenumber.
  type, public :: hyperdiffusion_config
    integer :: order = default_hyperdiffusion_order
    real(real64) :: timescale_seconds = default_hyperdiffusion_timescale_seconds
  end type hyperdiffusion_config

  !> Dry-atmosphere hyperdiffusion, whose e-folding times differ by prognostic variable.
  type, public :: dry_hyperdiffusion_config
    integer :: order = default_hyperdiffusion_order
    real(real64) :: vorticity_timescale_seconds = default_dry_vorticity_diffusion_seconds
    real(real64) :: divergence_timescale_seconds = default_dry_divergence_diffusion_seconds
    real(real64) :: temperature_timescale_seconds = default_dry_temperature_diffusion_seconds
    !> Specific humidity is damped like temperature (moist atmosphere only).
    real(real64) :: humidity_timescale_seconds = default_dry_humidity_diffusion_seconds
  end type dry_hyperdiffusion_config

  type, public :: model_numerics_config
    integer :: truncation = -1
    real(real64) :: time_step = 0.0_real64
    type(raw_filter_config) :: raw_filter
    type(hyperdiffusion_config) :: hyperdiffusion
  end type model_numerics_config

end module numerics_config
