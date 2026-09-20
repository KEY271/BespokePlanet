module spectral_hyperdiffusion
  use iso_fortran_env, only: real64
  implicit none
  private

  public :: apply_spectral_hyperdiffusion

contains

  !> Implicit hyperdiffusion damping.  Coefficients are supplied by the caller's
  !> numerics configuration; this utility owns no default values.
  subroutine apply_spectral_hyperdiffusion(truncation, interval, field, timescale_seconds, diffusion_order)
    integer, intent(in) :: truncation
    real(real64), intent(in) :: interval
    complex(real64), intent(inout) :: field(0:, 0:)
    real(real64), intent(in) :: timescale_seconds
    real(real64), intent(in) :: diffusion_order
    integer :: n, m
    real(real64) :: damping

    if (truncation < 1) error stop 'spectral hyperdiffusion requires positive T'
    if (interval < 0.0_real64) error stop 'spectral hyperdiffusion interval must be non-negative'
    if (timescale_seconds <= 0.0_real64) error stop 'spectral hyperdiffusion timescale must be positive'
    if (diffusion_order <= 0.0_real64) error stop 'spectral hyperdiffusion order must be positive'
    do m = 0, truncation
      do n = m, truncation
        damping = damping_rate(truncation, n, timescale_seconds, diffusion_order)
        field(n, m) = field(n, m)/(1.0_real64 + interval*damping)
      end do
    end do
  end subroutine apply_spectral_hyperdiffusion

  pure real(real64) function damping_rate(truncation, n, timescale_seconds, diffusion_order) result(rate)
    integer, intent(in) :: truncation, n
    real(real64), intent(in) :: timescale_seconds
    real(real64), intent(in) :: diffusion_order
    real(real64) :: ratio

    ratio = real(n*(n + 1), real64)/real(truncation*(truncation + 1), real64)
    rate = ratio**diffusion_order/timescale_seconds
  end function damping_rate

end module spectral_hyperdiffusion
