module spectral_hyperdiffusion
  use iso_fortran_env, only: real64
  implicit none
  private

  integer, parameter, public :: hyperdiffusion_order = 4
  real(real64), parameter, public :: hyperdiffusion_timescale_seconds = &
    6.0_real64*3600.0_real64

  public :: apply_spectral_hyperdiffusion

contains

  subroutine apply_spectral_hyperdiffusion(truncation, interval, field, timescale_seconds)
    integer, intent(in) :: truncation
    real(real64), intent(in) :: interval
    complex(real64), intent(inout) :: field(0:, 0:)
    real(real64), intent(in), optional :: timescale_seconds
    integer :: n, m
    real(real64) :: damping

    if (truncation < 1) error stop 'spectral hyperdiffusion requires positive T'
    if (interval < 0.0_real64) error stop 'spectral hyperdiffusion interval must be non-negative'
    do m = 0, truncation
      do n = m, truncation
        damping = damping_rate(truncation, n, timescale_seconds)
        field(n, m) = field(n, m)/(1.0_real64 + interval*damping)
      end do
    end do
  end subroutine apply_spectral_hyperdiffusion

  pure real(real64) function damping_rate(truncation, n, timescale_seconds) result(rate)
    integer, intent(in) :: truncation, n
    real(real64), intent(in), optional :: timescale_seconds
    real(real64) :: ratio, timescale

    timescale = hyperdiffusion_timescale_seconds
    if (present(timescale_seconds)) timescale = timescale_seconds
    if (timescale <= 0.0_real64) error stop 'spectral hyperdiffusion timescale must be positive'
    ratio = real(n*(n + 1), real64)/real(truncation*(truncation + 1), real64)
    rate = ratio**hyperdiffusion_order/timescale
  end function damping_rate

end module spectral_hyperdiffusion
