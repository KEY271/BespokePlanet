!> Bounded one-layer land-water reservoir (docs/tendency/bucket.md).
module land_bucket
  use iso_fortran_env, only: real64
  implicit none
  private
  public :: bucket_wetness, limit_bucket_evaporation, advance_bucket

contains

  pure real(real64) function bucket_wetness(water, capacity) result(wetness)
    real(real64), intent(in) :: water, capacity
    if (capacity <= 0.0_real64) then
      wetness = 0.0_real64
    else
      wetness = min(1.0_real64, max(0.0_real64, water/capacity))
    end if
  end function bucket_wetness

  pure real(real64) function limit_bucket_evaporation(potential_flux, water, capacity, interval) result(flux)
    real(real64), intent(in) :: potential_flux, water, capacity, interval
    flux = bucket_wetness(water, capacity)*potential_flux
    if (flux > 0.0_real64 .and. interval > 0.0_real64) then
      flux = min(flux, max(water, 0.0_real64)/interval)
    end if
  end function limit_bucket_evaporation

  pure subroutine advance_bucket(water, precipitation, evaporation, capacity, interval, tendency, runoff)
    real(real64), intent(in) :: water, precipitation, evaporation, capacity, interval
    real(real64), intent(out) :: tendency, runoff
    real(real64) :: bounded_water, candidate, next_water

    if (interval <= 0.0_real64 .or. capacity <= 0.0_real64) then
      tendency = 0.0_real64
      runoff = 0.0_real64
      return
    end if
    bounded_water = min(capacity, max(0.0_real64, water))
    candidate = bounded_water + interval*(precipitation - evaporation)
    runoff = max(candidate - capacity, 0.0_real64)/interval
    next_water = min(capacity, max(0.0_real64, candidate))
    tendency = (next_water - bounded_water)/interval
  end subroutine advance_bucket

end module land_bucket
