!> Diagnostic column cloud cover (docs/tendency/cloud.md).
!>
!> Clouds are not a prognostic field.  One effective cloud cover C in [0, 1] is
!> diagnosed per column from the largest layer relative humidity of the
!> provisional field left by the convective adjustments and the large-scale
!> condensation, and from the convective precipitation of the same step.  The
!> shortwave radiation reflects C alpha_c of the downward shortwave once; clouds
!> neither absorb shortwave nor enter the longwave.
module cloud_diagnostics
  use iso_fortran_env, only: real64
  use dry_physics_config, only: cloud_config
  use moist_thermodynamics, only: saturation_specific_humidity
  implicit none
  private
  public :: relative_humidity_cloud_cover, convective_cloud_cover, diagnose_cloud_cover

contains

  !> Slingo (1987) layer-cloud form applied to the column maximum relative humidity:
  !> [(RH_max - RH_c)/(1 - RH_c)]^2, zero below RH_c and one at saturation.
  pure real(real64) function relative_humidity_cloud_cover(config, maximum_relative_humidity) result(cover)
    type(cloud_config), intent(in) :: config
    real(real64), intent(in) :: maximum_relative_humidity
    real(real64) :: excess

    excess = (min(maximum_relative_humidity, 1.0_real64) - config%critical_relative_humidity)/ &
             (1.0_real64 - config%critical_relative_humidity)
    cover = max(excess, 0.0_real64)**2
  end function relative_humidity_cloud_cover

  !> Slingo (1987) convective cloud: c_0 + c_1 ln(P_conv/P_0), clipped to [0, C_max];
  !> zero without convective precipitation.  P_conv is in kg m^-2 s^-1.
  pure real(real64) function convective_cloud_cover(config, convective_precipitation) result(cover)
    type(cloud_config), intent(in) :: config
    real(real64), intent(in) :: convective_precipitation

    if (convective_precipitation <= 0.0_real64) then
      cover = 0.0_real64
      return
    end if
    cover = config%convective_intercept + &
            config%convective_slope*log(convective_precipitation/config%convective_reference_precipitation)
    cover = min(config%convective_maximum_cover, max(cover, 0.0_real64))
  end function convective_cloud_cover

  !> Effective cloud cover of one column: the larger of the relative-humidity
  !> and the convective cloud cover.  `humidity` is the provisional specific
  !> humidity; negative values count as zero.
  pure real(real64) function diagnose_cloud_cover(config, full_level_pressure, temperature, humidity, &
                                                  convective_precipitation) result(cover)
    type(cloud_config), intent(in) :: config
    real(real64), intent(in) :: full_level_pressure(:), temperature(:), humidity(:), convective_precipitation
    real(real64) :: maximum_relative_humidity
    integer :: k

    maximum_relative_humidity = 0.0_real64
    do k = 1, size(temperature)
      maximum_relative_humidity = max(maximum_relative_humidity, max(humidity(k), 0.0_real64)/ &
        saturation_specific_humidity(temperature(k), full_level_pressure(k)))
    end do
    cover = max(relative_humidity_cloud_cover(config, maximum_relative_humidity), &
                convective_cloud_cover(config, convective_precipitation))
  end function diagnose_cloud_cover

end module cloud_diagnostics
