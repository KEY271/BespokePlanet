!> Diagnostic column cloud cover (docs/tendency/cloud.md).
!>
!> Clouds are not a prognostic field.  One effective cloud cover C in [0, 1] is
!> diagnosed per column from the largest layer relative humidity of the
!> provisional field left by the convective adjustments and the large-scale
!> condensation, and from the convective precipitation of the same step.  The
!> grey shortwave radiation reflects C alpha_c of the downward shortwave once;
!> clouds neither absorb shortwave nor enter the grey longwave.
!>
!> The band radiation (docs/tendency/band-radiation.md) splits the column into a
!> clear, a large-scale-cloud and a convective-cloud sub-column.  The convective
!> cloud covers C^conv at the convection top; the large-scale cloud covers the
!> rest of C^RH at the uppermost level of the column-maximum relative humidity.
!> Where the two overlap the convective cloud, the higher one, is kept, so the
!> two fractions add up to the effective cover C = max(C^RH, C^conv).
module cloud_diagnostics
  use iso_fortran_env, only: real64
  use dry_physics_config, only: cloud_config
  use moist_thermodynamics, only: saturation_specific_humidity
  implicit none
  private
  public :: relative_humidity_cloud_cover, convective_cloud_cover, diagnose_cloud_cover, diagnose_cloud_layers

  !> Area fractions and levels of the two cloud sub-columns of one column.  A
  !> level is 0 when its fraction is 0; the clear fraction is the remainder.
  type, public :: cloud_layers
    real(real64) :: large_scale_fraction = 0.0_real64
    real(real64) :: convective_fraction = 0.0_real64
    integer :: large_scale_level = 0
    integer :: convective_level = 0
  end type cloud_layers

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

  !> The effective cover of diagnose_cloud_cover together with its split into a
  !> convective cloud (C^conv at the convection top) and a large-scale cloud
  !> (C^RH - C^conv at the uppermost level of the column-maximum relative
  !> humidity).  `convection_top` is the convection top of the same step, 0 when
  !> there was no precipitating convection; a convective cloud without a top
  !> falls back to the level of the relative-humidity maximum.
  pure subroutine diagnose_cloud_layers(config, full_level_pressure, temperature, humidity, &
                                        convective_precipitation, convection_top, cover, layers)
    type(cloud_config), intent(in) :: config
    real(real64), intent(in) :: full_level_pressure(:), temperature(:), humidity(:), convective_precipitation
    integer, intent(in) :: convection_top
    real(real64), intent(out) :: cover
    type(cloud_layers), intent(out) :: layers
    real(real64) :: relative_humidity(size(temperature)), maximum_relative_humidity
    real(real64) :: relative_humidity_cover, convective_cover
    integer :: k, maximum_level

    do k = 1, size(temperature)
      relative_humidity(k) = min(1.0_real64, max(humidity(k), 0.0_real64)/ &
        saturation_specific_humidity(temperature(k), full_level_pressure(k)))
    end do
    maximum_relative_humidity = maxval(relative_humidity)
    ! Levels run from the top down, so the first maximum is the uppermost one.
    maximum_level = 1
    do k = 1, size(temperature)
      if (relative_humidity(k) >= maximum_relative_humidity) then
        maximum_level = k
        exit
      end if
    end do
    relative_humidity_cover = relative_humidity_cloud_cover(config, maximum_relative_humidity)
    convective_cover = convective_cloud_cover(config, convective_precipitation)
    cover = max(relative_humidity_cover, convective_cover)
    layers = cloud_layers()
    if (convective_cover > 0.0_real64) then
      layers%convective_fraction = convective_cover
      layers%convective_level = convection_top
      if (convection_top < 1 .or. convection_top > size(temperature)) layers%convective_level = maximum_level
    end if
    if (relative_humidity_cover > convective_cover) then
      layers%large_scale_fraction = relative_humidity_cover - convective_cover
      layers%large_scale_level = maximum_level
    end if
  end subroutine diagnose_cloud_layers

end module cloud_diagnostics
