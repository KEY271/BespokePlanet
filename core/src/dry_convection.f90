!> Enthalpy-conserving dry convective adjustment
!> (docs/tendency/dry-convective-adjustment.md), with the moist-atmosphere
!> variant that judges stability by virtual potential temperature and mixes
!> specific humidity within each block
!> (docs/tendency/moist-convective-adjustment.md, "dry adjustment changes").
module dry_convection
  use iso_fortran_env, only: real64
  use dry_vertical_coordinate, only: dry_air_kappa, reference_surface_pressure
  use dry_physics_config, only: convection_config
  use moist_thermodynamics, only: virtual_temperature_coefficient
  implicit none
  private

  public :: dry_convective_adjustment_tendency

contains

  !> With `humidity` present the block potential temperature is the virtual one,
  !> theta_vc = (1 + delta_v Q/dP) H/W, and the humidity is mixed uniformly within
  !> each block; with q = 0 this reduces to the dry procedure.
  subroutine dry_convective_adjustment_tendency(config, pressure_half, temperature, temperature_tendency, &
                                                humidity, humidity_tendency)
    type(convection_config), intent(in) :: config
    real(real64), intent(in) :: pressure_half(0:), temperature(:)
    real(real64), intent(out) :: temperature_tendency(:)
    real(real64), intent(in), optional :: humidity(:)
    real(real64), intent(out), optional :: humidity_tendency(:)
    real(real64), allocatable :: exner(:), delta_p(:), clipped_humidity(:)
    real(real64), allocatable :: block_weight(:), block_enthalpy(:), block_water(:), block_mass(:)
    integer, allocatable :: block_top(:), block_bottom(:)
    real(real64) :: layer_log_pressure, alpha, merged_potential_temperature, merged_humidity
    integer :: number_of_levels, number_of_blocks, k, block
    logical :: moist

    number_of_levels = size(temperature)
    moist = present(humidity)
    if (moist .neqv. present(humidity_tendency)) then
      error stop 'dry convective adjustment humidity and its tendency must be supplied together'
    end if
    if (number_of_levels < 1 .or. size(pressure_half) /= number_of_levels + 1 .or. &
        size(temperature_tendency) /= number_of_levels) then
      error stop 'dry convective adjustment column has inconsistent vertical dimensions'
    end if
    if (moist) then
      if (size(humidity) /= number_of_levels .or. size(humidity_tendency) /= number_of_levels) then
        error stop 'dry convective adjustment humidity has inconsistent vertical dimensions'
      end if
    end if
    if (any(pressure_half <= 0.0_real64) .or. any(temperature <= 0.0_real64)) then
      error stop 'dry convective adjustment column contains a nonphysical state'
    end if

    allocate (exner(number_of_levels), delta_p(number_of_levels), clipped_humidity(number_of_levels))
    allocate (block_weight(number_of_levels), block_enthalpy(number_of_levels))
    allocate (block_water(number_of_levels), block_mass(number_of_levels))
    allocate (block_top(number_of_levels), block_bottom(number_of_levels))
    clipped_humidity = 0.0_real64
    if (moist) clipped_humidity = max(humidity, 0.0_real64)
    do k = 1, number_of_levels
      delta_p(k) = pressure_half(k) - pressure_half(k - 1)
      if (delta_p(k) <= 0.0_real64) then
        error stop 'dry convective adjustment pressures must increase downward'
      end if
      layer_log_pressure = log(pressure_half(k)/pressure_half(k - 1))
      alpha = 1.0_real64 - pressure_half(k - 1)*layer_log_pressure/delta_p(k)
      exner(k) = (pressure_half(k)*exp(-alpha)/reference_surface_pressure)**dry_air_kappa
    end do

    ! Pool adjacent unstable blocks while walking from the surface upward.  The
    ! stack is ordered from lower to upper atmosphere, so its last element is
    ! the newly added upper block.
    number_of_blocks = 0
    do k = number_of_levels, 1, -1
      number_of_blocks = number_of_blocks + 1
      block_top(number_of_blocks) = k
      block_bottom(number_of_blocks) = k
      block_weight(number_of_blocks) = exner(k)*delta_p(k)
      block_enthalpy(number_of_blocks) = temperature(k)*delta_p(k)
      block_water(number_of_blocks) = clipped_humidity(k)*delta_p(k)
      block_mass(number_of_blocks) = delta_p(k)
      do while (number_of_blocks >= 2)
        if (block_potential_temperature(number_of_blocks) >= block_potential_temperature(number_of_blocks - 1)) exit
        block_top(number_of_blocks - 1) = block_top(number_of_blocks)
        block_weight(number_of_blocks - 1) = block_weight(number_of_blocks - 1) + &
                                                block_weight(number_of_blocks)
        block_enthalpy(number_of_blocks - 1) = block_enthalpy(number_of_blocks - 1) + &
                                                  block_enthalpy(number_of_blocks)
        block_water(number_of_blocks - 1) = block_water(number_of_blocks - 1) + block_water(number_of_blocks)
        block_mass(number_of_blocks - 1) = block_mass(number_of_blocks - 1) + block_mass(number_of_blocks)
        number_of_blocks = number_of_blocks - 1
      end do
    end do

    temperature_tendency = 0.0_real64
    if (moist) humidity_tendency = 0.0_real64
    do block = 1, number_of_blocks
      if (block_top(block) == block_bottom(block)) cycle
      merged_potential_temperature = block_enthalpy(block)/block_weight(block)
      merged_humidity = block_water(block)/block_mass(block)
      do k = block_top(block), block_bottom(block)
        temperature_tendency(k) = -(temperature(k) - exner(k)*merged_potential_temperature)/ &
                                    config%adjustment_time
        if (moist) humidity_tendency(k) = -(clipped_humidity(k) - merged_humidity)/config%adjustment_time
      end do
    end do

  contains

    !> Block (virtual) potential temperature used for the stability comparison.
    real(real64) function block_potential_temperature(index) result(theta)
      integer, intent(in) :: index
      theta = block_enthalpy(index)/block_weight(index)
      if (moist) theta = (1.0_real64 + virtual_temperature_coefficient*block_water(index)/block_mass(index))*theta
    end function block_potential_temperature

  end subroutine dry_convective_adjustment_tendency

end module dry_convection
