module dry_convection
  use iso_fortran_env, only: real64
  use dry_vertical_coordinate, only: dry_air_kappa, reference_surface_pressure
  use dry_physics_config, only: convection_config
  implicit none
  private

  public :: dry_convective_adjustment_tendency

contains

  subroutine dry_convective_adjustment_tendency(config, pressure_half, temperature, temperature_tendency)
    type(convection_config), intent(in) :: config
    real(real64), intent(in) :: pressure_half(0:), temperature(:)
    real(real64), intent(out) :: temperature_tendency(:)
    real(real64), allocatable :: exner(:), delta_p(:)
    real(real64), allocatable :: block_weight(:), block_enthalpy(:)
    integer, allocatable :: block_top(:), block_bottom(:)
    real(real64) :: layer_log_pressure, alpha, merged_potential_temperature
    integer :: number_of_levels, number_of_blocks, k, block

    number_of_levels = size(temperature)
    if (number_of_levels < 1 .or. size(pressure_half) /= number_of_levels + 1 .or. &
        size(temperature_tendency) /= number_of_levels) then
      error stop 'dry convective adjustment column has inconsistent vertical dimensions'
    end if
    if (any(pressure_half <= 0.0_real64) .or. any(temperature <= 0.0_real64)) then
      error stop 'dry convective adjustment column contains a nonphysical state'
    end if

    allocate (exner(number_of_levels), delta_p(number_of_levels))
    allocate (block_weight(number_of_levels), block_enthalpy(number_of_levels))
    allocate (block_top(number_of_levels), block_bottom(number_of_levels))
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
      do while (number_of_blocks >= 2)
        if (block_enthalpy(number_of_blocks)/block_weight(number_of_blocks) >= &
            block_enthalpy(number_of_blocks - 1)/block_weight(number_of_blocks - 1)) exit
        block_top(number_of_blocks - 1) = block_top(number_of_blocks)
        block_weight(number_of_blocks - 1) = block_weight(number_of_blocks - 1) + &
                                                block_weight(number_of_blocks)
        block_enthalpy(number_of_blocks - 1) = block_enthalpy(number_of_blocks - 1) + &
                                                  block_enthalpy(number_of_blocks)
        number_of_blocks = number_of_blocks - 1
      end do
    end do

    temperature_tendency = 0.0_real64
    do block = 1, number_of_blocks
      if (block_top(block) == block_bottom(block)) cycle
      merged_potential_temperature = block_enthalpy(block)/block_weight(block)
      do k = block_top(block), block_bottom(block)
        temperature_tendency(k) = -(temperature(k) - exner(k)*merged_potential_temperature)/ &
                                    config%adjustment_time
      end do
    end do
  end subroutine dry_convective_adjustment_tendency

end module dry_convection
