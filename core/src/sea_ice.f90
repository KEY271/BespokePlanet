!> Local snow-free sea-ice thermodynamics (docs/tendency/sea-ice.md).
!> A and V are per ocean area; fluxes are per exposed tile area.
module sea_ice
  use iso_fortran_env, only: real64
  use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
  use dry_physics_config, only: sea_ice_config
  implicit none
  private
  public :: validate_sea_ice_config, solve_ice_surface, advance_sea_ice, equilibrate_sea_ice

  type, public :: sea_ice_budget
    real(real64) :: open_water_volume = 0.0_real64
    real(real64) :: surface_volume = 0.0_real64
    real(real64) :: basal_volume = 0.0_real64
    real(real64) :: excess_heat = 0.0_real64
    real(real64) :: energy_residual = 0.0_real64
  end type sea_ice_budget

  !> Maxima over active ocean cells; projection counts exclude rounding-only changes.
  type, public :: sea_ice_checks
    real(real64) :: energy_residual = 0.0_real64
    real(real64) :: surface_residual = 0.0_real64
    real(real64) :: projection_area = 0.0_real64
    real(real64) :: projection_volume = 0.0_real64
    real(real64) :: projection_temperature = 0.0_real64
    real(real64) :: projection_energy_residual = 0.0_real64
    integer :: projected_cells = 0
    integer :: checked_cells = 0
  contains
    procedure :: merge => merge_ice_checks
  end type sea_ice_checks

contains

  subroutine merge_ice_checks(this, other)
    class(sea_ice_checks), intent(inout) :: this
    type(sea_ice_checks), intent(in) :: other
    this%energy_residual = max(this%energy_residual, other%energy_residual)
    this%surface_residual = max(this%surface_residual, other%surface_residual)
    this%projection_area = max(this%projection_area, other%projection_area)
    this%projection_volume = max(this%projection_volume, other%projection_volume)
    this%projection_temperature = max(this%projection_temperature, other%projection_temperature)
    this%projection_energy_residual = max(this%projection_energy_residual, other%projection_energy_residual)
    this%projected_cells = this%projected_cells + other%projected_cells
    this%checked_cells = this%checked_cells + other%checked_cells
  end subroutine merge_ice_checks

  subroutine validate_sea_ice_config(config)
    type(sea_ice_config), intent(in) :: config
    real(real64) :: values(9)
    values = [config%freezing_temperature, config%melting_temperature, config%new_ice_thickness, &
              config%conductivity, config%density, config%latent_heat, config%albedo, &
              config%temperature_tolerance, config%flux_tolerance]
    if (.not. all(ieee_is_finite(values))) error stop 'nonfinite sea-ice configuration'
    if (any(values(1:6) <= 0.0_real64) .or. config%melting_temperature < config%freezing_temperature .or. &
        config%albedo < 0.0_real64 .or. config%albedo > 1.0_real64 .or. &
        any(values(8:9) <= 0.0_real64) .or. config%maximum_iterations < 1) then
      error stop 'invalid sea-ice configuration'
    end if
  end subroutine validate_sea_ice_config

  !> Solve h F(T) + k (Tf-T) = 0, or fix T at the surface melting point.
  !> absorbed_heat includes absorbed SW + downward LW; sensible = B*(T-Ta).
  subroutine solve_ice_surface(config, thickness, absorbed_heat, sigma, exchange, air_temperature, &
                               temperature, conductive_heat, melt_heat, residual, iterations)
    type(sea_ice_config), intent(in) :: config
    real(real64), intent(in) :: thickness, absorbed_heat, sigma, exchange, air_temperature
    real(real64), intent(out) :: temperature, conductive_heat, melt_heat
    real(real64), intent(out), optional :: residual
    integer, intent(out), optional :: iterations
    real(real64) :: lower, upper, value, flux, tolerance, roundoff
    integer :: iteration

    if (.not. all(ieee_is_finite([thickness, absorbed_heat, sigma, exchange, air_temperature])) .or. &
        thickness <= 0.0_real64 .or. sigma <= 0.0_real64 .or. exchange < 0.0_real64) then
      error stop 'invalid ice surface forcing'
    end if
    lower = 0.0_real64
    upper = config%melting_temperature
    temperature = upper
    flux = surface_flux(temperature)
    value = thickness*flux + config%conductivity*(config%freezing_temperature - temperature)
    if (value > 0.0_real64) then
      conductive_heat = config%conductivity*(config%freezing_temperature - temperature)/thickness
      melt_heat = flux + conductive_heat
      if (present(residual)) residual = 0.0_real64
      if (present(iterations)) iterations = 0
      return
    end if
    if (thickness*surface_flux(lower) + config%conductivity*config%freezing_temperature < 0.0_real64) then
      error stop 'ice surface root is not bracketed'
    end if
    do iteration = 1, config%maximum_iterations
      temperature = lower + 0.5_real64*(upper - lower)
      flux = surface_flux(temperature)
      value = thickness*flux + config%conductivity*(config%freezing_temperature - temperature)
      roundoff = 16.0_real64*epsilon(temperature)* &
        (abs(thickness*flux) + config%conductivity*config%freezing_temperature)
      tolerance = thickness*config%flux_tolerance + roundoff
      if (upper - lower <= config%temperature_tolerance .and. abs(value) <= tolerance) exit
      if (value > 0.0_real64) then
        lower = temperature
      else
        upper = temperature
      end if
    end do
    if (iteration > config%maximum_iterations .or. .not. ieee_is_finite(value)) then
      error stop 'ice surface temperature failed to converge'
    end if
    ! Avoid loss of significance in (Tf-T)/h for arbitrarily thin ice.
    conductive_heat = -flux
    melt_heat = 0.0_real64
    if (present(residual)) residual = value
    if (present(iterations)) iterations = iteration

  contains
    real(real64) function surface_flux(t) result(f)
      real(real64), intent(in) :: t
      f = absorbed_heat - sigma*t**4 - exchange*(t - air_temperature)
    end function surface_flux
  end subroutine solve_ice_surface

  subroutine advance_sea_ice(config, heat_capacity, interval, ocean_temperature, area, volume, &
                             open_water_flux, conductive_heat, melt_heat, new_temperature, new_area, new_volume, budget)
    type(sea_ice_config), intent(in) :: config
    real(real64), intent(in) :: heat_capacity, interval, ocean_temperature, area, volume
    real(real64), intent(in) :: open_water_flux, conductive_heat, melt_heat
    real(real64), intent(out) :: new_temperature, new_area, new_volume
    type(sea_ice_budget), intent(out) :: budget
    real(real64) :: latent, heat, freezing_heat, added_area, delta_volume, trial_volume, initial_energy

    if (.not. all(ieee_is_finite([heat_capacity, interval, ocean_temperature, area, volume, &
                                  open_water_flux, conductive_heat, melt_heat]))) error stop 'nonfinite sea-ice update'
    if (heat_capacity <= 0.0_real64 .or. interval <= 0.0_real64 .or. area < 0.0_real64 .or. &
        area > 1.0_real64 .or. volume < 0.0_real64 .or. melt_heat < 0.0_real64) error stop 'invalid sea-ice update'
    if ((area == 0.0_real64) .neqv. (volume == 0.0_real64)) error stop 'inconsistent sea-ice area and volume'
    latent = config%density*config%latent_heat
    budget = sea_ice_budget()
    heat = interval*(1.0_real64 - area)*open_water_flux
    ! Algebraically Co*(To* - Tf), without cancellation after the small warming step.
    freezing_heat = heat_capacity*(ocean_temperature - config%freezing_temperature) + heat
    new_temperature = ocean_temperature + heat/heat_capacity
    added_area = 0.0_real64
    if (freezing_heat < 0.0_real64 .or. volume > 0.0_real64) then
      budget%open_water_volume = -freezing_heat/latent
      new_temperature = config%freezing_temperature
      if (freezing_heat < 0.0_real64) then
        added_area = min(1.0_real64 - area, budget%open_water_volume/config%new_ice_thickness)
      end if
    end if
    if (volume > 0.0_real64) then
      budget%surface_volume = -interval*area*melt_heat/latent
      budget%basal_volume = interval*area*conductive_heat/latent
    end if
    delta_volume = budget%open_water_volume + budget%surface_volume + budget%basal_volume
    trial_volume = volume + delta_volume
    if (trial_volume <= 0.0_real64) then
      new_volume = 0.0_real64
      new_area = 0.0_real64
      budget%excess_heat = -latent*trial_volume
      new_temperature = new_temperature + budget%excess_heat/heat_capacity
    else
      new_volume = trial_volume
      new_area = area + added_area
      if (delta_volume < 0.0_real64 .and. volume > 0.0_real64) then
        new_area = new_area*sqrt(trial_volume/volume)
      end if
      new_area = min(1.0_real64, max(0.0_real64, new_area))
    end if
    initial_energy = heat_capacity*(ocean_temperature - config%freezing_temperature) - latent*volume
    budget%energy_residual = heat_capacity*(new_temperature - config%freezing_temperature) - &
      latent*new_volume - initial_energy - heat - interval*area*(melt_heat - conductive_heat)
  end subroutine advance_sea_ice

  !> Restore phase equilibrium after RAW without discarding sensible or latent heat.
  subroutine equilibrate_sea_ice(config, heat_capacity, temperature, area, volume, checks)
    type(sea_ice_config), intent(in) :: config
    real(real64), intent(in) :: heat_capacity
    real(real64), intent(inout) :: temperature, area, volume
    type(sea_ice_checks), intent(inout), optional :: checks
    real(real64) :: energy, latent, old_temperature, old_area, old_volume
    if (.not. all(ieee_is_finite([temperature, area, volume, heat_capacity]))) then
      error stop 'nonfinite sea-ice state'
    end if
    if (heat_capacity <= 0.0_real64) error stop 'invalid ocean heat capacity'
    latent = config%density*config%latent_heat
    old_temperature = temperature
    old_area = area
    old_volume = volume
    energy = heat_capacity*(temperature - config%freezing_temperature) - latent*volume
    if (energy >= 0.0_real64) then
      temperature = config%freezing_temperature + energy/heat_capacity
      area = 0.0_real64
      volume = 0.0_real64
    else
      temperature = config%freezing_temperature
      volume = -energy/latent
      area = min(1.0_real64, max(0.0_real64, area))
      if (area == 0.0_real64) area = min(1.0_real64, volume/config%new_ice_thickness)
    end if
    if (present(checks)) then
      checks%projection_area = max(checks%projection_area, abs(area - old_area))
      checks%projection_volume = max(checks%projection_volume, abs(volume - old_volume))
      checks%projection_temperature = max(checks%projection_temperature, abs(temperature - old_temperature))
      checks%projection_energy_residual = max(checks%projection_energy_residual, &
        abs(heat_capacity*(temperature - config%freezing_temperature) - latent*volume - energy))
      checks%checked_cells = checks%checked_cells + 1
      if (abs(area - old_area) > 32*epsilon(area)*max(1.0_real64, abs(old_area)) .or. &
          abs(volume - old_volume) > 32*epsilon(volume)*max(1.0_real64, abs(old_volume)) .or. &
          abs(temperature - old_temperature) > 32*epsilon(temperature)*max(1.0_real64, abs(old_temperature))) &
        checks%projected_cells = checks%projected_cells + 1
    end if
  end subroutine equilibrate_sea_ice
end module sea_ice
