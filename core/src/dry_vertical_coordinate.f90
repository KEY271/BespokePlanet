module dry_vertical_coordinate
  use iso_fortran_env, only: real64
  implicit none
  private

  integer, parameter, public :: default_number_of_levels = 12
  real(real64), parameter, public :: reference_surface_pressure = 1.0e5_real64
  real(real64), parameter, public :: dry_air_gas_constant = 287.0_real64
  real(real64), parameter, public :: dry_air_kappa = 2.0_real64/7.0_real64

  type, public :: hybrid_sigma_coordinate
    integer :: number_of_levels = 0
    real(real64), allocatable :: a_half(:), b_half(:), delta_b(:)
    real(real64), allocatable :: reference_p_half(:), reference_delta_p(:)
    real(real64), allocatable :: reference_l(:), reference_alpha(:)
    real(real64), allocatable :: reference_temperature(:), full_level_eta(:)
    real(real64), allocatable :: lambda(:), gamma(:), pressure_gradient_coefficient(:)
    real(real64), allocatable :: surface_pressure_geopotential_coefficient(:)
  contains
    procedure, public :: init => initialize_coordinate
    procedure, public :: init_default => initialize_default_coordinate
  end type hybrid_sigma_coordinate

  public :: jablonowski_mean_temperature

contains

  subroutine initialize_default_coordinate(this)
    class(hybrid_sigma_coordinate), intent(inout) :: this
    real(real64), parameter :: a(0:default_number_of_levels) = [ &
      100.0_real64, 300.0_real64, 1000.0_real64, 5000.0_real64, &
      10000.0_real64, 8000.0_real64, 8000.0_real64, 10000.0_real64, &
      12000.0_real64, 10000.0_real64, 7000.0_real64, 3000.0_real64, &
      0.0_real64]
    real(real64), parameter :: b(0:default_number_of_levels) = [ &
      0.0_real64, 0.0_real64, 0.0_real64, 0.0_real64, 0.0_real64, &
      0.1_real64, 0.2_real64, 0.3_real64, 0.4_real64, 0.55_real64, &
      0.7_real64, 0.85_real64, 1.0_real64]

    call this%init(a, b)
  end subroutine initialize_default_coordinate

  subroutine initialize_coordinate(this, a_half, b_half)
    class(hybrid_sigma_coordinate), intent(inout) :: this
    real(real64), intent(in) :: a_half(0:), b_half(0:)
    integer :: k, j, number_of_levels

    if (size(a_half) /= size(b_half) .or. size(a_half) < 2) then
      error stop 'hybrid sigma A and B arrays must have the same nontrivial size'
    end if
    number_of_levels = size(a_half) - 1
    if (abs(b_half(0)) > epsilon(1.0_real64) .or. &
        abs(b_half(number_of_levels) - 1.0_real64) > epsilon(1.0_real64)) then
      error stop 'hybrid sigma B must be zero at the top and one at the surface'
    end if

    this%number_of_levels = number_of_levels
    this%a_half = a_half
    this%b_half = b_half
    allocate (this%delta_b(number_of_levels))
    allocate (this%reference_p_half(0:number_of_levels))
    allocate (this%reference_delta_p(number_of_levels), this%reference_l(number_of_levels))
    allocate (this%reference_alpha(number_of_levels), this%reference_temperature(number_of_levels))
    allocate (this%full_level_eta(number_of_levels))
    allocate (this%lambda(number_of_levels), this%gamma(number_of_levels))
    allocate (this%pressure_gradient_coefficient(number_of_levels))
    allocate (this%surface_pressure_geopotential_coefficient(number_of_levels))

    this%reference_p_half = a_half + b_half*reference_surface_pressure
    do k = 1, number_of_levels
      this%delta_b(k) = b_half(k) - b_half(k - 1)
      this%reference_delta_p(k) = this%reference_p_half(k) - this%reference_p_half(k - 1)
      if (this%reference_delta_p(k) <= 0.0_real64) then
        error stop 'hybrid sigma reference pressures must increase downward'
      end if
      this%reference_l(k) = log(this%reference_p_half(k)/this%reference_p_half(k - 1))
      this%reference_alpha(k) = 1.0_real64 - &
        this%reference_p_half(k - 1)*this%reference_l(k)/this%reference_delta_p(k)
      ! eta is tied to the hybrid coefficients by eta_{k+1/2} = A_{k+1/2}/p0 + B_{k+1/2},
      ! so that p = eta*p0 wherever ps = p0.  Full levels are the arithmetic mean of
      ! the bounding half levels.
      this%full_level_eta(k) = 0.5_real64*(this%reference_p_half(k - 1) + &
        this%reference_p_half(k))/reference_surface_pressure
      this%reference_temperature(k) = jablonowski_mean_temperature(this%full_level_eta(k))
      this%lambda(k) = reference_surface_pressure*( &
        b_half(k)/this%reference_p_half(k) - b_half(k - 1)/this%reference_p_half(k - 1))
      this%gamma(k) = -reference_surface_pressure*b_half(k - 1)*this%reference_l(k)/ &
        this%reference_delta_p(k) + &
        this%reference_p_half(k - 1)*reference_surface_pressure*this%delta_b(k)* &
        this%reference_l(k)/this%reference_delta_p(k)**2 - &
        this%reference_p_half(k - 1)*this%lambda(k)/this%reference_delta_p(k)
      this%pressure_gradient_coefficient(k) = reference_surface_pressure/this%reference_delta_p(k)*( &
        b_half(k - 1)*this%reference_l(k) + this%reference_alpha(k)*this%delta_b(k))
    end do

    do k = 1, number_of_levels
      this%surface_pressure_geopotential_coefficient(k) = &
        this%reference_temperature(k)*(this%gamma(k) + this%pressure_gradient_coefficient(k))
      do j = k + 1, number_of_levels
        this%surface_pressure_geopotential_coefficient(k) = &
          this%surface_pressure_geopotential_coefficient(k) + &
          this%reference_temperature(j)*this%lambda(j)
      end do
    end do
  end subroutine initialize_coordinate

  pure real(real64) function jablonowski_mean_temperature(eta) result(temperature)
    real(real64), intent(in) :: eta
    real(real64), parameter :: base_temperature = 288.0_real64
    real(real64), parameter :: lapse_rate = 0.005_real64
    real(real64), parameter :: gravity = 9.80616_real64
    real(real64), parameter :: tropopause_eta = 0.2_real64
    real(real64), parameter :: stratospheric_adjustment = 4.8e5_real64

    if (eta <= 0.0_real64) error stop 'full-level eta must be positive'
    temperature = base_temperature*eta**(dry_air_gas_constant*lapse_rate/gravity)
    if (eta < tropopause_eta) then
      temperature = temperature + stratospheric_adjustment*(tropopause_eta - eta)**5
    end if
  end function jablonowski_mean_temperature

end module dry_vertical_coordinate
