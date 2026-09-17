module dry_gravity_wave
  use iso_fortran_env, only: real64
  use barotropic_vorticity, only: earth_radius
  use dry_vertical_coordinate, only: hybrid_sigma_coordinate, dry_air_gas_constant, dry_air_kappa, &
                                     reference_surface_pressure
  use lapack_interfaces, only: dgetrf, dgetri
  implicit none
  private

  real(real64), parameter, public :: dry_gravity_wave_implicitness = 0.5_real64

  type, public :: dry_gravity_wave_solver
    private
    integer :: truncation = -1
    integer :: number_of_levels = 0
    real(real64) :: time_step = 0.0_real64
    real(real64), allocatable :: operator_matrix(:, :, :)
    real(real64), allocatable :: inverse_matrix(:, :, :, :)
    real(real64) :: centered_intervals(3) = 0.0_real64
  contains
    procedure, public :: init => initialize_gravity_wave_solver
    procedure, public :: apply => apply_gravity_wave_operator
    procedure, public :: solve => solve_gravity_wave_step
  end type dry_gravity_wave_solver

contains

  subroutine initialize_gravity_wave_solver(this, truncation, time_step, coordinate, reference_temperature)
    class(dry_gravity_wave_solver), intent(inout) :: this
    integer, intent(in) :: truncation
    real(real64), intent(in) :: time_step
    type(hybrid_sigma_coordinate), intent(in) :: coordinate
    real(real64), intent(in) :: reference_temperature(:)
    real(real64), allocatable :: matrix(:, :), work(:)
    integer, allocatable :: pivots(:)
    integer :: dimension, n, interval_index, i, info

    if (truncation < 1) error stop 'dry gravity-wave truncation must be positive'
    if (time_step <= 0.0_real64) error stop 'dry gravity-wave time step must be positive'
    if (coordinate%number_of_levels < 1) error stop 'dry gravity-wave coordinate is not initialized'
    if (size(reference_temperature) /= coordinate%number_of_levels .or. &
        any(reference_temperature <= 0.0_real64)) then
      error stop 'dry gravity-wave reference temperature has an invalid shape or value'
    end if
    this%truncation = truncation
    this%number_of_levels = coordinate%number_of_levels
    this%time_step = time_step
    dimension = 2*this%number_of_levels + 1
    this%centered_intervals = [0.5_real64*time_step, time_step, 2.0_real64*time_step]
    allocate (this%operator_matrix(dimension, dimension, 0:truncation))
    allocate (this%inverse_matrix(dimension, dimension, 0:truncation, 3))
    allocate (matrix(dimension, dimension), pivots(dimension), work(max(1, 64*dimension)))

    do n = 0, truncation
      call build_operator_matrix(coordinate, reference_temperature, n, this%operator_matrix(:, :, n))
      do interval_index = 1, 3
        matrix = -this%centered_intervals(interval_index)*dry_gravity_wave_implicitness* &
                 this%operator_matrix(:, :, n)
        do i = 1, dimension
          matrix(i, i) = matrix(i, i) + 1.0_real64
        end do
        call dgetrf(dimension, dimension, matrix, dimension, pivots, info)
        if (info /= 0) error stop 'dry gravity-wave matrix factorization failed'
        call dgetri(dimension, matrix, dimension, pivots, work, size(work), info)
        if (info /= 0) error stop 'dry gravity-wave matrix inversion failed'
        this%inverse_matrix(:, :, n, interval_index) = matrix
      end do
    end do
  end subroutine initialize_gravity_wave_solver

  subroutine build_operator_matrix(coordinate, reference_temperature, degree, matrix)
    type(hybrid_sigma_coordinate), intent(in) :: coordinate
    real(real64), intent(in) :: reference_temperature(:)
    integer, intent(in) :: degree
    real(real64), intent(out) :: matrix(:, :)
    real(real64), allocatable :: delta_basis(:), cumulative(:)
    real(real64) :: laplacian_factor, mass_above, mass_below, vertical_term
    integer :: number_of_levels, k, j, delta_column, temperature_column

    number_of_levels = coordinate%number_of_levels
    matrix = 0.0_real64
    laplacian_factor = real(degree*(degree + 1), real64)/earth_radius**2

    do j = 1, number_of_levels
      matrix(1, 1 + j) = -coordinate%reference_delta_p(j)/reference_surface_pressure
    end do
    do k = 1, number_of_levels
      matrix(1 + k, 1) = laplacian_factor*dry_air_gas_constant* &
                         coordinate%surface_pressure_geopotential_coefficient(k)
      do j = k + 1, number_of_levels
        matrix(1 + k, 1 + number_of_levels + j) = &
          laplacian_factor*dry_air_gas_constant*coordinate%reference_l(j)
      end do
      matrix(1 + k, 1 + number_of_levels + k) = &
        matrix(1 + k, 1 + number_of_levels + k) + &
        laplacian_factor*dry_air_gas_constant*coordinate%reference_alpha(k)
    end do

    allocate (delta_basis(number_of_levels), cumulative(0:number_of_levels))
    do j = 1, number_of_levels
      delta_basis = 0.0_real64
      delta_basis(j) = 1.0_real64
      cumulative(0) = 0.0_real64
      do k = 1, number_of_levels
        cumulative(k) = cumulative(k - 1) + coordinate%reference_delta_p(k)*delta_basis(k)
      end do
      do k = 1, number_of_levels
        mass_above = coordinate%b_half(k)*cumulative(number_of_levels) - cumulative(k)
        mass_below = coordinate%b_half(k - 1)*cumulative(number_of_levels) - cumulative(k - 1)
        vertical_term = 0.0_real64
        if (k < number_of_levels) vertical_term = vertical_term + mass_above* &
          (reference_temperature(k + 1) - reference_temperature(k))
        if (k > 1) vertical_term = vertical_term + mass_below* &
          (reference_temperature(k) - reference_temperature(k - 1))
        vertical_term = -vertical_term/(2.0_real64*coordinate%reference_delta_p(k))
        vertical_term = vertical_term - dry_air_kappa*reference_temperature(k)* &
          (coordinate%reference_alpha(k)*delta_basis(k) + &
           coordinate%reference_l(k)*cumulative(k - 1)/coordinate%reference_delta_p(k))
        delta_column = 1 + j
        temperature_column = 1 + number_of_levels + k
        matrix(temperature_column, delta_column) = vertical_term
      end do
    end do
  end subroutine build_operator_matrix

  subroutine apply_gravity_wave_operator(this, surface_pressure, delta, temperature, &
                                         tendency_surface_pressure, tendency_delta, tendency_temperature)
    class(dry_gravity_wave_solver), intent(in) :: this
    complex(real64), intent(in) :: surface_pressure(0:, 0:), delta(0:, 0:, :), temperature(0:, 0:, :)
    complex(real64), allocatable, intent(out) :: tendency_surface_pressure(:, :)
    complex(real64), allocatable, intent(out) :: tendency_delta(:, :, :), tendency_temperature(:, :, :)
    complex(real64), allocatable :: state(:), tendency(:)
    integer :: dimension, n, m

    call check_initialized(this)
    dimension = 2*this%number_of_levels + 1
    allocate (tendency_surface_pressure(0:this%truncation + 1, 0:this%truncation))
    allocate (tendency_delta(0:this%truncation + 1, 0:this%truncation, this%number_of_levels))
    allocate (tendency_temperature(0:this%truncation + 1, 0:this%truncation, this%number_of_levels))
    allocate (state(dimension), tendency(dimension))
    tendency_surface_pressure = 0.0_real64
    tendency_delta = 0.0_real64
    tendency_temperature = 0.0_real64
    do m = 0, this%truncation
      do n = m, this%truncation
        call pack_state(this, n, m, surface_pressure, delta, temperature, state)
        tendency = matmul(this%operator_matrix(:, :, n), state)
        call unpack_state(this, n, m, tendency, tendency_surface_pressure, &
                          tendency_delta, tendency_temperature)
      end do
    end do
  end subroutine apply_gravity_wave_operator

  subroutine solve_gravity_wave_step(this, centered_interval, &
                                     previous_surface_pressure, previous_delta, previous_temperature, &
                                     current_surface_pressure, current_delta, current_temperature, &
                                     rhs_surface_pressure, rhs_delta, rhs_temperature, &
                                     next_surface_pressure, next_delta, next_temperature)
    class(dry_gravity_wave_solver), intent(in) :: this
    real(real64), intent(in) :: centered_interval
    complex(real64), intent(in) :: previous_surface_pressure(0:, 0:)
    complex(real64), intent(in) :: previous_delta(0:, 0:, :), previous_temperature(0:, 0:, :)
    complex(real64), intent(in) :: current_surface_pressure(0:, 0:)
    complex(real64), intent(in) :: current_delta(0:, 0:, :), current_temperature(0:, 0:, :)
    complex(real64), intent(in) :: rhs_surface_pressure(0:, 0:)
    complex(real64), intent(in) :: rhs_delta(0:, 0:, :), rhs_temperature(0:, 0:, :)
    complex(real64), allocatable, intent(out) :: next_surface_pressure(:, :)
    complex(real64), allocatable, intent(out) :: next_delta(:, :, :), next_temperature(:, :, :)
    complex(real64), allocatable :: previous_state(:), current_state(:), rhs(:), vector(:)
    integer :: dimension, interval_index, n, m

    call check_initialized(this)
    interval_index = nearest_interval(this, centered_interval)
    dimension = 2*this%number_of_levels + 1
    allocate (next_surface_pressure(0:this%truncation + 1, 0:this%truncation))
    allocate (next_delta(0:this%truncation + 1, 0:this%truncation, this%number_of_levels))
    allocate (next_temperature(0:this%truncation + 1, 0:this%truncation, this%number_of_levels))
    allocate (previous_state(dimension), current_state(dimension), rhs(dimension), vector(dimension))
    next_surface_pressure = 0.0_real64
    next_delta = 0.0_real64
    next_temperature = 0.0_real64

    do m = 0, this%truncation
      do n = m, this%truncation
        call pack_state(this, n, m, previous_surface_pressure, previous_delta, previous_temperature, &
                        previous_state)
        call pack_state(this, n, m, current_surface_pressure, current_delta, current_temperature, &
                        current_state)
        call pack_state(this, n, m, rhs_surface_pressure, rhs_delta, rhs_temperature, rhs)
        vector = previous_state + centered_interval*(rhs - &
          matmul(this%operator_matrix(:, :, n), current_state) + &
          (1.0_real64 - dry_gravity_wave_implicitness)* &
          matmul(this%operator_matrix(:, :, n), previous_state))
        vector = matmul(this%inverse_matrix(:, :, n, interval_index), vector)
        call unpack_state(this, n, m, vector, next_surface_pressure, next_delta, next_temperature)
      end do
    end do
  end subroutine solve_gravity_wave_step

  subroutine pack_state(this, n, m, surface_pressure, delta, temperature, vector)
    class(dry_gravity_wave_solver), intent(in) :: this
    integer, intent(in) :: n, m
    complex(real64), intent(in) :: surface_pressure(0:, 0:), delta(0:, 0:, :), temperature(0:, 0:, :)
    complex(real64), intent(out) :: vector(:)
    integer :: k

    vector(1) = surface_pressure(n, m)
    do k = 1, this%number_of_levels
      vector(1 + k) = delta(n, m, k)
      vector(1 + this%number_of_levels + k) = temperature(n, m, k)
    end do
  end subroutine pack_state

  subroutine unpack_state(this, n, m, vector, surface_pressure, delta, temperature)
    class(dry_gravity_wave_solver), intent(in) :: this
    integer, intent(in) :: n, m
    complex(real64), intent(in) :: vector(:)
    complex(real64), intent(inout) :: surface_pressure(0:, 0:), delta(0:, 0:, :), temperature(0:, 0:, :)
    integer :: k

    surface_pressure(n, m) = vector(1)
    do k = 1, this%number_of_levels
      delta(n, m, k) = vector(1 + k)
      temperature(n, m, k) = vector(1 + this%number_of_levels + k)
    end do
  end subroutine unpack_state

  integer function nearest_interval(this, centered_interval) result(index)
    class(dry_gravity_wave_solver), intent(in) :: this
    real(real64), intent(in) :: centered_interval
    real(real64) :: difference(3), tolerance

    difference = abs(this%centered_intervals - centered_interval)
    index = minloc(difference, dim=1)
    tolerance = 32.0_real64*epsilon(1.0_real64)*max(this%time_step, centered_interval)
    if (difference(index) > tolerance) then
      error stop 'dry gravity-wave solve requested an interval that was not precomputed'
    end if
  end function nearest_interval

  subroutine check_initialized(this)
    class(dry_gravity_wave_solver), intent(in) :: this
    if (this%truncation < 1 .or. .not. allocated(this%inverse_matrix)) then
      error stop 'dry gravity-wave solver is not initialized'
    end if
  end subroutine check_initialized

end module dry_gravity_wave
