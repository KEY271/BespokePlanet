!> Band radiation (docs/tendency/band-radiation.md).
!>
!> The longwave is split into five grey sub-bands (window, CO2 centre, CO2
!> wings, weak and strong water vapour) solved by the same non-scattering
!> two-stream as the grey scheme with the diffusivity factor D.  The shortwave
!> is split into four bands (UV, visible, weak and strong near-infrared water
!> vapour) of twelve g points in all.  The direct beam follows the magnified
!> path M(mu_0); light reflected by the surface, the Rayleigh layer or a cloud
!> goes up as diffuse light and is absorbed on the way.
!>
!> Every column is solved as a clear, a large-scale-cloud and a
!> convective-cloud sub-column weighted by their area fractions
!> (cloud_diagnostics).  A cloud sub-column holds one cloud: in the longwave the
!> transmittance of its layer is multiplied by 1 - epsilon, in the shortwave the
!> cloud reflects alpha of the direct beam at the top of its layer.  The clear
!> sub-column is always solved, so the clear-sky fluxes come without extra
!> passes.  With no cloud the clear sub-column has weight one and the result is
!> the cloud-free result bit for bit.
!>
!> The column is solved in two stages because the surface temperatures of the
!> tiles are found between them: the downward stage returns the optical
!> properties, the downward longwave and the whole shortwave (which needs only
!> the surface albedo), the upward stage takes the sub-band surface emission and
!> returns the upward longwave and the longwave heating.
!>
!> The Planck fractions f_b(T) come from a table of quintic Hermite pieces in T
!> on [100, 400] K (section 4.3), built once from the series by
!> prepare_band_planck_table.  The table belongs to the band edges and c_2 it
!> was built for; for other edges, outside its range or before it is built the
!> series is summed directly.  The table is written only by
!> prepare_band_planck_table, outside the parallel regions, and read by the
!> pure procedures.
module band_radiation
  use iso_fortran_env, only: real64
  use dry_physics_config, only: band_radiation_config
  use cloud_diagnostics, only: cloud_layers
  implicit none
  private

  integer, parameter, public :: band_longwave_subbands = 5
  integer, parameter, public :: band_shortwave_bands = 4
  integer, parameter, public :: band_shortwave_gpoints = 12
  integer, parameter, public :: band_window_subband = 1
  integer, parameter :: planck_series_terms = 40
  integer, parameter :: clear_subcolumn = 1, large_scale_subcolumn = 2, convective_subcolumn = 3
  real(real64), parameter :: pi = acos(-1.0_real64)
  integer :: series_index
  real(real64), parameter :: inverse_integers(planck_series_terms) = &
    [(1.0_real64/real(series_index, real64), series_index=1, planck_series_terms)]
  integer, parameter :: planck_edges = 8

  !> Planck-fraction table: on the interval j of [T_min + (j-1) h, T_min + j h]
  !> f_b = sum_n c(b, n, j) s^n with s = (T - T_min)/h - (j - 1), for the
  !> sub-bands 1-4; f_5 is the remainder as in the series.
  real(real64), parameter :: planck_table_minimum_temperature = 100.0_real64
  real(real64), parameter :: planck_table_maximum_temperature = 400.0_real64
  real(real64), parameter :: planck_table_spacing = 1.0_real64
  integer, parameter :: planck_table_intervals = 300
  integer, parameter :: planck_table_degree = 5
  logical :: planck_table_ready = .false.
  real(real64) :: planck_table_edges(planck_edges) = 0.0_real64
  real(real64) :: planck_table_constant = 0.0_real64
  real(real64) :: planck_table(band_longwave_subbands - 1, 0:planck_table_degree, planck_table_intervals) = 0.0_real64

  !> Column fluxes of one evaluation (W m^-2), area means over the sub-columns
  !> unless named clear.  The surface incident shortwave is what reaches the
  !> surface after the Rayleigh multiple reflection, of which a tile of albedo
  !> alpha absorbs 1 - alpha.
  type, public :: band_column_fluxes
    real(real64) :: incoming_shortwave = 0.0_real64
    real(real64) :: reflected_shortwave = 0.0_real64
    real(real64) :: clear_reflected_shortwave = 0.0_real64
    real(real64) :: surface_incident_shortwave = 0.0_real64
    real(real64) :: clear_surface_incident_shortwave = 0.0_real64
    real(real64) :: atmospheric_shortwave_absorption = 0.0_real64
    real(real64) :: surface_downward_longwave = 0.0_real64
    real(real64) :: surface_upward_longwave = 0.0_real64
    real(real64) :: outgoing_longwave = 0.0_real64
    real(real64) :: clear_outgoing_longwave = 0.0_real64
    real(real64) :: window_outgoing_longwave = 0.0_real64
  end type band_column_fluxes

  public :: band_planck_fractions, band_planck_fractions_series, band_surface_emission
  public :: prepare_band_planck_table
  public :: band_longwave_optics, band_longwave_downward, band_longwave_upward, band_shortwave
  public :: validate_band_radiation_config, validate_cloud_layers

contains

  !> Fraction of sigma T^4 emitted at wavenumbers above nu, as a function of
  !> x = c_2 nu/T: (15/pi^4) sum_m e^{-m x}/m (x^3 + 3x^2/m + 6x/m^2 + 6/m^3).
  pure real(real64) function planck_fraction_above(x) result(fraction)
    real(real64), intent(in) :: x
    real(real64) :: decay, term, series, inverse_m, x2, x3, addend
    integer :: m

    decay = exp(-x)
    x2 = x*x
    x3 = x2*x
    term = 1.0_real64
    series = 0.0_real64
    ! The terms fall off like e^{-m x}; stop once one no longer changes the sum.
    do m = 1, planck_series_terms
      term = term*decay
      inverse_m = inverse_integers(m)
      addend = term*inverse_m*(x3 + inverse_m*(3.0_real64*x2 + inverse_m*(6.0_real64*x + 6.0_real64*inverse_m)))
      series = series + addend
      if (addend <= epsilon(1.0_real64)*series) exit
    end do
    fraction = 15.0_real64/pi**4*series
  end function planck_fraction_above

  !> Sub-band sums of the fractions G(x_e) above the eight edges, or of their
  !> temperature derivatives, for the sub-bands 1-4.  Edges: 350, 500, 630, 700,
  !> 820, 1180, 1390, 1800 cm^-1.
  pure function subband_sums(above) result(fraction)
    real(real64), intent(in) :: above(planck_edges)
    real(real64) :: fraction(band_longwave_subbands - 1)

    fraction(1) = above(5) - above(6)
    fraction(2) = above(3) - above(4)
    fraction(3) = (above(2) - above(3)) + (above(4) - above(5))
    fraction(4) = (above(1) - above(2)) + (above(6) - above(7)) + above(8)
  end function subband_sums

  !> Fractions f_b(T) of sigma T^4 in the five longwave sub-bands, summed from
  !> the series.  The fifth, strong water vapour, is the remainder, so the
  !> fractions add up to one.
  pure function band_planck_fractions_series(band, temperature) result(fraction)
    type(band_radiation_config), intent(in) :: band
    real(real64), intent(in) :: temperature
    real(real64) :: fraction(band_longwave_subbands)
    real(real64) :: above(planck_edges)
    integer :: e

    do e = 1, planck_edges
      above(e) = planck_fraction_above(band%second_radiation_constant*band%longwave_edges(e)/temperature)
    end do
    fraction(1:4) = subband_sums(above)
    fraction(5) = 1.0_real64 - (fraction(1) + fraction(2) + fraction(3) + fraction(4))
  end function band_planck_fractions_series

  !> The table was built for the edges and c_2 of `band`.
  pure logical function planck_table_matches(band) result(matches)
    type(band_radiation_config), intent(in) :: band

    matches = planck_table_ready
    if (.not. matches) return
    matches = band%second_radiation_constant == planck_table_constant .and. &
              all(band%longwave_edges == planck_table_edges)
  end function planck_table_matches

  !> f_b(T) from the table; T must lie in its range.
  pure function planck_fractions_from_table(temperature) result(fraction)
    real(real64), intent(in) :: temperature
    real(real64) :: fraction(band_longwave_subbands)
    real(real64) :: position, piece(band_longwave_subbands - 1)
    integer :: j, n

    position = (temperature - planck_table_minimum_temperature)/planck_table_spacing
    j = min(int(position), planck_table_intervals - 1)
    position = position - real(j, real64)
    j = j + 1
    piece = planck_table(:, planck_table_degree, j)
    do n = planck_table_degree - 1, 0, -1
      piece = piece*position + planck_table(:, n, j)
    end do
    fraction(1:4) = piece
    fraction(5) = 1.0_real64 - (fraction(1) + fraction(2) + fraction(3) + fraction(4))
  end function planck_fractions_from_table

  !> f_b(T) from the table when `use_table` and T lies in its range, otherwise
  !> from the series.
  pure function planck_fractions(band, use_table, temperature) result(fraction)
    type(band_radiation_config), intent(in) :: band
    logical, intent(in) :: use_table
    real(real64), intent(in) :: temperature
    real(real64) :: fraction(band_longwave_subbands)

    if (use_table .and. temperature >= planck_table_minimum_temperature .and. &
        temperature <= planck_table_maximum_temperature) then
      fraction = planck_fractions_from_table(temperature)
    else
      fraction = band_planck_fractions_series(band, temperature)
    end if
  end function planck_fractions

  !> Fractions f_b(T) of sigma T^4 in the five longwave sub-bands, from the table
  !> where it applies (see the module header).
  pure function band_planck_fractions(band, temperature) result(fraction)
    type(band_radiation_config), intent(in) :: band
    real(real64), intent(in) :: temperature
    real(real64) :: fraction(band_longwave_subbands)

    fraction = planck_fractions(band, planck_table_matches(band), temperature)
  end function band_planck_fractions

  !> Builds the Planck-fraction table for the edges and c_2 of `band`; nothing
  !> to do when it is already built for them.  Each piece is the quintic that
  !> matches f_b, df_b/dT and d^2f_b/dT^2 at both ends, with
  !> dG/dT = (15/pi^4) x^4/((e^x - 1) T) and
  !> d^2G/dT^2 = -(15/pi^4) [5 x^4/(e^x - 1) - x^5 e^x/(e^x - 1)^2]/T^2.
  !> Must not be called inside a parallel region.
  subroutine prepare_band_planck_table(band)
    type(band_radiation_config), intent(in) :: band
    real(real64), dimension(band_longwave_subbands - 1, 0:planck_table_intervals) :: value, first, second
    real(real64), dimension(band_longwave_subbands - 1) :: f0, d0, e0, f1, d1, e1
    real(real64) :: first_above(planck_edges), second_above(planck_edges), fraction(band_longwave_subbands)
    real(real64) :: temperature, x, expm, h
    integer :: i, e

    if (planck_table_matches(band)) return
    planck_table_ready = .false.
    h = planck_table_spacing
    do i = 0, planck_table_intervals
      temperature = planck_table_minimum_temperature + h*real(i, real64)
      fraction = band_planck_fractions_series(band, temperature)
      value(:, i) = fraction(1:4)
      do e = 1, planck_edges
        x = band%second_radiation_constant*band%longwave_edges(e)/temperature
        expm = exp(x) - 1.0_real64
        first_above(e) = 15.0_real64/pi**4*x**4/(expm*temperature)
        second_above(e) = -15.0_real64/pi**4*(5.0_real64*x**4/expm - x**5*(expm + 1.0_real64)/expm**2)/ &
                          temperature**2
      end do
      first(:, i) = subband_sums(first_above)
      second(:, i) = subband_sums(second_above)
    end do
    do i = 1, planck_table_intervals
      f0 = value(:, i - 1)
      d0 = h*first(:, i - 1)
      e0 = h*h*second(:, i - 1)
      f1 = value(:, i)
      d1 = h*first(:, i)
      e1 = h*h*second(:, i)
      planck_table(:, 0, i) = f0
      planck_table(:, 1, i) = d0
      planck_table(:, 2, i) = 0.5_real64*e0
      planck_table(:, 3, i) = -10.0_real64*f0 - 6.0_real64*d0 - 1.5_real64*e0 + 10.0_real64*f1 - 4.0_real64*d1 + &
                              0.5_real64*e1
      planck_table(:, 4, i) = 15.0_real64*f0 + 8.0_real64*d0 + 1.5_real64*e0 - 15.0_real64*f1 + 7.0_real64*d1 - e1
      planck_table(:, 5, i) = -6.0_real64*f0 - 3.0_real64*d0 - 0.5_real64*e0 + 6.0_real64*f1 - 3.0_real64*d1 + &
                              0.5_real64*e1
    end do
    planck_table_edges = band%longwave_edges
    planck_table_constant = band%second_radiation_constant
    planck_table_ready = .true.
  end subroutine prepare_band_planck_table

  !> Black-body emission f_b(T) sigma T^4 of a surface in each sub-band.
  pure function band_surface_emission(band, sigma, temperature) result(emission)
    type(band_radiation_config), intent(in) :: band
    real(real64), intent(in) :: sigma, temperature
    real(real64) :: emission(band_longwave_subbands)

    emission = band_planck_fractions(band, temperature)*sigma*temperature**4
  end function band_surface_emission

  !> (p/p_0)^n without the power function for the common exponents 0 and 1;
  !> other exponents use exp(n log(p/p_0)) with `log_relative_pressure`, taken
  !> once per layer for all of them.
  pure real(real64) function pressure_scaling(relative_pressure, log_relative_pressure, exponent) result(scaling)
    real(real64), intent(in) :: relative_pressure, log_relative_pressure, exponent

    if (exponent == 0.0_real64) then
      scaling = 1.0_real64
    else if (exponent == 1.0_real64) then
      scaling = relative_pressure
    else
      scaling = exp(exponent*log_relative_pressure)
    end if
  end function pressure_scaling

  !> exp(-min(x, x_max)) of a nonnegative optical path x; exactly one for a
  !> transparent layer without calling the exponential.
  pure real(real64) function path_transmission(path, maximum_path) result(transmission)
    real(real64), intent(in) :: path, maximum_path

    if (path > 0.0_real64) then
      transmission = exp(-min(path, maximum_path))
    else
      transmission = 1.0_real64
    end if
  end function path_transmission

  !> Layer transmittances exp(-D dtau_b) and layer sources f_b(T) sigma T^4 of the
  !> clear sub-column.  `humidity` is the specific humidity of every layer
  !> (negative values count as zero) and `ozone_fraction` the fraction of the
  !> ozone column in every layer.
  pure subroutine band_longwave_optics(band, sigma, gravity, pressure_half, temperature, humidity, ozone_fraction, &
                                       transmission, source)
    type(band_radiation_config), intent(in) :: band
    real(real64), intent(in) :: sigma, gravity
    real(real64), intent(in) :: pressure_half(0:), temperature(:), humidity(:), ozone_fraction(:)
    real(real64), intent(out) :: transmission(:, :), source(:, :)
    real(real64) :: fraction(band_longwave_subbands), co2_factor(band_longwave_subbands)
    real(real64) :: pressure_thickness, mean_pressure, relative_pressure, water, water_path, co2_path, ozone_path
    real(real64) :: vapour, continuum_factor, emission, optical_depth, co2_ratio, co2_reference_mass
    real(real64) :: log_relative_pressure
    logical :: use_table, needs_log
    integer :: k, b

    co2_reference_mass = band%co2_reference_volume_mixing_ratio*band%co2_molar_mass_ratio
    co2_ratio = band%co2_volume_mixing_ratio/band%co2_reference_volume_mixing_ratio
    do b = 1, band_longwave_subbands
      if (co2_ratio > 0.0_real64) then
        co2_factor(b) = co2_ratio**band%longwave_co2_concentration_exponent(b)
      else
        co2_factor(b) = 0.0_real64
      end if
    end do
    use_table = planck_table_matches(band)
    needs_log = .false.
    do b = 1, band_longwave_subbands
      if (band%longwave_line(b) > 0.0_real64 .and. band%longwave_line_pressure_exponent(b) /= 0.0_real64 .and. &
          band%longwave_line_pressure_exponent(b) /= 1.0_real64) needs_log = .true.
      if (band%longwave_co2(b) > 0.0_real64 .and. band%longwave_co2_pressure_exponent(b) /= 0.0_real64 .and. &
          band%longwave_co2_pressure_exponent(b) /= 1.0_real64) needs_log = .true.
    end do
    log_relative_pressure = 0.0_real64
    do k = 1, size(temperature)
      pressure_thickness = pressure_half(k) - pressure_half(k - 1)
      mean_pressure = 0.5_real64*(pressure_half(k - 1) + pressure_half(k))
      relative_pressure = mean_pressure/band%reference_pressure
      if (needs_log) log_relative_pressure = log(relative_pressure)
      water = max(humidity(k), 0.0_real64)
      water_path = water*pressure_thickness/gravity
      co2_path = co2_reference_mass*pressure_thickness/gravity
      ozone_path = band%ozone_column*ozone_fraction(k)
      vapour = mean_pressure*water/(band%water_vapor_molar_mass_ratio*band%reference_pressure)
      continuum_factor = exp(band%self_continuum_temperature* &
        (1.0_real64/temperature(k) - 1.0_real64/band%self_continuum_reference_temperature))
      fraction = planck_fractions(band, use_table, temperature(k))
      emission = sigma*temperature(k)**4
      do b = 1, band_longwave_subbands
        optical_depth = band%longwave_ozone(b)*ozone_path
        if (band%longwave_line(b) > 0.0_real64) optical_depth = optical_depth + &
          band%longwave_line(b)*water_path* &
          pressure_scaling(relative_pressure, log_relative_pressure, band%longwave_line_pressure_exponent(b))
        if (band%longwave_self_continuum(b) > 0.0_real64) optical_depth = optical_depth + &
          band%longwave_self_continuum(b)*water_path*vapour*continuum_factor
        if (band%longwave_co2(b) > 0.0_real64) optical_depth = optical_depth + &
          band%longwave_co2(b)*co2_path*co2_factor(b)* &
          pressure_scaling(relative_pressure, log_relative_pressure, band%longwave_co2_pressure_exponent(b))
        transmission(b, k) = path_transmission(band%diffusivity*optical_depth, band%maximum_optical_path)
        source(b, k) = fraction(b)*emission
      end do
    end do
  end subroutine band_longwave_optics

  !> Sub-columns to solve: the clear one always, the cloud ones when their
  !> fraction is positive.  `weight` is the area fraction, `level` the cloud
  !> level (0 for clear), `emissivity` and `albedo` the cloud properties.
  pure subroutine band_subcolumns(band, clouds, count, kind, weight, level, emissivity, albedo)
    type(band_radiation_config), intent(in) :: band
    type(cloud_layers), intent(in) :: clouds
    integer, intent(out) :: count, kind(3), level(3)
    real(real64), intent(out) :: weight(3), emissivity(3), albedo(3)

    count = 1
    kind(1) = clear_subcolumn
    weight(1) = max(1.0_real64 - clouds%large_scale_fraction - clouds%convective_fraction, 0.0_real64)
    level(1) = 0
    emissivity(1) = 0.0_real64
    albedo(1) = 0.0_real64
    if (clouds%large_scale_fraction > 0.0_real64) then
      count = count + 1
      kind(count) = large_scale_subcolumn
      weight(count) = clouds%large_scale_fraction
      level(count) = clouds%large_scale_level
      emissivity(count) = band%large_scale_cloud_longwave_emissivity
      albedo(count) = band%large_scale_cloud_shortwave_albedo
    end if
    if (clouds%convective_fraction > 0.0_real64) then
      count = count + 1
      kind(count) = convective_subcolumn
      weight(count) = clouds%convective_fraction
      level(count) = clouds%convective_level
      emissivity(count) = band%convective_cloud_longwave_emissivity
      albedo(count) = band%convective_cloud_shortwave_albedo
    end if
  end subroutine band_subcolumns

  !> Area-mean downward longwave at every interface, summed over the sub-bands.
  !> No longwave enters at the top.
  pure subroutine band_longwave_downward(band, clouds, transmission, source, downward)
    type(band_radiation_config), intent(in) :: band
    type(cloud_layers), intent(in) :: clouds
    real(real64), intent(in) :: transmission(:, :), source(:, :)
    real(real64), intent(out) :: downward(0:)
    real(real64) :: subcolumn_downward(0:size(source, 2)), weight(3), emissivity(3), albedo(3), flux, layer
    integer :: count, kind(3), level(3), m, b, k, levels

    levels = size(source, 2)
    call band_subcolumns(band, clouds, count, kind, weight, level, emissivity, albedo)
    downward = 0.0_real64
    do m = 1, count
      if (weight(m) <= 0.0_real64) cycle
      subcolumn_downward = 0.0_real64
      do b = 1, band_longwave_subbands
        flux = 0.0_real64
        do k = 1, levels
          layer = transmission(b, k)
          if (k == level(m)) layer = (1.0_real64 - emissivity(m))*layer
          flux = layer*flux + (1.0_real64 - layer)*source(b, k)
          subcolumn_downward(k) = subcolumn_downward(k) + flux
        end do
      end do
      downward = downward + weight(m)*subcolumn_downward
    end do
  end subroutine band_longwave_downward

  !> Upward longwave from the sub-band surface emission, the longwave heating
  !> of the net flux against the area-mean `downward`, and the outgoing, clear
  !> and window fluxes at the top.
  pure subroutine band_longwave_upward(band, clouds, gravity, specific_heat, pressure_half, transmission, source, &
                                       downward, surface_emission, heating, fluxes)
    type(band_radiation_config), intent(in) :: band
    type(cloud_layers), intent(in) :: clouds
    real(real64), intent(in) :: gravity, specific_heat
    real(real64), intent(in) :: pressure_half(0:), transmission(:, :), source(:, :), downward(0:)
    real(real64), intent(in) :: surface_emission(:)
    real(real64), intent(out) :: heating(:)
    type(band_column_fluxes), intent(inout) :: fluxes
    real(real64) :: upward(0:size(source, 2)), subcolumn_upward(0:size(source, 2)), net(0:size(source, 2))
    real(real64) :: weight(3), emissivity(3), albedo(3), flux, layer, window
    integer :: count, kind(3), level(3), m, b, k, levels

    levels = size(source, 2)
    call band_subcolumns(band, clouds, count, kind, weight, level, emissivity, albedo)
    upward = 0.0_real64
    fluxes%window_outgoing_longwave = 0.0_real64
    do m = 1, count
      subcolumn_upward = 0.0_real64
      window = 0.0_real64
      do b = 1, band_longwave_subbands
        flux = surface_emission(b)
        subcolumn_upward(levels) = subcolumn_upward(levels) + flux
        do k = levels, 1, -1
          layer = transmission(b, k)
          if (k == level(m)) layer = (1.0_real64 - emissivity(m))*layer
          flux = layer*flux + (1.0_real64 - layer)*source(b, k)
          subcolumn_upward(k - 1) = subcolumn_upward(k - 1) + flux
        end do
        if (b == band_window_subband) window = flux
      end do
      if (kind(m) == clear_subcolumn) fluxes%clear_outgoing_longwave = subcolumn_upward(0)
      if (weight(m) <= 0.0_real64) cycle
      upward = upward + weight(m)*subcolumn_upward
      fluxes%window_outgoing_longwave = fluxes%window_outgoing_longwave + weight(m)*window
    end do
    net = upward - downward
    do k = 1, levels
      heating(k) = gravity/(specific_heat*(pressure_half(k) - pressure_half(k - 1)))*(net(k) - net(k - 1))
    end do
    fluxes%outgoing_longwave = upward(0)
    fluxes%surface_upward_longwave = upward(levels)
    fluxes%surface_downward_longwave = downward(levels)
  end subroutine band_longwave_upward

  !> Shortwave of the twelve g points in every sub-column.  `incoming` is
  !> S_0 max(0, mu_0); nothing is computed on the night side.  The surface is
  !> one reflector of area-mean albedo `surface_albedo` under a non-absorbing
  !> Rayleigh layer; the returned surface incident flux is the area mean over
  !> the sub-columns.
  pure subroutine band_shortwave(band, gravity, specific_heat, pressure_half, humidity, ozone_fraction, incoming, &
                                 cos_zenith, surface_albedo, clouds, heating, fluxes)
    type(band_radiation_config), intent(in) :: band
    real(real64), intent(in) :: gravity, specific_heat
    real(real64), intent(in) :: pressure_half(0:), humidity(:), ozone_fraction(:)
    real(real64), intent(in) :: incoming, cos_zenith, surface_albedo
    type(cloud_layers), intent(in) :: clouds
    real(real64), intent(out) :: heating(:)
    type(band_column_fluxes), intent(inout) :: fluxes
    real(real64), dimension(size(humidity)) :: water_path, ozone_path, optical_depth, direct, diffuse, absorbed
    real(real64) :: weight(3), emissivity(3), albedo(3)
    real(real64) :: magnification, normalization, top_flux, reflectance, diffuse_reflectance
    real(real64) :: flux, next_flux, cloud_reflected, incident, pressure_thickness, mean_pressure
    real(real64) :: relative_pressure, log_relative_pressure
    logical :: needs_log
    integer :: count, kind(3), level(3), m, j, b, k, levels

    levels = size(humidity)
    heating = 0.0_real64
    fluxes%incoming_shortwave = incoming
    fluxes%reflected_shortwave = 0.0_real64
    fluxes%clear_reflected_shortwave = 0.0_real64
    fluxes%surface_incident_shortwave = 0.0_real64
    fluxes%clear_surface_incident_shortwave = 0.0_real64
    fluxes%atmospheric_shortwave_absorption = 0.0_real64
    if (incoming <= 0.0_real64) return
    call band_subcolumns(band, clouds, count, kind, weight, level, emissivity, albedo)
    magnification = band%magnification_numerator/sqrt(band%magnification_quadratic*cos_zenith**2 + 1.0_real64)
    normalization = 0.0_real64
    do j = 1, band_shortwave_gpoints
      normalization = normalization + &
        band%shortwave_band_fraction(band%shortwave_gpoint_band(j))*band%shortwave_gpoint_weight(j)
    end do
    needs_log = band%shortwave_water_vapor_pressure_exponent /= 0.0_real64 .and. &
                band%shortwave_water_vapor_pressure_exponent /= 1.0_real64
    do k = 1, levels
      pressure_thickness = pressure_half(k) - pressure_half(k - 1)
      mean_pressure = 0.5_real64*(pressure_half(k - 1) + pressure_half(k))
      relative_pressure = mean_pressure/band%reference_pressure
      log_relative_pressure = 0.0_real64
      if (needs_log) log_relative_pressure = log(relative_pressure)
      water_path(k) = max(humidity(k), 0.0_real64)*pressure_thickness/gravity* &
        pressure_scaling(relative_pressure, log_relative_pressure, band%shortwave_water_vapor_pressure_exponent)
      ozone_path(k) = band%ozone_column*ozone_fraction(k)
    end do
    absorbed = 0.0_real64
    do j = 1, band_shortwave_gpoints
      b = band%shortwave_gpoint_band(j)
      top_flux = incoming*band%shortwave_band_fraction(b)*band%shortwave_gpoint_weight(j)/normalization
      if (band%shortwave_gpoint_water_vapor(j) > 0.0_real64 .or. band%shortwave_gpoint_ozone(j) > 0.0_real64) then
        optical_depth = band%shortwave_gpoint_water_vapor(j)*water_path + band%shortwave_gpoint_ozone(j)*ozone_path
        ! The ozone-only g points see no absorber below the ozone layer.
        do k = 1, levels
          direct(k) = path_transmission(magnification*optical_depth(k), band%maximum_optical_path)
          diffuse(k) = path_transmission(band%diffusivity*optical_depth(k), band%maximum_optical_path)
        end do
      else
        ! A window g point passes the atmosphere unabsorbed.
        direct = 1.0_real64
        diffuse = 1.0_real64
      end if
      reflectance = band%rayleigh_direct_amplitude(b)/(1.0_real64 + band%rayleigh_direct_slope(b)*cos_zenith)
      diffuse_reflectance = band%rayleigh_diffuse_reflectance(b)
      do m = 1, count
        flux = top_flux
        cloud_reflected = 0.0_real64
        do k = 1, levels
          if (k == level(m)) then
            cloud_reflected = albedo(m)*flux
            flux = flux - cloud_reflected
          end if
          next_flux = flux*direct(k)
          absorbed(k) = absorbed(k) + weight(m)*(flux - next_flux)
          flux = next_flux
        end do
        ! Rayleigh layer and surface with their multiple reflection.
        incident = (1.0_real64 - reflectance)*flux/(1.0_real64 - diffuse_reflectance*surface_albedo)
        flux = reflectance*flux + (1.0_real64 - diffuse_reflectance)*surface_albedo*incident
        do k = levels, 1, -1
          next_flux = flux*diffuse(k)
          absorbed(k) = absorbed(k) + weight(m)*(flux - next_flux)
          flux = next_flux
          if (k == level(m)) flux = flux + cloud_reflected
        end do
        if (kind(m) == clear_subcolumn) then
          fluxes%clear_reflected_shortwave = fluxes%clear_reflected_shortwave + flux
          fluxes%clear_surface_incident_shortwave = fluxes%clear_surface_incident_shortwave + incident
        end if
        fluxes%reflected_shortwave = fluxes%reflected_shortwave + weight(m)*flux
        fluxes%surface_incident_shortwave = fluxes%surface_incident_shortwave + weight(m)*incident
      end do
    end do
    do k = 1, levels
      heating(k) = gravity/(specific_heat*(pressure_half(k) - pressure_half(k - 1)))*absorbed(k)
    end do
    fluxes%atmospheric_shortwave_absorption = sum(absorbed)
  end subroutine band_shortwave

  !> Rejects coefficient sets that would break the energy closure.
  subroutine validate_band_radiation_config(band)
    type(band_radiation_config), intent(in) :: band
    real(real64) :: weight_sum
    integer :: b, j

    if (.not. (band%diffusivity > 0.0_real64) .or. band%co2_volume_mixing_ratio < 0.0_real64 .or. &
        .not. (band%co2_reference_volume_mixing_ratio > 0.0_real64) .or. band%ozone_column < 0.0_real64 .or. &
        .not. (band%reference_pressure > 0.0_real64) .or. .not. (band%maximum_optical_path > 0.0_real64)) &
      error stop 'invalid band radiation constants'
    if (any(band%longwave_edges(2:) <= band%longwave_edges(:size(band%longwave_edges) - 1)) .or. &
        band%longwave_edges(1) <= 0.0_real64) error stop 'band longwave edges must increase'
    if (any(band%longwave_line < 0.0_real64) .or. any(band%longwave_self_continuum < 0.0_real64) .or. &
        any(band%longwave_co2 < 0.0_real64) .or. any(band%longwave_ozone < 0.0_real64) .or. &
        any(band%shortwave_gpoint_water_vapor < 0.0_real64) .or. any(band%shortwave_gpoint_ozone < 0.0_real64)) &
      error stop 'band absorption coefficients must be nonnegative'
    if (abs(sum(band%shortwave_band_fraction) - 1.0_real64) > 1.0e-10_real64 .or. &
        any(band%shortwave_band_fraction < 0.0_real64)) error stop 'band solar fractions must add up to one'
    if (any(band%shortwave_gpoint_band < 1) .or. any(band%shortwave_gpoint_band > band_shortwave_bands)) &
      error stop 'band g point outside the shortwave bands'
    do b = 1, band_shortwave_bands
      weight_sum = 0.0_real64
      do j = 1, band_shortwave_gpoints
        if (band%shortwave_gpoint_band(j) == b) weight_sum = weight_sum + band%shortwave_gpoint_weight(j)
      end do
      if (abs(weight_sum - 1.0_real64) > 1.0e-10_real64) error stop 'band g point weights must add up to one'
    end do
    if (any(band%shortwave_gpoint_weight < 0.0_real64)) error stop 'band g point weights must be nonnegative'
    if (any(band%rayleigh_direct_amplitude < 0.0_real64) .or. any(band%rayleigh_direct_amplitude > 1.0_real64) .or. &
        any(band%rayleigh_direct_slope < 0.0_real64) .or. any(band%rayleigh_diffuse_reflectance < 0.0_real64) .or. &
        any(band%rayleigh_diffuse_reflectance >= 1.0_real64)) error stop 'invalid band Rayleigh reflectance'
    if (.not. (band%magnification_numerator > 0.0_real64) .or. band%magnification_quadratic < 0.0_real64) &
      error stop 'invalid band magnification'
    if (min(band%large_scale_cloud_shortwave_albedo, band%convective_cloud_shortwave_albedo, &
            band%large_scale_cloud_longwave_emissivity, band%convective_cloud_longwave_emissivity) < 0.0_real64 .or. &
        max(band%large_scale_cloud_shortwave_albedo, band%convective_cloud_shortwave_albedo, &
            band%large_scale_cloud_longwave_emissivity, band%convective_cloud_longwave_emissivity) > 1.0_real64) &
      error stop 'band cloud properties must lie in [0, 1]'
  end subroutine validate_band_radiation_config

  !> Rejects cloud sub-columns outside [0, 1] or without a level in the column.
  subroutine validate_cloud_layers(clouds, levels)
    type(cloud_layers), intent(in) :: clouds
    integer, intent(in) :: levels

    if (clouds%large_scale_fraction < 0.0_real64 .or. clouds%convective_fraction < 0.0_real64 .or. &
        clouds%large_scale_fraction + clouds%convective_fraction > 1.0_real64 + 1.0e-12_real64) &
      error stop 'cloud sub-column fractions are outside [0, 1]'
    if (clouds%large_scale_fraction > 0.0_real64 .and. &
        (clouds%large_scale_level < 1 .or. clouds%large_scale_level > levels)) &
      error stop 'large-scale cloud level is outside the column'
    if (clouds%convective_fraction > 0.0_real64 .and. &
        (clouds%convective_level < 1 .or. clouds%convective_level > levels)) &
      error stop 'convective cloud level is outside the column'
  end subroutine validate_cloud_layers

end module band_radiation
