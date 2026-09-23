!> Acceptance checks of the band radiation (docs/tendency/band-radiation.md, section 11).
!>
!> The regression values come from an independent Python implementation of the
!> document with the same constants, on the N = 12 levels of the model and the
!> global and tropical reference profiles of section 7.1 (mu_0 = 0.5, surface
!> albedo 0.3).
program check_band_radiation
  use iso_fortran_env, only: real64
  use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
  use dry_physics_config, only: radiation_config, band_radiation_config, radiation_scheme_band, cloud_config
  use cloud_diagnostics, only: cloud_layers, diagnose_cloud_layers, diagnose_cloud_cover
  use band_radiation, only: band_column_fluxes, band_longwave_subbands, band_planck_fractions, &
                            band_planck_fractions_series, prepare_band_planck_table, &
                            band_surface_emission, band_longwave_optics, band_longwave_downward, &
                            band_longwave_upward, band_shortwave, validate_band_radiation_config
  use dry_radiation, only: ozone_layer_fraction, ozone_layer_fractions, radiation_band_downward_column, radiation_band_upward_column, &
                           radiation_tendency, solar_zenith_cosine
  use moist_thermodynamics, only: saturation_specific_humidity
  implicit none
  integer, parameter :: levels = 12
  real(real64), parameter :: half_hpa(0:levels) = [1.0_real64, 3.0_real64, 10.0_real64, 50.0_real64, 100.0_real64, &
    180.0_real64, 280.0_real64, 400.0_real64, 520.0_real64, 650.0_real64, 770.0_real64, 880.0_real64, 1000.0_real64]
  real(real64), parameter :: global_temperature(levels) = [258.07516193921_real64, 246.2998922589474_real64, &
    231.02057724339736_real64, 221.86643925829367_real64, 215.63086961811536_real64, 217.77933890892183_real64, &
    234.58328666143404_real64, 248.46285541084768_real64, 260.0842933620071_real64, 269.84061277962604_real64, &
    277.6547526327356_real64, 284.6311301019913_real64]
  real(real64), parameter :: global_humidity(levels) = [3.0e-6_real64, 3.0e-6_real64, 3.0e-6_real64, 3.0e-6_real64, &
    3.0e-6_real64, 1.5622471533002694e-05_real64, 0.00010509504881026321_real64, 0.0004045113404181236_real64, &
    0.0011006971192993502_real64, 0.002363616694734935_real64, 0.004171461379578295_real64, &
    0.006723494653849436_real64]
  real(real64), parameter :: tropical_temperature(levels) = [258.5041972764614_real64, 238.95627773344185_real64, &
    213.5913537268977_real64, 198.39472913336937_real64, 206.41740593725655_real64, 226.8534780301269_real64, &
    244.35759027232712_real64, 258.8154743862997_real64, 270.9211389187574_real64, 281.0839716454438_real64, &
    289.22370065909956_real64, 296.4907605229076_real64]
  real(real64), parameter :: tropical_humidity(levels) = [3.0e-6_real64, 3.0e-6_real64, 3.0e-6_real64, 3.0e-6_real64, &
    3.3689143930282944e-06_real64, 4.476579548979525e-05_real64, 0.0002733098391367871_real64, &
    0.000983086962872639_real64, 0.0025450639484384707_real64, 0.0052628540732819334_real64, &
    0.00903286766845775_real64, 0.014223570443018415_real64]
  type(radiation_config) :: config
  real(real64) :: pressure_half(0:levels), ozone(levels)
  integer :: k

  config%scheme = radiation_scheme_band
  call validate_band_radiation_config(config%band)
  call check_planck_table()
  call check_ozone_layer_fractions()
  pressure_half = 100.0_real64*half_hpa
  do k = 1, levels
    ozone(k) = ozone_layer_fraction(config, pressure_half(k - 1), pressure_half(k))
  end do

  call check_planck_fractions()
  call check_reference_columns()
  call check_budgets()
  call check_isothermal_column()
  call check_transparent_column()
  call check_subcolumns()
  call check_cloud_layers()
  call check_radiation_tendency_budget()
  write (*, '(a)') 'check_band_radiation: all checks passed'

contains

  subroutine near(actual, target, tolerance, name)
    real(real64), intent(in) :: actual, target, tolerance
    character(*), intent(in) :: name
    if (.not. ieee_is_finite(actual) .or. abs(actual - target) > tolerance) then
      write (*, '(a,3es24.15)') name, actual, target, tolerance
      error stop 'band radiation assertion failed'
    end if
  end subroutine near

  !> Band radiation of one column at mu_0 = cos_zenith with a surface of one temperature.
  subroutine solve(band, temperature, humidity, surface_temperature, cos_zenith, albedo, clouds, fluxes, &
                   longwave_heating, shortwave_heating)
    type(band_radiation_config), intent(in) :: band
    real(real64), intent(in) :: temperature(:), humidity(:), surface_temperature, cos_zenith, albedo
    type(cloud_layers), intent(in) :: clouds
    type(band_column_fluxes), intent(out) :: fluxes
    real(real64), intent(out) :: longwave_heating(:), shortwave_heating(:)
    real(real64) :: transmission(band_longwave_subbands, size(temperature))
    real(real64) :: source(band_longwave_subbands, size(temperature)), downward(0:size(temperature))

    call band_longwave_optics(band, config%stefan_boltzmann_constant, config%gravity_acceleration, pressure_half, &
                              temperature, humidity, ozone, transmission, source)
    call band_longwave_downward(band, clouds, transmission, source, downward)
    call band_shortwave(band, config%gravity_acceleration, config%dry_air_specific_heat, pressure_half, humidity, &
                        ozone, config%solar_constant*max(cos_zenith, 0.0_real64), cos_zenith, albedo, clouds, &
                        shortwave_heating, fluxes)
    call band_longwave_upward(band, clouds, config%gravity_acceleration, config%dry_air_specific_heat, pressure_half, &
                              transmission, source, downward, &
                              band_surface_emission(band, config%stefan_boltzmann_constant, surface_temperature), &
                              longwave_heating, fluxes)
  end subroutine solve

  !> Column mass per unit area times c_p, for the energy of a heating profile.
  real(real64) function column_energy(heating) result(energy)
    real(real64), intent(in) :: heating(:)
    energy = sum(heating*(pressure_half(1:) - pressure_half(:levels - 1)))*config%dry_air_specific_heat/ &
      config%gravity_acceleration
  end function column_energy

  !> The table of the Planck fractions against the series (docs/tendency/band-radiation.md,
  !> 4.3).  Before the table is built and outside its range or for other edges
  !> the series itself is used.  Builds the table, so every later check goes
  !> through it.
  subroutine check_planck_table()
    type(band_radiation_config) :: shifted
    real(real64), parameter :: outside(4) = [99.9_real64, 60.0_real64, 400.1_real64, 450.0_real64]
    real(real64) :: table(band_longwave_subbands), series(band_longwave_subbands), temperature, error
    integer :: n

    if (any(band_planck_fractions(config%band, 250.3_real64) /= band_planck_fractions_series(config%band, 250.3_real64))) &
      error stop 'the Planck fractions before the table is built must be the series'
    call prepare_band_planck_table(config%band)
    call prepare_band_planck_table(config%band)
    error = 0.0_real64
    do n = 0, 300000
      temperature = 100.0_real64 + 300.0_real64*real(n, real64)/300000.0_real64 + 1.0e-7_real64*real(mod(n, 3), real64)
      temperature = min(temperature, 400.0_real64)
      table = band_planck_fractions(config%band, temperature)
      series = band_planck_fractions_series(config%band, temperature)
      error = max(error, maxval(abs(table - series)))
      if (abs(sum(table) - 1.0_real64) > 4.0e-16_real64 .or. any(table < 0.0_real64)) &
        error stop 'tabulated Planck fractions must be nonnegative and add up to one'
    end do
    call near(error, 0.0_real64, 1.0e-14_real64, 'Planck table against the series')
    do n = 1, 4
      temperature = outside(n)
      if (any(band_planck_fractions(config%band, temperature) /= band_planck_fractions_series(config%band, temperature))) &
        error stop 'the Planck fractions outside the table must be the series'
    end do
    shifted = config%band
    shifted%longwave_edges(1) = 360.0_real64
    if (any(band_planck_fractions(shifted, 250.3_real64) /= band_planck_fractions_series(shifted, 250.3_real64))) &
      error stop 'the Planck fractions for other edges must be the series'
  end subroutine check_planck_table

  !> The column ozone fractions equal the layer-by-layer ones bit for bit, for
  !> interfaces above, inside, on and below the ozone layer.
  subroutine check_ozone_layer_fractions()
    real(real64), parameter :: columns(0:6, 3) = reshape([ &
      0.5_real64, 1.0_real64, 3.0_real64, 10.0_real64, 50.0_real64, 100.0_real64, 1000.0_real64, &
      50.0_real64, 80.0_real64, 100.0_real64, 3000.0_real64, 9000.0_real64, 1.0e4_real64, 1.0e5_real64, &
      100.0_real64, 1.0e3_real64, 5.0e3_real64, 2.0e4_real64, 5.0e4_real64, 8.0e4_real64, 1.0e5_real64], [7, 3])
    real(real64) :: fraction(6)
    integer :: c, k

    do c = 1, 3
      call ozone_layer_fractions(config, columns(:, c), fraction)
      do k = 1, 6
        if (fraction(k) /= ozone_layer_fraction(config, columns(k - 1, c), columns(k, c))) &
          error stop 'column ozone fractions must equal the layer ones'
      end do
    end do
  end subroutine check_ozone_layer_fractions

  subroutine check_planck_fractions()
    real(real64) :: fraction(band_longwave_subbands)
    real(real64), parameter :: expected(band_longwave_subbands) = [0.2367121517344538_real64, &
      0.07381360761749178_real64, 0.2571566068966441_real64, 0.23874797252922997_real64, 0.19356966122218033_real64]
    real(real64) :: temperature
    integer :: b, n

    fraction = band_planck_fractions(config%band, 288.0_real64)
    do b = 1, band_longwave_subbands
      call near(fraction(b), expected(b), 1.0e-13_real64, 'Planck fraction at 288 K')
    end do
    do n = 0, 30
      temperature = 150.0_real64 + 8.0_real64*real(n, real64)
      fraction = band_planck_fractions(config%band, temperature)
      call near(sum(fraction), 1.0_real64, 4.0e-16_real64, 'Planck fractions add up to one')
      if (any(fraction < 0.0_real64)) error stop 'negative Planck fraction'
    end do
  end subroutine check_planck_fractions

  subroutine check_reference_columns()
    type(band_column_fluxes) :: fluxes
    real(real64) :: longwave(levels), shortwave(levels)
    type(cloud_layers) :: clouds
    real(real64), parameter :: tolerance = 1.0e-9_real64

    clouds = cloud_layers()
    call solve(config%band, global_temperature, global_humidity, 288.0_real64, 0.5_real64, 0.3_real64, clouds, &
               fluxes, longwave, shortwave)
    call near(fluxes%outgoing_longwave, 271.9718239461179_real64, tolerance, 'global OLR')
    call near(fluxes%clear_outgoing_longwave, 271.9718239461179_real64, tolerance, 'global clear OLR')
    call near(fluxes%window_outgoing_longwave, 83.18186547451538_real64, tolerance, 'global window OLR')
    call near(fluxes%surface_downward_longwave, 297.20104427123897_real64, tolerance, 'global downward longwave')
    call near(fluxes%reflected_shortwave, 181.67935150074848_real64, tolerance, 'global reflected shortwave')
    call near(fluxes%surface_incident_shortwave, 516.8527055620722_real64, tolerance, 'global surface shortwave')
    call near(fluxes%atmospheric_shortwave_absorption, 137.02375460580103_real64, tolerance, &
              'global atmospheric shortwave')

    clouds = cloud_layers(large_scale_fraction=0.45_real64, convective_fraction=0.2_real64, &
                          large_scale_level=10, convective_level=6)
    call solve(config%band, global_temperature, global_humidity, 288.0_real64, 0.5_real64, 0.3_real64, clouds, &
               fluxes, longwave, shortwave)
    call near(fluxes%outgoing_longwave, 232.61348758453158_real64, tolerance, 'cloudy global OLR')
    call near(fluxes%clear_outgoing_longwave, 271.9718239461179_real64, tolerance, 'cloudy global clear OLR')
    call near(fluxes%window_outgoing_longwave, 61.412770026770374_real64, tolerance, 'cloudy global window OLR')
    call near(fluxes%surface_downward_longwave, 324.0679448897703_real64, tolerance, 'cloudy global downward')
    call near(fluxes%reflected_shortwave, 300.4367636834485_real64, tolerance, 'cloudy global reflected')
    call near(fluxes%clear_reflected_shortwave, 181.67935150074848_real64, tolerance, 'cloudy clear reflected')
    call near(fluxes%surface_incident_shortwave, 372.39237435747305_real64, tolerance, 'cloudy surface shortwave')
    call near(fluxes%clear_surface_incident_shortwave, 516.8527055620722_real64, tolerance, 'cloudy clear surface')
    call near(fluxes%atmospheric_shortwave_absorption, 119.38857426632043_real64, tolerance, 'cloudy atmosphere')

    clouds = cloud_layers()
    call solve(config%band, tropical_temperature, tropical_humidity, 300.0_real64, 0.5_real64, 0.3_real64, clouds, &
               fluxes, longwave, shortwave)
    call near(fluxes%outgoing_longwave, 293.75909780320103_real64, tolerance, 'tropical OLR')
    call near(fluxes%window_outgoing_longwave, 96.33729842584695_real64, tolerance, 'tropical window OLR')
    call near(fluxes%surface_downward_longwave, 384.03566085836326_real64, tolerance, 'tropical downward longwave')
    call near(fluxes%reflected_shortwave, 172.93979329162164_real64, tolerance, 'tropical reflected shortwave')
    call near(fluxes%atmospheric_shortwave_absorption, 165.14105111034877_real64, tolerance, &
              'tropical atmospheric shortwave')
    clouds = cloud_layers(large_scale_fraction=0.45_real64, convective_fraction=0.2_real64, &
                          large_scale_level=10, convective_level=6)
    call solve(config%band, tropical_temperature, tropical_humidity, 300.0_real64, 0.5_real64, 0.3_real64, clouds, &
               fluxes, longwave, shortwave)
    call near(fluxes%outgoing_longwave, 254.34239555603608_real64, tolerance, 'cloudy tropical OLR')
    call near(fluxes%surface_downward_longwave, 403.42746358317606_real64, tolerance, 'cloudy tropical downward')
    call near(fluxes%reflected_shortwave, 290.79268333080773_real64, tolerance, 'cloudy tropical reflected')
    call near(fluxes%atmospheric_shortwave_absorption, 142.99431506081203_real64, tolerance, 'cloudy tropical atm')
  end subroutine check_reference_columns

  !> The longwave and shortwave column budgets close with and without clouds,
  !> at every sun angle including the grazing and the night side.
  subroutine check_budgets()
    type(band_column_fluxes) :: fluxes
    real(real64) :: longwave(levels), shortwave(levels), cos_zenith, albedo, surface_temperature, emission
    type(cloud_layers) :: clouds
    integer :: n, c

    do c = 1, 3
      select case (c)
      case (1)
        clouds = cloud_layers()
      case (2)
        clouds = cloud_layers(large_scale_fraction=1.0_real64, large_scale_level=12)
      case default
        clouds = cloud_layers(large_scale_fraction=0.3_real64, convective_fraction=0.5_real64, &
                              large_scale_level=9, convective_level=1)
      end select
      do n = 0, 6
        cos_zenith = -0.1_real64 + 0.18_real64*real(n, real64)
        albedo = 0.1_real64*real(n, real64)
        surface_temperature = 270.0_real64 + 5.0_real64*real(n, real64)
        call solve(config%band, tropical_temperature, tropical_humidity, surface_temperature, cos_zenith, albedo, &
                   clouds, fluxes, longwave, shortwave)
        call near(fluxes%incoming_shortwave, fluxes%atmospheric_shortwave_absorption + &
                  (1.0_real64 - albedo)*fluxes%surface_incident_shortwave + fluxes%reflected_shortwave, &
                  1.0e-10_real64, 'shortwave column budget')
        call near(column_energy(shortwave), fluxes%atmospheric_shortwave_absorption, 1.0e-10_real64, &
                  'shortwave heating is the absorption')
        emission = config%stefan_boltzmann_constant*surface_temperature**4
        call near(fluxes%surface_upward_longwave, emission, 1.0e-10_real64, 'surface emission sub-band sum')
        call near(column_energy(longwave), fluxes%surface_upward_longwave - fluxes%surface_downward_longwave - &
                  fluxes%outgoing_longwave, 1.0e-10_real64, 'longwave column budget')
        if (cos_zenith <= 0.0_real64 .and. (fluxes%reflected_shortwave /= 0.0_real64 .or. &
            any(shortwave /= 0.0_real64))) error stop 'shortwave on the night side'
      end do
    end do
  end subroutine check_budgets

  !> An isothermal column emits sigma T^4 upward at every interface, but it still
  !> cools: no longwave enters at the top, so the column loses the downward
  !> longwave it sends to the surface.
  subroutine check_isothermal_column()
    type(band_column_fluxes) :: fluxes
    real(real64) :: temperature(levels), longwave(levels), shortwave(levels)
    real(real64) :: transmission(band_longwave_subbands, levels), source(band_longwave_subbands, levels)
    real(real64) :: downward(0:levels), expected(0:levels), fraction(band_longwave_subbands), sigma_t4
    type(cloud_layers) :: clouds
    integer :: b

    temperature = 280.0_real64
    clouds = cloud_layers()
    sigma_t4 = config%stefan_boltzmann_constant*280.0_real64**4
    call solve(config%band, temperature, tropical_humidity, 280.0_real64, 0.0_real64, 0.3_real64, clouds, fluxes, &
               longwave, shortwave)
    call near(fluxes%outgoing_longwave, sigma_t4, 1.0e-10_real64, 'isothermal OLR is sigma T^4')
    call near(column_energy(longwave), -fluxes%surface_downward_longwave, 1.0e-10_real64, &
              'isothermal column cools by its downward longwave')
    if (.not. all(longwave < 0.0_real64)) error stop 'isothermal column must cool in every layer'
    call band_longwave_optics(config%band, config%stefan_boltzmann_constant, config%gravity_acceleration, &
                              pressure_half, temperature, tropical_humidity, ozone, transmission, source)
    call band_longwave_downward(config%band, clouds, transmission, source, downward)
    fraction = band_planck_fractions(config%band, 280.0_real64)
    expected = 0.0_real64
    do b = 1, band_longwave_subbands
      do k = 1, levels
        expected(k) = expected(k) + fraction(b)*sigma_t4*(1.0_real64 - product(transmission(b, 1:k)))
      end do
    end do
    do k = 0, levels
      call near(downward(k), expected(k), 1.0e-10_real64, 'isothermal downward longwave')
    end do
  end subroutine check_isothermal_column

  !> Without absorbers the longwave leaves unchanged and the shortwave is
  !> reflected only by the Rayleigh layer and the surface.
  subroutine check_transparent_column()
    type(band_radiation_config) :: band
    type(band_column_fluxes) :: fluxes
    real(real64) :: longwave(levels), shortwave(levels), incoming
    type(cloud_layers) :: clouds

    band = config%band
    band%co2_volume_mixing_ratio = 0.0_real64
    band%ozone_column = 0.0_real64
    clouds = cloud_layers()
    call solve(band, global_temperature, [(0.0_real64, k=1, levels)], 290.0_real64, 1.0_real64, 0.25_real64, clouds, &
               fluxes, longwave, shortwave)
    call near(fluxes%outgoing_longwave, config%stefan_boltzmann_constant*290.0_real64**4, 1.0e-10_real64, &
              'transparent OLR')
    call near(fluxes%surface_downward_longwave, 0.0_real64, 1.0e-12_real64, 'transparent downward longwave')
    call near(fluxes%atmospheric_shortwave_absorption, 0.0_real64, 1.0e-12_real64, 'transparent absorption')
    band%rayleigh_direct_amplitude = 0.0_real64
    band%rayleigh_diffuse_reflectance = 0.0_real64
    call solve(band, global_temperature, [(0.0_real64, k=1, levels)], 290.0_real64, 1.0_real64, 0.25_real64, clouds, &
               fluxes, longwave, shortwave)
    incoming = config%solar_constant
    call near(fluxes%surface_incident_shortwave, incoming, 1.0e-10_real64, 'transparent surface shortwave')
    call near(fluxes%reflected_shortwave, 0.25_real64*incoming, 1.0e-10_real64, 'transparent reflection')
  end subroutine check_transparent_column

  !> The clear sub-column alone reproduces the cloud-free result bit for bit, and
  !> two clouds of equal properties in the same layer act as one cloud of their
  !> summed cover.
  subroutine check_subcolumns()
    type(band_column_fluxes) :: clear, empty, split, single
    real(real64) :: lw_clear(levels), sw_clear(levels), lw_empty(levels), sw_empty(levels)
    real(real64) :: lw_split(levels), sw_split(levels), lw_single(levels), sw_single(levels)

    call solve(config%band, tropical_temperature, tropical_humidity, 300.0_real64, 0.6_real64, 0.2_real64, &
               cloud_layers(), clear, lw_clear, sw_clear)
    call solve(config%band, tropical_temperature, tropical_humidity, 300.0_real64, 0.6_real64, 0.2_real64, &
               cloud_layers(large_scale_level=4, convective_level=2), empty, lw_empty, sw_empty)
    if (any(lw_clear /= lw_empty) .or. any(sw_clear /= sw_empty) .or. &
        clear%outgoing_longwave /= empty%outgoing_longwave .or. &
        clear%reflected_shortwave /= empty%reflected_shortwave .or. &
        clear%surface_incident_shortwave /= empty%surface_incident_shortwave) &
      error stop 'zero cloud fractions must reproduce the clear column bit for bit'
    call solve(config%band, tropical_temperature, tropical_humidity, 300.0_real64, 0.6_real64, 0.2_real64, &
               cloud_layers(large_scale_fraction=0.3_real64, convective_fraction=0.2_real64, &
                            large_scale_level=8, convective_level=8), split, lw_split, sw_split)
    call solve(config%band, tropical_temperature, tropical_humidity, 300.0_real64, 0.6_real64, 0.2_real64, &
               cloud_layers(large_scale_fraction=0.5_real64, large_scale_level=8), single, lw_single, sw_single)
    call near(split%outgoing_longwave, single%outgoing_longwave, 1.0e-10_real64, 'split cloud OLR')
    call near(split%reflected_shortwave, single%reflected_shortwave, 1.0e-10_real64, 'split cloud reflection')
    call near(split%surface_downward_longwave, single%surface_downward_longwave, 1.0e-10_real64, 'split cloud DLR')
    call near(maxval(abs(lw_split - lw_single)), 0.0_real64, 1.0e-12_real64, 'split cloud longwave heating')
    call near(maxval(abs(sw_split - sw_single)), 0.0_real64, 1.0e-12_real64, 'split cloud shortwave heating')
    ! A cloud raises the downward longwave at the surface and the reflection, and lowers the OLR.
    if (.not. (single%outgoing_longwave < clear%outgoing_longwave .and. &
               single%reflected_shortwave > clear%reflected_shortwave .and. &
               single%surface_downward_longwave > clear%surface_downward_longwave)) &
      error stop 'cloud radiative effects have the wrong sign'
  end subroutine check_subcolumns

  !> The convective cloud takes C^conv at the convection top, the large-scale cloud
  !> the rest of C^RH at the uppermost saturated level, and their sum is the
  !> effective cover of the grey diagnosis.
  subroutine check_cloud_layers()
    type(cloud_config) :: clouds_config
    type(cloud_layers) :: layers
    real(real64) :: full_pressure(levels), humidity(levels), cover
    real(real64), parameter :: mm_per_day = 1.0_real64/86400.0_real64

    full_pressure = 0.5_real64*(pressure_half(1:) + pressure_half(:levels - 1))
    do k = 1, levels
      humidity(k) = 0.5_real64*saturation_specific_humidity(tropical_temperature(k), full_pressure(k))
    end do
    humidity(9) = saturation_specific_humidity(tropical_temperature(9), full_pressure(9))
    humidity(11) = saturation_specific_humidity(tropical_temperature(11), full_pressure(11))
    call diagnose_cloud_layers(clouds_config, full_pressure, tropical_temperature, humidity, 0.0_real64, 0, cover, &
                               layers)
    call near(cover, 1.0_real64, 1.0e-14_real64, 'saturated cover')
    if (layers%large_scale_level /= 9 .or. layers%convective_fraction /= 0.0_real64) &
      error stop 'large-scale cloud must sit at the uppermost saturated level'
    call diagnose_cloud_layers(clouds_config, full_pressure, tropical_temperature, humidity, &
                               10.0_real64*mm_per_day, 5, cover, layers)
    if (layers%convective_level /= 5 .or. layers%large_scale_level /= 9) error stop 'cloud levels'
    call near(layers%convective_fraction + layers%large_scale_fraction, cover, 1.0e-15_real64, 'cover split')
    call near(cover, diagnose_cloud_cover(clouds_config, full_pressure, tropical_temperature, humidity, &
                                          10.0_real64*mm_per_day), 0.0_real64, 'split cover equals grey cover')
    humidity = 0.5_real64*humidity
    call diagnose_cloud_layers(clouds_config, full_pressure, tropical_temperature, humidity, &
                               10.0_real64*mm_per_day, 5, cover, layers)
    if (layers%large_scale_fraction /= 0.0_real64 .or. layers%large_scale_level /= 0 .or. &
        layers%convective_level /= 5) error stop 'dry column keeps only the convective cloud'
  end subroutine check_cloud_layers

  !> The single-surface path closes the energy of atmosphere, ground and top of
  !> the atmosphere with the band scheme.
  subroutine check_radiation_tendency_budget()
    real(real64) :: tendency(levels), surface_tendency, deep_tendency, incoming, reflected, outgoing
    real(real64) :: stored, sin_latitude
    type(band_column_fluxes) :: fluxes

    sin_latitude = 0.3_real64
    call radiation_tendency(config, pressure_half, tropical_temperature, 299.0_real64, 297.0_real64, 3.0_real64, &
                            -2.0_real64, sin_latitude, 0.0_real64, 0.0_real64, tendency, surface_tendency, &
                            deep_tendency, incoming, reflected, outgoing, specific_humidity=tropical_humidity, &
                            clouds=cloud_layers(large_scale_fraction=0.4_real64, large_scale_level=10), &
                            band_fluxes=fluxes)
    if (.not. (incoming > 0.0_real64)) error stop 'the budget check needs daylight'
    call near(incoming, config%solar_constant*max(0.0_real64, &
              solar_zenith_cosine(config, sin_latitude, 0.0_real64, 0.0_real64)), 1.0e-12_real64, 'incoming')
    stored = column_energy(tendency) + config%surface_heat_capacity*surface_tendency + &
      config%deep_ground_heat_capacity*deep_tendency
    call near(stored, incoming - reflected - outgoing, 1.0e-9_real64, 'band single-surface energy budget')
    call near(reflected, fluxes%reflected_shortwave, 0.0_real64, 'reflected is the top-of-atmosphere flux')
  end subroutine check_radiation_tendency_budget

end program check_band_radiation
