!> Maps the command-line argument to the cases to run.  Everything a case
!> decides lives in its own module under cases/.
module case_runners
  use case_runtime, only: case_context
  use barotropic_case, only: run_barotropic_case
  use shallow_water_case, only: run_shallow_water_case
  use dry_case, only: run_dry_case
  use radiation_case, only: run_radiation_case
  implicit none
  private

  public :: run_requested_cases

contains

  subroutine run_requested_cases(requested_case)
    character(*), intent(in) :: requested_case
    type(case_context) :: context

    select case (trim(requested_case))
    case ('shallow-water', '--shallow-water', 'shallow_water', 'swe')
      call run_shallow_water_case(context, 'shallow_water_mountain', 1)
      call run_shallow_water_case(context, 'shallow_water_single_harmonic', 2)
    case ('barotropic', '--barotropic', 'barotropic-vorticity', 'bve')
      call run_barotropic_case(context, 'single_harmonic', 1)
      call run_barotropic_case(context, 'rossby_haurwitz_r4', 2)
      call run_barotropic_case(context, 'random_n8_n12_seed_20260913', 3)
    case ('dry', '--dry', 'dry-atmosphere')
      call run_dry_case(context, 'dry_jablonowski_williamson_steady', .false.)
      call run_dry_case(context, 'dry_jablonowski_williamson_perturbed', .true.)
    case ('held-suarez', '--held-suarez', 'held_suarez')
      call run_dry_case(context, 'dry_held_suarez', .false., held_suarez=.true.)
    case ('radiation', '--radiation', 'uniform-radiation', '--uniform-radiation', 'uniform_radiation')
      call run_radiation_case(context)
    case ('all')
      call run_dry_case(context, 'dry_jablonowski_williamson_steady', .false.)
      call run_dry_case(context, 'dry_jablonowski_williamson_perturbed', .true.)
      call run_shallow_water_case(context, 'shallow_water_mountain', 1)
      call run_shallow_water_case(context, 'shallow_water_single_harmonic', 2)
      call run_barotropic_case(context, 'single_harmonic', 1)
      call run_barotropic_case(context, 'rossby_haurwitz_r4', 2)
      call run_barotropic_case(context, 'random_n8_n12_seed_20260913', 3)
      call run_dry_case(context, 'dry_held_suarez', .false., held_suarez=.true.)
      call run_radiation_case(context)
    case ('--help', '-h', 'help')
      call print_usage()
    case default
      call print_usage()
      error stop 'unknown equation argument'
    end select
  end subroutine run_requested_cases

  subroutine print_usage()
    write (*, '(a)') 'Usage: core [shallow-water|barotropic|dry|held-suarez|radiation|all]'
    write (*, '(a)') '  shallow-water:           run the mountain and single-harmonic height cases'
    write (*, '(a)') '  barotropic:             run the three barotropic-vorticity cases'
    write (*, '(a)') '  dry (default):          run the 10-day Jablonowski-Williamson dry-atmosphere cases'
    write (*, '(a)') '                          (steady base state and localized wind perturbation)'
    write (*, '(a)') '  held-suarez:            run the 200-day forced dry-atmosphere case (output every 5 days)'
    write (*, '(a)') '  radiation:              run the 5-year diurnal/seasonal radiation case'
    write (*, '(a)') '  all:                    run every case'
  end subroutine print_usage

end module case_runners
