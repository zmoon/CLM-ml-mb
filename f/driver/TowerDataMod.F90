module TowerDataMod

  !-----------------------------------------------------------------------
  ! !DESCRIPTION:
  ! Parameters for flux tower sites
  !
  ! !USES:
  use shr_kind_mod, only : r8 => shr_kind_r8
  !
  ! !PUBLIC TYPES:
  implicit none
  save
  !-----------------------------------------------------------------------

  integer, parameter :: ntower = 17        ! Number of tower sites

  character(len=6) :: tower_id(ntower)     ! Flux tower sites
  real(r8) :: tower_lat(ntower)            ! Latitude (degrees)
  real(r8) :: tower_lon(ntower)            ! Longitude (degrees)
  integer :: tower_pft(ntower)             ! CLM PFT for each tower site
  character(len=15) :: tower_tex(ntower)   ! Soil texture class for each tower site
  integer :: tower_isoicol(ntower)         ! Soil color class for each tower site
  real(r8) :: tower_ht(ntower)             ! Flux tower height (m)
  real(r8) :: tower_canht(ntower)          ! Canopy height for each tower site (m)
  integer :: tower_time(ntower)            ! Time step of forcing data (minutes)

  data tower_id / 'US-Ha1', 'US-Ho1', 'US-MMS', 'US-UMB', 'US-Dk3', 'US-Me2', &
                  'US-Var', 'US-IB1', 'US-Ne3', 'US-ARM', 'US-Bo1', 'US-Dk1', 'US-Dk2', &
                  'BR-Ma2', 'BR-Ji2', 'BR-Sa1', 'BR-Sa3'/

  data tower_lat / 42.54_r8,  45.20_r8, 39.32_r8, 45.56_r8, 35.98_r8, 44.45_r8, &
                   38.41_r8,  41.86_r8, 41.18_r8, 36.61_r8, 40.01_r8, 35.97_r8, 35.97_r8, &
                   -2.61_r8, -10.08_r8, -2.86_r8, -3.02_r8/

  data tower_lon /  -72.17_r8,  -68.74_r8,  -86.41_r8,  -84.71_r8,  -79.09_r8, -121.56_r8, &
                   -120.95_r8,  -88.22_r8,  -96.44_r8,  -97.49_r8,  -88.29_r8,  -79.09_r8, -79.10_r8, &
                    -60.21_r8,  -61.93_r8,  -54.96_r8,  -54.97_r8/

  data tower_pft / 7, 2, 7, 7, 1, 2, 13, 15, 15, 15, 15, 13, 7, 4, 4, 4, 4/

  data tower_tex / 'loam', 'sandy loam', 'clay', 'sand', 'sandy loam', 'sandy loam', &
                   'silty loam','silty clay loam', 'clay loam', 'clay', &
                   'silty loam', 'sandy loam', 'sandy loam', &
                   'clay', 'loamy sand', 'clay', 'clay'/

  data tower_isoicol / 18, 16, 15, 17, 15, 20, 17, 15, 13, 13, 15, 15, 15, 16, 15, 14, 15/

  data tower_ht    / 30._r8, 29._r8, 48._r8,   46._r8, 22._r8, 32._r8, &
                     2.5_r8,  4._r8,  6._r8, -999._r8,  6._r8,  5._r8, 42._r8, &
                     50._r8, 60._r8, 63._r8,   64._r8/

  data tower_canht / 23._r8, 20._r8, 27._r8, 21._r8, 17._r8, 14._r8, &
                     0.6_r8, 0.9_r8, 0.9_r8, 0.5_r8, 0.9_r8, 0.5_r8, 25._r8, &
                     35._r8, 30._r8, 35._r8, 35._r8/

  data tower_time / 60, 30, 60, 60, 30, 30, 30, 30, 60, 30, 30, 30, 30, 60, 60, 60, 60/

  ! Specific year and month to process are controlled by namelist

  integer :: tower_num      ! Tower site index (maps to TowerDataMod arrays)
  integer :: tower_yrbeg    ! First year to process
  integer :: tower_yrend    ! Last year to process
  integer :: tower_month    ! Month to process

end module TowerDataMod
