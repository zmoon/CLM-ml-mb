module controlMod

  !-----------------------------------------------------------------------
  ! !DESCRIPTION:
  ! Initialize run control variables
  !
  ! !USES:
  use shr_kind_mod, only : r8 => shr_kind_r8
  use abortutils,   only : endrun
  !
  ! !PUBLIC TYPES:
  implicit none
  !
  ! !PUBLIC MEMBER FUNCTIONS:
  public :: control
  !-----------------------------------------------------------------------

contains

  !-----------------------------------------------------------------------
  subroutine control
    !
    ! !DESCRIPTION:
    ! Set default run control variables and then read these as namelist variables
    !
    ! !USES:
    use clm_varctl
    use TowerDataMod, only : ntower, tower_id, tower_num, tower_yrbeg, tower_yrend, tower_month
    ! !ARGUMENTS:
    implicit none
    !
    ! !LOCAL VARIABLES:
    character(len=6) :: tower            ! Flux tower site to process
    integer :: i                         ! Index

    logical :: use_hvap, use_init, run_spinup !!! Not needed for namelist !!!
    character(len=256) :: dirini              !!! Not needed for namelist !!!

    namelist /clm_inparm/ light, nsb, gstyp, use_colim, use_acclim, use_clm45kn, &
       use_tower, use_init, use_hvap, tower, tower_yrbeg, tower_yrend, tower_month, &
       diratm, dirclm, dirini, dirout, run_spinup, turb_type, fout_name_suffix
    !---------------------------------------------------------------------

    !---------------------------------------------------------------------
    ! Run control variables: set default values and then initialize via
    ! the [clm_inparm] namelist
    !---------------------------------------------------------------------

    ! Radiative transfer: 1 = Norman. 2 = Goudriaan. 3 = Two-stream

    light = 1

    ! # of equal-width sub-bands to divide each waveband into when doing canopy RT

    nsb = 1  ! results are equivalent to case of no sub-bands

    ! Stomatal conductance: Medlyn (0), Ball-Berry (1), or WUE optimization (2)

    gstyp = 2

    ! Use photosynthesis co-limitation

    use_colim = .false.

    ! Use photosynthesis acclimation

    use_acclim = .true.

    ! Use CLM leaf nitrogen decay coefficient, Kn

    use_clm45kn = .false.

    ! Use tower forcing height and canopy height

    use_tower = .true.

    ! Turbulence parameterization

    turb_type = 1

    ! Plant area density: No stem area (0) or use beta distribution (1) or uniform distribution (2)

    pad_type = 1

    ! True: no foliage clumping

    no_clumping = .true.

    ! True: reduce leaf heat capacity to near zero

    no_storage = .false.

    ! True: iterate leaf temperature during stomatal conductance calculation

    leaf_temp_iter = .false.

    ! True: ncan > ntop. False: ncan = ntop

    use_ncan = .true.

    ! True: allow for interception and evaporation

    use_h2ocan = .false.

    ! True: use within-canopy scalars. False: use values at canopy top

    use_scalars = .true.

    ! Flux tower site, year and month to process, and time step

    tower = ' '
    tower_yrbeg = 0
    tower_yrend = 0
    tower_month = 0

    ! Input and output directories

    diratm = ' '
    dirclm = ' '
    dirout = ' '

    ! Read values from namelist file

    print *, '  reading nml from stdin or waiting for nml to be input on cl ...'
    read (5, clm_inparm)  ! stdin

    ! Bit of validation
    print *, tower_yrbeg, tower_yrend
    if ( tower_yrbeg == 0 .or. tower_yrend == 0 ) stop 'beg/end years not set (possibly nml not directed to the program)'

    ! (Namelist overrides currently follow)

    light = 3                  ! Radiative transfer: 1 = Norman. 2 = Goudriaan. 3 = Two-stream
    gstyp = 2                  ! Stomatal conductance: Medlyn (0), Ball-Berry (1), or WUE optimization (2)
    pad_type = 1               ! Plant area density: No stem area (0) or use beta distribution (1) or uniform distribution (2)

    use_colim = .true.         ! True: use photosynthesis co-limitation
    use_acclim = .true.        ! True: use photosynthesis acclimation
    use_clm45kn = .false.      ! True: use CLM leaf nitrogen decay coefficient (Kn)

    no_clumping = .true.       ! True: no foliage clumping
    no_storage = .false.       ! True: reduce leaf heat capacity
    leaf_temp_iter = .false.   ! True: iterate leaf temperature during stomatal conductance calculation
    use_ncan = .true.          ! True: ncan > ntop. False: ncan = ntop
    use_h2ocan = .true.        ! True: allow for interception and evaporation
    use_scalars = .true.       ! True: use within-canopy scalars. False: use values at canopy top

!   turb_type = 0              ! No canopy profiles
!   turb_type = 1              ! H&F RSL (full RSL parameterization)
!   turb_type = 2              ! H&F but RSL neglected
!   turb_type = 3              ! H&F but RSL for above-canopy only
!   turb_type = 4              ! CLM MOST

    turb_type = 1              ! Turbulence parameterization

    select case (turb_type)
       case (4)
       use_ncan = .false.      ! CLM MOST is only for bulk surface layer
       case (0)
       leaf_temp_iter = .true. ! Must iterate leaf temperature
    end select

    use_tower = .true.

    ! Uncomment to match CLM4.5 simulations

!   gstyp = 1                  ! Stomatal conductance: Medlyn (0), Ball-Berry (1), or WUE optimization (2)
!   use_clm45kn = .true.       ! True: use CLM leaf nitrogen decay coefficient (Kn)
!   pad_type = 2               ! Plant area density: No stem area (0) or use beta distribution (1) or uniform distribution (2)
!   no_storage = .true.        ! True: reduce leaf heat capacity
!   turb_type = 4              ! CLM MOST

    ! Uncomment to incrementally add new features from CLM4.5 equivalent configuration

!   gstyp = 2
!   use_clm45kn = .false.
!   pad_type = 1
!   no_storage = .false.
!   turb_type = 3
!   turb_type = 1

    select case (turb_type)
       case (4)
       use_ncan = .false.      ! CLM MOST is only for bulk surface layer
    end select

    ! Input and output directories

    diratm = '../tower-forcing/'
    dirclm = '../clm4_5/'
    dirout = '../output/'

    !---------------------------------------------------------------------
    ! Match tower site to TowerDataMod arrays
    !---------------------------------------------------------------------

    tower_num = 0
    do i = 1, ntower
       if (tower == tower_id(i)) then
          tower_num = i
          exit
       else
          cycle
       end if
    end do

    if (tower_num == 0) then
       write (iulog,*) 'control error: tower site = ',tower, ' not found'
       call endrun()
    end if

  end subroutine control

end module controlMod
