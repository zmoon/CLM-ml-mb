module CLMml_driver

  !-----------------------------------------------------------------------
  ! !DESCRIPTION:
  ! Model driver
  !
  ! !USES:
  use shr_kind_mod, only : r8 => shr_kind_r8
  use abortutils,   only : endrun
  use decompMod ,   only : bounds_type
  !
  ! !PUBLIC TYPES:
  implicit none
  !
  ! !PUBLIC MEMBER FUNCTIONS:
  public :: CLMml_drv             ! Model driver
  !
  ! !PRIVATE MEMBER FUNCTIONS:
  private :: clm_subgrid          ! Patch and column mapping for CLM g/l/c/p hierarchy
  private :: init_acclim          ! Read tower meteorology data to get acclimation temperature
  private :: output               ! Write output files
  !-----------------------------------------------------------------------

contains

  !-----------------------------------------------------------------------
  subroutine CLMml_drv (bounds)
    !
    ! !DESCRIPTION:
    ! Model driver
    !
    ! !USES:
    use clm_varctl, only : iulog, ntim, diratm, dirclm, dirout, subdir
    use clm_time_manager
    use TowerDataMod, only : tower_id, tower_num, tower_yrbeg, tower_yrend, tower_month, tower_time
    use clm_varorb, only : eccen, mvelpp, lambm0, obliqr
    use shr_orb_mod, only : shr_orb_params
    use shr_file_mod, only : shr_file_getUnit, shr_file_freeUnit
    use controlMod, only : control
    use filterMod, only : setFilters, filter
    use TowerMetMod, only : TowerMet
    use CLMml_initializeMod, only : initialize
    use clm_driver, only : clm_drv
    use clm_instMod
    !
    ! !ARGUMENTS:
    implicit none
    type(bounds_type), intent(in) :: bounds
    !
    ! !LOCAL VARIABLES:
    integer  :: it                             ! Tower site index
    integer  :: iy                             ! Tower year index
    real(r8) :: obliq, mvelp                   ! Miscellaneous orbital parameters (not used)
    integer  :: istep                          ! Time stepping loop index
    integer  :: time_indx                      ! Time index from reference date (0Z January 1 of current year, when calday = 1.000)
    integer  :: yr                             ! Year (1900, ...)
    integer  :: mon                            ! Month (1, ..., 12)
    integer  :: day                            ! Day of month (1, ..., 31)
    integer  :: days_per_month                 ! Days in month
    integer  :: steps_per_day                  ! Number of time steps per day

    character(len=256) :: ext                  ! Local file name
    character(len=256) :: fin1, fin2           ! Full file name, including directory path
    character(len=256) :: fout1, fout2, fout3  ! Full file name, including directory path
    integer :: nout1, nout2, nout3             ! Fortran unit number
    !---------------------------------------------------------------------

    ! Initialize run control variables

    print *, "  initializing run control variables ..."
    call control

    ! Process the tower site (it) and year (iy)

    it = tower_num
    do iy = tower_yrbeg, tower_yrend

       !---------------------------------------------------------------
       ! Start date is 0Z July 1 of current year
       !
       ! start_date_ymd = year (yyyy), month (mm) and day (dd) of the
       ! simulation start date in yyyymmdd format (e.g., 19960701)
       !
       ! start_date_tod = time of day (UTC) of the simulation start date
       ! (seconds past 0Z)
       !---------------------------------------------------------------

       start_date_ymd = iy*10000 + tower_month*100 + 1
       start_date_tod = 0

       call get_curr_date (yr, mon, day, curr_date_tod)

      !  write (iulog,*) 'Processing: ',tower_id(it),yr,mon
       print *, '  processing: ',tower_id(it),yr,mon

       !---------------------------------------------------------------
       ! Number of time steps to execute
       !---------------------------------------------------------------

       ! Time step of forcing data (in seconds)

       dtstep = tower_time(it) * 60

       ! Number of days in the month

       if (isleap(yr, calkindflag)) then
          days_per_month = mdayleap(mon)
       else
          days_per_month = mday(mon)
       end if

       ! Number of time steps per day

       steps_per_day = 86400 / dtstep

       ! Number of time steps to execute

       ntim = steps_per_day * days_per_month

       !---------------------------------------------------------------
       ! Initialize variables
       !---------------------------------------------------------------

       print *, "  initializing variables ..."
       call initialize (it, iy, bounds%begp, bounds%endp, bounds%begc, bounds%endc, &
       soilstate_inst, waterstate_inst, temperature_inst, surfalb_inst, mlcanopy_inst)

       !---------------------------------------------------------------
       ! Build patch and column mapping for CLM g/l/c/p hierarchy. This
       ! code processes only one tower-year at a time. So there is one grid
       ! point. This grid point has one CLM soil column with one CLM patch.
       !---------------------------------------------------------------

       print *, "  building patch/column mapping for g/l/c/p ..."
       call clm_subgrid()

       ! Build CLM filters to process grid points

       print *, "  building CLM filters ..."
       call setFilters (filter)

       !---------------------------------------------------------------
       ! Read tower meteorology data once to get acclimation temperature
       !---------------------------------------------------------------

       print *, "  reading tower met data ..."
       call init_acclim (it, iy, mon, bounds%begp, bounds%endp, atm2lnd_inst, &
       temperature_inst, mlcanopy_inst)

       !---------------------------------------------------------------
       ! Orbital parameters for this year
       !---------------------------------------------------------------

       print *, "  calculating orbital params ..."
       call shr_orb_params (iy, eccen, obliq, mvelp, obliqr, lambm0, mvelpp)

       !---------------------------------------------------------------
       ! Open tower meteorology file (fin1) and CLM history file (fin2)
       ! to read forcing data
       !---------------------------------------------------------------

       write (ext,'(a6,"/",i4.4,"-",i2.2,".nc")') tower_id(it),iy,mon
       fin1 = diratm(1:len(trim(diratm)))//ext(1:len(trim(ext)))

       if (tower_id(it) == 'BR-Ma2' .or. tower_id(it) == 'BR-Ji2' .or. &
          tower_id(it) == 'BR-Sa1' .or. tower_id(it) == 'BR-Sa3') then
          write (ext,'("sp49wspin_",a6,"_I1PTCLM45.clm2.h1.",i4.4,".nc")') tower_id(it),iy
       else
          write (ext,'("lp67wspinPTCLM_",a6,"_I_2000_CLM45.clm2.h1.",i4.4,".nc")') tower_id(it),iy
       end if
       fin2 = dirclm(1:len(trim(dirclm)))//tower_id(it)//"/"//ext(1:len(trim(ext)))

       !---------------------------------------------------------------
       ! Open model output files
       !---------------------------------------------------------------

       print *, "  opening model output files ..."
       write (ext,'(a6,"_",i4.4,"-",i2.2,"_flux.out")') tower_id(it),iy,mon
       fout1 = dirout(1:len(trim(dirout)))//trim(subdir)//ext(1:len(trim(ext)))
       nout1 = shr_file_getUnit()
       open (unit=nout1, file=trim(fout1), action="write")

       write (ext,'(a6,"_",i4.4,"-",i2.2,"_aux.out")') tower_id(it),iy,mon
       fout2 = dirout(1:len(trim(dirout)))//trim(subdir)//ext(1:len(trim(ext)))
       nout2 = shr_file_getUnit()
       open (unit=nout2, file=trim(fout2), action="write")

       write (ext,'(a6,"_",i4.4,"-",i2.2,"_profile.out")') tower_id(it),iy,mon
       fout3 = dirout(1:len(trim(dirout)))//trim(subdir)//ext(1:len(trim(ext)))
       nout3 = shr_file_getUnit()
       open (unit=nout3, file=trim(fout3), action="write")

       !---------------------------------------------------------------
       ! Time stepping loop
       !---------------------------------------------------------------

       print *, "  starting the time loop ..."
       do istep = 1, ntim

          ! Get current date, time, and calendar day. These are for itim (at
          ! end of the time step). itim = time index from start date.
          ! curr_calday = current calendary day (1.000 on 0Z January 1 of
          ! current year)

          itim = istep
          call get_curr_date (yr, mon, day, curr_date_tod)
          call get_curr_time (curr_time_day, curr_time_sec)
          curr_calday = get_curr_calday()

          ! Time index from reference date (0Z January 1 of current year, when calday = 1.000)

          time_indx = (curr_calday - 1._r8) * 86400._r8 / dtstep

          ! Read tower meteorology for current time slice

          call TowerMet (it, itim, fin1, bounds%begp, bounds%endp, atm2lnd_inst, mlcanopy_inst)

          ! Calculate fluxes

          call clm_drv (bounds, istep, time_indx, it, fin2)

          ! Write output file

          call output (curr_calday, nout1, nout2, nout3, mlcanopy_inst, temperature_inst)

       end do

       !---------------------------------------------------------------
       ! Close files
       !---------------------------------------------------------------

       close (nout1)
       call shr_file_freeUnit (nout1)
       close (nout2)
       call shr_file_freeUnit (nout2)
       close (nout3)
       call shr_file_freeUnit (nout3)

    end do    ! End year loop

  end subroutine CLMml_drv

  !-----------------------------------------------------------------------
  subroutine clm_subgrid()
    !
    ! !DESCRIPTION:
    ! Build patch and column mapping for CLM grid cell (g), land unit (l),
    ! column (c), and patch (p) hierarchy. The code processes one point
    ! (one grid cell with one column and one patch).
    !
    ! !USES:
    use ColumnType, only : col
    use PatchType, only : patch
    !
    ! !LOCAL VARIABLES:
    integer :: c                        ! Column index for CLM g/l/c/p hierarchy
    integer :: p                        ! Patch index for CLM g/l/c/p hierarchy
    !---------------------------------------------------------------------

    associate ( &
    pcolumn     => patch%column    , &  ! Column of patch for CLM g/l/c/p hierarchy
    plandunit   => patch%landunit  , &  ! Land unit of corresponding patch for CLM g/l/c/p hierarchy
    pgridcell   => patch%gridcell  , &  ! Gridcell of corresponding patch for CLM g/l/c/p hierarchy
    pactive     => patch%active    , &  ! True => do computations on this patch
    pwtcol      => patch%wtcol     , &  ! Weight of patch relative to column
    clandunit   => col%landunit    , &  ! Land unit of corresponding column for CLM g/l/c/p hierarchy
    cgridcell   => col%gridcell    , &  ! Gridcell of corresponding column for CLM g/l/c/p hierarchy
    npatches    => col%npatches    , &  ! Number of patches on column
    patchi      => col%patchi        &  ! Beginning patch index for column
    )

    ! Column mapping for CLM g/l/c/p hierarchy

    c = 1
    clandunit(c) = c         ! Column's land unit
    cgridcell(c) = c         ! Column's grid cell
    npatches(c) = 1          ! One patch on column
    patchi(c) = 1            ! Index for first patch

    ! Patch mapping for CLM g/l/c/p hierarchy

    p = 1
    pcolumn(p) = p           ! Patch's column
    plandunit(p) = p         ! Patch's land unit
    pgridcell(p) = p         ! Patch's grid cell
    pactive(p) = .true.      ! True => do computations on this patch
    pwtcol(p) = 1._r8        ! Weight of patch relative to column

    end associate
  end subroutine clm_subgrid

  !-----------------------------------------------------------------------
  subroutine init_acclim (it, iy, mon, begp, endp, atm2lnd_inst, &
  temperature_inst, mlcanopy_inst)
    !
    ! !DESCRIPTION:
    ! Read tower meteorology data once to get acclimation temperature
    !
    ! !USES:
    use TowerDataMod, only : tower_id
    use clm_varctl, only : diratm, ntim
    use TowerMetMod, only : TowerMet
    use PatchType, only : patch
    use atm2lndType, only : atm2lnd_type
    use TemperatureType, only : temperature_type
    use CanopyFluxesMultilayerType, only : mlcanopy_type
    !
    ! !ARGUMENTS:
    implicit none
    integer, intent(in) :: it          ! Tower site index
    integer, intent(in) :: iy          ! Tower year index
    integer, intent(in) :: mon         ! Month (1, ..., 12)
    integer, intent(in) :: begp, endp  ! First and last patch
    type(atm2lnd_type), intent(inout) :: atm2lnd_inst
    type(temperature_type), intent(out) :: temperature_inst
    type(mlcanopy_type), intent(inout) :: mlcanopy_inst
    !
    ! !LOCAL VARIABLES:
    integer  :: p                      ! Patch index for CLM g/l/c/p hierarchy
    integer  :: c                      ! Column index for CLM g/l/c/p hierarchy
    integer  :: itim                   ! Time index

    character(len=256) :: ext          ! Local file name
    character(len=256) :: fin          ! Full file name, including directory path
    !---------------------------------------------------------------------

    associate ( &
    forc_t => atm2lnd_inst%forc_t_downscaled_col , & ! Atmospheric temperature (K)
    t10    => temperature_inst%t_a10_patch         & ! Average air temperature for acclimation (K)
    )

    ! File to read

    write (ext,'(a6,"/",i4.4,"-",i2.2,".nc")') tower_id(it),iy,mon
    fin = diratm(1:len(trim(diratm)))//ext(1:len(trim(ext)))

    ! Initialize accumulator to zero

    do p = begp, endp
       t10(p) = 0._r8
    end do

    ! Loop over all time slices to read tower data

    do itim = 1, ntim

       ! Read temperature for this time slice

       call TowerMet (it, itim, fin, begp, endp, atm2lnd_inst, mlcanopy_inst)

       ! Sum temperature

       do p = begp, endp
          c = patch%column(p)
          t10(p) = t10(p) + forc_t(c)
       end do

    end do

    ! Average temperature over all time slices

    do p = begp, endp
       t10(p) = t10(p) / float(ntim)
    end do

    end associate
  end subroutine init_acclim

  !-----------------------------------------------------------------------
  subroutine output (curr_calday, nout1, nout2, nout3, mlcanopy_inst, temperature_inst)
    !
    ! !DESCRIPTION:
    ! Write output
    !
    ! !USES:
    use clm_varpar, only : ivis, inir, isun, isha
    use clm_varcon, only : tfrz
    use ColumnType, only : col
    use TemperatureType, only : temperature_type
    use CanopyFluxesMultilayerType, only : mlcanopy_type
    !
    ! !ARGUMENTS:
    implicit none
    real(r8), intent(in) :: curr_calday  ! Current calendary day (1.000 on 0Z January 1 of current year)
    integer,  intent(in) :: nout1, nout2, nout3   ! Fortran unit number
    type(temperature_type), intent(in) :: temperature_inst
    type(mlcanopy_type), intent(in) :: mlcanopy_inst
    !
    ! !LOCAL VARIABLES:
    integer  :: ic                     ! Aboveground layer index
    integer  :: mid                    ! Mid-canopy layer index
    integer  :: p                      ! Patch index for CLM g/l/c/p hierarchy
    real(r8) :: swup                   ! Reflected solar radiation (W/m2)
    real(r8) :: missing_value          ! Missing value
    real(r8) :: zero_value             ! Zero
    !---------------------------------------------------------------------

    associate ( &
    ncan        => mlcanopy_inst%ncan          , &  ! Number of layers
    nbot        => mlcanopy_inst%nbot          , &  ! Index for bottom leaf layer
    ntop        => mlcanopy_inst%ntop          , &  ! Index for top leaf layer
    lai         => mlcanopy_inst%lai           , &  ! Leaf area index of canopy (m2/m2)
    sai         => mlcanopy_inst%sai           , &  ! Stem area index of canopy (m2/m2)
    dpai        => mlcanopy_inst%dpai          , &  ! Layer plant area index (m2/m2)
    sumpai      => mlcanopy_inst%sumpai        , &  ! Cumulative plant area index (m2/m2)
    zref        => mlcanopy_inst%zref          , &  ! Reference height (m)
    tref        => mlcanopy_inst%tref          , &  ! Air temperature at reference height (K)
    uref        => mlcanopy_inst%uref          , &  ! Wind speed at reference height (m/s)
    rhref       => mlcanopy_inst%rhref         , &  ! Relative humidity at reference height (%)
    pref        => mlcanopy_inst%pref          , &  ! Air pressure at reference height (Pa)
    co2ref      => mlcanopy_inst%co2ref        , &  ! Atmospheric CO2 at reference height (umol/mol)
    o2ref       => mlcanopy_inst%o2ref         , &  ! Atmospheric O2 at reference height (mmol/mol)
    solar_zen   => mlcanopy_inst%solar_zen     , &  ! Solar zenith angle (radians)
    swskyb      => mlcanopy_inst%swskyb        , &  ! Atmospheric direct beam solar radiation (W/m2)
    swskyd      => mlcanopy_inst%swskyd        , &  ! Atmospheric diffuse solar radiation (W/m2)
    irsky       => mlcanopy_inst%irsky         , &  ! Atmospheric longwave radiation (W/m2)
    qflx_rain   => mlcanopy_inst%qflx_rain     , &  ! Rainfall (mm H2O/s = kg H2O/m2/s)
    qflx_snow   => mlcanopy_inst%qflx_snow     , &  ! Snowfall (mm H2O/s = kg H2O/m2/s)
    eref        => mlcanopy_inst%eref          , &  ! Vapor pressure at reference height (Pa)
    qref        => mlcanopy_inst%qref          , &  ! Specific humidity at reference height (kg/kg)
    rhoair      => mlcanopy_inst%rhoair        , &  ! Air density at reference height (kg/m3)
    rhomol      => mlcanopy_inst%rhomol        , &  ! Molar density at reference height (mol/m3)
    mmair       => mlcanopy_inst%mmair         , &  ! Molecular mass of air at reference height (kg/mol)
    cpair       => mlcanopy_inst%cpair         , &  ! Specific heat of air at constant pressure, at reference height (J/mol/K)
    tacclim     => mlcanopy_inst%tacclim       , &  ! Average air temperature for acclimation (K)
    btran       => mlcanopy_inst%btran         , &  ! Ball-Berry soil wetness factor (-)
    psis        => mlcanopy_inst%psis          , &  ! Weighted soil water potential (MPa)
    rsoil       => mlcanopy_inst%rsoil         , &  ! Soil hydraulic resistance (MPa.s.m2/mmol H2O)
    taf         => mlcanopy_inst%taf           , &  ! Air temperature at canopy top (K)
    wind        => mlcanopy_inst%wind          , &  ! Wind speed profile (m/s)
    wind_most   => mlcanopy_inst%wind_most     , &  ! Wind speed profile from MOST (m/s)
    tair        => mlcanopy_inst%tair          , &  ! Air temperature profile (K)
    tair_most   => mlcanopy_inst%tair_most     , &  ! Air temperature profile from MOST (K)
    eair        => mlcanopy_inst%eair          , &  ! Vapor pressure profile (Pa)
    cair        => mlcanopy_inst%cair          , &  ! Atmospheric CO2 profile (umol/mol)
    tveg        => mlcanopy_inst%tveg          , &  ! Vegetation temperature profile (K)
    fracsun     => mlcanopy_inst%fracsun       , &  ! Sunlit fraction of canopy layer
    fracsha     => mlcanopy_inst%fracsha       , &  ! Shaded fraction of canopy layer
    irleaf      => mlcanopy_inst%irleaf        , &  ! Leaf absorbed longwave radiation for canopy layer (W/m2 leaf)
    sw_prof     => mlcanopy_inst%sw_prof       , &  ! Canopy layer absorbed solar radiation (W/m2)
    rn_prof     => mlcanopy_inst%rn_prof       , &  ! Canopy layer net radiation (W/m2)
    sh_prof     => mlcanopy_inst%sh_prof       , &  ! Canopy layer sensible heat flux (W/m2)
    lh_prof     => mlcanopy_inst%lh_prof       , &  ! Canopy layer latent heat flux (W/m2)
    et_prof     => mlcanopy_inst%et_prof       , &  ! Canopy layer water vapor flux (mol H2O/m2/s)
    fc_prof     => mlcanopy_inst%fc_prof       , &  ! Canopy layer CO2 flux (umol CO2/m2/s)
    ga_prof     => mlcanopy_inst%ga_prof       , &  ! Canopy layer aerodynamic conductance for scalars (mol/m2/s)
    zs          => mlcanopy_inst%zs            , &  ! Canopy height for scalar concentration and source (m)
    lwp         => mlcanopy_inst%lwp           , &  ! Leaf water potential of canopy layer (MPa)
    lsc         => mlcanopy_inst%lsc           , &  ! Leaf-specific conductance of canopy layer (mmol H2O/m2 leaf/s/MPa)
    h2ocan      => mlcanopy_inst%h2ocan        , &  ! Canopy layer intercepted water (kg H2O/m2)
    tleaf       => mlcanopy_inst%tleaf         , &  ! Leaf temperature (K)
    rnleaf      => mlcanopy_inst%rnleaf        , &  ! Leaf net radiation (W/m2 leaf)
    shleaf      => mlcanopy_inst%shleaf        , &  ! Leaf sensible heat flux (W/m2 leaf)
    lhleaf      => mlcanopy_inst%lhleaf        , &  ! Leaf latent heat flux (W/m2 leaf)
    swleaf      => mlcanopy_inst%swleaf        , &  ! Leaf absorbed solar radiation (W/m2 leaf)
    trleaf      => mlcanopy_inst%trleaf        , &  ! Leaf transpiration flux (mol H2O/m2 leaf/s)
    evleaf      => mlcanopy_inst%evleaf        , &  ! Leaf evaporation flux (mol H2O/m2 leaf/s)
    psil        => mlcanopy_inst%psil          , &  ! Leaf water potential (MPa)
    gbh         => mlcanopy_inst%gbh           , &  ! Leaf boundary layer conductance, heat (mol/m2 leaf/s)
    gbv         => mlcanopy_inst%gbv           , &  ! Leaf boundary layer conductance, H2O (mol H2O/m2 leaf/s)
    gbc         => mlcanopy_inst%gbc           , &  ! Leaf boundary layer conductance, CO2 (mol CO2/m2 leaf/s)
    apar        => mlcanopy_inst%apar          , &  ! Leaf absorbed PAR (umol photon/m2 leaf/s)
    ac          => mlcanopy_inst%ac            , &  ! Leaf rubisco-limited gross photosynthesis (umol CO2/m2 leaf/s)
    aj          => mlcanopy_inst%aj            , &  ! Leaf RuBP-limited gross photosynthesis (umol CO2/m2 leaf/s)
    ap          => mlcanopy_inst%ap            , &  ! Leaf product-limited(C3), CO2-limited(C4) gross photosynthesis (umol CO2/m2 leaf/s)
    ag          => mlcanopy_inst%ag            , &  ! Leaf gross photosynthesis (umol CO2/m2 leaf/s)
    an          => mlcanopy_inst%an            , &  ! Leaf net photosynthesis (umol CO2/m2 leaf/s)
    rd          => mlcanopy_inst%rd            , &  ! Leaf respiration rate (umol CO2/m2 leaf/s)
    ci          => mlcanopy_inst%ci            , &  ! Leaf intercellular CO2 (umol/mol)
    cs          => mlcanopy_inst%cs            , &  ! Leaf surface CO2 (umol/mol)
    hs          => mlcanopy_inst%hs            , &  ! Leaf fractional humidity at leaf surface (dimensionless)
    vpd         => mlcanopy_inst%vpd           , &  ! Leaf vapor pressure deficit (Pa)
    gs          => mlcanopy_inst%gs            , &  ! Leaf stomatal conductance (mol H2O/m2 leaf/s)
    swveg       => mlcanopy_inst%swveg         , &  ! Absorbed solar radiation, vegetation (W/m2)
    swvegsun    => mlcanopy_inst%swvegsun      , &  ! Absorbed solar radiation, sunlit canopy (W/m2)
    swvegsha    => mlcanopy_inst%swvegsha      , &  ! Absorbed solar radiation, shaded canopy (W/m2)
    irveg       => mlcanopy_inst%irveg         , &  ! Absorbed longwave radiation, vegetation (W/m2)
    irvegsun    => mlcanopy_inst%irvegsun      , &  ! Absorbed longwave radiation, sunlit canopy (W/m2)
    irvegsha    => mlcanopy_inst%irvegsha      , &  ! Absorbed longwave radiation, shaded canopy (W/m2)
    shveg       => mlcanopy_inst%shveg         , &  ! Sensible heat flux, vegetation (W/m2)
    shvegsun    => mlcanopy_inst%shvegsun      , &  ! Sensible heat flux, sunlit canopy (W/m2)
    shvegsha    => mlcanopy_inst%shvegsha      , &  ! Sensible heat flux, shaded canopy (W/m2)
    lhveg       => mlcanopy_inst%lhveg         , &  ! Latent heat flux, vegetation (W/m2)
    lhvegsun    => mlcanopy_inst%lhvegsun      , &  ! Latent heat flux, sunlit canopy (W/m2)
    lhvesha     => mlcanopy_inst%lhvegsha      , &  ! Latent heat flux, shaded canopy (W/m2)
    etveg       => mlcanopy_inst%etveg         , &  ! Water vapor flux, vegetation (mol H2O/m2/s)
    etvegsun    => mlcanopy_inst%etvegsun      , &  ! Water vapor flux, sunlit canopy (mol H2O/m2/s)
    etvesha     => mlcanopy_inst%etvegsha      , &  ! Water vapor flux, shaded canopy (mol H2O/m2/s)
    gppveg      => mlcanopy_inst%gppveg        , &  ! Gross primary production (umol CO2/m2/s)
    gppvegsun   => mlcanopy_inst%gppvegsun     , &  ! Gross primary production, sunlit canopy (umol CO2/m2/s)
    gppvegsha   => mlcanopy_inst%gppvegsha     , &  ! Gross primary production, shaded canopy (umol CO2/m2/s)
    albcan      => mlcanopy_inst%albcan        , &  ! Albedo above canopy
    ircan       => mlcanopy_inst%ircan         , &  ! Upward longwave radiation above canopy (W/m2)
    ustar       => mlcanopy_inst%ustar         , &  ! Friction velocity (m/s)
    obu         => mlcanopy_inst%obu           , &  ! Obukhov length (m)
    rnet        => mlcanopy_inst%rnet          , &  ! Net radiation (W/m2)
    stflx       => mlcanopy_inst%stflx         , &  ! Canopy storage heat flux (W/m2)
    shflx       => mlcanopy_inst%shflx         , &  ! Sensible heat flux (W/m2)
    lhflx       => mlcanopy_inst%lhflx         , &  ! Latent heat flux (W/m2)
    etflx       => mlcanopy_inst%etflx         , &  ! Water vapor flux (mol H2O/m2/s)
    fracminlwp  => mlcanopy_inst%fracminlwp    , &  ! Fraction of canopy with lwp < minlwp
    rnsoi       => mlcanopy_inst%rnsoi         , &  ! Net radiation, ground (W/m2)
    shsoi       => mlcanopy_inst%shsoi         , &  ! Sensible heat flux, ground (W/m2)
    lhsoi       => mlcanopy_inst%lhsoi         , &  ! Latent heat flux, ground (W/m2)
    gsoi        => mlcanopy_inst%gsoi          , &  ! Soil heat flux (W/m2)
    swsoi       => mlcanopy_inst%swsoi         , &  ! Absorbed solar radiation, ground (W/m2)
    irsoi       => mlcanopy_inst%irsoi         , &  ! Absorbed longwave radiation, ground (W/m2)
    etsoi       => mlcanopy_inst%etsoi         , &  ! Water vapor flux, ground (mol H2O/m2/s)
    dz          => col%dz                      , &  ! Soil layer thickness (m)
    z           => col%z                       , &  ! Soil layer depth (m)
    zi          => col%zi                        , &  ! Soil layer depth at layer interface (m)
    t_soisno    => temperature_inst%t_soisno_col , &  ! Soil temperature (K)
    t_grnd      => temperature_inst%t_grnd_col     &  ! Ground surface temperature (K)
    )

    missing_value = -999._r8
    zero_value = 0._r8

    p = 1
    swup = albcan(p,ivis)*(swskyb(p,ivis)+swskyd(p,ivis)) + albcan(p,inir)*(swskyb(p,inir)+swskyd(p,inir))
    write (nout1,'(13f10.3)') rnet(p), stflx(p), shflx(p), lhflx(p), gppveg(p), ustar(p), &
    swup, ircan(p), taf(p), gsoi(p), rnsoi(p), shsoi(p), lhsoi(p)

    mid = nbot(p) + (ntop(p)-nbot(p)+1)/2 - 1
    write (nout2,'(f10.4,5f10.3)') btran(p), lsc(p,ntop(p)), psis(p), lwp(p,ntop(p)), &
    lwp(p,mid), fracminlwp(p)

!   go to 100
    do ic = ntop(p), 0, -1
       if (ic >= nbot(p)) then ! Leaf layer
          write (nout3,'(f10.4,14f10.3)') curr_calday, zs(p,ic), dpai(p,ic), &
          rn_prof(p,ic), sh_prof(p,ic), lh_prof(p,ic), fc_prof(p,ic), &
          apar(p,ic,isun)*fracsun(p,ic)+apar(p,ic,isha)*fracsha(p,ic), &
          gs(p,ic,isun)*fracsun(p,ic)+gs(p,ic,isha)*fracsha(p,ic), lwp(p,ic), &
          tveg(p,ic,isun)*fracsun(p,ic)+tveg(p,ic,isha)*fracsha(p,ic), &
          wind(p,ic), tair(p,ic)-tref(p), eair(p,ic)/1000._r8, rhomol(p)/ga_prof(p,ic)
!         wind(p,ic), tair(p,ic)-tfrz, eair(p,ic)/1000._r8, rhomol(p)/ga_prof(p,ic)
       else ! Non-leaf layer or ground
          write (nout3,'(f10.4,14f10.3)') curr_calday, zs(p,ic), zero_value, &
          rn_prof(p,ic), sh_prof(p,ic), lh_prof(p,ic), fc_prof(p,ic), &
          missing_value, missing_value, missing_value, missing_value, &
          wind(p,ic), tair(p,ic)-tref(p), eair(p,ic)/1000._r8, rhomol(p)/ga_prof(p,ic)
!         wind(p,ic), tair(p,ic)-tfrz, eair(p,ic)/1000._r8, rhomol(p)/ga_prof(p,ic)
       end if
    end do

100 continue

    go to 200
    do ic = ncan(p), 1, -1
!      wind_most(p,ic) = missing_value
!      tair_most(p,ic) = missing_value
       if (dpai(p,ic) > 0._r8) then ! Leaf layer
          write (nout3,'(f10.4,14f10.3)') curr_calday, zs(p,ic), dpai(p,ic), &
          wind(p,ic), tair(p,ic), wind_most(p,ic), tair_most(p,ic), &
          tveg(p,ic,isun)*fracsun(p,ic)+tveg(p,ic,isha)*fracsha(p,ic)
       else
          write (nout3,'(f10.4,14f10.3)') curr_calday, zs(p,ic), zero_value, &
          wind(p,ic), tair(p,ic), wind_most(p,ic), tair_most(p,ic), &
          missing_value
       end if
    end do
200 continue

    end associate
  end subroutine output

end module CLMml_driver
