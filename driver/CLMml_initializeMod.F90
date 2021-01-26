module CLMml_initializeMod

  !-----------------------------------------------------------------------
  ! !DESCRIPTION:
  ! Performs model initialization
  !
  ! !USES:
  use shr_kind_mod, only : r8 => shr_kind_r8
  use abortutils, only : endrun, handle_err
  use ColumnType, only : col
  use PatchType, only : patch
  use SoilStateType, only : soilstate_type
  use WaterStateType, only : waterstate_type
  use TemperatureType, only : temperature_type
  use SurfaceAlbedoType, only : surfalb_type
  use CanopyFluxesMultilayerType, only : mlcanopy_type
  !
  ! !PUBLIC TYPES:
  implicit none
  include 'netcdf.inc'
  !
  ! !PUBLIC MEMBER FUNCTIONS:
  public :: initialize          ! Model initialization
  !
  ! !PRIVATE MEMBER FUNCTIONS:
  private :: iniTimeConst       ! Initialize time constant variables
  private :: TowerVeg           ! Read tower leaf area and define canopy layers
  private :: SoilInit           ! Initialize soil temperature and moisture profile
  !-----------------------------------------------------------------------

contains

  !-----------------------------------------------------------------------
  subroutine initialize (it, iy, begp, endp, begc, endc, soilstate_inst, &
  waterstate_inst, temperature_inst, surfalb_inst, mlcanopy_inst)
    !
    ! !DESCRIPTION:
    ! Initialization
    !
    ! !USES:
    use pftconMod, only : pftcon
    use clm_varctl, only : dirclm
    use TowerDataMod, only : tower_id
    use clm_time_manager
    !
    ! !ARGUMENTS:
    implicit none
    integer, intent(in) :: it                    ! Tower site index
    integer, intent(in) :: iy                    ! Tower year index
    integer, intent(in) :: begp, endp            ! First and last patch
    integer, intent(in) :: begc, endc            ! First and last column
    type(soilstate_type), intent(out) :: soilstate_inst
    type(waterstate_type), intent(out) :: waterstate_inst
    type(temperature_type), intent(out) :: temperature_inst
    type(surfalb_type), intent(out) :: surfalb_inst
    type(mlcanopy_type), intent(out) :: mlcanopy_inst
    !
    ! !LOCAL VARIABLES:
    integer  :: strt                             ! Time slice of data to retrieve
    character(len=256) :: ext                    ! Local file name
    character(len=256) :: fin                    ! Full file name, including directory path
    !---------------------------------------------------------------------

    ! Read list of PFTs and their parameter values

    call pftcon%Init()

    ! Initialize time constant variables

    call iniTimeConst (it, begc, endc, soilstate_inst, surfalb_inst)

    ! Time varying data are initialized from CLM history file. Define the
    ! file for this tower site/year and then choose first time slice to read.

    ! (a) Get file name for history file

    if (tower_id(it) == 'BR-Ma2' .or. tower_id(it) == 'BR-Ji2' .or. &
       tower_id(it) == 'BR-Sa1' .or. tower_id(it) == 'BR-Sa3') then
       write (ext,'("sp49wspin_",a6,"_I1PTCLM45.clm2.h1.",i4.4,".nc")') tower_id(it),iy
    else
       write (ext,'("lp67wspinPTCLM_",a6,"_I_2000_CLM45.clm2.h1.",i4.4,".nc")') tower_id(it),iy
    end if
    fin = dirclm(1:len(trim(dirclm)))//tower_id(it)//"/"//ext(1:len(trim(ext)))

    ! (b) For first timeslice (itim = 1) calculate time index from reference date
    ! (0Z January 1 of current year, when calday = 1.000)

    itim = 1
    curr_calday = get_curr_calday()
    strt = (curr_calday - 1._r8) * 86400._r8 / dtstep

    ! Read CLM history file to initialize canopy leaf area and define layers

    call TowerVeg (it, strt, fin, begp, endp, mlcanopy_inst)

    ! Read CLM history file to initialize soil temperature and moisture profile

    call SoilInit (strt, fin, begc, endc, soilstate_inst, waterstate_inst, temperature_inst)

  end subroutine initialize

  !-----------------------------------------------------------------------
  subroutine iniTimeConst (it, begc, endc, soilstate_inst, surfalb_inst)
    !
    ! !DESCRIPTION:
    ! Initialization
    !
    ! !USES:
    use clm_varpar, only : nlevgrnd, nlevsoi
    use clm_varcon, only : denh2o, mmh2o, grav
    use clm_varctl, only : iulog
    use TowerDataMod, only : tower_tex, tower_isoicol
    use SurfaceAlbedoMod, only: isoicol
    use SoilTexMod
    !
    ! !ARGUMENTS:
    implicit none
    integer, intent(in) :: it                    ! Tower site index
    integer, intent(in) :: begc, endc            ! First and last column
    type(soilstate_type), intent(out) :: soilstate_inst
    type(surfalb_type), intent(out) :: surfalb_inst
    !
    ! !LOCAL VARIABLES:
    integer :: c                                 ! Column index for CLM g/l/c/p hierarchy
    integer :: j                                 ! Soil layer index
    integer :: m, im
    real(r8) :: scalez = 0.025_r8                ! Soil layer thickness discretization (m)
    character(len=15) :: soiltyp                 ! Soil texture class
    !---------------------------------------------------------------------

    associate ( &
    dz        => col%dz                      , &  ! Soil layer thickness (m)
    z         => col%z                       , &  ! Soil layer depth (m)
    zi        => col%zi                      , &  ! Soil layer depth at layer interface (m)
    cellsand  => soilstate_inst%cellsand_col , &  ! Soil layer percent sand
    cellclay  => soilstate_inst%cellclay_col , &  ! Soil layer percent clay
    cellorg   => soilstate_inst%cellorg_col  , &  ! Soil layer organic fraction
    watsat    => soilstate_inst%watsat_col   , &  ! Soil layer volumetric water content at saturation (porosity)
    sucsat    => soilstate_inst%sucsat_col   , &  ! Soil layer suction (negative matric potential) at saturation (mm)
    hksat     => soilstate_inst%hksat_col    , &  ! Soil layer hydraulic conductivity at saturation (mm H2O/s)
    bsw       => soilstate_inst%bsw_col        &  ! Soil layer Clapp and Hornberger "b" parameter
    )

    do c = begc, endc

       ! Define CLM layer structure for soil

       do j = 1, nlevgrnd
          z(c,j) = scalez * (exp(0.5_r8*(j-0.5_r8)) - 1._r8)     ! Layer depths
       end do

       dz(c,1) = 0.5_r8 * (z(c,1) + z(c,2))                      ! Layer thickness
       do j = 2, nlevgrnd-1
          dz(c,j)= 0.5_r8 * (z(c,j+1) - z(c,j-1))
       end do
       dz(c,nlevgrnd) = z(c,nlevgrnd) - z(c,nlevgrnd-1)

       zi(c,0) = 0._r8
       do j = 1, nlevgrnd-1
          zi(c,j) = 0.5_r8 * (z(c,j) + z(c,j+1))                 ! Interface depths
       end do
       zi(c,nlevgrnd) = z(c,nlevgrnd) + 0.5_r8 * dz(c,nlevgrnd)

       ! Soil texture

       soiltyp = tower_tex(it)

       ! Organic fraction

       do j = 1, nlevsoi
          cellorg(c,j) = 0.02_r8
       end do

       ! Soil hydraulic properties: first find soil texture type

       m = 0
       do im = 1, ntex
          if (soiltyp == soil_tex(im)) then
             m = im
             exit
          else
             cycle
          end if
       end do
       if (m == 0) then
          write (iulog,*) 'iniTimeConst error: soil type = ',soiltyp, ' not found for c = ',c
          call endrun()
       end if

       ! Now assign texture based variables

       do j = 1, nlevsoi
          cellsand(c,j) = sand_tex(m) * 100._r8
          cellclay(c,j) = clay_tex(m) * 100._r8
       end do

       do j = 1, nlevgrnd
          watsat(c,j) = watsat_tex(m)
          sucsat(c,j) = -smpsat_tex(m)
          hksat(c,j) = hksat_tex(m) / 60._r8                                 ! mm/min -> mm/s
!         hksat(c,j) = hksat(c,j) * 1.e-03_r8 / (denh2o*grav*1.e-06_r8)      ! mm/s -> m/s -> m2/s/MPa
!         hksat(c,j) = hksat(c,j) * denh2o / mmh2o * 1000._r8                ! m2/s/MPa -> mmol/m/s/MPa
          bsw(c,j) = bsw_tex(m)
       end do

       ! Soil color

       isoicol(c) = tower_isoicol(it)

    end do

    end associate
  end subroutine iniTimeConst

  !-----------------------------------------------------------------------
  subroutine TowerVeg (it, strt, ncfilename, begp, endp, mlcanopy_inst)
    !
    ! !DESCRIPTION:
    ! Read CLM history file to initialize canopy leaf area and define layers
    !
    ! !USES:
    use clm_varctl, only : use_tower, pad_type, iulog, use_ncan
    use clm_varpar, only : mxpft, nlevcan
    use TowerDataMod, only : tower_pft, tower_canht, tower_id, tower_ht
    use MathToolsMod, only : beta_function, log_gamma_function
    !
    ! !ARGUMENTS:
    implicit none
    integer, intent(in) :: it                    ! Tower site index
    integer, intent(in) :: strt                  ! Current time slice of data to retrieve from CLM history file
    character(len=*), intent(in) :: ncfilename   ! CLM netcdf filename
    integer, intent(in) :: begp, endp            ! First and last patch
    type(mlcanopy_type), intent(out) :: mlcanopy_inst
    !
    ! !LOCAL VARIABLES:
    integer  :: ic                               ! Aboveground layer index
    integer  :: p                                ! Patch index for CLM g/l/c/p hierarchy
    integer  :: ncid                             ! netcdf file ID
    integer  :: status                           ! Function return status
    integer  :: varid                            ! netcdf variable id
    integer  :: start2(2), count2(2)             ! Start and count arrays for reading 2-D data from netcdf files
    integer  :: n_zref                           ! Number of above-canopy layers
    real(r8) :: dz_zref                          ! Atmospheric reference height - canopy height (m)
    real(r8) :: dht                              ! Height increment (m)
    real(r8) :: htop(0:mxpft)                    ! CLM canopy top height, by PFT (m)
    real(r8) :: elai_mod(1,1,1)                  ! Leaf area index (m2/m2)
    real(r8) :: esai_mod(1,1,1)                  ! Stem area index (m2/m2)
    real(r8) :: pbeta, qbeta                     ! Parameters for the beta function
    real(r8) :: beta                             ! Value of the beta function
    real(r8) :: zl, zu                           ! Bottom and top heights for a canopy layer (m)
    integer  :: num_int                          ! Number of layers for numerical integration of LAI between zl and zu
    real(r8) :: dz_int                           ! Height increment for numerical integration of LAI between zl and zu (m)
    integer  :: ic_int                           ! Do loop index for numerical integration of LAI between zl and zu
    real(r8) :: z_int                            ! Height for numerical integration of LAI between zl and zu (m)
    real(r8) :: zrel                             ! Height relative to canopy top (z_int/ztop)
    real(r8) :: beta_pdf                         ! Value for beta probibility density function
    real(r8) :: lad                              ! Layer leaf area density (m2/m3)
    real(r8) :: lai_err, sai_err, pai_err        ! Sum of leaf, stem, or plant area index over all layers (m2/m2)
    real(r8) :: lai_miss                         ! Missing leaf area after imposing dpai_min constraint on layers (m2/m2)
    real(r8) :: sai_miss                         ! Missing stem area after imposing dpai_min constraint on layers (m2/m2)
    real(r8), parameter :: dpai_min = 0.01_r8    ! Minimum plant area to be considered a vegetation layer (m2/m2)

    ! CLM top canopy height, by PFTs

    data (htop(i),i=1,mxpft) / 17._r8, 17._r8, 14._r8, 35._r8, 35._r8, 18._r8, 20._r8, 20._r8, &
                               0.5_r8, 0.5_r8, 0.5_r8, 0.5_r8, 0.5_r8, 0.5_r8, 0.5_r8, 0.5_r8 /
    !---------------------------------------------------------------------

    associate ( &
    zref         => mlcanopy_inst%zref          , &  ! Reference height (m)
    ztop         => mlcanopy_inst%ztop          , &  ! Canopy height (m)
    lai          => mlcanopy_inst%lai           , &  ! Leaf area index of canopy (m2/m2)
    sai          => mlcanopy_inst%sai           , &  ! Stem area index of canopy (m2/m2)
    root_biomass => mlcanopy_inst%root_biomass  , &  ! Fine root biomass (g biomass / m2)
    ncan         => mlcanopy_inst%ncan          , &  ! Number of layers
    nbot         => mlcanopy_inst%nbot          , &  ! Index for bottom leaf layer
    ntop         => mlcanopy_inst%ntop          , &  ! Index for top leaf layer
    dlai         => mlcanopy_inst%dlai          , &  ! Layer leaf area index (m2/m2)
    dsai         => mlcanopy_inst%dsai          , &  ! Layer stem area index (m2/m2)
    dpai         => mlcanopy_inst%dpai          , &  ! Layer plant area index (m2/m2)
    zs           => mlcanopy_inst%zs            , &  ! Canopy height for scalar concentration and source (m)
    zw           => mlcanopy_inst%zw              &  ! Canopy heights at layer interfaces (m)
    )

    !---------------------------------------------------------------------
    ! Read CLM history file
    !---------------------------------------------------------------------

    status = nf_open(ncfilename, nf_nowrite, ncid)
    if (status /= nf_noerr) call handle_err(status, ncfilename)

    ! Dimensions in FORTRAN are in column major order: the first array index varies the most rapidly.
    ! In NetCDF file the dimensions appear in the opposite order: lat, lon (2-D); time, lat, lon (3-D);
    ! time, levgrnd, lat, lon (4-D)

    start2 = (/ 1, strt /)
    count2 = (/ 1, 1 /)

    ! elai(nlndgrid, ntime): leaf area index (m2/m2)

    status = nf_inq_varid(ncid, "ELAI", varid)
    if (status /= nf_noerr) call handle_err(status, "ELAI")

    status = nf_get_vara_double(ncid, varid, start2, count2, elai_mod)
    if (status /= nf_noerr) call handle_err(status, "elai_mod")

    ! esai(nlndgrid, ntime): stem area index (m2/m2)

    status = nf_inq_varid(ncid, "ESAI", varid)
    if (status /= nf_noerr) call handle_err(status, "ESAI")

    status = nf_get_vara_double(ncid, varid, start2, count2, esai_mod)
    if (status /= nf_noerr) call handle_err(status, "esai_mod")

    if (pad_type == 0) esai_mod(1,1,1) = 0._r8

    status = nf_close(ncid)

    !---------------------------------------------------------------------
    ! Assign leaf area, stem area, PFT, and fine root biomass
    !---------------------------------------------------------------------

    do p = begp, endp
       lai(p) = elai_mod(1,1,1)
       sai(p) = esai_mod(1,1,1)
       patch%itype(p) = tower_pft(it)
       root_biomass(p) = 500._r8     ! Fine root biomass (g biomass / m2)
    end do

    !---------------------------------------------------------------------
    ! Define canopy layers
    !---------------------------------------------------------------------

    do p = begp, endp

       ! Top height of canopy

       ztop(p) = htop(patch%itype(p))

       ! Overwrite canopy height with tower data (if available)

       if (use_tower .and. tower_canht(it) /= -999._r8) then
          ztop(p) = tower_canht(it)
       end if

       ! Determine number of within-canopy layers by specifying height increment

       if (ztop(p) > 2._r8) then
          dht = 0.5_r8
       else
          dht = 0.1_r8
       end if

       ! Define ntop as the number of within-canopy layers and adjust dht if needed

       ntop(p) = nint(ztop(p) / dht)
       dht = ztop(p) / float(ntop(p))

       ! Calculate heights at layer interfaces (zw). These are the heights
       ! for the conductance between two scalar concentrations. They are
       ! defined for ic = 0 (ground) to ic = ntop (top of the canopy).

       ic = ntop(p)
       zw(p,ic) = ztop(p)
       do ic = ntop(p)-1, 0, -1
          zw(p,ic) = zw(p,ic+1) - dht
       end do

       if (zw(p,0) > 1.e-10_r8) then
          write (iulog,*) 'TowerVeg error: zw(p,0) improperly defined'
          call endrun()
       end if

       ! Now calculate the above-canopy layers and their heights

       if (use_ncan) then

          if (use_tower .and. nint(tower_ht(it)) /= -999) then
             zref(p) = tower_ht(it)
          else if (use_tower .and. tower_id(it) == 'US-ARM') then
             zref(p) = 30._r8
          else
             write (iulog,*) 'TowerVeg error: invalid tower height'
             call endrun()
          end if

          dz_zref = zref(p) - ztop(p)
          if (tower_id(it) == 'US-ARM') dht = 1._r8
          n_zref = nint(dz_zref / dht)
          dht = dz_zref / float(n_zref)
          ncan(p) = ntop(p) + n_zref

          if (ncan(p) > nlevcan) then
             write (iulog,*) 'TowerVeg error: ncan > nlevcan'
             call endrun()
          end if

          ic = ncan(p)
          zw(p,ic) = zref(p)
          do ic = ncan(p)-1, ntop(p)+1, -1
             zw(p,ic) = zw(p,ic+1) - dht
          end do

       else

          ncan(p) = ntop(p)

       end if

       ! Determine heights of the scalar concentration and scalar source
       ! (zs). These are physically centered between the conductance points
       ! (i.e., in the middle of the layer).

       zs(p,0) = 0._r8
       do ic = 1, ncan(p)
          zs(p,ic) = 0.5_r8 * (zw(p,ic) + zw(p,ic-1))
       end do

       ! Determine plant area index (PAI; lai+sai) increment for each layer.
       ! First calculate terms for the beta function.

       if (patch%itype(p) == 1) then
          pbeta = 11.5_r8    ! Pine (temperate NET)
          qbeta = 3.5_r8
       end if
       if (patch%itype(p) == 2 .or. patch%itype(p) == 7) then
          pbeta = 3.5_r8     ! BDT, spruce (boreal NET)
          qbeta = 2.0_r8
       end if
       if (patch%itype(p) == 13 .or. patch%itype(p) == 15) then
          pbeta = 2.5_r8     ! Grass, crop
          qbeta = 2.5_r8
       end if
       if (patch%itype(p) == 4) then
          pbeta = 3.5_r8     ! tropical BET
          qbeta = 2.0_r8
       end if
       if (tower_id(it) == 'US-Me2') then
          pbeta = 11.5_r8    ! Pine (temperate NET)
          qbeta = 3.5_r8
       end if

       beta = beta_function (pbeta, qbeta)

       ! Now calculate leaf+stem area at each height by numerically integrating
       ! the leaf area density distribution between the bottom and top
       ! heights for that layer and then adding the stem area distribution

       pai_err = 0._r8

       do ic = 1, ntop(p)
          zl = zw(p,ic-1)
          zu = zw(p,ic)

          dlai(p,ic) = 0._r8
          dsai(p,ic) = 0._r8

          ! Use beta distribution with numerical integration between zl and zu

          if (pad_type == 1) then
             num_int = 100                         ! 100 sublayers for numerical integration
             dz_int = (zu - zl) / float(num_int)   ! dz for numerical integration
             do ic_int = 1, num_int

                if (ic_int == 1) then
                   z_int = zl + 0.5_r8 * dz_int
                else
                   z_int = z_int + dz_int
                end if
                zrel = min(z_int/ztop(p), 1._r8)
                beta_pdf = (zrel**(pbeta-1._r8) * (1._r8 - zrel)**(qbeta-1._r8)) / beta

                ! Leaf area

                lad = (lai(p) / ztop(p)) * beta_pdf
                dlai(p,ic) = dlai(p,ic) + lad * dz_int

                ! Stem area

                lad = (sai(p) / ztop(p)) * beta_pdf
                dsai(p,ic) = dsai(p,ic) + lad * dz_int

             end do
          end if

          ! Use a uniform profile

          if (pad_type == 2) then
             dlai(p,ic) = (lai(p) / ztop(p)) * (zu - zl)
             dsai(p,ic) = (sai(p) / ztop(p)) * (zu - zl)
          end if

          ! Add stem area to leaf area and save this as the plant area for the canopy layer

          dpai(p,ic) = dlai(p,ic) + dsai(p,ic)

       end do

       ! Check to make sure sum of numerical PAI matches canopy PAI

       pai_err = sum(dpai(p,1:ntop(p)))
       if (abs(pai_err - (lai(p)+sai(p))) > 1.e-06_r8) then
          write (iulog,*) 'TowerVeg error: plant area does not match observations'
          call endrun()
       end if

       ! Check to make sure there are no non-plant layers at the top of the canopy.
       ! Non-plant layers must all be at the bottom of the canopy for the code to work correctly.

       if (dpai(p,ntop(p)) < dpai_min) then
          write (iulog,*) 'TowerVeg error: plant area at canopy top is < minimum'
          call endrun()
       end if

       ! Now determine the number of layers with plant area. Set the layers with
       ! plant area < dpai_min to zero and sum the "missing" leaf and stem area
       ! so that this can be distributed back across the vegetation layers.

       lai_miss = 0._r8; sai_miss = 0._r8
       do ic = 1, ntop(p)
          if (dpai(p,ic) < dpai_min) then
             lai_miss = lai_miss + dlai(p,ic)
             sai_miss = sai_miss + dsai(p,ic)
             dlai(p,ic) = 0._r8
             dsai(p,ic) = 0._r8
             dpai(p,ic) = 0._r8
          end if
       end do

       ! Distribute the missing leaf area across vegetation layers
       ! in proportion to the leaf area profile

       if (lai_miss > 0._r8) then
          lai_err = sum(dlai(p,1:ntop(p)))
          do ic = 1, ntop(p)
             dlai(p,ic) = dlai(p,ic) + lai_miss * (dlai(p,ic) / lai_err)
          end do
       end if

       ! Now do the same for stem area

       if (sai_miss > 0._r8) then
          sai_err = sum(dsai(p,1:ntop(p)))
          do ic = 1, ntop(p)
             dsai(p,ic) = dsai(p,ic) + sai_miss * (dsai(p,ic) / sai_err)
          end do
       end if

       ! Add leaf and stem area and save as plant area

       do ic = 1, ntop(p)
          dpai(p,ic) = dlai(p,ic) + dsai(p,ic)
       end do

       ! Find the lowest leaf layer

       nbot(p) = 0
       do ic = ntop(p), 1, -1
          if (dpai(p,ic) > 0._r8) nbot(p) = ic
       end do

       if (nbot(p) == 0) then
          write (iulog,*) 'TowerVeg error: nbot not defined'
          call endrun()
       end if

       ! Error check

       lai_err = sum(dlai(p,1:ntop(p)))
       sai_err = sum(dsai(p,1:ntop(p)))
       pai_err = sum(dpai(p,1:ntop(p)))

       if (abs(lai_err - lai(p)) > 1.e-06_r8) then
          write (iulog,*) 'TowerVeg error: leaf area does not match observations after redistribution'
          call endrun()
       end if

       if (abs(sai_err - sai(p)) > 1.e-06_r8) then
          write (iulog,*) 'TowerVeg error: stem area does not match observations after redistribution'
          call endrun()
       end if

       if (abs(pai_err - (lai(p)+sai(p))) > 1.e-06_r8) then
          write (iulog,*) 'TowerVeg error: plant area does not match observations after redistribution'
          call endrun()
       end if

       ! Zero out above-canopy layers

       do ic = ntop(p)+1, ncan(p)
          dlai(p,ic) = 0._r8
          dsai(p,ic) = 0._r8
          dpai(p,ic) = 0._r8
       end do

    end do

    end associate
  end subroutine TowerVeg

  !-----------------------------------------------------------------------
  subroutine SoilInit (strt, ncfilename, begc, endc, soilstate_inst, &
  waterstate_inst, temperature_inst)
    !
    ! !DESCRIPTION:
    ! Initialize soil temperature and soil moisture profile from CLM netcdf history file
    !
    ! !USES:
    use clm_varcon, only : denh2o
    use clm_varpar, only : nlevgrnd
    !
    ! !ARGUMENTS:
    implicit none
    integer, intent(in) :: strt                ! Current time slice of data to retrieve from CLM history file
    character(len=*), intent(in) :: ncfilename ! CLM netcdf filename
    integer, intent(in) :: begc, endc          ! First and last column
    type(soilstate_type), intent(inout) :: soilstate_inst
    type(waterstate_type), intent(out) :: waterstate_inst
    type(temperature_type), intent(out) :: temperature_inst
    !
    ! !LOCAL VARIABLES:
    integer  :: c                              ! Column index for CLM g/l/c/p hierarchy
    integer  :: j                              ! Soil layer index
    integer  :: ncid                           ! netcdf file ID
    integer  :: status                         ! Function return status
    integer  :: varid                          ! netcdf variable id
    integer  :: start3(3), count3(3)           ! Start and count arrays for reading 3-D data from netcdf files
    real(r8) :: tsoi_loc(1,1,nlevgrnd)         ! CLM: soil temperature (K)
    real(r8) :: h2osoi_loc(1,1,nlevgrnd)       ! CLM: volumetric soil moisture (m3/m3)
    real(r8) :: s                              ! Soil layer water content relative to saturation (fraction)
    !---------------------------------------------------------------------

    associate ( &
    dz          => col%dz                         , &  ! Soil layer thickness (m)
    watsat      => soilstate_inst%watsat_col      , &  ! Soil layer volumetric water content at saturation (porosity)
    sucsat      => soilstate_inst%sucsat_col      , &  ! Soil layer suction (negative matric potential) at saturation (mm)
    bsw         => soilstate_inst%bsw_col         , &  ! Soil layer Clapp and Hornberger "b" parameter
                                                       ! *** Output ***
    t_soisno    => temperature_inst%t_soisno_col  , &  ! Soil temperature (K)
    t_grnd      => temperature_inst%t_grnd_col    , &  ! Ground surface temperature (K)
    h2osoi_vol  => waterstate_inst%h2osoi_vol_col , &  ! Soil layer volumetric water content (m3/m3)
    h2osoi_ice  => waterstate_inst%h2osoi_ice_col , &  ! Soil layer ice lens (kg/m2)
    h2osoi_liq  => waterstate_inst%h2osoi_liq_col , &  ! Soil layer liquid water (kg/m2)
    smp_l       => soilstate_inst%smp_l_col         &  ! Soil layer matric potential (mm)
    )

    ! Open file

    status = nf_open(ncfilename, nf_nowrite, ncid)
    if (status /= nf_noerr) call handle_err(status, ncfilename)

    ! Dimensions in FORTRAN are in column major order: the first array index
    ! varies the most rapidly. In NetCDF file the dimensions appear in the
    ! opposite order: lat, lon (2-D); time, lat, lon (3-D); time, levgrnd, lat, lon (4-D)

    start3 = (/ 1,  1, strt /)
    count3 = (/ 1, nlevgrnd, 1 /)

    ! Read TSOI(nlndgrid, nlevgrnd, ntime): soil temperature

    status = nf_inq_varid(ncid, "TSOI", varid)
    if (status /= nf_noerr) call handle_err(status, "TSOI")

    status = nf_get_vara_double(ncid, varid, start3, count3, tsoi_loc)
    if (status /= nf_noerr) call handle_err(status, "tsoi_loc")

    ! Read H2OSOI(nlndgrid, nlevgrnd, ntime): volumetric soil water

    status = nf_inq_varid(ncid, "H2OSOI", varid)
    if (status /= nf_noerr) call handle_err(status, "H2OSOI")

    status = nf_get_vara_double(ncid, varid, start3, count3, h2osoi_loc)
    if (status /= nf_noerr) call handle_err(status, "h2osoi_loc")

    ! Close file

    status = nf_close(ncid)

    ! Copy data to model variables

    do c = begc, endc
       do j = 1, nlevgrnd
          t_soisno(c,j) = tsoi_loc(1,1,j)
          h2osoi_vol(c,j) = h2osoi_loc(1,1,j)
          h2osoi_liq(c,j) = h2osoi_vol(c,j) * dz(c,j) * denh2o
          h2osoi_ice(c,j) = 0._r8
          s = max(min(h2osoi_vol(c,j)/watsat(c,j), 1._r8), 0.01_r8)
          smp_l(c,j) = -sucsat(c,j) * s**(-bsw(c,j))
          smp_l(c,j) = max(smp_l(c,j), -1.e08_r8)
       end do
       t_grnd(c) = t_soisno(c,1)
    end do

    end associate
  end subroutine SoilInit

end module CLMml_initializeMod
