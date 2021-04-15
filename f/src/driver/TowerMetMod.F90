module TowerMetMod

  !-----------------------------------------------------------------------
  ! !DESCRIPTION:
  ! Read tower meteorology and assign this forcing to all CLM patches
  !
  ! !USES:
  use shr_kind_mod, only : r8 => shr_kind_r8
  use abortutils,   only : endrun, handle_err
  use PatchType, only : patch
  use atm2lndType, only: atm2lnd_type
  use CanopyFluxesMultilayerType, only : mlcanopy_type
  !
  ! !PUBLIC TYPES:
  implicit none
  include 'netcdf.inc'
  !
  ! !PUBLIC MEMBER FUNCTIONS:
  public :: TowerMet           ! Read atmospheric forcing
  !
  ! !PRIVATE MEMBER FUNCTIONS:
  private :: readTowerMet      ! Read variables from netcdf file
  !-----------------------------------------------------------------------

contains

  !-----------------------------------------------------------------------
  subroutine TowerMet (it, strt, ncfilename, begp, endp, atm2lnd_inst, mlcanopy_inst)
    !
    ! !DESCRIPTION:
    ! Read atmospheric input variables and assign this forcing to all CLM patches
    !
    ! !USES:
    use clm_varpar, only : ivis, inir
    use clm_varcon, only : sb, mmh2o, mmdry, pi => rpi
    use clm_varctl, only : use_tower, iulog
    use TowerDataMod, only : tower_ht, tower_lat, tower_lon
    use WaterVaporMod, only : SatVap
    use clm_varorb, only : eccen, mvelpp, lambm0, obliqr
    use shr_orb_mod, only : shr_orb_decl, shr_orb_cosz
    use clm_time_manager, only : get_prev_calday
    !
    ! !ARGUMENTS:
    implicit none
    integer, intent(in) :: it                    ! Tower site index
    integer, intent(in) :: strt                  ! Time slice of data to retrieve from tower meterology
    character(len=*), intent(in) :: ncfilename   ! Tower meteorology netcdf filename
    integer, intent(in) :: begp, endp            ! First and last patch
    type(atm2lnd_type), intent(out) :: atm2lnd_inst
    type(mlcanopy_type), intent(inout) :: mlcanopy_inst
    !
    ! !LOCAL VARIABLES:
    real(r8) :: eair                   ! Atmospheric vapor pressure (Pa)
    real(r8) :: esat                   ! Saturation vapor pressure (Pa)
    real(r8) :: desat                  ! Derivative of saturation vapor pressure (Pa/K)
    real(r8) :: emiss                  ! Atmospheric emissivity
    real(r8) :: fsds                   ! Solar radiation (W/m2)
    real(r8) :: fsds_vis               ! Solar radiation, visible waveband (W/m2)
    real(r8) :: fsds_nir               ! Solar radiation, near-infrared waveband (W/m2)
    real(r8) :: rvis                   ! Direct beam fraction of solar radiation
    real(r8) :: rnir                   ! Direct beam fraction of solar radiation
    integer  :: p                      ! Patch index for CLM g/l/c/p hierarchy
    integer  :: c                      ! Column index for CLM g/l/c/p hierarchy
    integer  :: g                      ! Gridcell index for CLM g/l/c/p hierarchy
    real(r8) :: declin                 ! Solar declination angle (radians)
    real(r8) :: eccf                   ! Earth orbit eccentricity factor
    real(r8) :: coszen                 ! Cosine solar zenith angle
    real(r8) :: lat, lon               ! Latitude and longitude (radians)
    real(r8) :: prev_calday            ! Calendar day at beginning of time step (itim-1)

    real(r8) :: zref_loc(1,1,1)        ! Tower: Reference height (m)
    real(r8) :: tref_loc(1,1,1)        ! Tower: Air temperature at reference height (K)
    real(r8) :: rhref_loc(1,1,1)       ! Tower: Relative humidity at reference height (%)
    real(r8) :: qref_loc(1,1,1)        ! Tower: Specific humidity at reference height (kg/kg)
    real(r8) :: uref_loc(1,1,1)        ! Tower: Wind speed at reference height (m/s)
    real(r8) :: fsds_loc(1,1,1)        ! Tower: Solar radiation (W/m2)
    real(r8) :: flds_loc(1,1,1)        ! Tower: Longwave radiation (W/m2)
    real(r8) :: pref_loc(1,1,1)        ! Tower: Air pressure at reference height (Pa)
    real(r8) :: prect_loc(1,1,1)       ! Tower: Precipitation (mm/s)

    real(r8), parameter :: a0 = 0.17639_r8
    real(r8), parameter :: a1 = 0.00380_r8
    real(r8), parameter :: a2 = -9.0039e-06_r8
    real(r8), parameter :: a3 = 8.1351e-09_r8

    real(r8), parameter :: b0 = 0.29548_r8
    real(r8), parameter :: b1 = 0.00504_r8
    real(r8), parameter :: b2 = -1.4957e-05_r8
    real(r8), parameter :: b3 = 1.4881e-08_r8
    !---------------------------------------------------------------------

    associate ( &
    ntop       => mlcanopy_inst%ntop                     , &  ! Index for top leaf layer
    zw         => mlcanopy_inst%zw                       , &  ! Canopy heights at layer interfaces (m)
                                                              ! *** Output ***
    solar_zen  => mlcanopy_inst%solar_zen                , &  ! Solar zenith angle (radians)
    forc_hgt   => atm2lnd_inst%forc_hgt_grc              , &  ! Atmospheric reference height (m)
    forc_u     => atm2lnd_inst%forc_u_grc                , &  ! Atmospheric wind speed in east direction (m/s)
    forc_v     => atm2lnd_inst%forc_v_grc                , &  ! Atmospheric wind speed in north direction (m/s)
    forc_rh    => atm2lnd_inst%forc_rh_grc               , &  ! Atmospheric relative humidity (%)
    forc_pco2  => atm2lnd_inst%forc_pco2_grc             , &  ! Atmospheric CO2 partial pressure (Pa)
    forc_po2   => atm2lnd_inst%forc_po2_grc              , &  ! Atmospheric O2 partial pressure (Pa)
    forc_solad => atm2lnd_inst%forc_solad_grc            , &  ! Atmospheric direct beam radiation (W/m2)
    forc_solai => atm2lnd_inst%forc_solai_grc            , &  ! Atmospheric diffuse radiation (W/m2)
    forc_t     => atm2lnd_inst%forc_t_downscaled_col     , &  ! Atmospheric temperature (K)
    forc_q     => atm2lnd_inst%forc_q_downscaled_col     , &  ! Atmospheric specific humidity (kg/kg)
    forc_pbot  => atm2lnd_inst%forc_pbot_downscaled_col  , &  ! Atmospheric pressure (Pa)
    forc_lwrad => atm2lnd_inst%forc_lwrad_downscaled_col , &  ! Atmospheric longwave radiation (W/m2)
    forc_rain  => atm2lnd_inst%forc_rain_downscaled_col  , &  ! Rainfall rate (mm/s)
    forc_snow  => atm2lnd_inst%forc_snow_downscaled_col    &  ! Snowfall rate (mm/s)
    )

    ! Read tower meteorology. These variables have the name var_loc

    call readTowerMet (ncfilename, strt, zref_loc, tref_loc, rhref_loc, qref_loc, uref_loc, &
    fsds_loc, flds_loc, pref_loc, prect_loc)

    ! Assign these variables to the appropriate CLM variables. All CLM patches use
    ! this meteorology

    do p = begp, endp
       c = patch%column(p)
       g = patch%gridcell(p)

       ! These CLM variables are defined at the grid cell

       forc_hgt(g) = zref_loc(1,1,1)
       forc_u(g) = uref_loc(1,1,1)
       forc_v(g) = 0._r8
       forc_rh(g) = rhref_loc(1,1,1)

       ! Incoming solar direct beam and diffuse radiation. The total solar
       ! radiation is split equally to the visible and near-infrared wavebands
       ! and is then split into direct beam and diffuse.

       fsds = max(fsds_loc(1,1,1), 0._r8)

       fsds_vis = 0.5_r8 * fsds
       rvis = a0 + fsds_vis*(a1 + fsds_vis*(a2 + fsds_vis*a3))
       rvis = max(0.01_r8, min(0.99_r8, rvis))

       fsds_nir = 0.5_r8 * fsds
       rnir = b0 + fsds_nir*(b1 + fsds_nir*(b2 + fsds_nir*b3))
       rnir = max(0.01_r8, min(0.99_r8, rnir))

       forc_solad(g,ivis) = fsds_vis * rvis
       forc_solai(g,ivis) = fsds_vis * (1._r8 - rvis)
       forc_solad(g,inir) = fsds_nir * rnir
       forc_solai(g,inir) = fsds_nir * (1._r8 - rnir)

       ! These CLM variables are defined at the column

       forc_t(c) = tref_loc(1,1,1)
       forc_q(c) = qref_loc(1,1,1)
       forc_pbot(c) = pref_loc(1,1,1)
       forc_lwrad(c) = flds_loc(1,1,1)
       forc_rain(c) = prect_loc(1,1,1)
       forc_snow(c) = 0._r8

    end do

    ! Adjustments and other forcing variables

    do p = begp, endp
       c = patch%column(p)
       g = patch%gridcell(p)

       ! Overwrite forcing height with tower data

       if (use_tower) forc_hgt(g) = tower_ht(it)

       ! Set forcing height to 30 m if tower forcing data has no height

       if (nint(forc_hgt(g)) == -999) forc_hgt(g) = 30._r8

       ! Set atmospheric pressure to surface value if missing

       if (nint(forc_pbot(c)) == -999) forc_pbot(c) = 101325._r8

       ! Reference height check

       if (forc_hgt(g) <= zw(p,ntop(p))) then
          write (iulog,*) 'TowerMet error: ',forc_hgt(g),zw(p,ntop(p))
          call endrun()
       end if

       ! CO2: umol/mol -> Pa and O2: mol/mol -> Pa

       forc_pco2(g) = 367._r8 / 1.e06_r8 * forc_pbot(c)
       forc_po2(g) = 0.209_r8 * forc_pbot(c)

       ! Relative humidity -> specific humidity
       ! Specific humidity -> relative humidity

       call SatVap (forc_t(c), esat, desat)
       if (nint(forc_rh(g)) /= -999) then
          eair = (forc_rh(g) / 100._r8) * esat
          forc_q(c) = mmh2o / mmdry * eair / (forc_pbot(c) - (1._r8 - mmh2o/mmdry) * eair)
       else if (nint(forc_q(c)) /= -999) then
          eair = forc_q(c) * forc_pbot(c) / (mmh2o / mmdry + (1._r8 - mmh2o / mmdry) * forc_q(c))
          forc_rh(g) = eair / esat * 100._r8
       else
          write (iulog,*) 'TowerMet error: rhref and qref not valid'
          call endrun()
       end if

       ! Calculate atmospheric longwave radiation if it is missing

       if (nint(forc_lwrad(c)) == -999) then
          emiss = 0.7_r8 + 5.95e-05_r8 * 0.01_r8 * eair * exp(1500._r8/forc_t(c))
          forc_lwrad(c) = emiss * sb * forc_t(c)**4
       end if

       ! Solar zenith angle. Need to subtract one time step because zenith angle is
       ! calculated for the beginning of the time step. So use calendar day at
       ! beginning of the time step (itim - 1)

       prev_calday = get_prev_calday()
       call shr_orb_decl (prev_calday, eccen, mvelpp, lambm0, obliqr, declin, eccf)
       lat = tower_lat(it) * pi / 180._r8
       lon = tower_lon(it) * pi / 180._r8
       coszen = shr_orb_cosz (prev_calday, lat, lon, declin)
       solar_zen(p) = acos(max(0.01_r8,coszen))

    end do

    end associate
  end subroutine TowerMet

  !-----------------------------------------------------------------------
  subroutine readTowerMet (ncfilename, strt, zbot, tbot, rhbot, qbot, ubot, &
  fsdsbot, fldsbot, pbot, prectbot)
    !
    ! !DESCRIPTION:
    ! Read variables from tower site atmospheric forcing netcdf files
    !
    ! !USES:
    !
    ! !ARGUMENTS:
    implicit none
    character(len=*), intent(in) :: ncfilename  ! netcdf filename
    integer,  intent(in)  :: strt               ! Time slice of data to retrieve
    real(r8), intent(out) :: zbot(1,1,1)        ! Reference height (m)
    real(r8), intent(out) :: tbot(1,1,1)        ! Air temperature at reference height (K)
    real(r8), intent(out) :: rhbot(1,1,1)       ! Relative humidity at reference height (%)
    real(r8), intent(out) :: qbot(1,1,1)        ! Specific humidity at reference height (kg/kg)
    real(r8), intent(out) :: ubot(1,1,1)        ! Wind speed at reference height (m/s)
    real(r8), intent(out) :: fsdsbot(1,1,1)     ! Solar radiation (W/m2)
    real(r8), intent(out) :: fldsbot(1,1,1)     ! Longwave radiation (W/m2)
    real(r8), intent(out) :: pbot(1,1,1)        ! Air pressure at reference height (Pa)
    real(r8), intent(out) :: prectbot(1,1,1)    ! Precipitation (mm/s)
    !
    ! !LOCAL VARIABLES:
    integer :: ncid                             ! netcdf file ID
    integer :: status                           ! function return status
    integer :: varid                            ! netcdf variable id
    integer :: start3(3), count3(3)             ! start and count arrays for reading 3-D data from netcdf files
    !---------------------------------------------------------------------

    status = nf_open(ncfilename, nf_nowrite, ncid)
    if (status /= nf_noerr) call handle_err(status, ncfilename)

    ! Dimensions in FORTRAN are in column major order: the first array index varies the most rapidly.
    ! In NetCDF file the dimensions appear in the opposite order: lat, lon (2-D); time, lat, lon (3-D);
    ! time, levgrnd, lat, lon (4-D)

    start3 = (/ 1, 1, strt /)
    count3 = (/ 1, 1, 1 /)

    ! flds(nlon, nlat, ntime): atmospheric longwave radiation (W/m2)
    status = nf_inq_varid(ncid, "FLDS", varid)
    if (status == nf_noerr) then
       status = nf_get_vara_double(ncid, varid, start3, count3, fldsbot)
       if (status /= nf_noerr) call handle_err(status, "fldsbot")
    else
       fldsbot(1,1,1) = -999._r8
    end if

    ! fsds(nlon, nlat, ntime): atmospheric incident solar radiation (W/m2)
    status = nf_inq_varid(ncid, "FSDS", varid)
    if (status /= nf_noerr) call handle_err(status, "FSDS")

    status = nf_get_vara_double(ncid, varid, start3, count3, fsdsbot)
    if (status /= nf_noerr) call handle_err(status, "fsdsbot")

    ! psrf(nlon, nlat, ntime): atmospheric pressure (Pa)
    status = nf_inq_varid(ncid, "PSRF", varid)
    if (status == nf_noerr) then
       status = nf_get_vara_double(ncid, varid, start3, count3, pbot)
       if (status /= nf_noerr) call handle_err(status, "pbot")
    else
       pbot(1,1,1) = -999._r8
    end if

    ! rh(nlon, nlat, ntime): atmospheric relative humidity (%)
    status = nf_inq_varid(ncid, "RH", varid)
    if (status == nf_noerr) then
       status = nf_get_vara_double(ncid, varid, start3, count3, rhbot)
       if (status /= nf_noerr) call handle_err(status, "rhbot")
    else
       rhbot(1,1,1) = -999._r8
    end if

    ! qbot(nlon, nlat, ntime): atmospheric specific humidity (kg/kg)
    status = nf_inq_varid(ncid, "QBOT", varid)
    if (status == nf_noerr) then
       status = nf_get_vara_double(ncid, varid, start3, count3, qbot)
       if (status /= nf_noerr) call handle_err(status, "qbot")
    else
       qbot(1,1,1) = -999._r8
    end if

    ! prectmms(nlon, nlat, ntime): precipitation (mm/sec)
    status = nf_inq_varid(ncid, "PRECTmms", varid)
    if (status /= nf_noerr) call handle_err(status, "PRECTmms")

    status = nf_get_vara_double(ncid, varid, start3, count3, prectbot)
    if (status /= nf_noerr) call handle_err(status, "prectbot")

    ! tbot(nlon, nlat, ntime): atmospheric air temperature (K)
    status = nf_inq_varid(ncid, "TBOT", varid)
    if (status /= nf_noerr) call handle_err(status, "TBOT")

    status = nf_get_vara_double(ncid, varid, start3, count3, tbot)
    if (status /= nf_noerr) call handle_err(status, "tbot")

    ! wind(nlon, nlat, ntime): atmospheric wind velocity magnitude (m/s)
    status = nf_inq_varid(ncid, "WIND", varid)
    if (status /= nf_noerr) call handle_err(status, "WIND")

    status = nf_get_vara_double(ncid, varid, start3, count3, ubot)
    if (status /= nf_noerr) call handle_err(status, "ubot")

    ! zbot(nlon, nlat, ntime): observational height (m)
    status = nf_inq_varid(ncid, "ZBOT", varid)
    if (status == nf_noerr) then
       status = nf_get_vara_double(ncid, varid, start3, count3, zbot)
       if (status /= nf_noerr) call handle_err(status, "zbot")
    else
       zbot(1,1,1) = -999._r8
    end if

    status = nf_close(ncid)

  end subroutine readTowerMet

end module TowerMetMod
