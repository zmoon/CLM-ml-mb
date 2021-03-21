module SoilTemperatureMod

  !-----------------------------------------------------------------------
  ! !DESCRIPTION:
  ! Calculates snow and soil temperatures including phase change
  !
  use shr_kind_mod, only : r8 => shr_kind_r8
  use abortutils, only : endrun
  use decompMod, only : bounds_type
  use EnergyFluxType, only : energyflux_type
  use WaterFluxType, only : waterflux_type
  use WaterStateType, only : waterstate_type
  use TemperatureType, only : temperature_type
  use SolarAbsorbedType, only : solarabs_type
  use CanopyStateType, only : canopystate_type
  use SoilStateType, only : soilstate_type
  use atm2lndType, only : atm2lnd_type
  use CanopyFluxesMultilayerType, only : mlcanopy_type
  !
  ! !PUBLIC TYPES:
  implicit none
  !
  ! !PUBLIC MEMBER FUNCTIONS:
  public :: SoilTemperature      ! CLM snow and soil temperatures including phase change
  public :: SoilTemperatureInit  ! Interface to CLM soil temperature calculation
  public :: SoilTemperatureAfter ! Interface to CLM soil temperature calculation
  public :: SoilThermProp        ! Thermal conductivity and heat capacity of snow/soil layers
  !
  ! !PRIVATE MEMBER FUNCTIONS:
  private :: PhaseChangeH2osfc   ! When surface water freezes move ice to bottom snow layer
  private :: PhaseChange_beta    ! Calculation of the phase change within snow and soil layers
  !-----------------------------------------------------------------------

contains

  !-----------------------------------------------------------------------
  subroutine SoilTemperatureAfter (num_nolakep, filter_nolakep, &
  energyflux_inst, waterflux_inst, mlcanopy_inst)
    !
    ! !DESCRIPTION:
    ! Interface to CLM soil temperature calculation. Copy CLM soil fluxes to
    ! multilayer canopy variables and update total fluxes (soil + veg).
    !
    ! !USES:
    use clm_varcon, only : mmh2o
    use clm_varctl, only : iulog
    use WaterVaporMod, only : LatVap
    use PatchType, only : patch
    !
    ! !ARGUMENTS:
    implicit none
    integer, intent(in) :: num_nolakep            ! Number of non-lake points in CLM patch filter
    integer, intent(in) :: filter_nolakep(:)      ! CLM patch filter for non-lake points
    type(energyflux_type) :: energyflux_inst
    type(waterflux_type) :: waterflux_inst
    type(mlcanopy_type) :: mlcanopy_inst
    !
    ! !LOCAL VARIABLES:
    integer  :: f                                 ! Filter index
    integer  :: p                                 ! Patch index for CLM g/l/c/p hierarchy
    integer  :: c                                 ! Column index for CLM g/l/c/p hierarchy
    real(r8) :: err                               ! Energy imbalance
    !-----------------------------------------------------------------------

    associate ( &
                                                                ! *** Input ***
    pcolumn        => patch%column                         , &  ! Column index of patch for CLM g/l/c/p hierarchy
    eflx_sh_grnd   => energyflux_inst%eflx_sh_grnd_patch   , &  ! CLM: Sensible heat flux from ground (W/m2)
    eflx_soil_grnd => energyflux_inst%eflx_soil_grnd_patch , &  ! CLM: Soil heat flux (W/m2)
    qflx_evap_soi  => waterflux_inst%qflx_evap_soi_patch   , &  ! CLM: Soil evaporation (kg H2O/m2/s)
    tref           => mlcanopy_inst%tref                   , &  ! Air temperature at reference height (K)
    shveg          => mlcanopy_inst%shveg                  , &  ! Sensible heat flux, vegetation (W/m2)
    lhveg          => mlcanopy_inst%lhveg                  , &  ! Latent heat flux, vegetation (W/m2)
    etveg          => mlcanopy_inst%etveg                  , &  ! Water vapor flux, vegetation (mol H2O/m2/s)
    rnet           => mlcanopy_inst%rnet                   , &  ! Net radiation (W/m2)
    stflx          => mlcanopy_inst%stflx                  , &  ! Canopy storage heat flux (W/m2)
    rnsoi          => mlcanopy_inst%rnsoi                  , &  ! Net radiation, ground (W/m2)
                                                                ! *** Output ***
    shsoi          => mlcanopy_inst%shsoi                  , &  ! Sensible heat flux, ground (W/m2)
    lhsoi          => mlcanopy_inst%lhsoi                  , &  ! Latent heat flux, ground (W/m2)
    etsoi          => mlcanopy_inst%etsoi                  , &  ! Water vapor flux, ground (mol H2O/m2/s)
    gsoi           => mlcanopy_inst%gsoi                   , &  ! Soil heat flux (W/m2)
    shflx          => mlcanopy_inst%shflx                  , &  ! Sensible heat flux (W/m2)
    lhflx          => mlcanopy_inst%lhflx                  , &  ! Latent heat flux (W/m2)
    etflx          => mlcanopy_inst%etflx                    &  ! Water vapor flux (mol H2O/m2/s)
    )

    do f = 1, num_nolakep
       p = filter_nolakep(f)
       c = pcolumn(p)

       ! Copy CLM fluxes to multilayer canopy variables

       shsoi(p) = eflx_sh_grnd(p)
       etsoi(p) = qflx_evap_soi(p) / mmh2o
       lhsoi(p) = etsoi(p) * LatVap(tref(p))
       gsoi(p) = eflx_soil_grnd(p)

       ! Total fluxes, including vegetation

       shflx(p) = shveg(p) + shsoi(p)
       lhflx(p) = lhveg(p) + lhsoi(p)
       etflx(p) = etveg(p) + etsoi(p)

       ! Error check

       err = rnet(p) - shflx(p) - lhflx(p) - gsoi(p) - stflx(p)
       if (abs(err) > 0.01_r8) then
          call endrun (msg=' ERROR: SoilTemperatureAfter: canopy fluxes failed to conserve energy')
!         write (iulog,'(i10,5f10.3)') p, rnet(p), shflx(p), lhflx(p), gsoi(p), stflx(p), err
       end if

       err = rnsoi(p) - shsoi(p) - lhsoi(p) - gsoi(p)
       if (abs(err) > 0.01_r8) then
          call endrun (msg=' ERROR: SoilTemperatureAfter: soil fluxes failed to conserve energy')
!         write (iulog,'(i10,5f10.3)') p, rnsoi(p), shsoi(p), lhsoi(p), gsoi(p), err
       end if

    end do

    end associate
  end subroutine SoilTemperatureAfter

  !-----------------------------------------------------------------------
  subroutine SoilTemperatureInit (num_nolakec, filter_nolakec, soilstate_inst, &
  waterstate_inst, energyflux_inst, temperature_inst, solarabs_inst, &
  atm2lnd_inst, mlcanopy_inst)
    !
    ! !DESCRIPTION:
    ! Interface to CLM soil temperature calculation
    !
    ! !USES:
    use clm_varcon, only : mmh2o
    use clm_varpar, only : ivis, inir, nlevgrnd, nlevsoi
    use WaterVaporMod, only : LatVap
    use ColumnType, only : col
    !
    ! !ARGUMENTS:
    implicit none
    integer, intent(in) :: num_nolakec                 ! Number of non-lake points in CLM column filter
    integer, intent(in) :: filter_nolakec(:)           ! CLM column filter for non-lake points
    type(energyflux_type) :: energyflux_inst
    type(waterstate_type) :: waterstate_inst
    type(temperature_type) :: temperature_inst
    type(solarabs_type) :: solarabs_inst
    type(soilstate_type) :: soilstate_inst
    type(atm2lnd_type) :: atm2lnd_inst
    type(mlcanopy_type) :: mlcanopy_inst
    !
    ! !LOCAL VARIABLES:
    integer  :: fc, c, p, j                            ! indices

    real(r8) :: bulk_dens_min                          ! Bulk density, mineral soil (kg/m3)
    real(r8) :: tkdry_min                              ! Dry thermal conductivity, mineral fraction (W/m/K)
    real(r8) :: quartz                                 ! Quartz fraction of soil (fraction)
    real(r8) :: tksol_other                            ! Thermal conductivity, other minerals (W/m/K)
    real(r8) :: tksol_min                              ! Thermal conductivity, minerals (W/m/K)
    real(r8) :: tkm                                    ! Mineral conductivity (W/m/K)

    real(r8), parameter :: tksol_quartz = 7.7_r8       ! Thermal conductivity, quartz (W/m/K)
    real(r8), parameter :: tksol_org = 0.25_r8         ! Thermal conductivity, organic material (W/m/K)
    real(r8), parameter :: tkdry_org = 0.05_r8         ! Thermal conductivity, dry organic material (W/m/K)
    real(r8), parameter :: cvsol = 1.926e06_r8         ! Heat capacity of soil solids (J/m3/K)
    real(r8), parameter :: cvorg = 2.5e06_r8           ! Heat capacity of organic material (J/m3/K)
    real(r8), parameter :: csol_bedrock = 2.0e6_r8     ! Vol. heat capacity of granite/sandstone (J/m3/K)
    !-----------------------------------------------------------------------

    associate ( &
    tref           => mlcanopy_inst%tref                     , &  ! Air temperature at reference height (K)
    swsoi          => mlcanopy_inst%swsoi                    , &  ! Absorbed solar radiation, ground (W/m2)
    irsky          => mlcanopy_inst%irsky                    , &  ! Atmospheric longwave radiation (W/m2)
    ircan          => mlcanopy_inst%ircan                    , &  ! Upward longwave radiation above canopy (W/m2)
    irsoi          => mlcanopy_inst%irsoi                    , &  ! Absorbed longwave radiation, ground (W/m2)
    emg            => temperature_inst%emg_col               , &  ! Ground (soil) emissivity
    cellsand       => soilstate_inst%cellsand_col            , &  ! Soil layer percent sand
    cellorg        => soilstate_inst%cellorg_col             , &  ! Soil layer organic fraction
    snl            => col%snl                                , &  ! Number of snow layers
    frac_sno_eff   => waterstate_inst%frac_sno_eff_col       , &  ! Effective fraction of ground covered by snow (0 to 1)
    frac_sno       => waterstate_inst%frac_sno_col           , &  ! Fraction of ground covered by snow (0 to 1)
    snow_depth     => waterstate_inst%snow_depth_col         , &  ! Snow height (m)
    h2osno         => waterstate_inst%h2osno_col             , &  ! Total snow water (kg/m2)
    sabg_snow      => solarabs_inst%sabg_snow_patch          , &  ! Solar radiation absorbed by snow (W/m2)
    frac_h2osfc    => waterstate_inst%frac_h2osfc_col        , &  ! Fraction of ground covered by surface water (0 to 1)
    h2osfc         => waterstate_inst%h2osfc_col             , &  ! Surface water (kg/m2)
    t_h2osfc       => temperature_inst%t_h2osfc_col          , &  ! Surface water temperature (K)
    t_h2osfc_bef   => temperature_inst%t_h2osfc_bef_col      , &  ! Saved surface water temperature (K)
    ulrad          => energyflux_inst%ulrad_patch            , &  ! Upward longwave radiation above the canopy (W/m2)
    sabg           => solarabs_inst%sabg_patch               , &  ! Solar radiation absorbed by ground (W/m2)
    sabg_soil      => solarabs_inst%sabg_soil_patch          , &  ! Solar radiation absorbed by soil (W/m2)
    sabg_lyr       => solarabs_inst%sabg_lyr_patch           , &  ! Absorbed radiation in each snow layer and top soil layer (W/m2)
    eflx_bot       => energyflux_inst%eflx_bot_col           , &  ! Heat flux from beneath column (W/m2) [+ = upward]
    htvp           => energyflux_inst%htvp_col               , &  ! Latent heat of vapor of water (or sublimation) (J/kg)
    dlrad          => energyflux_inst%dlrad_patch            , &  ! Downward longwave radiation below the canopy (W/m2)
    t_grnd         => temperature_inst%t_grnd_col            , &  ! Ground surface temperature (K)
    t_soisno       => temperature_inst%t_soisno_col          , &  ! Soil/snow temperature (K)
    tssbef         => temperature_inst%tssbef_col            , &  ! Soil/snow temperature before update
    watsat         => soilstate_inst%watsat_col              , &  ! Soil layer volumetric water content at saturation (porosity)
    tkmg           => soilstate_inst%tkmg_col                , &  ! Soil layer thermal conductivity, soil minerals  (W/m/K)
    tkdry          => soilstate_inst%tkdry_col               , &  ! Soil layer thermal conductivity, dry soil (W/m/K)
    csol           => soilstate_inst%csol_col                , &  ! Soil layer heat capacity, soil solids (J/m3/K)
    forc_lwrad     => atm2lnd_inst%forc_lwrad_downscaled_col   &  ! Downward infrared (longwave) radiation (W/m2)
    )

    ! Set CLM variables need for soil temperature: NO SNOW, NO SURFACE WATER

    do fc = 1, num_nolakec
       c = filter_nolakec(fc)
       p = c                       ! Multi-layer code has one patch one column

       ! Ground (soil) emissivity

       emg(c) = 0.96_r8

       ! Snow
       snl(c) = 0                  ! Number of snow layers
       frac_sno_eff(c) = 0._r8     ! Effective fraction of ground covered by snow (0 to 1)
       frac_sno(c) = 0._r8         ! Fraction of ground covered by snow (0 to 1)
       snow_depth(c) = 0._r8       ! Snow height (m)
       h2osno(c) = 0._r8           ! Total snow water (kg/m2)

       ! the code sets this equal to sabg even if there is no snow

       sabg_snow(p) = swsoi(p,ivis) + swsoi(p,inir)       ! Solar radiation absorbed by snow (W/m2)

       ! Surface water

       frac_h2osfc(c) = 0._r8      ! Fraction of ground covered by surface water (0 to 1)
       h2osfc(c) = 0._r8           ! Surface water (kg/m2)
       t_h2osfc(c) = t_grnd(c)     ! Surface water temperature (K)
       t_h2osfc_bef(c) = t_grnd(c) ! Saved surface water temperature (K)

       ! Radiation fluxes - change in longwave fluxes is excluded from soil temperature
       ! calculation, instead use net radiation as input (swsoi + irsoi)

!      forc_lwrad(c) = irsky(p)                  ! Downward longwave radiation (W/m2)
       ulrad(p) = ircan(p)                       ! Upward longwave radiation above the canopy (W/m2)
       dlrad(p) = irsoi(p)                       ! Longwave radiation absorbed by ground (W/m2)

       sabg(p) = swsoi(p,ivis) + swsoi(p,inir)   ! Solar radiation absorbed by ground (W/m2)
       sabg_soil(p) = sabg(p)                    ! Solar radiation absorbed by soil (W/m2)
       sabg_lyr(p,1) = sabg(p)                   ! Absorbed solar radiation for patch x lyr (W/m2)

       ! Miscellaneous

       eflx_bot(c) = 0._r8                       ! Heat flux from beneath column (W/m2) [+ = upward]
       htvp(c) = LatVap(tref(p)) / mmh2o         ! Latent heat of vapor of water (or sublimation) (J/kg)

       ! Soil thermal variables

       do j = 1, nlevgrnd

          tssbef(c,j) = t_soisno(c,j)            ! Soil/snow temperature before update

          ! Thermal conductivity, dry soil (W/m/K)

          if (j > nlevsoi) then
             bulk_dens_min = 2700._r8
             tkdry_min = (0.135_r8*bulk_dens_min + 64.7_r8) / (2700._r8 - 0.947_r8*bulk_dens_min)
             tkdry(c,j) = tkdry_min
          else
             bulk_dens_min = 2700._r8 * (1._r8 - watsat(c,j))
             tkdry_min = (0.135_r8*bulk_dens_min + 64.7_r8) / (2700._r8 - 0.947_r8*bulk_dens_min)
             tkdry(c,j) = (1._r8 - cellorg(c,j)) * tkdry_min + cellorg(c,j) * tkdry_org
          end if

          ! Soil solids thermal conductivity (W/m/K)

          if (j > nlevsoi) then
             tkmg(c,j) = 3._r8
          else
             quartz = cellsand(c,j) / 100._r8
             if (quartz > 0.2_r8) then
                tksol_other = 2._r8
             else
                tksol_other = 3._r8
             end if
             tksol_min = tksol_quartz**quartz * tksol_other**(1._r8 - quartz)
             tkm = (1._r8 - cellorg(c,j)) * tksol_min + cellorg(c,j) * tksol_org
             tkmg(c,j) = tkm ** (1._r8 - watsat(c,j))
          end if

          ! Heat capacity, soil solids (J/m3/K)

          if (j > nlevsoi) then
             csol(c,j) = csol_bedrock
          else
             csol(c,j) = cvsol * (1._r8 - cellorg(c,j)) + cvorg * cellorg(c,j)
          end if

       end do

    end do

    end associate
  end subroutine SoilTemperatureInit

  !-----------------------------------------------------------------------
  subroutine SoilTemperature (bounds, num_nolakec, filter_nolakec, num_nolakep, filter_nolakep, &
  xmf, fact, c_h2osfc, xmf_h2osfc, energyflux_inst, waterflux_inst, waterstate_inst, &
  temperature_inst, soilstate_inst, atm2lnd_inst, solarabs_inst, canopystate_inst)
    !
    ! !DESCRIPTION:
    ! Snow and soil temperatures including phase change
    ! o The volumetric heat capacity is calculated as a linear combination
    !   in terms of the volumetric fraction of the constituent phases.
    ! o The thermal conductivity of soil is computed from
    !   the algorithm of Johansen (as reported by Farouki 1981), and the
    !   conductivity of snow is from the formulation used in
    !   SNTHERM (Jordan 1991).
    ! o Boundary conditions:
    !   F = Rnet - Hg - LEg (top),  F= 0 (base of the soil column).
    ! o Soil / snow temperature is predicted from heat conduction
    !   in 10 soil layers and up to 5 snow layers.
    !   The thermal conductivities at the interfaces between two
    !   neighboring layers (j, j+1) are derived from an assumption that
    !   the flux across the interface is equal to that from the node j
    !   to the interface and the flux from the interface to the node j+1.
    !   The equation is solved using the Crank-Nicholson method and
    !   results in a tridiagonal system equation.
    !
    ! !USES:
    use clm_time_manager, only : get_step_size
    use clm_varctl, only : iulog
    use clm_varcon, only : sb, capr, cnfac, hvap, hsub, denh2o, cpliq
    use clm_varpar, only : nlevsno, nlevgrnd, nlevsoi,max_patch_per_col
    use BandDiagonalMod, only : BandDiagonal
    use PatchType, only : patch
    use ColumnType, only : col
    !
    ! !ARGUMENTS:
    implicit none
    type(bounds_type), intent(in) :: bounds                     ! bounds
    integer , intent(in)  :: num_nolakep                        ! Number of non-lake points in CLM pft filter
    integer , intent(in)  :: filter_nolakep(:)                  ! CLM pft filter for non-lake points
    integer , intent(in)  :: num_nolakec                        ! number of column non-lake points in column filter
    integer , intent(in)  :: filter_nolakec(:)                  ! column filter for non-lake points
    real(r8), intent(out) :: xmf( bounds%begc: )                ! total latent heat of phase change of ground water [col]
    real(r8), intent(out) :: fact( bounds%begc: , -nlevsno+1: ) ! used in computing tridiagonal matrix [col, lev]
    real(r8), intent(out) :: xmf_h2osfc( bounds%begc: )         ! latent heat of phase change of surface water [col]
    real(r8), intent(out) :: c_h2osfc( bounds%begc: )           ! heat capacity of surface water [col]
    type(energyflux_type) :: energyflux_inst
    type(waterflux_type) :: waterflux_inst
    type(waterstate_type) :: waterstate_inst
    type(temperature_type) :: temperature_inst
    type(soilstate_type) :: soilstate_inst
    type(atm2lnd_type) :: atm2lnd_inst
    type(solarabs_type) :: solarabs_inst
    type(canopystate_type) :: canopystate_inst
    !
    ! !LOCAL VARIABLES:
    real(r8) :: scalez = 0.025_r8                                  ! Soil layer thickness discretization (m)
    integer  :: j,c,p,l,g,pi,jj,kk                                 ! indices
    integer  :: fc                                                 ! lake filtered column indices
    integer  :: fp                                                 ! lake filtered pft indices
    integer  :: jtop(bounds%begc:bounds%endc)                      ! top level at each column
    real(r8) :: dtime                                              ! land model time step (sec)
    real(r8) :: at (bounds%begc:bounds%endc,-nlevsno+1:nlevgrnd)   ! "a" vector for tridiagonal matrix
    real(r8) :: bt (bounds%begc:bounds%endc,-nlevsno+1:nlevgrnd)   ! "b" vector for tridiagonal matrix
    real(r8) :: ct (bounds%begc:bounds%endc,-nlevsno+1:nlevgrnd)   ! "c" vector for tridiagonal matrix
    real(r8) :: rt (bounds%begc:bounds%endc,-nlevsno+1:nlevgrnd)   ! "r" vector for tridiagonal solution
    real(r8) :: cv (bounds%begc:bounds%endc,-nlevsno+1:nlevgrnd)   ! heat capacity [J/(m2 K)]
    real(r8) :: tk (bounds%begc:bounds%endc,-nlevsno+1:nlevgrnd)   ! thermal conductivity [W/(m K)]
    real(r8) :: fn (bounds%begc:bounds%endc,-nlevsno+1:nlevgrnd)   ! heat diffusion through the layer interface [W/m2]
    real(r8) :: fn1(bounds%begc:bounds%endc,-nlevsno+1:nlevgrnd)   ! heat diffusion through the layer interface [W/m2]
    real(r8) :: dzm                                                ! used in computing tridiagonal matrix
    real(r8) :: dzp                                                ! used in computing tridiagonal matrix
    real(r8) :: hs(bounds%begc:bounds%endc)                        ! net energy flux into the surface (w/m2)
    real(r8) :: lwrad_emit(bounds%begc:bounds%endc)                ! emitted longwave radiation
    real(r8) :: dlwrad_emit(bounds%begc:bounds%endc)               ! time derivative of emitted longwave radiation
    integer  :: lyr_top                                            ! index of top layer of snowpack (-4 to 0) [idx]
    real(r8) :: sabg_lyr_col(bounds%begc:bounds%endc,-nlevsno+1:1) ! absorbed solar radiation (col,lyr) [W/m2]
    real(r8) :: eflx_gnet_top                                      ! net energy flux into surface layer, pft-level [W/m2]
    real(r8) :: hs_top(bounds%begc:bounds%endc)                    ! net energy flux into surface layer (col) [W/m2]
    real(r8) :: fn_h2osfc(bounds%begc:bounds%endc)
    real(r8) :: fn1_h2osfc(bounds%begc:bounds%endc)
    real(r8) :: dz_h2osfc(bounds%begc:bounds%endc)

    integer, parameter :: nband=5
    real(r8) :: bmatrix(bounds%begc:bounds%endc,nband,-nlevsno:nlevgrnd)
    real(r8) :: tvector(bounds%begc:bounds%endc,-nlevsno:nlevgrnd)
    real(r8) :: rvector(bounds%begc:bounds%endc,-nlevsno:nlevgrnd)
    real(r8) :: tk_h2osfc(bounds%begc:bounds%endc)
    real(r8) :: dhsdT(bounds%begc:bounds%endc)
    real(r8) :: hs_snow(bounds%begc:bounds%endc)
    real(r8) :: hs_soil(bounds%begc:bounds%endc)
    real(r8) :: hs_top_snow(bounds%begc:bounds%endc)
    real(r8) :: hs_top_soil(bounds%begc:bounds%endc)
    real(r8) :: hs_h2osfc(bounds%begc:bounds%endc)
    real(r8) :: lwrad_emit_snow(bounds%begc:bounds%endc)
    real(r8) :: lwrad_emit_soil(bounds%begc:bounds%endc)
    real(r8) :: lwrad_emit_h2osfc(bounds%begc:bounds%endc)
    real(r8) :: eflx_gnet_snow
    real(r8) :: eflx_gnet_soil
    real(r8) :: eflx_gnet_h2osfc
    integer  :: jbot(bounds%begc:bounds%endc)                      ! bottom level at each column

    real(r8) :: t_grnd0(bounds%begc:bounds%endc)                   ! t_grnd of previous time step
    real(r8) :: tinc(bounds%begc:bounds%endc)                      ! temperature difference of two time step
    real(r8) :: egsmax(bounds%begc:bounds%endc)                    ! max. evaporation which soil can provide at one time step
    real(r8) :: topsoil_evap_tot(bounds%begc:bounds%endc)          ! column-level total evaporation from top soil layer
    real(r8) :: sumwt(bounds%begc:bounds%endc)                     ! temporary
    real(r8) :: egirat(bounds%begc:bounds%endc)                    ! ratio of topsoil_evap_tot : egsmax
    real(r8) :: save_qflx_evap_soi                                 ! temporary storage for qflx_evap_soi
    real(r8) :: lw_grnd                                            ! temporary longwave radiation from ground
    !-----------------------------------------------------------------------

   associate( &
   forc_lwrad     => atm2lnd_inst%forc_lwrad_downscaled_col , & ! Input: [real(r8) (:)]  downward infrared (longwave) radiation (W/m**2)
   frac_sno_eff   => waterstate_inst%frac_sno_eff_col       , & ! Input:  [real(r8) (:)]  eff. fraction of ground covered by snow (0 to 1)
   frac_sno       => waterstate_inst%frac_sno_col           , & ! Input:  [real(r8) (:)]  fraction of ground covered by snow (0 to 1)
   frac_h2osfc    => waterstate_inst%frac_h2osfc_col        , & ! Input:  [real(r8) (:)]  fraction of ground covered by surface water (0 to 1)
   h2osfc         => waterstate_inst%h2osfc_col             , & ! Input:  [real(r8) (:)]  surface water (mm)
   t_h2osfc       => temperature_inst%t_h2osfc_col          , & ! Input:  [real(r8) (:)]  surface water temperature
   t_h2osfc_bef   => temperature_inst%t_h2osfc_bef_col      , & ! Input:  [real(r8) (:)]  saved surface water temperature
   snow_depth     => waterstate_inst%snow_depth_col         , & ! Input:  [real(r8) (:)]  snow height (m)
   eflx_sh_snow   => energyflux_inst%eflx_sh_snow_patch     , & ! Input:  [real(r8) (:)]  sensible heat flux from snow (W/m**2)
   eflx_sh_soil   => energyflux_inst%eflx_sh_soil_patch     , & ! Input:  [real(r8) (:)]  sensible heat flux from soil (W/m**2)
   eflx_sh_h2osfc => energyflux_inst%eflx_sh_h2osfc_patch   , & ! Input:  [real(r8) (:)]  sensible heat flux from surface water (W/m**2)
   qflx_ev_snow   => waterflux_inst%qflx_ev_snow_patch      , & ! Input:  [real(r8) (:)]  evaporation flux from snow (W/m**2)
   qflx_ev_soil   => waterflux_inst%qflx_ev_soil_patch      , & ! Input:  [real(r8) (:)]  evaporation flux from soil (kg/m**2/s)
   qflx_ev_h2osfc => waterflux_inst%qflx_ev_h2osfc_patch    , & ! Input:  [real(r8) (:)]  evaporation flux from h2osfc (W/m**2)
   sabg_soil      => solarabs_inst%sabg_soil_patch          , & ! Input:  [real(r8) (:)]  solar radiation absorbed by soil (W/m**2)
   sabg_snow      => solarabs_inst%sabg_snow_patch          , & ! Input:  [real(r8) (:)]  solar radiation absorbed by snow (W/m**2)
   sabg_chk       => solarabs_inst%sabg_chk_patch           , & ! Input:  [real(r8) (:)]  sum of soil/snow using current fsno, for balance check
   sabg           => solarabs_inst%sabg_patch               , & ! Input:  [real(r8) (:)]  solar radiation absorbed by ground (W/m**2)
   sabg_lyr       => solarabs_inst%sabg_lyr_patch           , & ! Input:  [real(r8) (:,:)]  absorbed solar radiation (pft,lyr) [W/m2]
   clandunit      => col%landunit                           , & ! Input:  [integer (:)]  column's landunit
   npfts          => col%npatches                           , & ! Input:  [integer (:)]  column's number of patches
   pfti           => col%patchi                             , & ! Input:  [integer (:)]  column's beginning patch index
   snl            => col%snl                                , & ! Input:  [integer (:)]  number of snow layers
   htvp           => energyflux_inst%htvp_col               , & ! Input:  [real(r8) (:)]  latent heat of vapor of water (or sublimation) [j/kg]
   emg            => temperature_inst%emg_col               , & ! Input:  [real(r8) (:)]  ground emissivity
   t_grnd         => temperature_inst%t_grnd_col            , & ! Input:  [real(r8) (:)]  ground surface temperature [K]
   hc_soi         => temperature_inst%hc_soi_col            , & ! Input:  [real(r8) (:)]  soil heat content (MJ/m2)
   hc_soisno      => temperature_inst%hc_soisno_col         , & ! Input:  [real(r8) (:)]  soil plus snow plus lake heat content (MJ/m2)
   eflx_fgr12     => energyflux_inst%eflx_fgr12_col         , & ! Input:  [real(r8) (:)]  heat flux between soil layer 1 and 2 (W/m2)
   eflx_fgr       => energyflux_inst%eflx_fgr_col           , & ! Input:  [real(r8) (:,:)]  (rural) soil downward heat flux (W/m2) (1:nlevgrnd)
   zi             => col%zi                                 , & ! Input:  [real(r8) (:,:)]  interface level below a "z" level (m)
   dz             => col%dz                                 , & ! Input:  [real(r8) (:,:)]  layer depth (m)
   z              => col%z                                  , & ! Input:  [real(r8) (:,:)]  layer thickness (m)
   t_soisno       => temperature_inst%t_soisno_col          , & ! Input:  [real(r8) (:,:)]  soil temperature (Kelvin)
   eflx_bot       => energyflux_inst%eflx_bot_col           , & ! Input:  [real(r8) (:)]  heat flux from beneath column (W/m**2) [+ = upward]
   pactive        => patch%active                           , & ! Input:  [logical (:)]  true=>do computations on this pft
   pgridcell      => patch%gridcell                         , & ! Input:  [integer (:)]  pft's gridcell index
   plandunit      => patch%landunit                         , & ! Input:  [integer (:)]  pft's landunit index
   pcolumn        => patch%column                           , & ! Input:  [integer (:)]  pft's column index
   pwtcol         => patch%wtcol                            , & ! Input:  [real(r8) (:)]  weight of pft relative to column
   frac_veg_nosno => canopystate_inst%frac_veg_nosno_patch  , & ! Input:  [integer (:)]  fraction of vegetation not covered by snow (0 OR 1)
   cgrnd          => energyflux_inst%cgrnd_patch            , & ! Input:  [real(r8) (:)]  deriv. of soil energy flux wrt to soil temp [w/m2/k]
   cgrnds         => energyflux_inst%cgrnds_patch           , & ! Input:  [real(r8) (:)]  deriv, of soil sensible heat flux wrt soil temp [w/m2/k]
   cgrndl         => energyflux_inst%cgrndl_patch           , & ! Input:  [real(r8) (:)]  deriv of soil latent heat flux wrt soil temp [kg/m**2/k]
   dlrad          => energyflux_inst%dlrad_patch            , & ! Input:  [real(r8) (:)]  downward longwave radiation blow the canopy [W/m2]
   eflx_sh_grnd   => energyflux_inst%eflx_sh_grnd_patch     , & ! Input:  [real(r8) (:)]  sensible heat flux from ground (W/m**2)
   qflx_evap_soi  => waterflux_inst%qflx_evap_soi_patch     , & ! Input:  [real(r8) (:)]  soil evaporation (mm H2O/s)
   h2osno         => waterstate_inst%h2osno_col             , & ! Input:  [real(r8) (:)]  total snow water (col) [kg/m2]
   tssbef         => temperature_inst%tssbef_col            , & ! Input:  [real(r8) (:,:)]  soil/snow temperature before update
   h2osoi_liq     => waterstate_inst%h2osoi_liq_col         , & ! Input:  [real(r8) (:,:)]  liquid water (col,lyr) [kg/m2]
   h2osoi_ice     => waterstate_inst%h2osoi_ice_col         , & ! Input:  [real(r8) (:,:)]  ice content (col,lyr) [kg/m2]
   watsat         => soilstate_inst%watsat_col              , & ! Input:  [real(r8) (:,:)]  volumetric soil water at saturation (porosity)
   tkmg           => soilstate_inst%tkmg_col                , & ! Input:  [real(r8) (:,:)]  thermal conductivity, soil minerals  [W/m-K]
   tkdry          => soilstate_inst%tkdry_col               , & ! Input:  [real(r8) (:,:)]  thermal conductivity, dry soil (W/m/Kelvin)
   csol           => soilstate_inst%csol_col                , & ! Input:  [real(r8) (:,:)]  heat capacity, soil solids (J/m**3/Kelvin)
   bsw            => soilstate_inst%bsw_col                 , & ! Input:  [real(r8) (:,:)] Clapp and Hornberger "b"
   sucsat         => soilstate_inst%sucsat_col              , & ! Input:  [real(r8) (:,:)] minimum soil suction (mm)
   ulrad          => energyflux_inst%ulrad_patch            , & ! Input:  [real(r8) (:)]  upward longwave radiation above the canopy [W/m2]
   eflx_gnet      => energyflux_inst%eflx_gnet_patch        , & ! Output: [real(r8) (:)]  net ground heat flux into the surface (W/m**2)
   dgnetdT        => energyflux_inst%dgnetdT_patch          , & ! Output: [real(r8) (:)]  temperature derivative of ground net heat flux
   eflx_soil_grnd => energyflux_inst%eflx_soil_grnd_patch   , & ! Output: [real(r8) (:)]  soil heat flux (W/m**2) [+ = into soil]
   eflx_lwrad_out => energyflux_inst%eflx_lwrad_out_patch   , & ! Output: [real(r8) (:)]  emitted infrared (longwave) radiation (W/m**2)
   eflx_lwrad_net => energyflux_inst%eflx_lwrad_net_patch   , & ! Output: [real(r8) (:)]  net infrared (longwave) rad (W/m**2) [+ = to atm]
   begc           => bounds%begc                            , &
   endc           => bounds%endc                              &
   )

    ! Get step size

    dtime = get_step_size()

    ! Compute ground surface and soil temperatures

    ! Thermal conductivity and Heat capacity

    tk_h2osfc(begc:endc) = 1.e36_r8
    call SoilThermProp(bounds, num_nolakec, filter_nolakec, &
    tk(begc:endc, :), cv(begc:endc, :), tk_h2osfc(begc:endc), &
    waterstate_inst, temperature_inst, soilstate_inst)

    ! Net ground heat flux into the surface and its temperature derivative
    ! Added a pfts loop here to get the average of hs and dhsdT over
    ! all PFTs on the column. Precalculate the terms that do not depend on PFT.

    do fc = 1,num_nolakec
       c = filter_nolakec(fc)
       lwrad_emit(c)  =    emg(c) * sb * t_grnd(c)**4
       dlwrad_emit(c) = 4._r8*emg(c) * sb * t_grnd(c)**3

       ! fractionate lwrad_emit; balanced in CanopyFluxes & Biogeophysics2
       lwrad_emit_snow(c)    =    emg(c) * sb * t_soisno(c,snl(c)+1)**4
       lwrad_emit_soil(c)    =    emg(c) * sb * t_soisno(c,1)**4
       lwrad_emit_h2osfc(c)  =    emg(c) * sb * t_h2osfc(c)**4
    end do

    hs_snow(begc:endc)   = 0._r8
    hs_soil(begc:endc)   = 0._r8
    hs_h2osfc(begc:endc) = 0._r8
    hs(begc:endc)        = 0._r8
    dhsdT(begc:endc)     = 0._r8
    do pi = 1,max_patch_per_col
       do fc = 1,num_nolakec
          c = filter_nolakec(fc)
          if ( pi <= npfts(c) ) then
             p = pfti(c) + pi - 1
             l = plandunit(p)
             g = pgridcell(p)

             if (pactive(p)) then
                eflx_gnet(p) = sabg(p) + dlrad(p) &
!                              + (1-frac_veg_nosno(p))*emg(c)*forc_lwrad(c) - lwrad_emit(c) &
                               - (eflx_sh_grnd(p)+qflx_evap_soi(p)*htvp(c))
                ! save sabg for balancecheck, in case frac_sno is set to zero later
                sabg_chk(p) = frac_sno_eff(c) * sabg_snow(p) + (1._r8 - frac_sno_eff(c) ) * sabg_soil(p)

                eflx_gnet_snow = sabg_snow(p) + dlrad(p) &
!                    + (1-frac_veg_nosno(p))*emg(c)*forc_lwrad(c) - lwrad_emit_snow(c) &
                     - (eflx_sh_snow(p)+qflx_ev_snow(p)*htvp(c))

                eflx_gnet_soil = sabg_soil(p) + dlrad(p) &
!                    + (1-frac_veg_nosno(p))*emg(c)*forc_lwrad(c) - lwrad_emit_soil(c) &
                     - (eflx_sh_soil(p)+qflx_ev_soil(p)*htvp(c))

                eflx_gnet_h2osfc = sabg_soil(p) + dlrad(p) &
!                    + (1-frac_veg_nosno(p))*emg(c)*forc_lwrad(c) - lwrad_emit_h2osfc(c) &
                     - (eflx_sh_h2osfc(p)+qflx_ev_h2osfc(p)*htvp(c))
             end if
!            dgnetdT(p) = - cgrnd(p) - dlwrad_emit(c)
             dgnetdT(p) = - cgrnd(p)
             hs(c) = hs(c) + eflx_gnet(p) * pwtcol(p)
             dhsdT(c) = dhsdT(c) + dgnetdT(p) * pwtcol(p)
             ! separate surface fluxes for soil/snow
             hs_snow(c) = hs_snow(c) + eflx_gnet_snow * pwtcol(p)
             hs_soil(c) = hs_soil(c) + eflx_gnet_soil * pwtcol(p)
             hs_h2osfc(c) = hs_h2osfc(c) + eflx_gnet_h2osfc * pwtcol(p)
          end if
       end do
    end do

    !       Additional calculations with SNICAR:
    !       Set up tridiagonal matrix in a new manner. There is now
    !       absorbed solar radiation in each snow layer, instead of
    !       only the surface. Following the current implementation,
    !       absorbed solar flux should be: S + ((delS/delT)*dT),
    !       where S is absorbed radiation, and T is temperature. Now,
    !       assume delS/delT is zero, then it is OK to just add S
    !       to each layer

    ! Initialize:
    sabg_lyr_col(begc:endc,-nlevsno+1:1) = 0._r8
    hs_top(begc:endc)      = 0._r8
    hs_top_snow(begc:endc) = 0._r8
    hs_top_soil(begc:endc) = 0._r8

    do pi = 1,max_patch_per_col
       do fc = 1,num_nolakec
          c = filter_nolakec(fc)
          lyr_top = snl(c) + 1
          if ( pi <= npfts(c) ) then
             p = pfti(c) + pi - 1
             if (pactive(p)) then
                g = pgridcell(p)
                l = plandunit(p)
                eflx_gnet_top = sabg_lyr(p,lyr_top) + dlrad(p) &
!                    + (1-frac_veg_nosno(p))*emg(c)*forc_lwrad(c) - lwrad_emit(c) &
                     - (eflx_sh_grnd(p)+qflx_evap_soi(p)*htvp(c))

                hs_top(c) = hs_top(c) + eflx_gnet_top*pwtcol(p)

                eflx_gnet_snow = sabg_lyr(p,lyr_top) + dlrad(p) &
!                    + (1-frac_veg_nosno(p))*emg(c)*forc_lwrad(c) - lwrad_emit_snow(c) &
                     - (eflx_sh_snow(p)+qflx_ev_snow(p)*htvp(c))

                eflx_gnet_soil = sabg_lyr(p,lyr_top) + dlrad(p) &
!                    + (1-frac_veg_nosno(p))*emg(c)*forc_lwrad(c) - lwrad_emit_soil(c) &
                     - (eflx_sh_soil(p)+qflx_ev_soil(p)*htvp(c))

                hs_top_snow(c) = hs_top_snow(c) + eflx_gnet_snow*pwtcol(p)
                hs_top_soil(c) = hs_top_soil(c) + eflx_gnet_soil*pwtcol(p)

                do j = lyr_top,1,1
                   sabg_lyr_col(c,j) = sabg_lyr_col(c,j) + sabg_lyr(p,j) * pwtcol(p)
                enddo
             endif
          endif
       enddo
    enddo

    ! Determine heat diffusion through the layer interface and factor used in computing
    ! tridiagonal matrix and set up vector r and vectors a, b, c that define tridiagonal
    ! matrix and solve system

    do j = -nlevsno+1,nlevgrnd
       do fc = 1,num_nolakec
          c = filter_nolakec(fc)
          l = clandunit(c)
             if (j >= snl(c)+1) then
                if (j == snl(c)+1) then
                   fact(c,j) = dtime/cv(c,j) * dz(c,j) / (0.5_r8*(z(c,j)-zi(c,j-1)+capr*(z(c,j+1)-zi(c,j-1))))
                   fn(c,j) = tk(c,j)*(t_soisno(c,j+1)-t_soisno(c,j))/(z(c,j+1)-z(c,j))
                else if (j <= nlevgrnd-1) then
                   fact(c,j) = dtime/cv(c,j)
                   fn(c,j) = tk(c,j)*(t_soisno(c,j+1)-t_soisno(c,j))/(z(c,j+1)-z(c,j))
                   dzm     = (z(c,j)-z(c,j-1))
                else if (j == nlevgrnd) then
                   fact(c,j) = dtime/cv(c,j)
                   fn(c,j) = eflx_bot(c)
                end if
             end if
       end do
    end do

    do j = -nlevsno+1,nlevgrnd
       do fc = 1,num_nolakec
          c = filter_nolakec(fc)
             if (j >= snl(c)+1) then
                if (j == snl(c)+1) then
                   dzp     = z(c,j+1)-z(c,j)
                   at(c,j) = 0._r8
                   bt(c,j) = 1+(1._r8-cnfac)*fact(c,j)*tk(c,j)/dzp-fact(c,j)*dhsdT(c)
                   ct(c,j) =  -(1._r8-cnfac)*fact(c,j)*tk(c,j)/dzp
                   rt(c,j) = t_soisno(c,j) +  fact(c,j)*( hs_top_snow(c) &
                        - dhsdT(c)*t_soisno(c,j) + cnfac*fn(c,j) )
                else if (j == 1) then
                   ! this is the snow/soil interface layer
                   dzm     = (z(c,j)-z(c,j-1))
                   dzp     = (z(c,j+1)-z(c,j))

                   at(c,j) =   - frac_sno_eff(c) * (1._r8-cnfac) * fact(c,j) &
                        * tk(c,j-1)/dzm

                   bt(c,j) = 1._r8 + (1._r8-cnfac)*fact(c,j)*(tk(c,j)/dzp &
                        + frac_sno_eff(c) * tk(c,j-1)/dzm) &
                        - (1._r8 - frac_sno_eff(c))*fact(c,j)*dhsdT(c)

                   ct(c,j) = - (1._r8-cnfac)*fact(c,j)*tk(c,j)/dzp

                   rt(c,j) = t_soisno(c,j) + fact(c,j) &
                        *((1._r8-frac_sno_eff(c))*(hs_soil(c) - dhsdT(c)*t_soisno(c,j)) &
                        + cnfac*(fn(c,j) - frac_sno_eff(c) * fn(c,j-1)))

                   rt(c,j) = rt(c,j) +  frac_sno_eff(c)*fact(c,j)*sabg_lyr_col(c,j)
                else if (j <= nlevgrnd-1) then
                   dzm     = (z(c,j)-z(c,j-1))
                   dzp     = (z(c,j+1)-z(c,j))
                   at(c,j) =   - (1._r8-cnfac)*fact(c,j)* tk(c,j-1)/dzm
                   bt(c,j) = 1._r8+ (1._r8-cnfac)*fact(c,j)*(tk(c,j)/dzp + tk(c,j-1)/dzm)
                   ct(c,j) =   - (1._r8-cnfac)*fact(c,j)* tk(c,j)/dzp

                   rt(c,j) = t_soisno(c,j) + cnfac*fact(c,j)*( fn(c,j) - fn(c,j-1) )
                   if (j < 1) rt(c,j) = rt(c,j) + fact(c,j)*sabg_lyr_col(c,j)

                else if (j == nlevgrnd) then
                   dzm     = (z(c,j)-z(c,j-1))
                   at(c,j) =   - (1._r8-cnfac)*fact(c,j)*tk(c,j-1)/dzm
                   bt(c,j) = 1._r8+ (1._r8-cnfac)*fact(c,j)*tk(c,j-1)/dzm
                   ct(c,j) = 0._r8
                   rt(c,j) = t_soisno(c,j) - cnfac*fact(c,j)*fn(c,j-1) + fact(c,j)*fn(c,j)
                end if
             end if
       enddo
    end do

    ! compute thermal properties of h2osfc
    do fc = 1,num_nolakec
       c = filter_nolakec(fc)
       dz_h2osfc(c)=max(1.0e-6_r8,1.0e-3*h2osfc(c))
       c_h2osfc(c)=cpliq*denh2o*dz_h2osfc(c) !"areametric" heat capacity [J/K/m^2]

    enddo

    ! set up compact matrix for band diagonal solver, requires additional
    !     sub/super diagonals (1 each), and one additional row for t_h2osfc
    jtop = -9999
    do fc = 1,num_nolakec
       c = filter_nolakec(fc)
       jtop(c) = snl(c)
       ! compute jbot
       jbot(c) = nlevgrnd
    end do

    ! initialize matrices for BandDiagonal
    bmatrix(begc:endc, :, :)=0.0
    tvector(begc:endc, :) = 1.e36_r8
    rvector(begc:endc, :) = 1.e36_r8

    ! the solution will be organized as (snow:h2osfc:soil) to minimize
    !     bandwidth; this requires a 5-element band instead of 3
    do fc = 1,num_nolakec
       c = filter_nolakec(fc)
!=========================================================================
       do j = snl(c)+1, 0
          ! snow layers; bottom layer will have one offset coefficient
          bmatrix(c,2,j-1)=ct(c,j)
          bmatrix(c,3,j-1)=bt(c,j)
          bmatrix(c,4,j-1)=at(c,j)

          rvector(c,j-1)=rt(c,j)
          tvector(c,j-1)=t_soisno(c,j)
       end do
       ! bottom snow layer has super coef shifted to 2nd super diagonal
       bmatrix(c,2,-1)=0.0
       bmatrix(c,1,-1)=ct(c,0) !flux to top soil layer
!=========================================================================
       ! surface water layer has two coefficients
       dzm=(0.5*dz_h2osfc(c)+z(c,1))

       bmatrix(c,2,0)= -(1._r8-cnfac)*(dtime/c_h2osfc(c))*tk_h2osfc(c)/dzm !flux to top soil layer
       bmatrix(c,3,0)= 1+(1._r8-cnfac)*(dtime/c_h2osfc(c)) &
            *tk_h2osfc(c)/dzm -(dtime/c_h2osfc(c))*dhsdT(c) !interaction from atm

       fn_h2osfc(c)=tk_h2osfc(c)*(t_soisno(c,1)-t_h2osfc(c))/dzm
       rvector(c,0)= t_h2osfc(c) +  (dtime/c_h2osfc(c)) &
            *( hs_h2osfc(c) - dhsdT(c)*t_h2osfc(c) + cnfac*fn_h2osfc(c) )!rhs for h2osfc

       tvector(c,0)=t_h2osfc(c)

!=========================================================================
       ! soil layers; top layer will have one offset and one extra coefficient
       bmatrix(c,2,1:nlevgrnd)=ct(c,1:nlevgrnd)
       bmatrix(c,3,1:nlevgrnd)=bt(c,1:nlevgrnd)
       bmatrix(c,4,1:nlevgrnd)=at(c,1:nlevgrnd)
       ! top soil layer has sub coef shifted to 2nd super diagonal
       if ( frac_h2osfc(c) /= 0.0_r8 )then
          bmatrix(c,4,1)=  - frac_h2osfc(c) * (1._r8-cnfac) * fact(c,1) &
               * tk_h2osfc(c)/dzm !flux from h2osfc
       else
          bmatrix(c,4,1)= 0.0_r8
       end if
       bmatrix(c,5,1)=at(c,1)
       ! diagonal element correction for presence of h2osfc
       if ( frac_h2osfc(c) /= 0.0_r8 )then
          bmatrix(c,3,1)=bmatrix(c,3,1)+ frac_h2osfc(c) &
               *((1._r8-cnfac)*fact(c,1)*tk_h2osfc(c)/dzm + fact(c,1)*dhsdT(c))
       end if

       rvector(c,1:nlevgrnd)=rt(c,1:nlevgrnd)
       ! rhs correction for h2osfc
       if ( frac_h2osfc(c) /= 0.0_r8 )then
          rvector(c,1)=rvector(c,1) &
               -frac_h2osfc(c)*fact(c,1)*((hs_soil(c) - dhsdT(c)*t_soisno(c,1)) &
               +cnfac*fn_h2osfc(c))
       end if

       tvector(c,1:nlevgrnd)=t_soisno(c,1:nlevgrnd)

    enddo

    call BandDiagonal(bounds, -nlevsno, nlevgrnd, jtop(begc:endc), jbot(begc:endc), &
    num_nolakec, filter_nolakec, nband, bmatrix(begc:endc, :, :), &
    rvector(begc:endc, :), tvector(begc:endc, :))

    ! return temperatures to original array
    do fc = 1,num_nolakec
       c = filter_nolakec(fc)
       do j = snl(c)+1, 0
          t_soisno(c,j) = tvector(c,j-1) !snow layers
       end do
       t_soisno(c,1:nlevgrnd)   = tvector(c,1:nlevgrnd)  !soil layers

       if(frac_h2osfc(c) == 0._r8) then
          t_h2osfc(c)=t_soisno(c,1)
       else
          t_h2osfc(c)              = tvector(c,0)           !surface water
       endif
    enddo

    ! Melting or Freezing

    do j = -nlevsno+1,nlevgrnd
       do fc = 1,num_nolakec
          c = filter_nolakec(fc)
          l = clandunit(c)
             if (j >= snl(c)+1) then
                if (j <= nlevgrnd-1) then
                   fn1(c,j) = tk(c,j)*(t_soisno(c,j+1)-t_soisno(c,j))/(z(c,j+1)-z(c,j))
                else if (j == nlevgrnd) then
                   fn1(c,j) = 0._r8
                end if
             end if
       end do
    end do

    ! compute terms needed for phase change of h2osfc
    do fc = 1,num_nolakec
       c = filter_nolakec(fc)
       dzp=(0.5*dz_h2osfc(c)+z(c,1))
       fn1_h2osfc(c)=tk_h2osfc(c)*(t_soisno(c,1)-t_h2osfc(c))/dzp
    enddo

    xmf_h2osfc=0.
    ! compute phase change of h2osfc
    call PhaseChangeH2osfc (bounds, num_nolakec, filter_nolakec, &
    fact(bounds%begc:bounds%endc, :), &
    dhsdT(bounds%begc:bounds%endc), &
    c_h2osfc(bounds%begc:bounds%endc), &
    xmf_h2osfc(bounds%begc:bounds%endc), &
    waterflux_inst, waterstate_inst, temperature_inst)

    call PhaseChange_beta (bounds, num_nolakec, filter_nolakec, &
    fact(bounds%begc:bounds%endc, :), &
    dhsdT(bounds%begc:bounds%endc), &
    xmf(bounds%begc:bounds%endc), &
    energyflux_inst, waterflux_inst, waterstate_inst, temperature_inst, soilstate_inst)

    do fc = 1,num_nolakec
       c = filter_nolakec(fc)
       ! this expression will (should) work whether there is snow or not
       if (snl(c) < 0) then
          if(frac_h2osfc(c) /= 0._r8) then
              t_grnd(c) = frac_sno_eff(c) * t_soisno(c,snl(c)+1) &
                   + (1.0_r8 - frac_sno_eff(c) - frac_h2osfc(c)) * t_soisno(c,1) &
                   + frac_h2osfc(c) * t_h2osfc(c)
          else
              t_grnd(c) = frac_sno_eff(c) * t_soisno(c,snl(c)+1) &
                   + (1.0_r8 - frac_sno_eff(c)) * t_soisno(c,1)
          end if

       else
          if(frac_h2osfc(c) /= 0._r8) then
             t_grnd(c) = (1 - frac_h2osfc(c)) * t_soisno(c,1) + frac_h2osfc(c) * t_h2osfc(c)
          else
             t_grnd(c) = t_soisno(c,1)
          end if
       endif
    end do

    ! Initialize soil heat content
    do fc = 1,num_nolakec
       c = filter_nolakec(fc)
       l = clandunit(c)
       hc_soisno(c) = 0._r8
       hc_soi(c)    = 0._r8
       eflx_fgr12(c)= 0._r8
    end do

    ! Calculate soil heat content and soil plus snow heat content
    do j = -nlevsno+1,nlevgrnd
       do fc = 1,num_nolakec
          c = filter_nolakec(fc)
          l = clandunit(c)

          if (j == 1) then ! this only needs to be done once
             eflx_fgr12(c) = -cnfac*fn(c,1) - (1._r8-cnfac)*fn1(c,1)
          end if
          if (j > 0 .and. j < nlevgrnd) then
             eflx_fgr(c,j) = -cnfac*fn(c,j) - (1._r8-cnfac)*fn1(c,j)
          else if (j == nlevgrnd) then
             eflx_fgr(c,j) = 0._r8
          end if

          if (j >= snl(c)+1) then
             hc_soisno(c) = hc_soisno(c) + cv(c,j)*t_soisno(c,j) / 1.e6_r8
          endif
          if (j >= 1) then
             hc_soi(c) = hc_soi(c) + cv(c,j)*t_soisno(c,j) / 1.e6_r8
          end if
       end do
    end do

    ! The rest of the code below is from Biogeophysics2Mod.F90 after
    ! SoilTemperature has been called

    do fc = 1,num_nolakec
       c = filter_nolakec(fc)
       j = snl(c)+1

       ! Calculate difference in soil temperature from last time step, for
       ! flux corrections

       if (snl(c) < 0) then
          t_grnd0(c) = frac_sno_eff(c) * tssbef(c,snl(c)+1) &
               + (1 - frac_sno_eff(c) - frac_h2osfc(c)) * tssbef(c,1) &
               + frac_h2osfc(c) * t_h2osfc_bef(c)
       else
          t_grnd0(c) = (1 - frac_h2osfc(c)) * tssbef(c,1) + frac_h2osfc(c) * t_h2osfc_bef(c)
       endif

       tinc(c) = t_grnd(c) - t_grnd0(c)

       ! Determine ratio of topsoil_evap_tot

       egsmax(c) = (h2osoi_ice(c,j)+h2osoi_liq(c,j)) / dtime

       ! added to trap very small negative soil water,ice

       if (egsmax(c) < 0._r8) then
          egsmax(c) = 0._r8
       end if
    end do

    ! A preliminary pft loop to determine if corrections are required for
    ! excess evaporation from the top soil layer... Includes new logic
    ! to distribute the corrections between pfts on the basis of their
    ! evaporative demands.
    ! egirat holds the ratio of demand to availability if demand is
    ! greater than availability, or 1.0 otherwise.
    ! Correct fluxes to present soil temperature

    do fp = 1,num_nolakep
       p = filter_nolakep(fp)
       c = pcolumn(p)
       eflx_sh_grnd(p) = eflx_sh_grnd(p) + tinc(c)*cgrnds(p)
       qflx_evap_soi(p) = qflx_evap_soi(p) + tinc(c)*cgrndl(p)

       qflx_ev_snow(p) = qflx_ev_snow(p) + tinc(c)*cgrndl(p)
       qflx_ev_soil(p) = qflx_ev_soil(p) + tinc(c)*cgrndl(p)
       qflx_ev_h2osfc(p) = qflx_ev_h2osfc(p) + tinc(c)*cgrndl(p)
    end do

    ! Set the column-average qflx_evap_soi as the weighted average over all pfts
    ! but only count the pfts that are evaporating

    do fc = 1,num_nolakec
       c = filter_nolakec(fc)
       topsoil_evap_tot(c) = 0._r8
       sumwt(c) = 0._r8
    end do

    do pi = 1,max_patch_per_col
       do fc = 1,num_nolakec
          c = filter_nolakec(fc)
          if ( pi <= npfts(c) ) then
             p = pfti(c) + pi - 1
             if (pactive(p)) then
                topsoil_evap_tot(c) = topsoil_evap_tot(c) + qflx_evap_soi(p) * pwtcol(p)
             end if
          end if
       end do
    end do

    ! Calculate ratio for rescaling pft-level fluxes to meet availability

    do fc = 1,num_nolakec
       c = filter_nolakec(fc)
       if (topsoil_evap_tot(c) > egsmax(c)) then
          egirat(c) = (egsmax(c)/topsoil_evap_tot(c))
       else
          egirat(c) = 1.0_r8
       end if
    end do

    do fp = 1,num_nolakep
       p = filter_nolakep(fp)
       c = pcolumn(p)
       l = plandunit(p)
       g = pgridcell(p)
       j = snl(c)+1

       ! Correct soil fluxes for possible evaporation in excess of top layer water
       ! excess energy is added to the sensible heat flux from soil

       if (egirat(c) < 1.0_r8) then
          save_qflx_evap_soi = qflx_evap_soi(p)
          qflx_evap_soi(p) = qflx_evap_soi(p) * egirat(c)
          eflx_sh_grnd(p) = eflx_sh_grnd(p) + (save_qflx_evap_soi - qflx_evap_soi(p))*htvp(c)
          qflx_ev_snow(p) = qflx_ev_snow(p) * egirat(c)
          qflx_ev_soil(p) = qflx_ev_soil(p) * egirat(c)
          qflx_ev_h2osfc(p) = qflx_ev_h2osfc(p) * egirat(c)
       end if

       ! Ground heat flux

       lw_grnd=(frac_sno_eff(c)*tssbef(c,snl(c)+1)**4 &
            +(1._r8-frac_sno_eff(c)-frac_h2osfc(c))*tssbef(c,1)**4 &
            +frac_h2osfc(c)*t_h2osfc_bef(c)**4)

       eflx_soil_grnd(p) = ((1._r8- frac_sno_eff(c))*sabg_soil(p) + frac_sno_eff(c)*sabg_snow(p)) + dlrad(p) &
!           + (1-frac_veg_nosno(p))*emg(c)*forc_lwrad(c) &
!           - emg(c)*sb*lw_grnd - emg(c)*sb*t_grnd0(c)**3*(4._r8*tinc(c)) &
            - (eflx_sh_grnd(p)+qflx_evap_soi(p)*htvp(c))

    end do

    ! Soil Energy balance check: errsoi is not in CLM5

!   do fp = 1,num_nolakep
!      p = filter_nolakep(fp)
!      c = pcolumn(p)
!      errsoi_pft(p) = eflx_soil_grnd(p) - xmf(c) - xmf_h2osfc(c) &
!           - frac_h2osfc(c)*(t_h2osfc(c)-t_h2osfc_bef(c)) &
!           *(c_h2osfc(c)/dtime)
!   end do
!   do j = -nlevsno+1,nlevgrnd
!      do fp = 1,num_nolakep
!         p = filter_nolakep(fp)
!         c = pcolumn(p)
!
!         ! area weight heat absorbed by snow layers
!         if (j >= snl(c)+1 .and. j < 1) errsoi_pft(p) = errsoi_pft(p) &
!              - frac_sno_eff(c)*(t_soisno(c,j)-tssbef(c,j))/fact(c,j)
!         if (j >= 1) errsoi_pft(p) = errsoi_pft(p) &
!              - (t_soisno(c,j)-tssbef(c,j))/fact(c,j)
!      end do
!   end do

    ! Outgoing long-wave radiation from vegetation + ground
    ! For conservation we put the increase of ground longwave to outgoing

    do fp = 1,num_nolakep
       p = filter_nolakep(fp)
       c = pcolumn(p)
       l = plandunit(p)
       g = pgridcell(p)
       j = snl(c)+1

       lw_grnd=(frac_sno_eff(c)*tssbef(c,snl(c)+1)**4 &
            +(1._r8-frac_sno_eff(c)-frac_h2osfc(c))*tssbef(c,1)**4 &
            +frac_h2osfc(c)*t_h2osfc_bef(c)**4)

!      eflx_lwrad_out(p) = ulrad(p) &
!           + (1-frac_veg_nosno(p))*(1.-emg(c))*forc_lwrad(c) &
!           + (1-frac_veg_nosno(p))*emg(c)*sb*lw_grnd &
!           + 4._r8*emg(c)*sb*t_grnd0(c)**3*tinc(c)
       eflx_lwrad_out(p) = ulrad(p)

       eflx_lwrad_net(p) = eflx_lwrad_out(p) - forc_lwrad(c)
    end do

    end associate
   end subroutine SoilTemperature

   !-----------------------------------------------------------------------
   subroutine SoilThermProp (bounds, num_nolakec, filter_nolakec, tk, cv, tk_h2osfc, &
waterstate_inst, temperature_inst, soilstate_inst)
     !
     ! !DESCRIPTION:
     ! Calculation of thermal conductivities and heat capacities of
     ! snow/soil layers
     ! (1) The volumetric heat capacity is calculated as a linear combination
     !     in terms of the volumetric fraction of the constituent phases.
     !
     ! (2) The thermal conductivity of soil is computed from the algorithm of
     !     Johansen (as reported by Farouki 1981), and of snow is from the
     !     formulation used in SNTHERM (Jordan 1991).
     ! The thermal conductivities at the interfaces between two neighboring
     ! layers (j, j+1) are derived from an assumption that the flux across
     ! the interface is equal to that from the node j to the interface and the
     ! flux from the interface to the node j+1.
     !
     ! !USES:
     use clm_varcon, only : denh2o, denice, tfrz, tkwat, tkice, tkair, &
                            cpice,  cpliq, thk_bedrock
     use clm_varpar, only : nlevsno, nlevgrnd, nlevsoi
     use clm_varctl, only : iulog
     use ColumnType, only : col
     !
     ! !ARGUMENTS:
     implicit none
     type(bounds_type), intent(in) :: bounds                   ! bounds
     integer , intent(in)  :: num_nolakec                      ! number of column non-lake points in column filter
     integer , intent(in)  :: filter_nolakec(:)                ! column filter for non-lake points
     real(r8), intent(out) :: cv( bounds%begc: , -nlevsno+1: ) ! heat capacity [J/(m2 K)] [col, lev]
     real(r8), intent(out) :: tk( bounds%begc: , -nlevsno+1: ) ! thermal conductivity at the layer interface [W/(m K)] [col, lev]
     real(r8), intent(out) :: tk_h2osfc( bounds%begc: )        ! thermal conductivity of h2osfc [W/(m K)] [col]
     type(waterstate_type) :: waterstate_inst
     type(temperature_type) :: temperature_inst
     type(soilstate_type) :: soilstate_inst
     !
     ! !LOCAL VARIABLES:
     integer  :: l,c,j                     ! indices
     integer  :: fc                        ! lake filtered column indices
     real(r8) :: dksat                     ! thermal conductivity for saturated soil (j/(k s m))
     real(r8) :: dke                       ! kersten number
     real(r8) :: fl                        ! volume fraction of liquid or unfrozen water to total water
     real(r8) :: satw                      ! relative total water content of soil.
     real(r8) :: zh2osfc
     !-----------------------------------------------------------------------

   associate( &
   h2osfc      =>  waterstate_inst%h2osfc_col       , & ! Input:  [real(r8) (:)]  surface (mm H2O)
   frac_sno    =>  waterstate_inst%frac_sno_eff_col , & ! Input:  [real(r8) (:)]  fractional snow covered area
   clandunit   =>  col%landunit                     , & ! Input:  [integer (:)]  column's landunit
   snl         =>  col%snl                          , & ! Input:  [integer (:)]  number of snow layers
   h2osno      =>  waterstate_inst%h2osno_col       , & ! Input:  [real(r8) (:)]  snow water (mm H2O)
   watsat      =>  soilstate_inst%watsat_col        , & ! Input:  [real(r8) (:,:)]  volumetric soil water at saturation (porosity)
   tkmg        =>  soilstate_inst%tkmg_col          , & ! Input:  [real(r8) (:,:)]  thermal conductivity, soil minerals  [W/m-K]
   tkdry       =>  soilstate_inst%tkdry_col         , & ! Input:  [real(r8) (:,:)]  thermal conductivity, dry soil (W/m/Kelvin)
   csol        =>  soilstate_inst%csol_col          , & ! Input:  [real(r8) (:,:)]  heat capacity, soil solids (J/m**3/Kelvin)
   thk         =>  soilstate_inst%thk_col           , & ! Output: [real(r8) (:,:)]  thermal conductivity of each layer  [W/m-K] (-nlevsno+1:nlevgrnd)
   dz          =>  col%dz                           , & ! Input:  [real(r8) (:,:)]  layer depth (m)
   zi          =>  col%zi                           , & ! Input:  [real(r8) (:,:)]  interface level below a "z" level (m)
   z           =>  col%z                            , & ! Input:  [real(r8) (:,:)]  layer thickness (m)
   t_soisno    =>  temperature_inst%t_soisno_col    , & ! Input:  [real(r8) (:,:)]  soil temperature (Kelvin)
   h2osoi_liq  =>  waterstate_inst%h2osoi_liq_col   , & ! Input:  [real(r8) (:,:)]  liquid water (kg/m2)
   h2osoi_ice  =>  waterstate_inst%h2osoi_ice_col   , & ! Input:  [real(r8) (:,:)]  ice lens (kg/m2)
   bw          =>  waterstate_inst%bw_col             & ! Output: [real(r8) (:,:)]  partial density of water in the snow pack (ice + liquid) [kg/m3] (-nlevsno+1:0)
   )

    ! Thermal conductivity of soil from Farouki (1981)
    ! Urban values are from Masson et al. 2002, Evaluation of the Town Energy Balance (TEB)
    ! scheme with direct measurements from dry districts in two cities, J. Appl. Meteorol.,
    ! 41, 1011-1026.

    do j = -nlevsno+1,nlevgrnd
       do fc = 1, num_nolakec
          c = filter_nolakec(fc)

          ! Only examine levels from 1->nlevgrnd
          if (j >= 1) then
             l = clandunit(c)
                satw = (h2osoi_liq(c,j)/denh2o + h2osoi_ice(c,j)/denice)/(dz(c,j)*watsat(c,j))
                satw = min(1._r8, satw)
                if (satw > .1e-6_r8) then
                   if (t_soisno(c,j) >= tfrz) then       ! Unfrozen soil
                      dke = max(0._r8, log10(satw) + 1.0_r8)
                   else                               ! Frozen soil
                      dke = satw
                   end if
                   fl = (h2osoi_liq(c,j)/(denh2o*dz(c,j))) / (h2osoi_liq(c,j)/(denh2o*dz(c,j)) + &
                                                              h2osoi_ice(c,j)/(denice*dz(c,j)))
                   dksat = tkmg(c,j)*tkwat**(fl*watsat(c,j))*tkice**((1._r8-fl)*watsat(c,j))
                   thk(c,j) = dke*dksat + (1._r8-dke)*tkdry(c,j)
                else
                   thk(c,j) = tkdry(c,j)
                endif
                if (j > nlevsoi) thk(c,j) = thk_bedrock
          endif

          ! Thermal conductivity of snow, which from Jordan (1991) pp. 18
          ! Only examine levels from snl(c)+1 -> 0 where snl(c) < 1
          if (snl(c)+1 < 1 .AND. (j >= snl(c)+1) .AND. (j <= 0)) then
             bw(c,j) = (h2osoi_ice(c,j)+h2osoi_liq(c,j))/(frac_sno(c)*dz(c,j))
             thk(c,j) = tkair + (7.75e-5_r8 *bw(c,j) + 1.105e-6_r8*bw(c,j)*bw(c,j))*(tkice-tkair)
          end if

       end do
    end do

    ! Thermal conductivity at the layer interface

    do j = -nlevsno+1,nlevgrnd
       do fc = 1,num_nolakec
          c = filter_nolakec(fc)
            if (j >= snl(c)+1 .AND. j <= nlevgrnd-1) then
               tk(c,j) = thk(c,j)*thk(c,j+1)*(z(c,j+1)-z(c,j)) &
                         /(thk(c,j)*(z(c,j+1)-zi(c,j))+thk(c,j+1)*(zi(c,j)-z(c,j)))
            else if (j == nlevgrnd) then
               tk(c,j) = 0._r8
            end if
       end do
    end do

    ! calculate thermal conductivity of h2osfc
    do fc = 1, num_nolakec
       c = filter_nolakec(fc)
       zh2osfc=1.0e-3*(0.5*h2osfc(c)) !convert to [m] from [mm]
       tk_h2osfc(c)= tkwat*thk(c,1)*(z(c,1)+zh2osfc) &
                       /(tkwat*z(c,1)+thk(c,1)*zh2osfc)
    enddo

    ! Soil heat capacity, from de Vires (1963)
    ! Urban values are from Masson et al. 2002, Evaluation of the Town Energy Balance (TEB)
    ! scheme with direct measurements from dry districts in two cities, J. Appl. Meteorol.,
    ! 41, 1011-1026.

    do j = 1, nlevgrnd
       do fc = 1,num_nolakec
          c = filter_nolakec(fc)
          l = clandunit(c)
          cv(c,j) = csol(c,j)*(1-watsat(c,j))*dz(c,j) + (h2osoi_ice(c,j)*cpice + h2osoi_liq(c,j)*cpliq)
          if (j == 1) then
             if (snl(c)+1 == 1 .AND. h2osno(c) > 0._r8) then
                cv(c,j) = cv(c,j) + cpice*h2osno(c)
             end if
          end if
       enddo
    end do

    ! Snow heat capacity

    do j = -nlevsno+1,0
       do fc = 1,num_nolakec
          c = filter_nolakec(fc)
          if (snl(c)+1 < 1 .and. j >= snl(c)+1) then
             cv(c,j) = cpliq*h2osoi_liq(c,j) + cpice*h2osoi_ice(c,j)
          end if
       end do
    end do

    end associate
   end subroutine SoilThermProp

   !-----------------------------------------------------------------------
   subroutine PhaseChangeH2osfc (bounds, num_nolakec, filter_nolakec, fact, &
   dhsdT, c_h2osfc, xmf_h2osfc, waterflux_inst, waterstate_inst, temperature_inst)
     !
     ! !DESCRIPTION:
     ! Only freezing is considered.  When water freezes, move ice to bottom snow layer.
     !
     ! !USES:
     use clm_time_manager, only : get_step_size
     use clm_varcon  , only : tfrz, hfus, grav,denice,cnfac,cpice
     use clm_varpar  , only : nlevsno, nlevgrnd
     use clm_varctl  , only : iulog
     use ColumnType, only : col
     !
     ! !ARGUMENTS:
     implicit none
     type(bounds_type), intent(in) :: bounds                         ! bounds
     integer , intent(in)    :: num_nolakec                          ! number of column non-lake points in column filter
     integer , intent(in)    :: filter_nolakec(:)                    ! column filter for non-lake points
     real(r8), intent(inout) :: fact  ( bounds%begc: , -nlevsno+1: ) ! temporary [col, lev]
     real(r8), intent(in)    :: dhsdT ( bounds%begc: )               ! temperature derivative of "hs" [col]
     real(r8), intent(in)    :: c_h2osfc( bounds%begc: )             ! heat capacity of surface water [col]
     real(r8), intent(out)   :: xmf_h2osfc( bounds%begc: )           ! latent heat of phase change of surface water [col]
     type(waterflux_type) :: waterflux_inst
     type(waterstate_type) :: waterstate_inst
     type(temperature_type) :: temperature_inst
     !
     ! !LOCAL VARIABLES:
     integer  :: j,c,g                              !do loop index
     integer  :: fc                                 !lake filtered column indices
     real(r8) :: dtime                              !land model time step (sec)
     real(r8) :: heatr                              !energy residual or loss after melting or freezing
     real(r8) :: temp1                              !temporary variables [kg/m2]
     real(r8) :: hm(bounds%begc:bounds%endc)                        !energy residual [W/m2]
     real(r8) :: xm(bounds%begc:bounds%endc)                        !melting or freezing within a time step [kg/m2]
     real(r8) :: tinc                               !t(n+1)-t(n) (K)
     real(r8) :: smp                                !frozen water potential (mm)
     real(r8) :: rho_avg
     real(r8) :: z_avg
     real(r8) :: dcv(bounds%begc:bounds%endc)
     real(r8) :: t_h2osfc_new
     real(r8) :: c1
     real(r8) :: c2
     real(r8) :: h_excess
     real(r8) :: c_h2osfc_ice
     !-----------------------------------------------------------------------

   associate( &
   frac_sno           =>  waterstate_inst%frac_sno_eff_col       , & ! Input:  [real(r8) (:)]  fraction of ground covered by snow (0 to 1)
   frac_h2osfc        =>  waterstate_inst%frac_h2osfc_col        , & ! Input:  [real(r8) (:)]  fraction of ground covered by surface water (0 to 1)
   t_h2osfc           =>  temperature_inst%t_h2osfc_col          , & ! Input:  [real(r8) (:)]  surface water temperature
   t_h2osfc_bef       =>  temperature_inst%t_h2osfc_bef_col      , & ! Input:  [real(r8) (:)]  saved surface water temperature
   h2osfc             =>  waterstate_inst%h2osfc_col             , & ! Input:  [real(r8) (:)]  surface water (mm)
   int_snow           =>  waterstate_inst%int_snow_col           , & ! Input:  [real(r8) (:)]  integrated snowfall [mm]
   qflx_h2osfc_to_ice =>  waterflux_inst%qflx_h2osfc_to_ice_col  , & ! Input:  [real(r8) (:)]  conversion of h2osfc to ice
   snl                =>  col%snl                                , & ! Input:  [integer (:)]  number of snow layers
   h2osno             =>  waterstate_inst%h2osno_col             , & ! Input:  [real(r8) (:)]  snow water (mm H2O)
   snow_depth         =>  waterstate_inst%snow_depth_col         , & ! Input:  [real(r8) (:)] snow height (m)
   h2osoi_ice         =>  waterstate_inst%h2osoi_ice_col         , & ! Input:  [real(r8) (:,:)] ice lens (kg/m2) (new)
   t_soisno           =>  temperature_inst%t_soisno_col          , & ! Input:  [real(r8) (:,:)] soil temperature (Kelvin)
   dz                 =>  col%dz                                   & ! Input:  [real(r8) (:,:)] layer thickness (m)
   )

    ! Get step size

    dtime = get_step_size()

    ! Initialization

    do fc = 1,num_nolakec
       c = filter_nolakec(fc)
       xmf_h2osfc(c) = 0._r8
       hm(c) = 0._r8
       xm(c) = 0._r8
       qflx_h2osfc_to_ice(c) = 0._r8
    end do

    ! Freezing identification
    do fc = 1,num_nolakec
       c = filter_nolakec(fc)
       ! If liquid exists below melt point, freeze some to ice.
       if ( frac_h2osfc(c) > 0._r8 .AND. t_h2osfc(c) <= tfrz) then
          tinc = t_h2osfc(c)-tfrz
          t_h2osfc(c) = tfrz
          ! energy absorbed beyond freezing temperature
          hm(c) = dhsdT(c)*tinc - tinc*c_h2osfc(c)/dtime

          ! mass of water converted from liquid to ice
          xm(c) = hm(c)*dtime/hfus
          temp1 = h2osfc(c) - xm(c)

          ! compute change in cv due to additional ice
          dcv(c)=cpice*min(xm(c),h2osfc(c))

          z_avg=frac_sno(c)*snow_depth(c)
          if (z_avg > 0._r8) then
             rho_avg=min(800._r8,h2osno(c)/z_avg)
          else
             rho_avg=200._r8
          endif
!=====================  xm < h2osfc  ====================================
          if(temp1 >= 0._r8) then ! add some frozen water to snow column
             ! add ice to snow column
             h2osno(c) = h2osno(c) + xm(c)
             int_snow(c) = int_snow(c) + xm(c)

             if(snl(c) < 0) h2osoi_ice(c,0) = h2osoi_ice(c,0) + xm(c)

             ! remove ice from h2osfc
             h2osfc(c) = h2osfc(c) - xm(c)

             xmf_h2osfc(c) = -frac_h2osfc(c)*hm(c)

             qflx_h2osfc_to_ice(c) = xm(c)/dtime

             ! update snow depth
             if (frac_sno(c) > 0 .and. snl(c) < 0) then
                snow_depth(c)=h2osno(c)/(rho_avg*frac_sno(c))
             else
                snow_depth(c)=h2osno(c)/denice
             endif
!=========================  xm > h2osfc  =============================
          else !all h2osfc converted to ice, apply residual heat to top soil layer

             rho_avg=(h2osno(c)*rho_avg + h2osfc(c)*denice)/(h2osno(c) + h2osfc(c))
             h2osno(c) = h2osno(c) + h2osfc(c)
             int_snow(c) = int_snow(c) + h2osfc(c)

             qflx_h2osfc_to_ice(c) = h2osfc(c)/dtime

             ! excess energy is used to cool ice layer
             if(snl(c) < 0) h2osoi_ice(c,0) = h2osoi_ice(c,0) + h2osfc(c)

             ! compute heat capacity of frozen h2osfc layer
             c_h2osfc_ice=cpice*denice*(1.0e-3*h2osfc(c)) !h2osfc in [m]

             ! cool frozen h2osfc layer with extra heat
             t_h2osfc_new = t_h2osfc(c) - temp1*hfus/(dtime*dhsdT(c) - c_h2osfc_ice)

             ! next, determine equilibrium temperature of combined ice/snow layer
             xmf_h2osfc(c) = -frac_h2osfc(c)*hm(c)
             if (snl(c) == 0) then
                t_soisno(c,0) = t_h2osfc_new
             else if (snl(c) == -1) then
                c1=frac_sno(c)/fact(c,0) - dhsdT(c)*dtime
                if ( frac_h2osfc(c) /= 0.0_r8 )then
                   c2=frac_h2osfc(c)*(c_h2osfc_ice/dtime)
                else
                   c2=0.0_r8
                end if
                ! account for the change in t_soisno(c,0) via xmf_h2osfc(c)
                xmf_h2osfc(c) = xmf_h2osfc(c) + frac_sno(c)*t_soisno(c,0)/fact(c,0)
                t_soisno(c,0) = (c1*t_soisno(c,0)+ c2*t_h2osfc_new) &
                                   /(c1 + c2)
                xmf_h2osfc(c) = xmf_h2osfc(c) - frac_sno(c)*t_soisno(c,0)/fact(c,0)

             else
                c1=frac_sno(c)/fact(c,0)
                if ( frac_h2osfc(c) /= 0.0_r8 )then
                   c2=frac_h2osfc(c)*(c_h2osfc_ice/dtime)
                else
                   c2=0.0_r8
                end if
                xmf_h2osfc(c) = xmf_h2osfc(c) + c1*t_soisno(c,0)
                t_soisno(c,0) = (c1*t_soisno(c,0)+ c2*t_h2osfc_new) &
                               /(c1 + c2)
                xmf_h2osfc(c) = xmf_h2osfc(c) - c1*t_soisno(c,0)
             endif

             ! set h2osfc to zero (all liquid converted to ice)
             h2osfc(c) = 0._r8

             ! update snow depth
             if (frac_sno(c) > 0 .and. snl(c) < 0) then
                snow_depth(c)=h2osno(c)/(rho_avg*frac_sno(c))
             else
                snow_depth(c)=h2osno(c)/denice
             endif

          endif
       endif
    enddo
    end associate
   end subroutine PhaseChangeH2osfc

   !-----------------------------------------------------------------------
   subroutine PhaseChange_beta (bounds, num_nolakec, filter_nolakec, fact, dhsdT, xmf, &
   energyflux_inst, waterflux_inst, waterstate_inst, temperature_inst, soilstate_inst)
     !
     ! !DESCRIPTION:
     ! Calculation of the phase change within snow and soil layers:
     ! (1) Check the conditions for which the phase change may take place,
     !     i.e., the layer temperature is great than the freezing point
     !     and the ice mass is not equal to zero (i.e. melting),
     !     or the layer temperature is less than the freezing point
     !     and the liquid water mass is greater than the allowable supercooled
     !     liquid water calculated from freezing point depression (i.e. freezing).
     ! (2) Assess the rate of phase change from the energy excess (or deficit)
     !     after setting the layer temperature to freezing point.
     ! (3) Re-adjust the ice and liquid mass, and the layer temperature
     !
     ! !USES:
     use clm_time_manager, only : get_step_size
     use clm_varcon, only : tfrz, hfus, grav
     use clm_varpar, only : nlevsno, nlevgrnd
     use clm_varctl, only : iulog
     use ColumnType, only : col
     !
     ! !ARGUMENTS:
     implicit none
     type(bounds_type), intent(in) :: bounds                      ! bounds
     integer , intent(in) :: num_nolakec                          ! number of column non-lake points in column filter
     integer , intent(in) :: filter_nolakec(:)                    ! column filter for non-lake points
     real(r8), intent(in) :: fact  ( bounds%begc: , -nlevsno+1: ) ! temporary [col, lev]
     real(r8), intent(in) :: dhsdT ( bounds%begc: )               ! temperature derivative of "hs" [col]
     real(r8), intent(out):: xmf   ( bounds%begc: )               ! total latent heat of phase change [col]
     type(energyflux_type) :: energyflux_inst
     type(waterflux_type) :: waterflux_inst
     type(waterstate_type) :: waterstate_inst
     type(temperature_type) :: temperature_inst
     type(soilstate_type) :: soilstate_inst
     !
     ! !LOCAL VARIABLES:
     integer  :: j,c,g,l                            !do loop index
     integer  :: fc                                 !lake filtered column indices
     real(r8) :: dtime                              !land model time step (sec)
     real(r8) :: heatr                              !energy residual or loss after melting or freezing
     real(r8) :: temp1                              !temporary variables [kg/m2]
     real(r8) :: hm(bounds%begc:bounds%endc,-nlevsno+1:nlevgrnd)    !energy residual [W/m2]
     real(r8) :: xm(bounds%begc:bounds%endc,-nlevsno+1:nlevgrnd)    !melting or freezing within a time step [kg/m2]
     real(r8) :: wmass0(bounds%begc:bounds%endc,-nlevsno+1:nlevgrnd)!initial mass of ice and liquid (kg/m2)
     real(r8) :: wice0 (bounds%begc:bounds%endc,-nlevsno+1:nlevgrnd)!initial mass of ice (kg/m2)
     real(r8) :: wliq0 (bounds%begc:bounds%endc,-nlevsno+1:nlevgrnd)!initial mass of liquid (kg/m2)
     real(r8) :: supercool(bounds%begc:bounds%endc,nlevgrnd)        !supercooled water in soil (kg/m2)
     real(r8) :: propor                             !proportionality constant (-)
     real(r8) :: tinc(bounds%begc:bounds%endc,-nlevsno+1:nlevgrnd)  !t(n+1)-t(n) (K)
     real(r8) :: smp                                !frozen water potential (mm)
     !-----------------------------------------------------------------------

   associate( &
   qflx_snow_drain  =>  waterflux_inst%qflx_snow_drain_col   , & ! Input:  [real(r8) (:)]  drainage from snow pack
   frac_sno_eff     =>  waterstate_inst%frac_sno_eff_col     , & ! Input:  [real(r8) (:)]  eff. fraction of ground covered by snow (0 to 1)
   frac_sno         =>  waterstate_inst%frac_sno_col         , & ! Input:  [real(r8) (:)]  fraction of ground covered by snow (0 to 1)
   frac_h2osfc      =>  waterstate_inst%frac_h2osfc_col      , & ! Input:  [real(r8) (:)]  fraction of ground covered by surface water (0 to 1)
   snl              =>  col%snl                              , & ! Input:  [integer (:)]  number of snow layers
   h2osno           =>  waterstate_inst%h2osno_col           , & ! Input:  [real(r8) (:)]  snow water (mm H2O)
   snow_depth       =>  waterstate_inst%snow_depth_col       , & ! Input:  [real(r8) (:)]  snow height (m)
   qflx_snomelt     =>  waterflux_inst%qflx_snomelt_col      , & ! Output: [real(r8) (:)]  snow melt (mm H2O /s)
   eflx_snomelt     =>  energyflux_inst%eflx_snomelt_col     , & ! Output: [real(r8) (:)] snow melt heat flux (W/m**2)
   h2osoi_liq       =>  waterstate_inst%h2osoi_liq_col       , & ! Input:  [real(r8) (:,:)] liquid water (kg/m2) (new)
   h2osoi_ice       =>  waterstate_inst%h2osoi_ice_col       , & ! Input:  [real(r8) (:,:)] ice lens (kg/m2) (new)
   imelt            =>  temperature_inst%imelt_col           , & ! Output: [integer (:,:)] flag for melting (=1), freezing (=2), Not=0 (new)
   t_soisno         =>  temperature_inst%t_soisno_col        , & ! Input:  [real(r8) (:,:)] soil temperature (Kelvin)
   bsw              =>  soilstate_inst%bsw_col               , & ! Input:  [real(r8) (:,:)] Clapp and Hornberger "b"
   sucsat           =>  soilstate_inst%sucsat_col            , & ! Input:  [real(r8) (:,:)] minimum soil suction (mm)
   watsat           =>  soilstate_inst%watsat_col            , & ! Input:  [real(r8) (:,:)] volumetric soil water at saturation (porosity)
   dz               =>  col%dz                               , & ! Input:  [real(r8) (:,:)] layer thickness (m)
   qflx_snofrz_lyr  =>  waterflux_inst%qflx_snofrz_lyr_col   , & ! Output: [real(r8) (:,:)] snow freezing rate (positive definite) (col,lyr) [kg m-2 s-1]
   qflx_snofrz_col  =>  waterflux_inst%qflx_snofrz_col       , & ! Output: [real(r8) (:)] column-integrated snow freezing rate (positive definite) [kg m-2 s-1]
   qflx_glcice      =>  waterflux_inst%qflx_glcice_col       , & ! Output: [real(r8) (:)] flux of new glacier ice (mm H2O/s) [+ = ice grows]
   qflx_glcice_melt =>  waterflux_inst%qflx_glcice_melt_col    & ! Output: [real(r8) (:)] ice melt (positive definite) (mm H2O/s)
   )

    ! Get step size

    dtime = get_step_size()

    ! Initialization

    do fc = 1,num_nolakec
       c = filter_nolakec(fc)
       l = col%landunit(c)

       qflx_snomelt(c) = 0._r8
       xmf(c) = 0._r8
       qflx_snofrz_lyr(c,-nlevsno+1:0) = 0._r8
       qflx_snofrz_col(c) = 0._r8
       qflx_snow_drain(c) = 0._r8
    end do

    do j = -nlevsno+1,nlevgrnd       ! all layers
       do fc = 1,num_nolakec
          c = filter_nolakec(fc)
          if (j >= snl(c)+1) then

             ! Initialization
             imelt(c,j) = 0
             hm(c,j) = 0._r8
             xm(c,j) = 0._r8
             wice0(c,j) = h2osoi_ice(c,j)
             wliq0(c,j) = h2osoi_liq(c,j)
             wmass0(c,j) = h2osoi_ice(c,j) + h2osoi_liq(c,j)
          endif   ! end of snow layer if-block
       end do   ! end of column-loop
    enddo   ! end of level-loop

!--  snow layers  ---------------------------------------------------
    do j = -nlevsno+1,0
       do fc = 1,num_nolakec
          c = filter_nolakec(fc)
          if (j >= snl(c)+1) then

             ! Melting identification
             ! If ice exists above melt point, melt some to liquid.
             if (h2osoi_ice(c,j) > 0._r8 .AND. t_soisno(c,j) > tfrz) then
                imelt(c,j) = 1
!                tinc(c,j) = t_soisno(c,j) - tfrz
                tinc(c,j) = tfrz - t_soisno(c,j)
                t_soisno(c,j) = tfrz
             endif

             ! Freezing identification
             ! If liquid exists below melt point, freeze some to ice.
             if (h2osoi_liq(c,j) > 0._r8 .AND. t_soisno(c,j) < tfrz) then
                imelt(c,j) = 2
!                tinc(c,j) = t_soisno(c,j) - tfrz
                tinc(c,j) = tfrz - t_soisno(c,j)
                t_soisno(c,j) = tfrz
             endif
          endif   ! end of snow layer if-block
       end do   ! end of column-loop
    enddo   ! end of level-loop

!-- soil layers   ---------------------------------------------------
    do j = 1,nlevgrnd
       do fc = 1,num_nolakec
          c = filter_nolakec(fc)
          l = col%landunit(c)
          supercool(c,j) = 0.0_r8

          if (h2osoi_ice(c,j) > 0. .AND. t_soisno(c,j) > tfrz) then
             imelt(c,j) = 1
!             tinc(c,j) = t_soisno(c,j) - tfrz
             tinc(c,j) = tfrz - t_soisno(c,j)
             t_soisno(c,j) = tfrz
          endif

          ! from Zhao (1997) and Koren (1999)
          supercool(c,j) = 0.0_r8
          if(t_soisno(c,j) < tfrz) then
             smp = hfus*(tfrz-t_soisno(c,j))/(grav*t_soisno(c,j)) * 1000._r8  !(mm)
             supercool(c,j) = watsat(c,j)*(smp/sucsat(c,j))**(-1._r8/bsw(c,j))
             supercool(c,j) = supercool(c,j)*dz(c,j)*1000._r8       ! (mm)
          endif

          if (h2osoi_liq(c,j) > supercool(c,j) .AND. t_soisno(c,j) < tfrz) then
             imelt(c,j) = 2
!             tinc(c,j) = t_soisno(c,j) - tfrz
             tinc(c,j) = tfrz - t_soisno(c,j)
             t_soisno(c,j) = tfrz
          endif

          ! If snow exists, but its thickness is less than the critical value (0.01 m)
          if (snl(c)+1 == 1 .AND. h2osno(c) > 0._r8 .AND. j == 1) then
             if (t_soisno(c,j) > tfrz) then
                imelt(c,j) = 1
!                tincc,j) = t_soisno(c,j) - tfrz
                tinc(c,j) = tfrz - t_soisno(c,j)
                t_soisno(c,j) = tfrz
             endif
          endif

       end do
    enddo


    do j = -nlevsno+1,nlevgrnd       ! all layers
       do fc = 1,num_nolakec
          c = filter_nolakec(fc)

             if (j >= snl(c)+1) then

                ! Calculate the energy surplus and loss for melting and freezing
                if (imelt(c,j) > 0) then

                   ! added unique cases for this calculation,
                   ! to account for absorbed solar radiation in each layer

                   !==================================================================
                   if (j == snl(c)+1) then ! top layer
                      hm(c,j) = dhsdT(c)*tinc(c,j) - tinc(c,j)/fact(c,j)

                      if ( j==1 .and. frac_h2osfc(c) /= 0.0_r8 ) then
                         hm(c,j) = hm(c,j) - frac_h2osfc(c)*(dhsdT(c)*tinc(c,j))
                      end if
                   else if (j == 1) then
                      hm(c,j) = (1.0_r8 - frac_sno_eff(c) - frac_h2osfc(c)) &
                           *dhsdT(c)*tinc(c,j) - tinc(c,j)/fact(c,j)
                   else ! non-interfacial snow/soil layers
                      hm(c,j) = - tinc(c,j)/fact(c,j)
                   endif
                endif

                ! These two errors were checked carefully (Y. Dai).  They result from the
                ! computed error of "Tridiagonal-Matrix" in subroutine "thermal".
                if (imelt(c,j) == 1 .AND. hm(c,j) < 0._r8) then
                   hm(c,j) = 0._r8
                   imelt(c,j) = 0
                endif
                if (imelt(c,j) == 2 .AND. hm(c,j) > 0._r8) then
                   hm(c,j) = 0._r8
                   imelt(c,j) = 0
                endif

                ! The rate of melting and freezing

                if (imelt(c,j) > 0 .and. abs(hm(c,j)) > 0._r8) then
                   xm(c,j) = hm(c,j)*dtime/hfus                           ! kg/m2

                   ! If snow exists, but its thickness is less than the critical value
                   ! (1 cm). Note: more work is needed to determine how to tune the
                   ! snow depth for this case
                   if (j == 1) then
                      if (snl(c)+1 == 1 .AND. h2osno(c) > 0._r8 .AND. xm(c,j) > 0._r8) then
                         temp1 = h2osno(c)                           ! kg/m2
                         h2osno(c) = max(0._r8,temp1-xm(c,j))
                         propor = h2osno(c)/temp1
                         snow_depth(c) = propor * snow_depth(c)
                         heatr = hm(c,j) - hfus*(temp1-h2osno(c))/dtime   ! W/m2
                         if (heatr > 0._r8) then
                            xm(c,j) = heatr*dtime/hfus                    ! kg/m2
                            hm(c,j) = heatr                               ! W/m2
                         else
                            xm(c,j) = 0._r8
                            hm(c,j) = 0._r8
                         endif
                         qflx_snomelt(c) = max(0._r8,(temp1-h2osno(c)))/dtime   ! kg/(m2 s)
                         xmf(c) = hfus*qflx_snomelt(c)
                         qflx_snow_drain(c) = qflx_snomelt(c)
                      endif
                   endif

                   heatr = 0._r8
                   if (xm(c,j) > 0._r8) then
                      h2osoi_ice(c,j) = max(0._r8, wice0(c,j)-xm(c,j))
                      heatr = hm(c,j) - hfus*(wice0(c,j)-h2osoi_ice(c,j))/dtime
                   else if (xm(c,j) < 0._r8) then
                      if (j <= 0) then
                         h2osoi_ice(c,j) = min(wmass0(c,j), wice0(c,j)-xm(c,j))  ! snow
                      else
                         if (wmass0(c,j) < supercool(c,j)) then
                            h2osoi_ice(c,j) = 0._r8
                         else
                            h2osoi_ice(c,j) = min(wmass0(c,j) - supercool(c,j),wice0(c,j)-xm(c,j))
                         endif
                      endif
                      heatr = hm(c,j) - hfus*(wice0(c,j)-h2osoi_ice(c,j))/dtime
                   endif

                   h2osoi_liq(c,j) = max(0._r8,wmass0(c,j)-h2osoi_ice(c,j))

                   if (abs(heatr) > 0._r8) then
                      if (j == snl(c)+1) then

                         if(j==1) then
                            t_soisno(c,j) = t_soisno(c,j) + fact(c,j)*heatr &
                                 /(1._r8-(1.0_r8 - frac_h2osfc(c))*fact(c,j)*dhsdT(c))
                         else
                            t_soisno(c,j) = t_soisno(c,j) + fact(c,j)*heatr &
                                 /(1._r8-fact(c,j)*dhsdT(c))
                         endif

                      else if (j == 1) then

                         t_soisno(c,j) = t_soisno(c,j) + fact(c,j)*heatr &
                              /(1._r8-(1.0_r8 - frac_sno_eff(c) - frac_h2osfc(c))*fact(c,j)*dhsdT(c))
                      else
                         t_soisno(c,j) = t_soisno(c,j) + fact(c,j)*heatr
                      endif

                      if (j <= 0) then    ! snow
                         if (h2osoi_liq(c,j)*h2osoi_ice(c,j)>0._r8) t_soisno(c,j) = tfrz
                      end if
                   endif  ! end of heatr > 0 if-block

                   if (j >= 1) then
                      xmf(c) = xmf(c) + hfus*(wice0(c,j)-h2osoi_ice(c,j))/dtime
                   else
                      xmf(c) = xmf(c) + frac_sno_eff(c)*hfus*(wice0(c,j)-h2osoi_ice(c,j))/dtime
                   endif

                   if (imelt(c,j) == 1 .AND. j < 1) then
                      qflx_snomelt(c) = qflx_snomelt(c) + max(0._r8,(wice0(c,j)-h2osoi_ice(c,j)))/dtime


                   endif

                   ! layer freezing mass flux (positive):
                   if (imelt(c,j) == 2 .AND. j < 1) then
                      qflx_snofrz_lyr(c,j) = max(0._r8,(h2osoi_ice(c,j)-wice0(c,j)))/dtime
                   endif

                endif

             endif   ! end of snow layer if-block

       end do   ! end of column-loop
    enddo   ! end of level-loop

    ! Needed for history file output

    do fc = 1,num_nolakec
       c = filter_nolakec(fc)
       eflx_snomelt(c) = qflx_snomelt(c) * hfus
       l = col%landunit(c)
    end do

    do j = -nlevsno+1,0
       do fc = 1,num_nolakec
          c = filter_nolakec(fc)
          qflx_snofrz_col(c) = qflx_snofrz_col(c) + qflx_snofrz_lyr(c,j)
       end do
    end do

    end associate
   end subroutine PhaseChange_beta

end module SoilTemperatureMod
