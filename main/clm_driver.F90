module clm_driver

  !-----------------------------------------------------------------------
  ! !DESCRIPTION:
  ! Main CLM model driver to calculate fluxes
  !
  ! !USES:
  use shr_kind_mod, only : r8 => shr_kind_r8
  use abortutils, only : endrun
  use decompMod, only : bounds_type
  !
  ! !PUBLIC TYPES:
  implicit none
  !
  ! !PUBLIC MEMBER FUNCTIONS:
  public :: clm_drv
  !-----------------------------------------------------------------------

contains

  !-----------------------------------------------------------------------
  subroutine clm_drv (bounds, istep, time_indx, it, fin2)
    !
    ! !DESCRIPTION:
    ! Main CLM model driver to calculate fluxes
    !
    ! !USES:
    use clm_varctl, only : iulog
    use clm_varpar, only : nlevgrnd, nlevsoi, nlevsno, isun, isha
    use clm_varcon, only : mmh2o, mmdry
    use TowerDataMod, only : tower_id
    use PatchType, only : patch
    use ColumnType, only : col
    use pftconMod, only : pftcon
    use filterMod, only : filter, setExposedvegpFilter
    use clm_instMod
    use clmDataMod, only : clmData
    use SurfaceAlbedoMod, only : SoilAlbedo
    use SurfaceResistanceMod, only : calc_soilevap_resis
    use SoilTemperatureMod, only : SoilTemperature, SoilTemperatureInit, SoilTemperatureAfter
    use CanopyFluxesMultilayerMod, only : CanopyFluxesMultilayer
    use CanopyTurbulenceMod, only : LookupPsihatINI
    !
    ! !ARGUMENTS:
    implicit none
    type(bounds_type), intent(in) :: bounds     ! CLM bounds
    integer, intent(in) :: istep                ! Time stepping loop index
    integer, intent(in) :: time_indx            ! Time index from reference date (0Z January 1 of current year, when calday = 1.000
    integer, intent(in) :: it                   ! Tower site index
    character(len=256) :: fin2                  ! File name
    !
    ! !LOCAL VARIABLES:
    integer  :: ic                               ! Aboveground layer index
    integer  :: j                                ! Soil layer index
    integer  :: f                                ! Filter index
    integer  :: p                                ! Patch index for CLM g/l/c/p hierarchy
    integer  :: c                                ! Column index for CLM g/l/c/p hierarchy
    integer  :: g                                ! Grid cell index for CLM g/l/c/p hierarchy

    real(r8) :: totpai                           ! Canopy lai+sai for error check (m2/m2)
    real(r8) :: lai_bef(bounds%begp:bounds%endp) ! Leaf area index of canopy, from previous timestep (m2/m2)
    real(r8) :: sai_bef(bounds%begp:bounds%endp) ! Stem area index of canopy, from previous timestep (m2/m2)

    ! For CLM soil temperature
    real(r8) :: xmf(bounds%begc:bounds%endc)                       ! Dummy output
    real(r8) :: fact(bounds%begc:bounds%endc, -nlevsno+1:nlevgrnd) ! Dummy output
    real(r8) :: c_h2osfc(bounds%begc:bounds%endc)                  ! Dummy output
    real(r8) :: xmf_h2osfc(bounds%begc:bounds%endc)                ! Dummy output
    !---------------------------------------------------------------------

    associate ( &
    ncan           => mlcanopy_inst%ncan                     , &  ! Number of layers
    lai            => mlcanopy_inst%lai                      , &  ! Leaf area index of canopy (m2/m2)
    sai            => mlcanopy_inst%sai                      , &  ! Stem area index of canopy (m2/m2)
    dlai           => mlcanopy_inst%dlai                     , &  ! Layer leaf area index (m2/m2)
    dsai           => mlcanopy_inst%dsai                     , &  ! Layer stem area index (m2/m2)
    dpai           => mlcanopy_inst%dpai                     , &  ! Layer plant area index (m2/m2)
!   swskyb         => mlcanopy_inst%swskyb                   , &  ! Atmospheric direct beam solar radiation (W/m2)
!   swskyd         => mlcanopy_inst%swskyd                   , &  ! Atmospheric diffuse solar radiation (W/m2)
    wind           => mlcanopy_inst%wind                     , &  ! Wind speed profile (m/s)
    tair           => mlcanopy_inst%tair                     , &  ! Air temperature profile (K)
    eair           => mlcanopy_inst%eair                     , &  ! Vapor pressure profile (Pa)
    cair           => mlcanopy_inst%cair                     , &  ! Atmospheric CO2 profile (umol/mol)
    tveg           => mlcanopy_inst%tveg                     , &  ! Vegetation temperature profile (K)
    tleaf          => mlcanopy_inst%tleaf                    , &  ! Leaf temperature (K)
    h2ocan         => mlcanopy_inst%h2ocan                   , &  ! Canopy layer intercepted water (kg H2O/m2)
    taf            => mlcanopy_inst%taf                      , &  ! Air temperature at canopy top (K)
    eaf            => mlcanopy_inst%eaf                      , &  ! Vapor pressure at canopy top (Pa)
    qaf            => mlcanopy_inst%qaf                      , &  ! Specific humidity at canopy top (kg/kg)
    ga_prof        => mlcanopy_inst%ga_prof                  , &  ! Canopy layer aerodynamic conductance for scalars (mol/m2/s)
    tg             => mlcanopy_inst%tg                       , &  ! Soil surface temperature (K)
    zref_old       => mlcanopy_inst%zref_old                 , &  ! Reference height for previous timestep (m)
    zi             => col%zi                                 , &  ! CLM: Soil layer depth at layer interface (m)
    roota_par      => pftcon%roota_par                       , &  ! CLM: Rooting distribution parameter (1/m)
    rootb_par      => pftcon%rootb_par                       , &  ! CLM: Rooting distribution parameter (1/m)
    forc_u         => atm2lnd_inst%forc_u_grc                , &  ! CLM: Atmospheric wind speed in east direction (m/s)
    forc_v         => atm2lnd_inst%forc_v_grc                , &  ! CLM: Atmospheric wind speed in north direction (m/s)
    forc_pco2      => atm2lnd_inst%forc_pco2_grc             , &  ! CLM: Atmospheric CO2 partial pressure (Pa)
    forc_t         => atm2lnd_inst%forc_t_downscaled_col     , &  ! CLM: Atmospheric temperature (K)
    forc_q         => atm2lnd_inst%forc_q_downscaled_col     , &  ! CLM: Atmospheric specific humidity (kg/kg)
    forc_pbot      => atm2lnd_inst%forc_pbot_downscaled_col  , &  ! CLM: Atmospheric pressure (Pa)
    forc_hgt       => atm2lnd_inst%forc_hgt_grc              , &  ! CLM: Atmospheric reference height (m)
!   forc_solad     => atm2lnd_inst%forc_solad_grc            , &  ! CLM: Atmospheric direct beam radiation (W/m2)
!   forc_solai     => atm2lnd_inst%forc_solai_grc            , &  ! CLM: Atmospheric diffuse radiation (W/m2)
    rootfr         => soilstate_inst%rootfr_patch            , &  ! CLM: Fraction of roots in each soil layer
    frac_veg_nosno => canopystate_inst%frac_veg_nosno_patch  , &  ! CLM: Fraction of vegetation not covered by snow (0 or 1)
    elai           => canopystate_inst%elai_patch            , &  ! CLM: Leaf area index of canopy (m2/m2)
    esai           => canopystate_inst%esai_patch              &  ! CLM: Stem area index of canopy (m2/m2)
    )

    ! Read CLM data for current time slice

    call clmData (time_indx, fin2, bounds%begp, bounds%endp, bounds%begc, bounds%endc, &
    soilstate_inst, waterstate_inst, energyflux_inst, canopystate_inst)

    ! Set CLM frac_veg_nosno

    do p = bounds%begp, bounds%endp
       frac_veg_nosno(p) = 1
    end do

    ! Now set CLM frac_veg_nosno filter (filter%exposedvegp)

    call setExposedvegpFilter (filter, frac_veg_nosno)

    ! Initialize scalar profiles and leaf and ground temperatures

    if (istep == 1) then
       do f = 1, filter%num_exposedvegp
          p = filter%exposedvegp(f)
          c = patch%column(p)
          g = patch%gridcell(p)
          do ic = 0, ncan(p)
             wind(p,ic) = max (1.0_r8, sqrt(forc_u(g)*forc_u(g)+forc_v(g)*forc_v(g)))
             tair(p,ic) = forc_t(c)
             eair(p,ic) = forc_q(c) * forc_pbot(c) / (mmh2o/mmdry + (1._r8-mmh2o/mmdry) * forc_q(c))
             cair(p,ic) = forc_pco2(g) / forc_pbot(c) * 1.e06_r8
             tveg(p,ic,isun) = forc_t(c)
             tveg(p,ic,isha) = forc_t(c)
          end do
          do ic = 1, ncan(p)
             tleaf(p,ic,isun) = forc_t(c)
             tleaf(p,ic,isha) = forc_t(c)
             h2ocan(p,ic) = 0._r8
          end do
          taf(p) = forc_t(c)
          eaf(p) = forc_q(c) * forc_pbot(c) / (mmh2o/mmdry + (1._r8-mmh2o/mmdry) * forc_q(c))
          qaf(p) = forc_q(c)
          tg(p) = tair(p,0)
          ga_prof(p,0) = 1._r8 / 10._r8 * 42.3_r8
          zref_old(p) = forc_hgt(g)
       end do

       ! Initialize the look-up tables needed to calculate the RSL psihat functions

       call LookupPsihatINI

    end if

    ! Update leaf and stem area for current values (CLM: elai, esai)

    do f = 1, filter%num_exposedvegp
       p = filter%exposedvegp(f)

       ! Save values for previous time step

       lai_bef(p) = lai(p)
       sai_bef(p) = sai(p)

       ! Values for current time step

       lai(p) = elai(p)
       sai(p) = esai(p)

       ! Update profile for current value

       totpai = 0._r8
       do ic = 1, ncan(p)
          if (lai_bef(p) > 0._r8) dlai(p,ic) = dlai(p,ic) * lai(p) / lai_bef(p)
          if (sai_bef(p) > 0._r8) dsai(p,ic) = dsai(p,ic) * sai(p) / sai_bef(p)
          dpai(p,ic) = dlai(p,ic) + dsai(p,ic)
          totpai = totpai + dpai(p,ic)
       end do

       if (abs(totpai - (lai(p)+sai(p))) > 1.e-06_r8) then
          write (iulog,*) 'clm_drv error: Plant area index not updated correctly'
          write (6,*) lai(p)+sai(p), totpai
          call endrun()
       end if

    end do

    ! CLM root fraction

    do f = 1, filter%num_exposedvegp
       p = filter%exposedvegp(f)
       c = patch%column(p)

       rootfr(p,:) = 0._r8
       do j = 1, nlevsoi-1
          rootfr(p,j) = 0.5_r8 * &
          (exp(-roota_par(patch%itype(p))*zi(c,j-1)) + exp(-rootb_par(patch%itype(p))*zi(c,j-1)) - &
           exp(-roota_par(patch%itype(p))*zi(c,j  )) - exp(-rootb_par(patch%itype(p))*zi(c,j  )))
       end do

       j = nlevsoi
       rootfr(p,j) = 0.5_r8 * (exp(-roota_par(patch%itype(p))*zi(c,j-1)) + &
       exp(-rootb_par(patch%itype(p))*zi(c,j-1)))

    end do

    ! Soil albedo - CLM processes non-urban columns

    call SoilAlbedo (bounds, filter%num_nourbanc, filter%nourbanc, waterstate_inst, surfalb_inst)

    ! Solar radiation transfer through the canopy - CLM processes non-urban patches

!   do f = 1, filter%num_exposedvegp
!      p = filter%exposedvegp(f)
!      g = patch%gridcell(p)
!      swskyb(p,ivis) = forc_solad(g,ivis)
!      swskyb(p,inir) = forc_solad(g,inir)
!      swskyd(p,ivis) = forc_solai(g,ivis)
!      swskyd(p,inir) = forc_solai(g,inir)
!   end do

    ! Moisture stress/resistance for soil evaporation - CLM processes non-lake columns

    call calc_soilevap_resis (bounds, filter%num_nolakec, filter%nolakec, &
    soilstate_inst, waterstate_inst, temperature_inst)

    ! Need to initialize some variables to calculate thermal conductivity

    call SoilTemperatureInit (filter%num_nolakec, filter%nolakec, soilstate_inst, &
    waterstate_inst, energyflux_inst, temperature_inst, solarabs_inst, &
    atm2lnd_inst, mlcanopy_inst)

    ! Canopy and soil fluxes

    call CanopyFluxesMultilayer (bounds, filter%num_exposedvegp, filter%exposedvegp, &
    filter%num_nolakec, filter%nolakec, atm2lnd_inst, temperature_inst, waterstate_inst, &
    waterflux_inst, energyflux_inst, frictionvel_inst, soilstate_inst, surfalb_inst, &
    mlcanopy_inst)

    ! Update soil temperature

    call SoilTemperatureInit (filter%num_nolakec, filter%nolakec, soilstate_inst, &
    waterstate_inst, energyflux_inst, temperature_inst, solarabs_inst, &
    atm2lnd_inst, mlcanopy_inst)

    call SoilTemperature (bounds, filter%num_nolakec, filter%nolakec, &
    filter%num_nolakep, filter%nolakep, xmf(bounds%begc:bounds%endc), &
    fact(bounds%begc:bounds%endc, :), c_h2osfc(bounds%begc:bounds%endc), &
    xmf_h2osfc(bounds%begc:bounds%endc), energyflux_inst, waterflux_inst, waterstate_inst, &
    temperature_inst, soilstate_inst, atm2lnd_inst, solarabs_inst, canopystate_inst)

!   call SoilTemperatureAfter (filter%num_nolakep, filter%nolakep, energyflux_inst, &
!   waterflux_inst, mlcanopy_inst)

    end associate
  end subroutine clm_drv

end module clm_driver
