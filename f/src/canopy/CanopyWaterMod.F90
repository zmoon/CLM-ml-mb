module CanopyWaterMod

  !-----------------------------------------------------------------------
  ! !DESCRIPTION:
  ! Update canopy water
  !
  ! !USES:
  use shr_kind_mod, only : r8 => shr_kind_r8
  use abortutils, only : endrun
  use CanopyFluxesMultilayerType, only : mlcanopy_type
  !
  ! !PUBLIC TYPES:
  implicit none
  save
  !
  ! !PUBLIC MEMBER FUNCTIONS:
  public :: CanopyInterception     ! Interception and throughfall
  public :: CanopyEvaporation      ! Update canopy intercepted water for evaporation
  !-----------------------------------------------------------------------

contains

  !-----------------------------------------------------------------------
  subroutine CanopyInterception (num_exposedvegp, filter_exposedvegp, mlcanopy_inst)
    !
    ! !DESCRIPTION:
    ! Interception and throughfall
    !
    ! !USES:
    use clm_varctl, only : dtime_sub, use_h2ocan
    !
    ! !ARGUMENTS:
    integer, intent(in) :: num_exposedvegp       ! Number of non-snow-covered veg points in CLM patch filter
    integer, intent(in) :: filter_exposedvegp(:) ! CLM patch filter for non-snow-covered vegetation
    type(mlcanopy_type), intent(inout) :: mlcanopy_inst
    !
    ! !LOCAL VARIABLES:
    integer  :: f                                ! Filter index
    integer  :: p                                ! Patch index for CLM g/l/c/p hierarchy
    integer  :: ic                               ! Aboveground layer index
    integer  :: n                                ! Number of leaf layers
    real(r8) :: dtime                            ! Model time step (s)
    real(r8) :: dewmx                            ! Maximum allowed interception (kg H2O/m2 leaf)
    real(r8) :: maximum_leaf_wetted_fraction     ! Maximum fraction of leaf that may be wet
    real(r8) :: fracrain                         ! Fraction of precipitation that is rain
    real(r8) :: fracsnow                         ! Fraction of precipitation that is snow
    real(r8) :: interception_fraction            ! Fraction of intercepted precipitation
    real(r8) :: fpi                              ! Fraction of precipitation intercepted
    real(r8) :: qflx_through_rain                ! Rain precipitation direct through canopy (kg H2O/m2/s)
    real(r8) :: qflx_through_snow                ! Snow precipitation direct through canopy (kg H2O/m2/s)
    real(r8) :: qflx_candrip                     ! Flux of water falling off canopy (kg H2O/m2/s)
    real(r8) :: h2ocanmx                         ! Maximum allowed water on canopy layer (kg H2O/m2)
    real(r8) :: xrun                             ! Excess water that exceeds the leaf capacity (kg H2O/m2/s)
    real(r8) :: qflx_prec_grnd_rain              ! Total rain throughfall onto ground (kg H2O/m2/s)
    real(r8) :: qflx_prec_grnd_snow              ! Total snow throughfall onto ground (kg H2O/m2/s)
    !---------------------------------------------------------------------

    associate ( &
                                                         ! *** Input ***
    lai            => mlcanopy_inst%lai             , &  ! Leaf area index of canopy (m2/m2)
    sai            => mlcanopy_inst%sai             , &  ! Stem area index of canopy (m2/m2)
    ncan           => mlcanopy_inst%ncan            , &  ! Number of aboveground layers
    dlai           => mlcanopy_inst%dlai            , &  ! Layer leaf area index (m2/m2)
    dpai           => mlcanopy_inst%dpai            , &  ! Layer plant area index (m2/m2)
    qflx_rain      => mlcanopy_inst%qflx_rain       , &  ! Rainfall (mm H2O/s = kg H2O/m2/s)
    qflx_snow      => mlcanopy_inst%qflx_snow       , &  ! Snowfall (mm H2O/s = kg H2O/m2/s)
                                                         ! *** Input/Output ***
    h2ocan         => mlcanopy_inst%h2ocan          , &  ! Canopy layer intercepted water (kg H2O/m2)
                                                         ! *** Output ***
    fwet           => mlcanopy_inst%fwet            , &  ! Fraction of plant area index that is wet
    fdry           => mlcanopy_inst%fdry            , &  ! Fraction of plant area index that is green and dry
    qflx_prec_intr => mlcanopy_inst%qflx_prec_intr    &  ! Intercepted precipitation (kg H2O/m2/s)
    )

    ! Time step

    dtime = dtime_sub

    ! Maximum allowed interception (kg H2O/m2 leaf)

    dewmx = 0.1_r8

    do f = 1, num_exposedvegp
       p = filter_exposedvegp(f)

       ! Fraction of precipitation that is rain and snow

       if ((qflx_snow(p) + qflx_rain(p)) > 0._r8) then
          fracrain = qflx_rain(p) / (qflx_snow(p) + qflx_rain(p))
          fracsnow = qflx_snow(p) / (qflx_snow(p) + qflx_rain(p))
       else
          fracrain = 0._r8
          fracsnow = 0._r8
       end if

       ! Fraction of precipitation that is intercepted

!      fpi = 0.25_r8 * (1._r8 - exp(-0.5_r8*(lai(p) + sai(p))))   ! CLM4.5
       interception_fraction = 1._r8
       fpi = interception_fraction * tanh(lai(p) + sai(p))        ! CLM5

       ! Direct throughfall

       qflx_through_rain = qflx_rain(p) * (1._r8 - fpi)
       qflx_through_snow = qflx_snow(p) * (1._r8 - fpi)

       ! Intercepted precipitation

       qflx_prec_intr(p) = (qflx_snow(p) + qflx_rain(p)) * fpi

       ! Find number of layers with lai+sai

       n = 0
       do ic = 1, ncan(p)
          if (dpai(p,ic) > 0._r8) n = n + 1
       end do

       ! Loop through layers for water balance calculation

       qflx_candrip = 0._r8
       do ic = 1, ncan(p)

          if (dpai(p,ic) > 0._r8) then ! leaf layer

             ! Maximum external water held in layer

             h2ocanmx = dewmx * dpai(p,ic)

             ! Water storage of intercepted precipitation. Intercepted water
             ! is applied equally to all layers.

             if (use_h2ocan) then
                h2ocan(p,ic) = h2ocan(p,ic) + qflx_prec_intr(p) * dtime / float(n)
             else
                h2ocan(p,ic) = 0._r8
             end if

             ! Excess water that exceeds the maximum capacity. If xrun > 0
             ! then h2ocan is set to h2ocanmx and excess water is added to
             ! throughfall.

             xrun = (h2ocan(p,ic) - h2ocanmx) / dtime
             if (xrun > 0._r8) then
                qflx_candrip = qflx_candrip + xrun
                h2ocan(p,ic) = h2ocanmx
             end if

             ! Wetted fraction of canopy

             maximum_leaf_wetted_fraction = 0.05_r8
             fwet(p,ic) = max((h2ocan(p,ic)/h2ocanmx),0._r8)**0.67_r8
             fwet(p,ic) = min (fwet(p,ic), maximum_leaf_wetted_fraction)

             ! Fraction of canopy that is green and dry

             fdry(p,ic) = (1._r8 - fwet(p,ic)) * (dlai(p,ic) / dpai(p,ic))

          else ! non-leaf layer

             h2ocan(p,ic) = 0._r8
             fwet(p,ic) = 0._r8
             fdry(p,ic) = 0._r8

          end if

       end do

       ! Total throughfall onto ground

       qflx_prec_grnd_rain = qflx_through_rain + qflx_candrip * fracrain
       qflx_prec_grnd_snow = qflx_through_snow + qflx_candrip * fracsnow

    end do

    end associate
  end subroutine CanopyInterception

  !-----------------------------------------------------------------------
  subroutine CanopyEvaporation (p, mlcanopy_inst)
    !
    ! !DESCRIPTION:
    ! Update canopy intercepted water for evaporation and dew
    !
    ! !USES:
    use clm_varctl, only : dtime_sub, use_h2ocan
    use clm_varpar, only : isun, isha
    use clm_varcon, only : mmh2o
    !
    ! !ARGUMENTS:
    integer, intent(in) :: p                        ! Patch index for CLM g/l/c/p hierarchy
    type(mlcanopy_type), intent(inout) :: mlcanopy_inst
    !
    ! !LOCAL VARIABLES:
    integer  :: ic                                  ! Aboveground layer index
    real(r8) :: dtime                               ! Model time step (s)
    real(r8) :: dew                                 ! Water (kg H2O/m2)
    !---------------------------------------------------------------------

    associate ( &
                                                    ! *** Input ***
    ncan           => mlcanopy_inst%ncan       , &  ! Number of aboveground layers
    dpai           => mlcanopy_inst%dpai       , &  ! Layer plant area index (m2/m2)
    trleaf         => mlcanopy_inst%trleaf     , &  ! Leaf transpiration flux (mol H2O/m2 leaf/s)
    evleaf         => mlcanopy_inst%evleaf     , &  ! Leaf evaporation flux (mol H2O/m2 leaf/s)
    fracsun        => mlcanopy_inst%fracsun    , &  ! Sunlit fraction of canopy layer
    fracsha        => mlcanopy_inst%fracsha    , &  ! Shaded fraction of canopy layer
                                                    ! *** Input/Output ***
    h2ocan         => mlcanopy_inst%h2ocan       &  ! Canopy layer intercepted water (kg H2O/m2)
    )

    dtime = dtime_sub

    do ic = 1, ncan(p)

       if (dpai(p,ic) > 0._r8) then

          ! Add dew, from both evaporation and transpiration

          dew = (evleaf(p,ic,isun) + trleaf(p,ic,isun)) * fracsun(p,ic) * dpai(p,ic) * mmh2o * dtime
          if (dew < 0._r8) then
             h2ocan(p,ic) = h2ocan(p,ic) - dew
          end if

          dew = (evleaf(p,ic,isha) + trleaf(p,ic,isha)) * fracsha(p,ic) * dpai(p,ic) * mmh2o * dtime
          if (dew < 0._r8) then
             h2ocan(p,ic) = h2ocan(p,ic) - dew
          end if

          ! Evaporate intercepted water

          if (evleaf(p,ic,isun) > 0._r8) then
             h2ocan(p,ic) = h2ocan(p,ic) - evleaf(p,ic,isun) * fracsun(p,ic) * dpai(p,ic) * mmh2o * dtime
          end if

          if (evleaf(p,ic,isha) > 0._r8) then
             h2ocan(p,ic) = h2ocan(p,ic) - evleaf(p,ic,isha) * fracsha(p,ic) * dpai(p,ic) * mmh2o * dtime
          end if

          ! Do not allow h2ocan to go negative

!         h2ocan(p,ic) = max (0._r8, h2ocan(p,ic))

          ! Exclude evaporation of intercepted water

          if (.not. use_h2ocan) then
             h2ocan(p,ic) = 0._r8
          end if

       end if

    end do

    end associate
  end subroutine CanopyEvaporation

end module CanopyWaterMod
