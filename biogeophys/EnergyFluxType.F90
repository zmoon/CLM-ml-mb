module EnergyFluxType

  !-----------------------------------------------------------------------
  ! !DESCRIPTION:
  ! Energy flux variables
  !
  ! !USES:
  use shr_kind_mod, only : r8 => shr_kind_r8
  use clm_varpar, only : nlevgrnd, nlevsno
  use clm_varcon, only : ispval, nan => spval
  use decompMod, only : bounds_type
  !
  ! !PUBLIC TYPES:
  implicit none
  save
  private
  !
  !PUBLIC DATA TYPES:
  type, public :: energyflux_type

    ! Fluxes
    real(r8), pointer :: eflx_sh_grnd_patch   (:)   ! Sensible heat flux from ground (W/m2)
    real(r8), pointer :: eflx_sh_snow_patch   (:)   ! Sensible heat flux from snow (W/m2)
    real(r8), pointer :: eflx_sh_soil_patch   (:)   ! Sensible heat flux from soil (W/m2)
    real(r8), pointer :: eflx_sh_h2osfc_patch (:)   ! Sensible heat flux from surface water (W/m2)
    real(r8), pointer :: eflx_sh_tot_patch    (:)   ! Total sensible heat flux (W/m2)
    real(r8), pointer :: eflx_lh_tot_patch    (:)   ! Total latent heat flux (W/m2)
    real(r8), pointer :: eflx_soil_grnd_patch (:)   ! Soil heat flux (W/m2)
    real(r8), pointer :: eflx_lwrad_out_patch (:)   ! Emitted longwave radiation (W/m2)
    real(r8), pointer :: eflx_lwrad_net_patch (:)   ! Net longwave radiation (W/m2) [+ = to atm]
    real(r8), pointer :: eflx_gnet_patch      (:)   ! Net heat flux into ground (W/m2)
    real(r8), pointer :: eflx_snomelt_col     (:)   ! Snow melt heat flux (W/m2)
    real(r8), pointer :: eflx_bot_col         (:)   ! Heat flux from beneath the soil or ice column (W/m2)
    real(r8), pointer :: eflx_fgr12_col       (:)   ! Ground heat flux between soil layers 1 and 2 (W/m2)
    real(r8), pointer :: eflx_fgr_col         (:,:) ! (rural) soil downward heat flux (W/m2) [for nlevgrnd layers]

    ! Derivatives of energy fluxes
    real(r8), pointer :: dgnetdT_patch        (:)   ! Temperature derivative ground net heat flux (W/m2/K)
    real(r8), pointer :: cgrnd_patch          (:)   ! Deriv. of soil energy flux wrt to soil temp (W/m2/K)
    real(r8), pointer :: cgrndl_patch         (:)   ! Deriv. of soil evaporation flux wrt soil temp (kg/m2/s/K)
    real(r8), pointer :: cgrnds_patch         (:)   ! Deriv. of soil sensible heat flux wrt soil temp (W/m2/K)

    ! Canopy radiation
    real(r8), pointer :: dlrad_patch          (:)   ! Downward longwave radiation below the canopy (W/m2)
    real(r8), pointer :: ulrad_patch          (:)   ! Upward longwave radiation above the canopy (W/m2)

    ! Wind Stress
    real(r8), pointer :: taux_patch           (:)   ! Wind (shear) stress: east-west (kg/m/s2)
    real(r8), pointer :: tauy_patch           (:)   ! Wind (shear) stress: north-south (kg/m/s2)

    ! Transpiration
    real(r8), pointer :: btran_patch          (:)   ! Transpiration wetness factor (0 to 1)

    ! Latent heat
    real(r8), pointer :: htvp_col             (:)   ! Latent heat of vapor of water (or sublimation) (J/kg)

  contains

    procedure, public  :: Init
    procedure, private :: InitAllocate

  end type energyflux_type
  !-----------------------------------------------------------------------

contains

  !-----------------------------------------------------------------------
  subroutine Init (this, bounds)

    class(energyflux_type) :: this
    type(bounds_type), intent(in) :: bounds

    call this%InitAllocate (bounds)

  end subroutine Init

  !-----------------------------------------------------------------------
  subroutine InitAllocate (this, bounds)
    !
    ! !DESCRIPTION:
    ! Initialize module data structure
    !
    ! !ARGUMENTS:
    class(energyflux_type) :: this
    type(bounds_type), intent(in) :: bounds
    !
    ! !LOCAL VARIABLES:
    integer :: begp, endp   ! Patch indices
    integer :: begc, endc   ! Column indices
    !---------------------------------------------------------------------

    begp = bounds%begp ; endp = bounds%endp
    begc = bounds%begc ; endc = bounds%endc

    allocate (this%eflx_sh_grnd_patch   (begp:endp))            ; this%eflx_sh_grnd_patch   (:)   = nan
    allocate (this%eflx_sh_snow_patch   (begp:endp))            ; this%eflx_sh_snow_patch   (:)   = nan
    allocate (this%eflx_sh_soil_patch   (begp:endp))            ; this%eflx_sh_soil_patch   (:)   = nan
    allocate (this%eflx_sh_h2osfc_patch (begp:endp))            ; this%eflx_sh_h2osfc_patch (:)   = nan
    allocate (this%eflx_sh_tot_patch    (begp:endp))            ; this%eflx_sh_tot_patch    (:)   = nan
    allocate (this%eflx_lh_tot_patch    (begp:endp))            ; this%eflx_lh_tot_patch    (:)   = nan
    allocate (this%eflx_soil_grnd_patch (begp:endp))            ; this%eflx_soil_grnd_patch (:)   = nan
    allocate (this%eflx_lwrad_out_patch (begp:endp))            ; this%eflx_lwrad_out_patch (:)   = nan
    allocate (this%eflx_lwrad_net_patch (begp:endp))            ; this%eflx_lwrad_net_patch (:)   = nan
    allocate (this%eflx_gnet_patch      (begp:endp))            ; this%eflx_gnet_patch      (:)   = nan
    allocate (this%eflx_snomelt_col     (begc:endc))            ; this%eflx_snomelt_col     (:)   = nan
    allocate (this%eflx_bot_col         (begc:endc))            ; this%eflx_bot_col         (:)   = nan
    allocate (this%eflx_fgr12_col       (begc:endc))            ; this%eflx_fgr12_col       (:)   = nan
    allocate (this%eflx_fgr_col         (begc:endc,1:nlevgrnd)) ; this%eflx_fgr_col         (:,:) = nan
    allocate (this%dgnetdT_patch        (begp:endp))            ; this%dgnetdT_patch        (:)   = nan
    allocate (this%cgrnd_patch          (begp:endp))            ; this%cgrnd_patch          (:)   = nan
    allocate (this%cgrndl_patch         (begp:endp))            ; this%cgrndl_patch         (:)   = nan
    allocate (this%cgrnds_patch         (begp:endp))            ; this%cgrnds_patch         (:)   = nan
    allocate (this%dlrad_patch          (begp:endp))            ; this%dlrad_patch          (:)   = nan
    allocate (this%ulrad_patch          (begp:endp))            ; this%ulrad_patch          (:)   = nan
    allocate (this%taux_patch           (begp:endp))            ; this%taux_patch           (:)   = nan
    allocate (this%tauy_patch           (begp:endp))            ; this%tauy_patch           (:)   = nan
    allocate (this%btran_patch          (begp:endp))            ; this%btran_patch          (:)   = nan
    allocate (this%htvp_col             (begc:endc))            ; this%htvp_col             (:)   = nan

  end subroutine InitAllocate

end module EnergyFluxType
