module WaterStateType

  !-----------------------------------------------------------------------
  ! !DESCRIPTION:
  ! Water state variables
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

  type, public :: waterstate_type

    real(r8), pointer :: int_snow_col     (:)    ! Integrated snowfall (mm H2O)
    real(r8), pointer :: bw_col           (:,:)  ! Partial density of water in the snow pack (ice + liquid) (kg/m3)
    real(r8), pointer :: h2osno_col       (:)    ! Total snow water (kg/m2)
    real(r8), pointer :: h2osoi_liq_col   (:,:)  ! Soil layer liquid water (kg/m2) [for -nlevsno+1:nlevgrnd layers]
    real(r8), pointer :: h2osoi_ice_col   (:,:)  ! Soil layer ice lens (kg/m2) [for -nlevsno+1:nlevgrnd layers]
    real(r8), pointer :: h2osoi_vol_col   (:,:)  ! Soil layer volumetric water content (m3/m3) [for nlevgrnd layers]
    real(r8), pointer :: h2osfc_col       (:)    ! Surface water (kg/m2)
    real(r8), pointer :: snow_depth_col   (:)    ! Snow height (m)
    real(r8), pointer :: frac_sno_eff_col (:)    ! Effective fraction of ground covered by snow (0 to 1)
    real(r8), pointer :: frac_sno_col     (:)    ! Fraction of ground covered by snow (0 to 1)
    real(r8), pointer :: frac_h2osfc_col  (:)    ! Fraction of ground covered by surface water (0 to 1)
    real(r8), pointer :: q_ref2m_patch    (:)    ! 2-m height surface specific humidity (kg/kg)

  contains

    procedure, public  :: Init
    procedure, private :: InitAllocate

  end type waterstate_type
  !-----------------------------------------------------------------------

contains

  !-----------------------------------------------------------------------
  subroutine Init (this, bounds)

    class(waterstate_type) :: this
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
    class(waterstate_type) :: this
    type(bounds_type), intent(in) :: bounds
    !
    ! !LOCAL VARIABLES:
    integer :: begp, endp   ! Patch indices
    integer :: begc, endc   ! Column indices
    !---------------------------------------------------------------------

    begp = bounds%begp ; endp = bounds%endp
    begc = bounds%begc ; endc = bounds%endc

    allocate (this%int_snow_col     (begc:endc))                     ; this%int_snow_col     (:)   = nan
    allocate (this%bw_col           (begc:endc,-nlevsno+1:0))        ; this%bw_col           (:,:) = nan
    allocate (this%h2osno_col       (begc:endc))                     ; this%h2osno_col       (:)   = nan
    allocate (this%h2osoi_liq_col   (begc:endc,-nlevsno+1:nlevgrnd)) ; this%h2osoi_liq_col   (:,:) = nan
    allocate (this%h2osoi_ice_col   (begc:endc,-nlevsno+1:nlevgrnd)) ; this%h2osoi_ice_col   (:,:) = nan
    allocate (this%h2osoi_vol_col   (begc:endc,1:nlevgrnd))          ; this%h2osoi_vol_col   (:,:) = nan
    allocate (this%h2osfc_col       (begc:endc))                     ; this%h2osfc_col       (:)   = nan
    allocate (this%snow_depth_col   (begc:endc))                     ; this%snow_depth_col   (:)   = nan
    allocate (this%frac_sno_eff_col (begc:endc))                     ; this%frac_sno_eff_col (:)   = nan
    allocate (this%frac_sno_col     (begc:endc))                     ; this%frac_sno_col     (:)   = nan
    allocate (this%frac_h2osfc_col  (begc:endc))                     ; this%frac_h2osfc_col  (:)   = nan
    allocate (this%q_ref2m_patch    (begp:endp))                     ; this%q_ref2m_patch    (:)   = nan

  end subroutine InitAllocate

end module WaterStateType
