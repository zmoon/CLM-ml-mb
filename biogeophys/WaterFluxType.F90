module WaterFluxType

  !-----------------------------------------------------------------------
  ! !DESCRIPTION:
  ! Water flux variables
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

  type, public :: waterflux_type

    real(r8), pointer :: qflx_h2osfc_to_ice_col (:)   ! Conversion of h2osfc to ice
    real(r8), pointer :: qflx_snomelt_col       (:)   ! Snow melt (mm H2O /s)
    real(r8), pointer :: qflx_snow_drain_col    (:)   ! Drainage from snow pack
    real(r8), pointer :: qflx_snofrz_lyr_col    (:,:) ! Snow freezing rate (positive definite) (kg/m2/s) [for nlevsno layers]
    real(r8), pointer :: qflx_snofrz_col        (:)   ! Column-integrated snow freezing rate (positive definite) (kg/m2/s)
    real(r8), pointer :: qflx_glcice_col        (:)   ! Net flux of new glacial ice (growth - melt) (mm H2O/s), passed to GLC
    real(r8), pointer :: qflx_glcice_melt_col   (:)   ! Ice melt (positive definite) (mm H2O/s)
    real(r8), pointer :: qflx_ev_h2osfc_patch   (:)   ! Evaporation flux from h2osfc (kg/m2/s)
    real(r8), pointer :: qflx_ev_soil_patch     (:)   ! Evaporation flux from soil (kg/m2/s)
    real(r8), pointer :: qflx_ev_snow_patch     (:)   ! Evaporation flux from snow (kg/m2/s)
    real(r8), pointer :: qflx_evap_soi_patch    (:)   ! Soil evaporation (kg/m2/s)
    real(r8), pointer :: qflx_evap_tot_patch    (:)   ! qflx_evap_soi + qflx_evap_veg + qflx_tran_veg

  contains

    procedure, public  :: Init
    procedure, private :: InitAllocate

  end type waterflux_type
  !-----------------------------------------------------------------------

contains

  !-----------------------------------------------------------------------
  subroutine Init (this, bounds)

    class(waterflux_type) :: this
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
    class(waterflux_type) :: this
    type(bounds_type), intent(in) :: bounds
    !
    ! !LOCAL VARIABLES:
    integer :: begp, endp   ! Patch indices
    integer :: begc, endc   ! Column indices
    !---------------------------------------------------------------------

    begp = bounds%begp ; endp = bounds%endp
    begc = bounds%begc ; endc = bounds%endc

    allocate (this%qflx_h2osfc_to_ice_col (begc:endc))              ; this%qflx_h2osfc_to_ice_col (:)   = nan
    allocate (this%qflx_snomelt_col       (begc:endc))              ; this%qflx_snomelt_col       (:)   = nan
    allocate (this%qflx_snow_drain_col    (begc:endc))              ; this%qflx_snow_drain_col    (:)   = nan
    allocate (this%qflx_snofrz_lyr_col    (begc:endc,-nlevsno+1:0)) ; this%qflx_snofrz_lyr_col    (:,:) = nan
    allocate (this%qflx_snofrz_col        (begc:endc))              ; this%qflx_snofrz_col        (:)   = nan
    allocate (this%qflx_glcice_col        (begc:endc))              ; this%qflx_glcice_col        (:)   = nan
    allocate (this%qflx_glcice_melt_col   (begc:endc))              ; this%qflx_glcice_melt_col   (:)   = nan
    allocate (this%qflx_ev_h2osfc_patch   (begp:endp))              ; this%qflx_ev_h2osfc_patch   (:)   = nan
    allocate (this%qflx_ev_soil_patch     (begp:endp))              ; this%qflx_ev_soil_patch     (:)   = nan
    allocate (this%qflx_ev_snow_patch     (begp:endp))              ; this%qflx_ev_snow_patch     (:)   = nan
    allocate (this%qflx_evap_soi_patch    (begp:endp))              ; this%qflx_evap_soi_patch    (:)   = nan
    allocate (this%qflx_evap_tot_patch    (begp:endp))              ; this%qflx_evap_tot_patch    (:)   = nan

  end subroutine InitAllocate

end module WaterFluxType
