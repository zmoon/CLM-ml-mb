module TemperatureType

  !-----------------------------------------------------------------------
  ! !DESCRIPTION:
  ! Temperature variables
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

  type, public :: temperature_type

    real(r8), pointer :: t_h2osfc_col     (:)    ! Surface water temperature (K)
    real(r8), pointer :: t_h2osfc_bef_col (:)    ! Surface water temperature from time-step before (K)
    real(r8), pointer :: t_grnd_col       (:)    ! Ground surface temperature (K)
    real(r8), pointer :: t_soisno_col     (:,:)  ! Soil/snow temperature (K) [for -nlevsno+1:nlevgrnd layers]
    real(r8), pointer :: tssbef_col       (:,:)  ! Soil/snow temperature before update (K) [for -nlevsno+1:nlevgrnd layers]
    real(r8), pointer :: hc_soi_col       (:)    ! Soil heat content (MJ/m2)
    real(r8), pointer :: hc_soisno_col    (:)    ! Soil plus snow heat content (MJ/m2)
    integer , pointer :: imelt_col        (:,:)  ! Flag for melting (=1), freezing (=2), Not=0 [for -nlevsno+1:nlevgrnd layers]
    real(r8), pointer :: emg_col          (:)    ! Ground emissivity
    real(r8), pointer :: t_a10_patch      (:)    ! 10-day running mean of the 2-m temperature (K)
    real(r8), pointer :: t_ref2m_patch    (:)    ! 2-m height surface air temperature (K)

  contains

    procedure, public  :: Init
    procedure, private :: InitAllocate

  end type temperature_type
  !-----------------------------------------------------------------------

contains

  !-----------------------------------------------------------------------
  subroutine Init (this, bounds)

    class(temperature_type) :: this
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
    class(temperature_type) :: this
    type(bounds_type), intent(in) :: bounds
    !
    ! !LOCAL VARIABLES:
    integer :: begp, endp   ! Patch indices
    integer :: begc, endc   ! Column indices
    !---------------------------------------------------------------------

    begp = bounds%begp ; endp = bounds%endp
    begc = bounds%begc ; endc = bounds%endc

    allocate (this%t_h2osfc_col     (begc:endc))                     ; this%t_h2osfc_col     (:)   = nan
    allocate (this%t_h2osfc_bef_col (begc:endc))                     ; this%t_h2osfc_bef_col (:)   = nan
    allocate (this%t_grnd_col       (begc:endc))                     ; this%t_grnd_col       (:)   = nan
    allocate (this%t_soisno_col     (begc:endc,-nlevsno+1:nlevgrnd)) ; this%t_soisno_col     (:,:) = nan
    allocate (this%tssbef_col       (begc:endc,-nlevsno+1:nlevgrnd)) ; this%tssbef_col       (:,:) = nan
    allocate (this%hc_soi_col       (begc:endc))                     ; this%hc_soi_col       (:)   = nan
    allocate (this%hc_soisno_col    (begc:endc))                     ; this%hc_soisno_col    (:)   = nan
    allocate (this%imelt_col        (begc:endc,-nlevsno+1:nlevgrnd)) ; this%imelt_col        (:,:) = ispval
    allocate (this%emg_col          (begc:endc))                     ; this%emg_col          (:)   = nan
    allocate (this%t_a10_patch      (begp:endp))                     ; this%t_a10_patch      (:)   = nan
    allocate (this%t_ref2m_patch    (begp:endp))                     ; this%t_ref2m_patch    (:)   = nan

  end subroutine InitAllocate

end module TemperatureType
