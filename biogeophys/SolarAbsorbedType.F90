module SolarAbsorbedType

  !-----------------------------------------------------------------------
  ! !DESCRIPTION:
  ! Solar absorbed variables
  !
  ! !USES:
  use shr_kind_mod, only : r8 => shr_kind_r8
  use clm_varpar, only : nlevsno
  use clm_varcon, only : ispval, nan => spval
  use decompMod, only : bounds_type
  !
  ! !PUBLIC TYPES:
  implicit none
  save
  private
  !
  !PUBLIC DATA TYPES:

  type, public :: solarabs_type

    real(r8), pointer :: sabg_patch      (:)    ! Solar radiation absorbed by ground (W/m2)
    real(r8), pointer :: sabg_snow_patch (:)    ! Solar radiation absorbed by snow (W/m2)
    real(r8), pointer :: sabg_soil_patch (:)    ! Solar radiation absorbed by soil (W/m2)
    real(r8), pointer :: sabg_chk_patch  (:)    ! fsno weighted sum (needed by balancecheck, because fsno changes midway)
    real(r8), pointer :: sabg_lyr_patch  (:,:)  ! Absorbed radiation in each snow layer and top soil layer (W/m2)

  contains

    procedure, public  :: Init
    procedure, private :: InitAllocate

  end type solarabs_type
  !-----------------------------------------------------------------------

contains

  !-----------------------------------------------------------------------
  subroutine Init (this, bounds)

    class(solarabs_type) :: this
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
    class(solarabs_type) :: this
    type(bounds_type), intent(in) :: bounds
    !
    ! !LOCAL VARIABLES:
    integer :: begp, endp   ! Patch indices
    !---------------------------------------------------------------------

    begp = bounds%begp
    endp = bounds%endp

    allocate (this%sabg_patch      (begp:endp))              ; this%sabg_patch      (:)   = nan
    allocate (this%sabg_snow_patch (begp:endp))              ; this%sabg_snow_patch (:)   = nan
    allocate (this%sabg_soil_patch (begp:endp))              ; this%sabg_soil_patch (:)   = nan
    allocate (this%sabg_chk_patch  (begp:endp))              ; this%sabg_chk_patch  (:)   = nan
    allocate (this%sabg_lyr_patch  (begp:endp,-nlevsno+1:1)) ; this%sabg_lyr_patch  (:,:) = nan

  end subroutine InitAllocate

end module SolarAbsorbedType
