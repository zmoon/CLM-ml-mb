module CanopyStateType

  !-----------------------------------------------------------------------
  ! !DESCRIPTION:
  ! Canopy state variables
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
  type, public :: canopystate_type

     integer , pointer :: frac_veg_nosno_patch (:)   ! Fraction of vegetation not covered by snow (0 or 1)
     real(r8), pointer :: elai_patch           (:)   ! Canopy one-sided leaf area index with burying by snow
     real(r8), pointer :: esai_patch           (:)   ! Canopy one-sided stem area index with burying by snow

  contains

    procedure, public  :: Init
    procedure, private :: InitAllocate

  end type canopystate_type
  !-----------------------------------------------------------------------

contains

  !-----------------------------------------------------------------------
  subroutine Init (this, bounds)

    class(canopystate_type) :: this
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
    class(canopystate_type) :: this
    type(bounds_type), intent(in) :: bounds
    !
    ! !LOCAL VARIABLES:
    integer :: begp, endp   ! Patch indices
    !---------------------------------------------------------------------

    begp = bounds%begp
    endp = bounds%endp

    allocate (this%frac_veg_nosno_patch (begp:endp))  ; this%frac_veg_nosno_patch (:) = nan
    allocate (this%elai_patch           (begp:endp))  ; this%elai_patch           (:) = nan
    allocate (this%esai_patch           (begp:endp))  ; this%esai_patch           (:) = nan

  end subroutine InitAllocate

end module CanopyStateType
