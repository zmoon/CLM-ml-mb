module PatchType

  !-----------------------------------------------------------------------
  ! !DESCRIPTION:
  ! Patch data type
  !
  ! !USES:
  use shr_kind_mod, only : r8 => shr_kind_r8
  use clm_varcon, only : ispval, nan => spval
  !
  ! !PUBLIC TYPES:
  implicit none
  save
  private
  !
  !PUBLIC DATA TYPES:
  type, public :: patch_type

    integer,  pointer :: column   (:)    ! Column of corresponding patch for CLM g/l/c/p hierarchy
    integer,  pointer :: landunit (:)    ! Landunit of corresponding patch for CLM g/l/c/p hierarchy
    integer,  pointer :: gridcell (:)    ! Gridcell of corresponding patch for CLM g/l/c/p hierarchy
    logical,  pointer :: active   (:)    ! True => do computations on this patch
    real(r8), pointer :: wtcol    (:)    ! Weight of patch relative to column
    integer , pointer :: itype    (:)    ! Vegetation type

  contains

    procedure, public  :: Init

  end type patch_type
  type(patch_type), public, target :: patch
  !-----------------------------------------------------------------------

contains

  !-----------------------------------------------------------------------
  subroutine Init (this, begp, endp)
    !
    ! !DESCRIPTION:
    ! Initialize module data structure
    !
    ! !ARGUMENTS:
    class(patch_type) :: this
    integer, intent(in) :: begp, endp    ! Patch indices
    !---------------------------------------------------------------------

    allocate (this%column   (begp:endp))   ; this%column   (:)   = ispval
    allocate (this%landunit (begp:endp))   ; this%landunit (:)   = ispval
    allocate (this%gridcell (begp:endp))   ; this%gridcell (:)   = ispval
    allocate (this%active   (begp:endp))   ; this%active   (:)   = .false.
    allocate (this%wtcol    (begp:endp))   ; this%wtcol    (:)   = nan
    allocate (this%itype    (begp:endp))   ; this%itype    (:)   = ispval

  end subroutine Init

end module PatchType
