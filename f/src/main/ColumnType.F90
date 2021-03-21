module ColumnType

  !-----------------------------------------------------------------------
  ! !DESCRIPTION:
  ! Column data type
  !
  ! !USES:
  use shr_kind_mod, only : r8 => shr_kind_r8
  use clm_varpar, only : nlevgrnd, nlevsno
  use clm_varcon, only : ispval, nan => spval
  !
  ! !PUBLIC TYPES:
  implicit none
  save
  private
  !
  !PUBLIC DATA TYPES:
  type, public :: column_type

    integer,  pointer :: landunit (:)     ! Landunit of corresponding column for CLM g/l/c/p hierarchy
    integer,  pointer :: gridcell (:)     ! Gridcell of corresponding column for CLM g/l/c/p hierarchy
    integer,  pointer :: patchi   (:)     ! Beginning patches index for column
    integer,  pointer :: npatches (:)     ! Number of patches on column
    integer,  pointer :: snl      (:)     ! Number of snow layers
    real(r8), pointer :: dz       (:,:)   ! Soil layer thickness (m) [for -nlevsno+1:nlevgrnd layers]
    real(r8), pointer :: z        (:,:)   ! Soil layer depth (m) [for -nlevsno+1:nlevgrnd layers]
    real(r8), pointer :: zi       (:,:)   ! Soil layer depth at layer interface (m) [for -nlevsno:nlevgrnd layers]

  contains

    procedure, public  :: Init

  end type column_type
  type(column_type), public, target :: col
  !-----------------------------------------------------------------------

contains

  !-----------------------------------------------------------------------
  subroutine Init (this, begc, endc)
    !
    ! !DESCRIPTION:
    ! Initialize module data structure
    !
    ! !ARGUMENTS:
    class(column_type) :: this
    integer, intent(in) :: begc, endc     ! Column indices
    !---------------------------------------------------------------------

    allocate (this%landunit (begc:endc))                     ; this%landunit (:)   = ispval
    allocate (this%gridcell (begc:endc))                     ; this%gridcell (:)   = ispval
    allocate (this%patchi   (begc:endc))                     ; this%patchi   (:)   = ispval
    allocate (this%npatches (begc:endc))                     ; this%npatches (:)   = ispval
    allocate (this%snl      (begc:endc))                     ; this%snl      (:)   = ispval
    allocate (this%dz       (begc:endc,-nlevsno+1:nlevgrnd)) ; this%dz       (:,:) = nan
    allocate (this%z        (begc:endc,-nlevsno+1:nlevgrnd)) ; this%z        (:,:) = nan
    allocate (this%zi       (begc:endc,-nlevsno+0:nlevgrnd)) ; this%zi       (:,:) = nan

  end subroutine Init

end module ColumnType
