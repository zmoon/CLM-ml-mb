module SurfaceAlbedoType

  !-----------------------------------------------------------------------
  ! !DESCRIPTION:
  ! Surface albedo variables
  !
  ! !USES:
  use shr_kind_mod, only : r8 => shr_kind_r8
  use clm_varpar, only : numrad
  use clm_varcon, only : ispval, nan => spval
  use decompMod, only : bounds_type
  !
  ! !PUBLIC TYPES:
  implicit none
  save
  private
  !
  !PUBLIC DATA TYPES:

  type, public :: surfalb_type

    real(r8), pointer :: albgrd_col  (:,:)   ! Direct beam albedo of ground (soil) [for numrad wavebands]
    real(r8), pointer :: albgri_col  (:,:)   ! Diffuse albedo of ground (soil) [for numrad wavebands]

  contains

    procedure, public  :: Init
    procedure, private :: InitAllocate

  end type surfalb_type
  !-----------------------------------------------------------------------

contains

  !-----------------------------------------------------------------------
  subroutine Init (this, bounds)

    class(surfalb_type) :: this
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
    class(surfalb_type) :: this
    type(bounds_type), intent(in) :: bounds
    !
    ! !LOCAL VARIABLES:
    integer :: begc, endc   ! Column indices
    !---------------------------------------------------------------------

    begc = bounds%begc
    endc = bounds%endc

    allocate (this%albgrd_col (begc:endc,1:numrad))  ; this%albgrd_col  (:,:) = nan
    allocate (this%albgri_col (begc:endc,1:numrad))  ; this%albgri_col  (:,:) = nan

  end subroutine InitAllocate

end module SurfaceAlbedoType
