module FrictionVelocityMod

  !-----------------------------------------------------------------------
  ! !DESCRIPTION:
  ! Friction velocity variables
  !
  ! !USES:
  use shr_kind_mod, only : r8 => shr_kind_r8
  use clm_varcon, only : ispval, nan => spval
  use decompMod, only : bounds_type
  !
  ! !PUBLIC TYPES:
  implicit none
  save
  !
  !PUBLIC DATA TYPES:

  type, public :: frictionvel_type

    real(r8), pointer, public :: u10_clm_patch  (:)   ! 10-m wind (m/s)
    real(r8), pointer, public :: fv_patch       (:)   ! Friction velocity (m/s)
    real(r8), pointer, public :: ram1_patch     (:)   ! Aerodynamical resistance (s/m)

  contains

    procedure, public  :: Init
    procedure, private :: InitAllocate

  end type frictionvel_type
  !-----------------------------------------------------------------------

contains

  !-----------------------------------------------------------------------
  subroutine Init (this, bounds)

    class(frictionvel_type) :: this
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
    class(frictionvel_type) :: this
    type(bounds_type), intent(in) :: bounds
    !
    ! !LOCAL VARIABLES:
    integer :: begp, endp   ! Patch indices
    !---------------------------------------------------------------------

    begp = bounds%begp
    endp = bounds%endp

    allocate (this%u10_clm_patch (begp:endp))  ; this%u10_clm_patch (:) = nan
    allocate (this%fv_patch      (begp:endp))  ; this%fv_patch      (:) = nan
    allocate (this%ram1_patch    (begp:endp))  ; this%ram1_patch    (:) = nan

  end subroutine InitAllocate

end module FrictionVelocityMod
