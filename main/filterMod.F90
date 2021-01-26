module filterMod

  !-----------------------------------------------------------------------
  ! !DESCRIPTION:
  ! Module of filters used for processing CLM columns and patches
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

  type, public :: clumpfilter

    integer :: num_exposedvegp                ! Number of points with non-zero frac_veg_nosno in CLM patch filter
    integer, pointer :: exposedvegp(:)        ! CLM patch filter for non-zero frac_veg_nosno

    integer :: num_noexposedvegp              ! Number of points with frac_veg_nosno = 0 in CLM patch filter
    integer, pointer :: noexposedvegp(:)      ! CLM patch filter for frac_veg_nosno = 0

    integer :: num_nolakep                    ! Number of non-lake points in CLM patch filter
    integer, pointer :: nolakep(:)            ! CLM patch filter for non-lake points

    integer :: num_nourbanp                   ! Number of non-urban points in CLM patch filter
    integer, pointer :: nourbanp(:)           ! CLM patch filter for non-urban points

    integer :: num_nolakeurbanp               ! Number of patcges in non-lake, non-urban filter
    integer, pointer :: nolakeurbanp(:)       ! CLM patch filter for non-lake, non-urban points

    integer :: num_nolakec                    ! Number of non-lake points in CLM column filter
    integer, pointer :: nolakec(:)            ! CLM column filter for non-lake points

    integer :: num_nourbanc                   ! Number of non-urban points in CLM column filter
    integer, pointer :: nourbanc(:)           ! CLM column filter for non-urban points

  end type clumpfilter
  type(clumpfilter), public, target :: filter

  ! !PUBLIC MEMBER FUNCTIONS:
  public :: allocFilters                      ! Initialize data structure
  public :: setFilters                        ! Set filters
  public :: setExposedvegpFilter              ! Set the exposedvegp and noexposedvegp filters
  !-----------------------------------------------------------------------

contains

  !-----------------------------------------------------------------------
  subroutine allocFilters (filter, begp, endp, begc, endc)
    !
    ! !DESCRIPTION:
    ! Initialize module data structure
    !
    ! !ARGUMENTS:
    type(clumpfilter), intent(inout) :: filter
    integer, intent(in) :: begp, endp            ! Patch indices
    integer, intent(in) :: begc, endc            ! Column indices
    !---------------------------------------------------------------------

    allocate (filter%exposedvegp   (begp:endp))
    allocate (filter%noexposedvegp (begp:endp))
    allocate (filter%nolakep       (begp:endp))
    allocate (filter%nourbanp      (begp:endp))
    allocate (filter%nolakeurbanp  (begp:endp))
    allocate (filter%nolakec       (begc:endc))
    allocate (filter%nourbanc      (begc:endc))

  end subroutine allocFilters

  !-----------------------------------------------------------------------
  subroutine setFilters (filter)
    !
    ! !DESCRIPTION:
    ! Set CLM filters
    !
    ! !ARGUMENTS:
    type(clumpfilter), intent(inout) :: filter
    !---------------------------------------------------------------------

    filter%num_nolakep      = 1  ;  filter%nolakep(:)      = 1
    filter%num_nourbanp     = 1  ;  filter%nourbanp(:)     = 1
    filter%num_nolakeurbanp = 1  ;  filter%nolakeurbanp(:) = 1
    filter%num_nolakec      = 1  ;  filter%nolakec(:)      = 1
    filter%num_nourbanc     = 1  ;  filter%nourbanc(:)     = 1

  end subroutine setFilters

  !-----------------------------------------------------------------------
  subroutine setExposedvegpFilter (filter, frac_veg_nosno)
    !
    ! !DESCRIPTION:
    ! Set the exposedvegp and noexposedvegp patch filters.
    ! The exposedvegp filter includes patches for which frac_veg_nosno > 0.
    ! The noexposedvegp filter includes patches for which frac_veg_nosno = 0.
    ! However, neither filter includes urban or lake points.
    !
    ! !ARGUMENTS:
    type(clumpfilter), intent(inout) :: filter
    integer, intent(in) :: frac_veg_nosno(:)     ! Fraction of vegetation not covered by snow
    !
    ! !LOCAL VARIABLES:
    integer :: fp       ! Filter index
    integer :: p        ! Patch index
    integer :: fe, fn   ! Filter counts
    !---------------------------------------------------------------------

    fe = 0
    fn = 0
    do fp = 1, filter%num_nolakeurbanp
       p = filter%nolakeurbanp(fp)
       if (frac_veg_nosno(p) > 0) then
          fe = fe + 1
          filter%exposedvegp(fe) = p
       else
          fn = fn + 1
          filter%noexposedvegp(fn) = p
       end if
    end do
    filter%num_exposedvegp = fe
    filter%num_noexposedvegp = fn

  end subroutine setExposedvegpFilter

end module filterMod
