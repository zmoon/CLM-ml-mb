module ncdio_pio

  !-----------------------------------------------------------------------
  ! !DESCRIPTION:
  ! Generic interfaces to write fields to netcdf files for CLM
  !
  ! !PUBLIC TYPES:
  implicit none
  private
  save
  !
  ! !PUBLIC MEMBER FUNCTIONS:
  !
  public :: ncd_defvar
  public :: ncd_io
  public :: ncd_inqvdlen

  integer, public :: ncd_double
  integer, public :: ncd_int

  type, public :: file_desc_t
     integer :: ncid
  end type file_desc_t
  !-----------------------------------------------------------------------

contains

  subroutine ncd_defvar
  end subroutine ncd_defvar

  subroutine ncd_io
  end subroutine ncd_io

  subroutine ncd_inqvdlen
  end subroutine ncd_inqvdlen

end module ncdio_pio
