module clm_varpar

  !-----------------------------------------------------------------------
  ! !DESCRIPTION:
  ! Module containing various model parameters
  !-----------------------------------------------------------------------

  use shr_kind_mod, only : r8 => shr_kind_r8
  implicit none

  !-----------------------------------------------------------------------
  ! Existing CLM parameters: set here to specific values to test canopy code
  !-----------------------------------------------------------------------

  integer, parameter :: numrad = 2            ! Number of radiation wavebands
  integer, parameter :: nlevcan = 100         ! Number of layers in canopy model
  integer, parameter :: nlevsoi = 10          ! Number of hydrologically active soil layers
  integer, parameter :: nlevgrnd = 15         ! Number of ground layers (including hydrologically inactive)
  integer, parameter :: nlevsno  = 5          ! Maximum number of snow layer
  integer, parameter :: mxpft = 16            ! Maximum number of plant functional types
  integer, parameter :: ivis = 1              ! Visible waveband index
  integer, parameter :: inir = 2              ! Near-infrared waveband index
  integer, parameter :: max_patch_per_col = 1 ! Maximum number of patches per column

  !-----------------------------------------------------------------------
  ! New parameters for multilayer canopy
  !-----------------------------------------------------------------------

  integer, parameter :: nleaf = 2             ! Number of leaf types (sunlit and shaded)
  integer, parameter :: isun = 1              ! Sunlit leaf index
  integer, parameter :: isha = 2              ! Shaded leaf index

end module clm_varpar
