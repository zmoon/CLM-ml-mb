module MultibandSolarMod

  !-----------------------------------------------------------------------
  ! !DESCRIPTION:
  ! Distribute the two bands visible/PAR and NIR into multiple, using reference 
  ! leaf and solar spectra. After applying canopy radiative transfer, integrate 
  ! to obtain the needed PAR and NIR values.
  !
  ! !USES:
  use shr_kind_mod, only : r8 => shr_kind_r8
  use abortutils, only : endrun
  use decompMod, only : bounds_type
  use pftconMod, only : pftcon
  use PatchType, only : patch
  use SurfaceAlbedoType, only : surfalb_type
  use CanopyFluxesMultilayerType, only : mlcanopy_type
  !
  ! !PUBLIC TYPES:
  implicit none
  !
  ! !PUBLIC MEMBER FUNCTIONS:
  ! public :: SolarRadiation           ! Main driver for radiative transfer
  !
  ! !PRIVATE MEMBER FUNCTIONS:
  ! private :: NormanRadiation         ! Norman radiative transfer
  !-----------------------------------------------------------------------

contains

  ! load leaf fn


  ! load solar fn
  ! Planck as another option?


  ! load soil fn


  ! TUV smear fn


  ! main distribute rad subroutine
  ! 1. Load reference spectra
  ! 2. Smear spectra to n number of bands (or specific specified bands?)
  ! 3. Check integral (solar) or integrated average (optical props) for consistency with original 
  ! 4. Modify necessary variables in 
  !   - mlcanopy_inst (swskyb, swskyd, ...)
  !   - surfalb_inst (albgrd_col, albgri_col, ...)
  !   - clm_varpar (numrad, ...)
  !   Or instead create an input type that has all of the necessary things


  ! Integration, for application after CRT
  ! * reset modified variables in structs to their original values?
  ! * replace CRT outputs by the PAR and NIR integrated ones


end module MultibandSolarMod
