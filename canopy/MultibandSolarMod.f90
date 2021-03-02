module MultibandSolarMod

  !-----------------------------------------------------------------------
  ! !DESCRIPTION:
  ! Distribute the two bands visible/PAR and NIR into multiple, using reference 
  ! leaf and solar spectra. After applying canopy radiative transfer, integrate 
  ! to obtain the needed PAR and NIR values.
  !
  ! !USES:
  use shr_kind_mod, only : r8 => shr_kind_r8
  ! use abortutils, only : endrun
  ! use decompMod, only : bounds_type
  ! use pftconMod, only : pftcon
  ! use PatchType, only : patch
  ! use SurfaceAlbedoType, only : surfalb_type
  ! use CanopyFluxesMultilayerType, only : mlcanopy_type
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


  ! Number of wavelengths in each of the spectrum shape definitions
  integer, parameter :: nwl_leaf = 2101, nwl_solar = 122, nwl_soil = 2101


contains

  !> Load the sample leaf spectrum
  subroutine load_leaf_spectrum(wl, rl, tl)
    integer, parameter :: n = nwl_leaf  ! number of values/lines
    real(r8), dimension(n), intent(out) :: wl, rl, tl  ! wavelength, leaf reflectance, leaf transmittance
    character(len=*), parameter :: fp = "./PROSPECT_sample.txt"  ! file path
    integer :: i

    open(unit=10, file=fp)
      ! no header on this file
      do i = 1, n
        read(10, *) wl(i), rl(i), tl(i)
      end do
    close(10)
    wl = wl / 1000  ! nm -> um
  end subroutine load_leaf_spectrum

  !> Load the sample solar spectrum
  subroutine load_solar_spectrum(wl, si_dr, si_df)
    integer, parameter :: n = nwl_solar  ! number of values/lines
    real(r8), dimension(n), intent(out) :: wl, si_dr, si_df  ! wavelength, downwelling *spectral* direct irradiance, " " diffuse
    character(len=*), parameter :: fp = "./SPCTRAL2_xls_default-spectrum.csv"
    integer :: i

    open(unit=10, file=fp)
      read(10, *)  ! one header line
      do i = 1, n
        read(10, *) wl(i), si_dr(i), si_dr(i)
      end do
    close(10)
  end subroutine load_solar_spectrum

  ! TODO: Planck as another option for solar shape?

  !> Load the sample soil spectrum
  subroutine load_soil_spectrum(wl, rs)
    integer, parameter :: n = nwl_soil  ! number of values/lines
    real(r8), dimension(n), intent(out) :: wl, rs  ! wavelength, soil reflectivity
    real(r8), dimension(n) :: rs_dry, rs_wet
    character(len=*), parameter :: fp = "./PROSAIL_sample-soil.txt"  ! file path
    real(r8), parameter :: f_wet = 0
    integer :: i

    !> Wavelength isn't in the file, but it is 400:1:2500 in nm just like the PROSPECT one
    wl = [ (real(i, r8), i = 400, 2500, 1) ] / 1000

    open(unit=10, file=fp)
      ! no header on this file
      do i = 1, n
        read(10, *) rs_dry(i), rs_wet(i)
      end do
    close(10)

    !> One soil reflectivity from combining dry and wet
    rs = f_wet * rs_wet + (1 - f_wet) * rs_dry 
  end subroutine load_soil_spectrum


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
