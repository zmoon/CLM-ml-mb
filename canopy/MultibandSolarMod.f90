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

  ! Note: the paths will need to modified, `./` -> `../canopy/` to work in the real model

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
        read(10, *) wl(i), si_dr(i), si_df(i)
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

  !> Smear y(x) values over the single bin defined by [xgl, xgu]
  !> Using cumulative extrapolative trapezoidal integration, based on TUV's smear
  pure real(r8) function smear1(x, y, xgl, xgu) result(yg) 
    real(r8), dimension(:), intent(in) :: x, y
    real(r8), intent(in) :: xgl, xgu
    real(r8) :: area, a1, a2, slope, b1, b2
    integer :: k, n
    n = size(x)

    area = 0
    do k = 1, n - 1
      if ( x(k+1) < xgl .or. x(k) > xgu ) cycle  ! outside window, go to next
      a1 = max(x(k), xgl)
      a2 = min(x(k+1), xgu)
      slope = (y(k+1) - y(k)) / (x(k+1) - x(k))
      b1 = y(k) + slope * (a1 - x(k))
      b2 = y(k) + slope * (a2 - x(k))
      area = area + (a2 - a1) * (b2 + b1) / 2
    end do
    yg = area / (xgu - xgl)
  end function smear1

  !> Smear y(x) values to x-bins (bin edges!) using smear1
  pure function smear(x, y, bins) result(ynew)
    real(r8), dimension(:), intent(in) :: x, y
    real(r8), dimension(:), intent(in) :: bins
    real(r8), dimension(:), allocatable :: ynew
    integer :: i, n
    n = size(bins)  ! number of bins!

    allocate(ynew(n-1))
    do i = 1, n - 1
      ynew(i) = smear1(x, y, bins(i), bins(i+1))
    end do
  end function smear

  !> Trapezoidal integral of y(x)
  pure real(r8) function trapz(x, y) result(res)
    real(r8), dimension(:), intent(in) :: x, y
    integer :: i, n
    n = size(x)

    res = 0  ! don't really need loop to do this
    do i = 1, n - 1
      res = res + (y(i) + y(i+1))/2 * (x(i+1) - x(i))
    end do
  end function trapz

  ! main distribute rad subroutine
  ! 1. Load reference spectra
  ! 2. Smear spectra to n number of bands (or specific specified bands?)
  ! 3. Check integral (solar) or integrated average (optical props) for consistency with original 
  ! 4. Modify necessary variables in 
  !   - mlcanopy_inst (swskyb, swskyd, ...)
  !   - surfalb_inst (albgrd_col, albgri_col, ...)
  !   - clm_varpar (numrad, ...)
  !   Or instead create an input type that has all of the necessary things
  subroutine distribute_rad(  &
    wlbi, rli, tli, rsi, idri, idfi,  &  ! TODO: clean up by making some types?
    wle,  &
    wl, dwl, rl, tl, rs, idr, idf  &
  )
    real(r8), intent(in) :: wlbi(2)  ! wl bounds for the single band for the input values
    real(r8), intent(in) :: rli, tli, rsi, idri, idfi  ! input values (single band)
    real(r8), dimension(:), intent(in) :: wle  ! wavelength bounds for new grid
    real(r8), dimension(:), intent(out) :: wl, dwl
    real(r8), dimension(:), intent(out) :: rl, tl, rs, idr, idf  ! new values (spectral)

    ! Reference spectra (later could be an input)
    real(r8), dimension(nwl_leaf) :: wl0_leafsoil, rl0, tl0, rs0
    real(r8), dimension(nwl_solar) :: wl0_solar, sidr0, sidf0

    ! Local variables
    integer :: n, nbins
    real(r8), dimension(:), allocatable :: sidr, sidf

    ! New wavelength grid
    n = ubound(wle, dim=1)  ! number of bin edges (nbins + 1)
    if ( wle(1) /= wle(1) .or. wle(n) /= wlbi(2) ) stop 'wle edges should match orig band'
    dwl = wle(2:n) - wle(1:n-1)  ! new wavelength band widths
    wl = wl(1:n-1) + dwl  ! wavelength band centers (bin midpoints)
    nbins = n - 1
    allocate(sidr(n-1), sidf(n-1))

    ! Load and smear optical properties
    call load_leaf_spectrum(wl0_leafsoil, rl0, tl0)
    rl = smear(wl0_leafsoil, rl0, wle)
    tl = smear(wl0_leafsoil, tl0, wle)
    call load_soil_spectrum(wl0_leafsoil, rs0)
    rs = smear(wl0_leafsoil, rs0, wle)

    ! Correct optical property spectra based on the input values
    ! We want the bandwidth-weighted average to match
    rl = rl * rli / (sum(rl * dwl) / (wle(n) - wle(1)))
    tl = tl * tli / (sum(tl * dwl) / (wle(n) - wle(1)))
    rs = rs * rsi / (sum(rs * dwl) / (wle(n) - wle(1)))

    ! Load and smear irradiances
    call load_solar_spectrum(wl0_solar, sidr0, sidf0)
    sidr = smear(wl0_solar, sidr0, wle)
    sidf = smear(wl0_solar, sidf0, wle)
    idr = sidr * dwl  ! W m-2 um-1 -> W m-2
    idf = sidf * dwl

    ! Correct irradiances based on the input values
    ! We want the integral (sum of in-band irradiances) to match
    ! Since we have already multiplied by `dwl`, the integral is just the sum
    idr = idr * idri / sum(idr)
    idf = idf * idfi / sum(idf)

  end subroutine distribute_rad


  ! Integration, for application after CRT
  ! * reset modified variables in structs to their original values?
  ! * replace CRT outputs by the PAR and NIR integrated ones


end module MultibandSolarMod
