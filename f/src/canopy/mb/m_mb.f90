!> This module provides routines for re-binning spectra
!> and creating a set of sub-band values from a single broadband value for
!> * leaf reflectance/transmittance
!> * soil reflectance (albedo)
!> * direct/diffuse irradiance
module m_mb
  use m_mb_data, only: rk
  implicit none

  ! Number of wavelengths in each of the spectrum shape definitions
  ! Only needed when loading from file
  integer, parameter :: nwl_leaf = 2101, nwl_solar = 122, nwl_soil = 2101

  ! Relative path to the location of the data files (so it matters where the code is run from)
  character(len=*), private, parameter :: &
    refspecbasepath = "../canopy/mb/"  ! when running CLM-ml or running the Ninja test

  ! Wavelength bounds for some common bands (um)
  real(rk), dimension(2), parameter :: &
    wlb_par = [0.4_rk, 0.7_rk], &
    wlb_nir = [0.7_rk, 2.5_rk]  ! no use going further in wavelength since don't have leaf/soil data

  ! Constants for Planck
  real(rk), private, parameter :: &
    h = 6.62607e-34_rk, &  ! Planck constant (J s)
    c = 2.99792e8_rk, &  ! speed of light (m s-1)
    k_B = 1.38065e-23_rk  ! Boltzmann constant (J K-1)


contains


  !> Print constants to check precision/accuracy
  subroutine print_constants()
    print *, 'h, c, k_B', h, c, k_B
  end subroutine print_constants


  !> Compute Planck radiance L(T, wl)
  pure elemental real(rk) function l_wl_planck(T_K, wl_um) result(l)
    real(rk), intent(in) :: T_K, wl_um
    real(rk) :: wl
    wl = wl_um * 1.e-6_rk  ! -> m
    l = (2 * h * c**2) / ( &
      wl**5 * (exp(h * c / (wl * k_B * T_K)) - 1) &
    )
  end function l_wl_planck


  !> Estimate the definite integral of Planck radiance L(T, wl) over a region
  !> Currently using trapz but might be better to define and use a simps routine
  pure elemental real(rk) function l_wl_plank_integ(T_K, wla_um, wlb_um) result(res)
    real(rk), intent(in) :: T_K, wla_um, wlb_um
    integer, parameter :: n = 200  ! discrete points in the trapz
    real(rk) :: wl_um(n), l(n)
    integer :: i
    do i = 1, n  ! linspace
      wl_um(i) = wla_um + (wlb_um - wla_um) * (i - 1) / (n - 1)
    end do
    l = l_wl_planck(T_K, wl_um)
    res = trapz(wl_um * 1.e-6_rk, l)
  end function l_wl_plank_integ


  !> Generate equal-width wavelength sub-bands for a certain spectral region
  !> Returns band edges, so the result will have size *n+1, not n*
  pure function wle_equal_width(wlb, n) result(wle)
    real(rk), intent(in) :: wlb(2)  ! wavelength bounds of the region
    integer, intent(in) :: n  ! desired number of bands
    real(rk) :: wle(n+1)

    integer :: i
    real(rk) :: dwl
    dwl = (wlb(2) - wlb(1)) / n

    wle(1) = wlb(1)
    do i = 2, n
      wle(i) = wle(i-1) + dwl
    end do
    wle(n+1) = wlb(2)
  end function wle_equal_width


  !> Load the sample leaf spectrum
  subroutine load_leaf_spectrum(wl, rl, tl)
    integer, parameter :: n = nwl_leaf  ! number of values/lines
    real(rk), dimension(n), intent(out) :: wl, rl, tl  ! wavelength, leaf reflectance, leaf transmittance
    character(len=*), parameter :: fp = refspecbasepath // "PROSPECT_sample.txt"  ! file path
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
    real(rk), dimension(n), intent(out) :: wl, si_dr, si_df  ! wavelength, downwelling *spectral* direct irradiance, " " diffuse
    character(len=*), parameter :: fp = refspecbasepath // "SPCTRAL2_xls_default-spectrum.csv"
    integer :: i

    open(unit=10, file=fp)
      read(10, *)  ! one header line
      do i = 1, n
        read(10, *) wl(i), si_dr(i), si_df(i)
      end do
    close(10)
  end subroutine load_solar_spectrum


  !> Load the sample soil spectrum
  subroutine load_soil_spectrum(wl, rs)
    integer, parameter :: n = nwl_soil  ! number of values/lines
    real(rk), dimension(n), intent(out) :: wl, rs  ! wavelength, soil reflectivity
    real(rk), dimension(n) :: rs_dry, rs_wet
    character(len=*), parameter :: fp = refspecbasepath // "PROSAIL_sample-soil.txt"  ! file path
    real(rk), parameter :: f_wet = 0
    integer :: i

    !> Wavelength isn't in the file, but it is 400:1:2500 in nm just like the PROSPECT one
    wl = [ (real(i, rk), i = 400, 2500, 1) ] / 1000

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
  pure real(rk) function smear1(x, y, xgl, xgu) result(yg)
    real(rk), dimension(:), intent(in) :: x, y
    real(rk), intent(in) :: xgl, xgu
    real(rk) :: area, a1, a2, slope, b1, b2
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
    real(rk), dimension(:), intent(in) :: x, y
    real(rk), dimension(:), intent(in) :: bins
    real(rk), dimension(:), allocatable :: ynew
    integer :: i, n
    n = size(bins)

    allocate(ynew(n-1))
    do i = 1, n - 1
      ynew(i) = smear1(x, y, bins(i), bins(i+1))
    end do
  end function smear


  !> AOP smear
  pure function smear_aop200_planck6000(x, y, bins) result(ynew)
    real(rk), dimension(:), intent(in) :: x, y
      !! `x` must be wavelength in microns!
    real(rk), dimension(:), intent(in) :: bins
    real(rk), dimension(:), allocatable :: ynew

    integer :: i, n
    integer, parameter :: n2 = 200
      !! number of sub-bins within each bin
      ! TODO: compute for each bin based on bin size?
    real(rk) :: yi(n2), wi(n2), binsi(n2+1)

    n = size(bins)  ! number of bin edges

    allocate(ynew(n-1))
    do i = 1, n - 1

      ! Smear to sub-bins within bin
      binsi = wle_equal_width([bins(i), bins(i+1)], n2)
      yi = smear(x, y, binsi)

      ! Construct weights using Planck
      wi = l_wl_plank_integ(6000._rk, binsi(1:n2), binsi(2:n2+1))

      ! Current bin value as weighted average of the sub-bins
      ynew(i) = sum(yi * wi) / sum(wi)

    end do

  end function smear_aop200_planck6000


  !> Trapezoidal integral of y(x)
  pure real(rk) function trapz(x, y) result(res)
    real(rk), dimension(:), intent(in) :: x, y
    integer :: i, n
    n = size(x)

    res = 0  ! don't really need loop to do this
    do i = 1, n - 1
      res = res + (y(i) + y(i+1))/2 * (x(i+1) - x(i))
    end do
  end function trapz


  !> Distribute single rl, tl, rs, idr, and idf reference values into sub-bands
  !> defined by `wle`. Example/reference high-res spectra give the shapes.
  !> Note that this subroutine doesn't currently have the `weight` option that `distribute` has.
  subroutine distribute_rad(  &
    wlbi, rli, tli, rsi, idri, idfi,  &  ! TODO: clean up by making some types?
    wle,  &
    wl, dwl, rl, tl, rs, idr, idf  &
  )
    use m_mb_data  ! reference spectra
    real(rk), intent(in) :: wlbi(2)  ! wl bounds for the single band for the input values
    real(rk), intent(in) :: rli, tli, rsi, idri, idfi  ! input values (single band)
    real(rk), dimension(:), intent(in) :: wle  ! wavelength bounds for new grid
    real(rk), dimension(:), intent(out) :: wl, dwl
    real(rk), dimension(:), intent(out) :: rl, tl, rs, idr, idf  ! new values (spectral)

    ! Local variables
    integer :: n, nbins
    real(rk), dimension(:), allocatable :: sidr, sidf

    ! Comment `use` above and uncomment below to load reference spectra from text files instead
    ! real(rk), dimension(nwl_leaf) :: wl0_leafsoil, rl0, tl0, rs0
    ! real(rk), dimension(nwl_solar) :: wl0_solar, sidr0, sidf0
    ! call load_leaf_spectrum(wl0_leafsoil, rl0, tl0)
    ! call load_soil_spectrum(wl0_leafsoil, rs0)
    ! call load_solar_spectrum(wl0_solar, sidr0, sidf0)

    ! New wavelength grid
    n = ubound(wle, dim=1)  ! number of bin edges (nbins + 1)
    if ( wle(1) /= wle(1) .or. wle(n) /= wlbi(2) ) stop 'wle edges should match orig band'
    dwl = wle(2:n) - wle(1:n-1)  ! new wavelength band widths
    wl = wl(1:n-1) + dwl  ! wavelength band centers (bin midpoints)
    nbins = n - 1
    allocate(sidr(n-1), sidf(n-1))

    ! Smear optical properties
    rl = smear(wl0_leafsoil, rl0, wle)
    tl = smear(wl0_leafsoil, tl0, wle)
    rs = smear(wl0_leafsoil, rs0, wle)

    ! Correct optical property spectra based on the input values
    ! We want the bandwidth-weighted average to match
    rl = rl * rli / (sum(rl * dwl) / (wle(n) - wle(1)))
    tl = tl * tli / (sum(tl * dwl) / (wle(n) - wle(1)))
    rs = rs * rsi / (sum(rs * dwl) / (wle(n) - wle(1)))

    ! Smear irradiances
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


  !> Distribute one value `yi` in band `wlbi` into sub-bands defined by edges `wle`
  function distribute(wlbi, yi, wle, which, weight, weight_method) result(y)
    use m_mb_data  ! reference spectra
    real(rk), intent(in) :: wlbi(2), yi, wle(:)
    character(len=*), intent(in) :: which
    logical, intent(in), optional :: weight
    integer, intent(in), optional :: weight_method
      !! 1 - one correction factor for all bins
      !! 2 - AOP smearing method (all bins corrected individually)
    real(rk), allocatable :: y(:)

    real(rk), allocatable :: wl0(:), y0(:), dwl(:), w(:)
    integer :: n, iwm

    n = ubound(wle, dim=1)
    if ( wle(1) /= wle(1) .or. wle(n) /= wlbi(2) ) stop 'wle edges should match orig band'
    allocate(y(n-1), dwl(n-1), w(n-1))
    dwl = wle(2:n) - wle(1:n-1)

    iwm = 1  ! default weight method
    if (present(weight_method)) then
      iwm = weight_method
    end if

    w = 1  ! default weight for each bin for the overall average
    if ( present(weight) ) then
      if ( weight ) w = l_wl_plank_integ(6000._rk, wle(1:n-1), wle(2:n))
    end if

    ! Load reference wl positions and validate `which`
    select case (which)
      case ('rl', 'tl', 'rs')
        wl0 = wl0_leafsoil
      case ('idr', 'idf')
        wl0 = wl0_solar
      case default
        stop "invalid `which`. Valid options are: 'rl', 'tl'"
    end select

    ! Load reference spectrum
    select case (which)
      case ('rl')
        y0 = rl0
      case ('tl')
        y0 = tl0
      case ('rs')
        y0 = rs0
      case ('idr')
        y0 = sidr0
      case ('idf')
        y0 = sidf0
    end select

    ! Smear reference spectrum to the desired bins
    y = smear(wl0, y0, wle)

    ! Correct based on input value `yi`
    which_var: select case (which)
      case ('idr', 'idf')
        y = y * dwl  ! W m-2 um-1 -> W m-2
        y = y * yi / sum(y)

      case default  ! ('rl', 'tl', 'rs')

        w = w * dwl  ! add waveband width as weight

        aop_weight_method: select case (iwm)
          case (1)
            ! Weighted average correction factor
            ! Weighting by waveband width and, if activated, solar intensity (6000 K Planck)
            continue  ! do nothing

          case (2)
            ! Use AOP smear instead of standard smear for y
            y = smear_aop200_planck6000(wl0, y0, wle)

          case default
            stop 'invalid `weight_method`'

        end select aop_weight_method

        ! Correct using weighted average
        y = y * yi / sum(y * w) * sum(w)

    end select which_var

  end function distribute


  ! Integration, for application after CRT
  ! * reset modified variables in structs to their original values?
  ! * replace CRT outputs by the PAR and NIR integrated ones


end module m_mb
