!> Test out the functionality/API of the multi-band mod
!> In Bash, use this one-liner in this directory to run:
!> gfortran -c -ffree-line-length-none MultibandSolarDataMod.f90 && gfortran -c -Wall MultibandSolarMod.f90 -I../exe && gfortran -Wall test_MultibandSolarMod.f90 -I../exe MultibandSolarMod.o MultibandSolarDataMod.o && ./a.out
program test
  use m_mb_data, only: rk
  use m_mb
  implicit none

  real(rk), parameter :: pi = 4*atan(1._rk)

  !> Initial data, loaded from the files
  real(rk) :: &
    wl0_leaf(nwl_leaf), rl0(nwl_leaf), tl0(nwl_leaf), &
    wl0_solar(nwl_solar), si_dr(nwl_solar), si_df(nwl_solar), &
    wl0_soil(nwl_soil), rs0(nwl_soil)

  !> Smeared
  real(rk) :: &
    x1(11), y1(11), bins1(6), ynew1(5), dbins1(5), &
    wle1(4), dwl1(3), rl1(3)

  !> Distribute rad
  real(rk) :: wlbi(2), wlbi2(2), rli, tli, rsi, idri, idfi, wle(4), wle2(2), wle3(4)
  real(rk), dimension(3) :: wl, dwl, rl, tl, rs, idr, idf
  real(rk), dimension(1) :: wl2, dwl2, rl2, tl2, rs2, idr2, idf2

  integer, parameter :: nwl_print = 3

  integer :: i, n

  character(len=40) :: fmt1, fmt2

  !> Planck radiance
  print *, '!> Planck radiance L'
  call print_constants()
  print *, 'L(6000 K, 1 um):', l_wl_planck(6000._rk, 1._rk) / 1.e9 / 1.e4
  print *, '  should be: ~ 1.191 W/(sr m2)/m'
  print *, 'L(6000 K, [0.5, 1.0, 1.5] um):', l_wl_planck(6000.0_rk, [0.5_rk, 1._rk, 1.5_rk]) / 1.e9 / 1.e4

  !> Planck radiance definite integrals
  print *
  print *, '!> Planck radiance definite integrals'
  print *, 'L(6000 K) 0.5:1.5:', l_wl_plank_integ(6000._rk, 0.5_rk, 1.5_rk)
  print *, 'L(6000 K) 0.3:5:', l_wl_plank_integ(6000._rk, 0.3_rk, 5._rk)
  print *, 'L(6000 K) 0.1:10:', l_wl_plank_integ(6000._rk, 0.1_rk, 10._rk)
  print *, '0.5:1.5 frac of total (0.1:10):', l_wl_plank_integ(6000._rk, 0.5_rk, 1.5_rk)/l_wl_plank_integ(6000._rk, 0.1_rk, 10._rk)
  print *, '0.5:1.5 frac of total (S-B):', pi*l_wl_plank_integ(6000._rk, 0.5_rk, 1.5_rk)/(5.6704e-8_rk * 6000._rk**4)
  print *, '  should be ~ 0.6168'
  print *, 'multiple bounds (0.4:0.5, 0.5:0.6):', l_wl_plank_integ(6000._rk, [0.4_rk, 0.5_rk], [0.5_rk, 0.6_rk])

  !> Check that leaf was loaded correctly
  call load_leaf_spectrum(wl0_leaf, rl0, tl0)
  print *
  print *, "!> leaf spectrum: wavelength (um), leaf element reflectance, transmittance"
  do i = 1, nwl_print
    print *, i, wl0_leaf(i), rl0(i), tl0(i)
  end do
  print *, "..."
  do i = nwl_leaf - nwl_print + 1, nwl_leaf
    print *, i, wl0_leaf(i), rl0(i), tl0(i)
  end do

  !> Check that solar was loaded correctly
  call load_solar_spectrum(wl0_solar, si_dr, si_df)
  print *
  print *, "!> solar spectrum: wavelength (um), downwelling spectral direct, diffuse"
  do i = 1, nwl_print
    print *, i, wl0_solar(i), si_dr(i), si_df(i)
  end do
  print *, "..."
  do i = nwl_solar - nwl_print + 1, nwl_solar
    print *, i, wl0_solar(i), si_dr(i), si_df(i)
  end do

  !> Check that soil was loaded correctly
  call load_soil_spectrum(wl0_soil, rs0)
  print *
  print *, "!> soil spectrum: wavelength (um), soil reflectivity"
  do i = 1, nwl_print
    print *, i, wl0_soil(i), rs0(i)
  end do
  print *, "..."
  do i = nwl_soil - nwl_print + 1, nwl_soil
    print *, i, wl0_soil(i), rs0(i)
  end do

  !> Basic smearing test
  x1 = [ (real(i, rk), i = 0, 10) ]
  y1 = 25 - (x1 - 5)**2  ! inverted parabola, true integral is 500/3 ~ 166.67
  bins1 = [0._rk, 2.5_rk, 7._rk, 8._rk, 9._rk, 10._rk]
  n = ubound(bins1, dim=1)
  dbins1 = bins1(2:n) - bins1(1:n-1)
  ynew1 = smear(x1, y1, bins1)
  print *
  print *, "!> basic smearing"
  ! print *, "orig:", y1
  print *, "orig trapz:", trapz(x1, y1)
  ! print *, "bins:", bins1
  ! print *, "new:", ynew1
  print *, "new integ.:", sum(ynew1 * dbins1)
  print *, 'smear of ones:', smear(x1, [(1._rk, i = 1, size(x1))], bins1)
  if ( sum(ynew1 * dbins1) /= trapz(x1, y1) ) stop "smear no good"  ! TODO: function for this check

  !> Smear leaf
  wle1 = [0.4_rk, 0.7_rk, 1.0_rk, 2.5_rk]  ! note: same endpoints as orig
  n = ubound(wle1, dim=1)
  dwl1 = wle1(2:n) - wle1(1:n-1)
  rl1 = smear(wl0_leaf, rl0, wle1)
  print *
  print *, "!> smearing leaf reflectance"
  print *, "orig trapz:", trapz(wl0_leaf, rl0)
  print *, "new integ.:", sum(rl1 * dwl1)
  print *, "orig avg.:", trapz(wl0_leaf, rl0) / (wle1(n) - wle1(1))
  print *, "new avg.:", sum(rl1 * dwl1) / (wle1(n) - wle1(1))

  !> Smear soil

  !> Full distribute-rad routine
  wlbi = [0.4_rk, 0.7_rk]  ! (um) PAR band edges
  ! wlbi = [0.7_rk, 2.5_rk]  ! NIR
  rli = 0.2_rk
  tli = 0.15_rk
  rsi = 0.25_rk
  idri = 500
  idfi = 100
  wle = [0.4_rk, 0.5_rk, 0.6_rk, 0.7_rk]  ! new band edges for PAR region
  ! wle = [0.7, 1.0, 1.5, 2.5]  ! NIR
  call distribute_rad(wlbi, rli, tli, rsi, idri, idfi, wle, wl, dwl, rl, tl, rs, idr, idf)
  print *
  print *, "!> distribute rad -- to 3 bands"
  fmt1 = "(1x, a4, 2x, g8.3, a, 2x, 3(g10.3))"
  print fmt1, 'rl', rli, '->', rl
  print fmt1, 'tl', tli, '->', tl
  print fmt1, 'rs', rsi, '->', rs
  print fmt1, 'idr', idri, '->', idr
  print fmt1, 'idf', idfi, '->', idf

  !> Distribute one spectrum at once
  print *
  print *, "!> distribute one spectrum at once (should be same result as above)"
  print fmt1, 'rl', rli, '->', distribute(wlbi, rli, wle, 'rl')
  print fmt1, 'rlw', rli, '->', distribute(wlbi, rli, wle, 'rl', weight=.true.)
  print fmt1, 'rlw2', rli, '->', distribute(wlbi, rli, wle, 'rl', weight=.true., weight_method=2)
  print fmt1, 'tl', tli, '->', distribute(wlbi, tli, wle, 'tl')
  print fmt1, 'tlw', tli, '->', distribute(wlbi, tli, wle, 'tl', weight=.true.)
  print fmt1, 'rs', rsi, '->', distribute(wlbi, rsi, wle, 'rs')
  print fmt1, 'idr', idri, '->', distribute(wlbi, idri, wle, 'idr')
  print fmt1, 'idf', idfi, '->', distribute(wlbi, idfi, wle, 'idf')

  !> Smear rl in NIR
  wlbi2 = [0.8_rk, 2.5_rk]
  wle3 = [0.8_rk, 1.1_rk, 1.8_rk, 2.5_rk]
  print *
  print *, "!> distribute rl to the NIR"
  print fmt1, "rl", rli, '->', distribute(wlbi2, rli, wle3, 'rl')
  print fmt1, "rlw1", rli, '->', distribute(wlbi2, rli, wle3, 'rl', weight=.true.)
  print fmt1, "rlw2", rli, '->', distribute(wlbi2, rli, wle3, 'rl', weight=.true., weight_method=2)

  !> Check that if we have only one bin we get the same as what we put in
  wle2 = wlbi
  call distribute_rad(wlbi, rli, tli, rsi, idri, idfi, wle2, wl2, dwl2, rl2, tl2, rs2, idr2, idf2)
  print *
  print *, "!> distribute rad -- one band only (results should be same as input, though size-1 array)"
  print *, "      orig                        new                         diff"
  print *, 'rl', rli, rl2, rl2 - rli
  print *, 'tl', tli, tl2, tl2 - tli
  print *, 'rs', rsi, rs2, rs2 - rsi
  print *, 'idr', idri, idr2, idr2 - idri
  print *, 'idf', idfi, idf2, idf2 - idfi

  !> Evenly-spaced wavelength grid generation
  print *
  print *, "!> evenly spaced grid generation?"
  fmt2 = "(2(g8.3), a, 4x, 11(g8.3))"
  print fmt2, wle2, "->", wle_equal_width(wle2, 10)

end program test
