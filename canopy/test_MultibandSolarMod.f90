!> Test out the functionality of the multi-band mod
!> In Bash, use this one-liner in this directory to run:
!> gfortran -c -Wall MultibandSolarMod.f90 -I../exe && gfortran -Wall test_MultibandSolarMod.f90 -I../exe MultibandSolarMod.o && ./a.out
program test
  use shr_kind_mod, only : r8 => shr_kind_r8
  use MultibandSolarMod
  implicit none

  !> Initial data, loaded from the files
  real(r8) :: &
    wl0_leaf(nwl_leaf), rl0(nwl_leaf), tl0(nwl_leaf), &
    wl0_solar(nwl_solar), si_dr(nwl_solar), si_df(nwl_solar), &
    wl0_soil(nwl_soil), rs0(nwl_soil)

  !> Smeared
  real(r8) :: & 
    x1(11), y1(11), bins1(6), ynew1(5), dbins1(5), &
    wle1(4), dwl1(3), rl1(3)

  integer, parameter :: nwl_print = 3

  integer :: i, n

  !> Check that leaf was loaded correctly
  call load_leaf_spectrum(wl0_leaf, rl0, tl0)
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
  x1 = [ (real(i, r8), i = 0, 10) ]
  y1 = 25 - (x1 - 5)**2  ! inverted parabola, true integral is 500/3 ~ 166.67
  bins1 = [0., 2.5, 7., 8., 9., 10.]
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
  if ( sum(ynew1 * dbins1) /= trapz(x1, y1) ) stop "smear no good"  ! TODO: function for this check

  !> Smear leaf
  wle1 = [0.4, 0.7, 1.0, 2.5]  ! note: same endpoints as orig
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

end program test
