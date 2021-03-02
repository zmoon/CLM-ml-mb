!> Test out the functionality of the multi-band mod
!> In Bash, use this one-liner in this directory to run:
!> gfortran -c MultibandSolarMod.f90 -I../exe && gfortran test_MultibandSolarMod.f90 -I../exe MultibandSolarMod.o && ./a.out
program test
  use shr_kind_mod, only : r8 => shr_kind_r8
  use MultibandSolarMod
  implicit none

  real(r8) :: &
    wl0_leaf(nwl_leaf), rl0(nwl_leaf), tl0(nwl_leaf), &
    wl0_solar(nwl_solar), si_dr(nwl_solar), si_df(nwl_solar), &
    wl0_soil(nwl_soil), rs0(nwl_soil)

  integer, parameter :: nwl_print = 3

  integer :: i

  !> Check that leaf was loaded correctly
  call load_leaf_spectrum(wl0_leaf, rl0, tl0)
  print *, "leaf spectrum: wavelength (um), refl., trans."
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
  print *, "solar spectrum: wavelength (um), downwelling spectral direct, diffuse"
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
  print *, "soil spectrum: wavelength (um), soil reflectivity"
  do i = 1, nwl_print
    print *, i, wl0_soil(i), rs0(i)
  end do
  print *, "..."
  do i = nwl_soil - nwl_print + 1, nwl_soil
    print *, i, wl0_soil(i), rs0(i)
  end do

end program test
