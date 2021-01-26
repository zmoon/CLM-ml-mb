module clm_varcon

  !-----------------------------------------------------------------------
  ! !DESCRIPTION:
  ! Module containing various model constants
  !-----------------------------------------------------------------------

  use shr_kind_mod, only : r8 => shr_kind_r8
  implicit none

  !-----------------------------------------------------------------------
  ! Mathematical constants used in CLM
  !-----------------------------------------------------------------------

  real(r8), parameter :: rpi = 3.141592654_r8          ! pi

  !-----------------------------------------------------------------------
  ! Physical constants used in CLM
  !-----------------------------------------------------------------------

  real(r8), parameter :: tfrz = 273.15_r8              ! Freezing point of water (K)
  real(r8), parameter :: sb = 5.67e-08_r8              ! Stefan-Boltzmann constant (W/m2/K4)
  real(r8), parameter :: grav = 9.80665_r8             ! Gravitational acceleration (m/s2)
  real(r8), parameter :: vkc = 0.4_r8                  ! von Karman constant
  real(r8), parameter :: denh2o = 1000._r8             ! Density of liquid water (kg/m3)
  real(r8), parameter :: denice =  917._r8             ! Density of ice (kg/m3)
  real(r8), parameter :: tkwat = 0.57_r8               ! Thermal conductivity of water (W/m/K)
  real(r8), parameter :: tkice = 2.29_r8               ! Thermal conductivity of ice (W/m/K)
  real(r8), parameter :: tkair = 0.023_r8              ! Thermal conductivity of air (W/m/K)
  real(r8), parameter :: hfus = 0.3337e6_r8            ! Latent heat of fusion for water at 0 C (J/kg)
  real(r8), parameter :: hvap = 2.501e6_r8             ! Latent heat of evaporation (J/kg)
  real(r8), parameter :: hsub = hfus + hvap            ! Latent heat of sublimation (J/kg)
  real(r8), parameter :: rgas = 8314.46_r8             ! Universal gas constant (J/K/kmol)
  real(r8), parameter :: cpice = 2.11727e3_r8          ! Specific heat of ice (J/kg/K)
  real(r8), parameter :: cpliq = 4.188e3_r8            ! Specific heat of water (J/kg/K)
  real(r8), parameter :: cnfac  = 0.5_r8               ! Crank-Nicolson factor (between 0 and 1)
  real(r8), parameter :: capr   = 0.34_r8              ! Tuning factor to turn first soil layer temp. into surface temp.
  real(r8), parameter :: thk_bedrock = 3.0_r8          ! Thermal conductivity of saturated granitic rock (W/m/K)

  !-----------------------------------------------------------------------
  ! Special value flags used in CLM
  !-----------------------------------------------------------------------

  real(r8), parameter ::  spval = 1.e36_r8             ! Special value for real data
  integer , parameter :: ispval = -9999                ! Special value for integer data

  !-----------------------------------------------------------------------
  ! New constants for multilayer canopy
  !-----------------------------------------------------------------------

  ! Physical constants for multilayer canopy

  real(r8), parameter :: rgasc = rgas * 1.e-03_r8      ! Universal gas constant (J/K/mol)
  real(r8), parameter :: mmdry = 28.97_r8 * 1.e-03_r8  ! Molecular mass of dry air (kg/mol)
  real(r8), parameter :: mmh2o = 18.02_r8 * 1.e-03_r8  ! Molecular mass of water vapor (kg/mol)
  real(r8), parameter :: cpd = 1005._r8                ! Specific heat of dry air at constant pressure (J/kg/K)
  real(r8), parameter :: cpw = 1846._r8                ! Specific heat of water vapor at constant pressure (J/kg/K)
  real(r8), parameter :: visc0 = 13.3e-06_r8           ! Kinematic viscosity at 0C and 1013.25 hPa (m2/s)
  real(r8), parameter :: dh0 = 18.9e-06_r8             ! Molecular diffusivity (heat) at 0C and 1013.25 hPa (m2/s)
  real(r8), parameter :: dv0 = 21.8e-06_r8             ! Molecular diffusivity (H2O) at 0C and 1013.25 hPa (m2/s)
  real(r8), parameter :: dc0 = 13.8e-06_r8             ! Molecular diffusivity (CO2) at 0C and 1013.25 hPa (m2/s)
  real(r8), parameter :: cpbio = cpliq / 3._r8         ! Specific heat of dry biomass (J/kg/K)

  ! Constants used in the RSL psihat look-up tables

  integer, parameter :: nZ = 276, nL = 41 ! Dimensions of RSL psihat look-up tables
  real(r8) :: zdtgridM(nZ,1)              ! Grid of zdt on which psihat is given for momentum
  real(r8) :: dtLgridM(1,nL)              ! Grid of dtL on which psihat is given for momentum
  real(r8) :: psigridM(nZ,nL)             ! Grid of psihat values for momentum
  real(r8) :: zdtgridH(nZ,1)              ! Grid of zdt on which psihat is given for heat
  real(r8) :: dtLgridH(1,nL)              ! Grid of dtL on which psihat is given for heat
  real(r8) :: psigridH(nZ,nL)             ! Grid of psihat values for heat

end module clm_varcon
