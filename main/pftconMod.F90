module pftconMod

  !-----------------------------------------------------------------------
  ! !DESCRIPTION:
  ! Module containing and initializing constants for plant functional types
  !
  ! !USES:
  use shr_kind_mod, only : r8 => shr_kind_r8
  use clm_varpar, only : mxpft, numrad, ivis, inir
  use clm_varctl, only : gstyp, no_clumping
  !
  ! !PUBLIC TYPES:
  implicit none
  !
  ! PFT constants
  !
  type, public :: pftcon_type

    ! Existing CLM constants
    real(r8), allocatable :: xl           (:)   ! Departure of leaf angle from spherical orientation (-)
    real(r8), allocatable :: clump_fac    (:)   ! Foliage clumping index (-)
    real(r8), allocatable :: rhol         (:,:) ! Leaf reflectance (-) [for numrad wavebands]
    real(r8), allocatable :: taul         (:,:) ! Leaf transmittance (-) [for numrad wavebands]
    real(r8), allocatable :: rhos         (:,:) ! Stem reflectance (-) [for numrad wavebands]
    real(r8), allocatable :: taus         (:,:) ! Stem transmittance (-) [for numrad wavebands]
    real(r8), allocatable :: c3psn        (:)   ! Photosynthetic pathway: 1. = c3 plant and 0. = c4 plant
    real(r8), allocatable :: dleaf        (:)   ! Leaf dimension (m)
    real(r8), allocatable :: slatop       (:)   ! Specific leaf area at top of canopy (m2/gC)
    real(r8), allocatable :: z0mr         (:)   ! Ratio of momentum roughness length to canopy top height (-)
    real(r8), allocatable :: displar      (:)   ! Ratio of displacement height to canopy top height (-)
    real(r8), allocatable :: roota_par    (:)   ! CLM rooting distribution parameter (1/m)
    real(r8), allocatable :: rootb_par    (:)   ! CLM rooting distribution parameter (1/m)

    ! New constants for multilayer canopy
    real(r8), allocatable :: vcmaxpft     (:)   ! Maximum carboxylation rate at 25C (umol/m2/s)
    real(r8), allocatable :: minlwp       (:)   ! Minimum leaf water potential (MPa)
    real(r8), allocatable :: gplant       (:)   ! Stem (xylem-to-leaf) hydraulic conductance (mmol H2O/m2 leaf area/s/Mpa)
    real(r8), allocatable :: capac        (:)   ! Plant capacitance (mmol H2O/m2 leaf area/MPa)
    real(r8), allocatable :: iota         (:)   ! Stomatal water-use efficiency (umol CO2/ mol H2O)
    real(r8), allocatable :: root_radius  (:)   ! Fine root radius (m)
    real(r8), allocatable :: root_density (:)   ! Fine root density (g biomass / m3 root)
    real(r8), allocatable :: root_resist  (:)   ! Hydraulic resistivity of root tissue (MPa.s.g/mmol H2O)
    real(r8), allocatable :: g0opt        (:)   ! Ball-Berry minimum leaf conductance, unstressed (mol H2O/m2/s)
    real(r8), allocatable :: g1opt        (:)   ! Ball-Berry slope of conductance-photosynthesis relationship, unstressed
    real(r8), allocatable :: emleaf       (:)   ! Leaf emissivity (-)

  contains

    procedure, public  :: Init
    procedure, private :: InitAllocate
    procedure, private :: InitRead

  end type pftcon_type

  type(pftcon_type), public :: pftcon
  !-----------------------------------------------------------------------

contains

  !-----------------------------------------------------------------------
  subroutine Init (this)

    class(pftcon_type) :: this

    call this%InitAllocate()
    call this%InitRead()

  end subroutine Init

  !-----------------------------------------------------------------------
  subroutine InitAllocate (this)
    !
    ! !DESCRIPTION:
    ! Allocate memory for pft data structure
    !
    ! !ARGUMENTS:
    class(pftcon_type) :: this
    !---------------------------------------------------------------------

    ! Existing CLM constants
    allocate (this%xl           (0:mxpft))
    allocate (this%clump_fac    (0:mxpft))
    allocate (this%rhol         (0:mxpft,numrad))
    allocate (this%taul         (0:mxpft,numrad))
    allocate (this%rhos         (0:mxpft,numrad))
    allocate (this%taus         (0:mxpft,numrad))
    allocate (this%c3psn        (0:mxpft))
    allocate (this%dleaf        (0:mxpft))
    allocate (this%slatop       (0:mxpft))
    allocate (this%z0mr         (0:mxpft))
    allocate (this%displar      (0:mxpft))
    allocate (this%roota_par    (0:mxpft))
    allocate (this%rootb_par    (0:mxpft))

    ! New constants for multilayer canopy
    allocate (this%vcmaxpft     (0:mxpft))
    allocate (this%minlwp       (0:mxpft))
    allocate (this%gplant       (0:mxpft))
    allocate (this%capac        (0:mxpft))
    allocate (this%iota         (0:mxpft))
    allocate (this%root_radius  (0:mxpft))
    allocate (this%root_density (0:mxpft))
    allocate (this%root_resist  (0:mxpft))
    allocate (this%g0opt        (0:mxpft))
    allocate (this%g1opt        (0:mxpft))
    allocate (this%emleaf       (0:mxpft))

  end subroutine InitAllocate

  !-----------------------------------------------------------------------
  subroutine InitRead (this)
    !
    ! !DESCRIPTION:
    ! Read and initialize vegetation (PFT) constants
    !
    ! !ARGUMENTS:
    class(pftcon_type) :: this
    !---------------------------------------------------------------------

    ! PFTs
    ! -------------------------------
    !  1 = Needleleaf evergreen tree - temperate
    !  2 = Needleleaf evergreen tree - boreal
    !  3 = Needleleaf deciduous tree - boreal
    !  4 = Broadleaf evergreen tree  - tropical
    !  5 = Broadleaf evergreen tree  - temperate
    !  6 = Broadleaf deciduous tree  - tropical
    !  7 = Broadleaf deciduous tree  - temperate
    !  8 = Broadleaf deciduous tree  - boreal
    !  9 = Broadleaf evergreen shrub - temperate
    ! 10 = Broadleaf deciduous shrub - temperate
    ! 11 = Broadleaf deciduous shrub - boreal
    ! 12 = C3 grass - arctic
    ! 13 = C3 grass
    ! 14 = C4 grass
    ! 15 = Crop
    ! 16 = Crop
    ! -------------------------------

    ! Leaf angle

    this%xl( 1: 3) = 0.01_r8
    this%xl( 4: 5) = 0.10_r8
    this%xl( 6: 6) = 0.01_r8
    this%xl( 7: 8) = 0.25_r8
    this%xl( 9: 9) = 0.01_r8
    this%xl(10:11) = 0.25_r8
    this%xl(12:16) = -0.30_r8

    ! Foliage clumping index

    this%clump_fac( 1: 2) = 0.74_r8
    this%clump_fac( 3: 3) = 0.78_r8
    this%clump_fac( 4: 5) = 0.66_r8
    this%clump_fac( 6: 8) = 0.70_r8
    this%clump_fac( 9:11) = 0.75_r8
    this%clump_fac(12:16) = 0.75_r8

    if (no_clumping) then
       this%clump_fac(1:16) = 1._r8
    end if

    ! Leaf reflectance and transmittance: visible and near-infrared

    this%rhol( 1: 3,ivis) = 0.07_r8
    this%rhol( 4: 8,ivis) = 0.10_r8
    this%rhol( 9: 9,ivis) = 0.07_r8
    this%rhol(10:11,ivis) = 0.10_r8
    this%rhol(12:16,ivis) = 0.11_r8

    this%rhol( 1: 3,inir) = 0.35_r8
    this%rhol( 4: 8,inir) = 0.45_r8
    this%rhol( 9: 9,inir) = 0.35_r8
    this%rhol(10:11,inir) = 0.45_r8
    this%rhol(12:16,inir) = 0.35_r8

    this%taul(:,ivis) = 0.05_r8

    this%taul( 1: 3,inir) = 0.10_r8
    this%taul( 4: 8,inir) = 0.25_r8
    this%taul( 9: 9,inir) = 0.10_r8
    this%taul(10:11,inir) = 0.25_r8
    this%taul(12:16,inir) = 0.34_r8

    this%rhos( 1:11,ivis) = 0.16_r8
    this%rhos(12:16,ivis) = 0.31_r8

    this%rhos( 1:11,inir) = 0.39_r8
    this%rhos(12:16,inir) = 0.53_r8

    this%taus( 1:11,ivis) = 0.001_r8
    this%taus(12:16,ivis) = 0.12_r8

    this%taus( 1:11,inir) = 0.001_r8
    this%taus(12:16,inir) = 0.25_r8

    ! Photosynthetic pathway: 1. = c3 plant and 0. = c4 plant

    this%c3psn(:)  = 1._r8
    this%c3psn(14) = 0._r8

    ! Leaf dimension

    this%dleaf(:) = 0.04_r8

    ! Specific leaf area at top of canopy, projected area basis (m2/gC)

    this%slatop( 1) = 0.010_r8
    this%slatop( 2) = 0.008_r8
    this%slatop( 3) = 0.024_r8
    this%slatop( 4) = 0.012_r8
    this%slatop( 5) = 0.012_r8
    this%slatop( 6) = 0.030_r8
    this%slatop( 7) = 0.030_r8
    this%slatop( 8) = 0.030_r8
    this%slatop( 9) = 0.012_r8
    this%slatop(10) = 0.030_r8
    this%slatop(11) = 0.030_r8
    this%slatop(12) = 0.030_r8
    this%slatop(13) = 0.030_r8
    this%slatop(14) = 0.030_r8
    this%slatop(15) = 0.030_r8
    this%slatop(16) = 0.030_r8

    ! Ratio of momentum roughness length to canopy top height

    this%z0mr = (/ 0._r8, 0.055_r8, 0.055_r8, 0.055_r8, 0.075_r8, 0.075_r8, 0.055_r8, 0.055_r8, 0.055_r8, &
                          0.120_r8, 0.120_r8, 0.120_r8, 0.120_r8, 0.120_r8, 0.120_r8, 0.120_r8, 0.120_r8 /)

    ! Ratio of displacement height to canopy top height

    this%displar = (/ 0._r8, 0.67_r8, 0.67_r8, 0.67_r8, 0.67_r8, 0.67_r8, 0.67_r8, 0.67_r8, 0.67_r8, &
                             0.68_r8, 0.68_r8, 0.68_r8, 0.68_r8, 0.68_r8, 0.68_r8, 0.68_r8, 0.68_r8 /)

    ! CLM rooting distribution parameters (1/m)

    this%roota_par( 1: 5) = 7._r8
    this%roota_par( 6: 8) = 6._r8
    this%roota_par( 9:11) = 7._r8
    this%roota_par(12:14) = 11._r8
    this%roota_par(15:16) = 6._r8

    this%rootb_par( 1: 3) = 2._r8
    this%rootb_par( 4: 5) = 1._r8
    this%rootb_par( 6: 8) = 2._r8
    this%rootb_par( 9:11) = 1.5_r8
    this%rootb_par(12:14) = 2._r8
    this%rootb_par(15:16) = 3._r8

    ! vcmax

    this%vcmaxpft( 1) = 62.5_r8
    this%vcmaxpft( 2) = 62.5_r8
    this%vcmaxpft( 3) = 39.1_r8
    this%vcmaxpft( 4) = 41.0_r8
    this%vcmaxpft( 5) = 61.4_r8
    this%vcmaxpft( 6) = 41.0_r8
    this%vcmaxpft( 7) = 57.7_r8
    this%vcmaxpft( 8) = 57.7_r8
    this%vcmaxpft( 9) = 61.7_r8
    this%vcmaxpft(10) = 54.0_r8
    this%vcmaxpft(11) = 54.0_r8
    this%vcmaxpft(12) = 78.2_r8
    this%vcmaxpft(13) = 78.2_r8
    this%vcmaxpft(14) = 51.6_r8
    this%vcmaxpft(15) = 100.7_r8
    this%vcmaxpft(16) = 100.7_r8

    ! Plant hydraulics

    this%minlwp(:) = -2._r8
    this%gplant(:) = 4._r8
    this%capac( 1:11) = 2500.0_r8
    this%capac(12:16) = 500.0_r8

    ! Root hydraulics

    this%root_radius(:) = 0.29e-03_r8
    this%root_density(:) = 0.31e06_r8
    this%root_resist(:) = 25._r8
!   this%root_resist(1:3) = 75._r8

    ! Stomatal optimization

    this%iota(:) = 750._r8
    this%iota(2:3) = 1500._r8
    this%iota(4) = 500._r8
!   this%iota(15:16) = 500._r8

    ! Ball-Berry stomatal conductance parameters

    if (gstyp == 1) then
       this%g0opt(:)  = 0.01_r8
       this%g0opt(14) = 0.04_r8
       this%g1opt(:)  = 9._r8
       this%g1opt(14) = 4._r8
    end if

    ! Medlyn stomatal conductance parameters

    if (gstyp == 0) then
       this%g0opt(:)  = 0.0001_r8
       this%g1opt( 1) = 2.35_r8
       this%g1opt( 2) = 2.35_r8
       this%g1opt( 3) = 2.35_r8
       this%g1opt( 4) = 4.12_r8
       this%g1opt( 5) = 4.12_r8
       this%g1opt( 6) = 4.45_r8
       this%g1opt( 7) = 4.45_r8
       this%g1opt( 8) = 4.45_r8
       this%g1opt( 9) = 4.70_r8
       this%g1opt(10) = 4.70_r8
       this%g1opt(11) = 4.70_r8
       this%g1opt(12) = 2.22_r8
       this%g1opt(13) = 5.25_r8
       this%g1opt(14) = 1.62_r8
       this%g1opt(15) = 5.79_r8
       this%g1opt(16) = 5.79_r8
    end if

    ! Leaf emissivity

    this%emleaf(:) = 0.98_r8

  end subroutine InitRead

end module pftconMod
