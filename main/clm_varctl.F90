module clm_varctl

  !-----------------------------------------------------------------------
  ! !DESCRIPTION:
  ! Module containing run control variables
  !-----------------------------------------------------------------------

  use shr_kind_mod, only : r8 => shr_kind_r8
  implicit none

  integer :: light                  ! Radiative transfer: 1 = Norman. 2 = Goudriaan. 3 = Two-stream
  integer :: gstyp                  ! Stomatal conductance: Medlyn (0), Ball-Berry (1), or WUE optimization (2)
  integer :: turb_type              ! Turbulence parameterization
  integer :: pad_type               ! Plant area density: No stem area (0) or use beta distribution (1) or uniform distribution (2)
  logical :: use_colim              ! True: use photosynthesis co-limitation
  logical :: use_acclim             ! True: use photosynthesis acclimation
  logical :: use_clm45kn            ! True: use CLM leaf nitrogen decay coefficient (Kn)
  logical :: no_clumping            ! True: no foliage clumping
  logical :: no_storage             ! True: reduce leaf heat capacity
  logical :: leaf_temp_iter         ! True: iterate leaf temperature during stomatal conductance calculation
  logical :: use_ncan               ! True: ncan > ntop. False: ncan = ntop
  logical :: use_tower              ! True: use tower forcing height and canopy height
  logical :: use_h2ocan             ! True: allow for interception and evaporation
  logical :: use_scalars            ! True: use within-canopy scalars. False: use values at canopy top

  character(len=256) :: diratm      ! Tower input file directory path
  character(len=256) :: dirclm      ! CLM input file directory path
  character(len=256) :: dirout      ! Model output file directory path

  integer, save :: iulog = 6        ! "stdout" log file unit number

  integer :: ntim                   ! Number of time slices to process

  real(r8) :: dtime_sub             ! Model sub-timestep (s)

end module clm_varctl
