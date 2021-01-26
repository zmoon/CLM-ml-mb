module clmDataMod

  !-----------------------------------------------------------------------
  ! !DESCRIPTION:
  ! Read CLM forcing data for tower site
  !
  ! !USES:
  use shr_kind_mod, only : r8 => shr_kind_r8
  use abortutils,   only : endrun, handle_err
  use ColumnType, only : col
  use SoilStateType, only : soilstate_type
  use WaterStateType, only : waterstate_type
  use EnergyFluxType, only : energyflux_type
  use CanopyStateType, only : canopystate_type
  !
  ! !PUBLIC TYPES:
  implicit none
  include 'netcdf.inc'
  !
  ! !PUBLIC MEMBER FUNCTIONS:
  public :: clmData               ! Read CLM forcing data
  !
  ! !PRIVATE MEMBER FUNCTIONS:
  private :: readCLMdata          ! Read single-level variables from CLM netcdf history file
  private :: readCLMhyd           ! Read volumetric soil water from CLM netcdf history file
  !-----------------------------------------------------------------------

contains

  !-----------------------------------------------------------------------
  subroutine clmData (strt, ncfilename, begp, endp, begc, endc, &
  soilstate_inst, waterstate_inst, energyflux_inst, canopystate_inst)
    !
    ! !DESCRIPTION:
    ! Read variables from CLM netcdf history file
    !
    ! !USES:
    use clm_varcon, only : denh2o
    use clm_varpar, only : nlevgrnd 
    !
    ! !ARGUMENTS:
    implicit none
    integer, intent(in) :: strt                  ! Current time slice of data to retrieve from CLM history file
    character(len=*), intent(in) :: ncfilename   ! CLM netcdf filename
    integer, intent(in) :: begp, endp            ! First and last patch
    integer, intent(in) :: begc, endc            ! First and last column
    type(soilstate_type), intent(inout) :: soilstate_inst
    type(waterstate_type) , intent(out) :: waterstate_inst
    type(energyflux_type), intent(out)  :: energyflux_inst
    type(canopystate_type) , intent(out)  :: canopystate_inst
    !
    ! !LOCAL VARIABLES:
    integer  :: p                                ! Patch index for CLM g/l/c/p hierarchy
    integer  :: c                                ! Column index for CLM g/l/c/p hierarchy
    integer  :: j                                ! Soil layer index
    real(r8) :: s                                ! Soil layer water content relative to saturation (fraction)

    real(r8) :: btran_loc(1,1,1)                 ! CLM: Ball-Berry soil wetness factor (-)
    real(r8) :: coszen_loc(1,1,1)                ! CLM: Cosine solar zenith angle
    real(r8) :: elai_loc(1,1,1)                  ! CLM: Leaf area index (m2/m2)
    real(r8) :: esai_loc(1,1,1)                  ! CLM: Stem area index (m2/m2)
    real(r8) :: h2osoi_loc(1,1,nlevgrnd)         ! CLM: Volumetric soil moisture (m3/m3)
    !---------------------------------------------------------------------

    associate ( &
                                                      ! *** Input ***
    dz         => col%dz                         , &  ! Soil layer thickness (m)
    watsat     => soilstate_inst%watsat_col      , &  ! Soil layer volumetric water content at saturation (porosity)
    sucsat     => soilstate_inst%sucsat_col      , &  ! Soil layer suction (negative matric potential) at saturation (mm)
    bsw        => soilstate_inst%bsw_col         , &  ! Soil layer Clapp and Hornberger "b" parameter
                                                      ! *** Output ***
    btran      => energyflux_inst%btran_patch    , &  ! Ball-Berry soil wetness factor (-)
    elai       => canopystate_inst%elai_patch    , &  ! Leaf area index of canopy (m2/m2)
    esai       => canopystate_inst%esai_patch    , &  ! Stem area index of canopy (m2/m2)
    h2osoi_vol => waterstate_inst%h2osoi_vol_col , &  ! Soil layer volumetric water content (m3/m3)
    h2osoi_liq => waterstate_inst%h2osoi_liq_col , &  ! Soil layer liquid water (kg/m2)
    h2osoi_ice => waterstate_inst%h2osoi_ice_col , &  ! Soil layer ice lens (kg/m2)
    smp_l      => soilstate_inst%smp_l_col         &  ! Soil layer matric potential (mm)
    )

    ! Read CLM data for current time step. CLM coszen is for the next model time step,
    ! so must read value at previous time step (i.e., use strt - 1) to get correct coszen

    call readCLMdata (ncfilename, strt, btran_loc, coszen_loc, elai_loc, esai_loc)

    ! All patches (p) get these values

    do p = begp, endp
       btran(p) = btran_loc(1,1,1)
       elai(p) = elai_loc(1,1,1)
       esai(p) = esai_loc(1,1,1)
    end do

    ! Read CLM volumetric soil water for current time step

    call readCLMhyd (ncfilename, strt, h2osoi_loc)

    ! All columns get this soil water

    do j = 1, nlevgrnd
       do c = begc, endc
          h2osoi_vol(c,j) = h2osoi_loc(1,1,j)
          h2osoi_liq(c,j) = h2osoi_vol(c,j) * dz(c,j) * denh2o
          h2osoi_ice(c,j) = 0._r8
          s = max(min(h2osoi_vol(c,j)/watsat(c,j), 1._r8), 0.01_r8)
          smp_l(c,j) = -sucsat(c,j) * s**(-bsw(c,j))
          smp_l(c,j) = max(smp_l(c,j), -1.e08_r8)
       end do
    end do

    end associate
  end subroutine clmData

  !-----------------------------------------------------------------------
  subroutine readCLMdata (ncfilename, strt, btran_mod, coszen_mod, elai_mod, esai_mod)
    !
    ! !DESCRIPTION:
    ! Read single-level variables from CLM netcdf history file
    !
    ! !USES:
    use clm_varctl, only : pad_type 
    !
    ! !ARGUMENTS:
    implicit none
    character(len=*), intent(in) :: ncfilename     ! netcdf filename
    integer,  intent(in)  :: strt                  ! Time slice of data to retrieve
    real(r8), intent(out) :: btran_mod(1,1,1)      ! Ball-Berry soil wetness factor (-)
    real(r8), intent(out) :: coszen_mod(1,1,1)     ! Cosine solar zenith angle
    real(r8), intent(out) :: elai_mod(1,1,1)       ! Leaf area index (m2/m2)
    real(r8), intent(out) :: esai_mod(1,1,1)       ! Stem area index (m2/m2)
    !
    ! !LOCAL VARIABLES:
    integer :: ncid                                ! netcdf file ID
    integer :: status                              ! Function return status
    integer :: varid                               ! netcdf variable id
    integer :: start2(2), count2(2)                ! Start and count arrays for reading 2-D data from netcdf files
    !-----------------------------------------------------------------------

    status = nf_open(ncfilename, nf_nowrite, ncid)
    if (status /= nf_noerr) call handle_err(status, ncfilename)

    ! Dimensions in FORTRAN are in column major order: the first array index varies the most rapidly.
    ! In NetCDF file the dimensions appear in the opposite order: lat, lon (2-D); time, lat, lon (3-D);
    ! time, levgrnd, lat, lon (4-D)
 
    start2 = (/ 1, strt /)
    count2 = (/ 1, 1 /)

    ! btran(nlndgrid, ntime): transpiration beta factor (-)

    status = nf_inq_varid(ncid, "BTRAN", varid)
    if (status /= nf_noerr) call handle_err(status, "BTRAN")

    status = nf_get_vara_double(ncid, varid, start2, count2, btran_mod)
    if (status /= nf_noerr) call handle_err(status, "btran_mod")

    ! coszen(nlndgrid, ntime): cosine of solar zenith angle

    status = nf_inq_varid(ncid, "COSZEN", varid)
    if (status /= nf_noerr) call handle_err(status, "COSZEN")

    status = nf_get_vara_double(ncid, varid, start2, count2, coszen_mod)
    if (status /= nf_noerr) call handle_err(status, "coszen_mod")

    ! elai(nlndgrid, ntime): leaf area index (m2/m2)

    status = nf_inq_varid(ncid, "ELAI", varid)
    if (status /= nf_noerr) call handle_err(status, "ELAI")

    status = nf_get_vara_double(ncid, varid, start2, count2, elai_mod)
    if (status /= nf_noerr) call handle_err(status, "elai_mod")

    ! esai(nlndgrid, ntime): stem area index (m2/m2)

    status = nf_inq_varid(ncid, "ESAI", varid)
    if (status /= nf_noerr) call handle_err(status, "ESAI")

    status = nf_get_vara_double(ncid, varid, start2, count2, esai_mod)
    if (status /= nf_noerr) call handle_err(status, "esai_mod")

    if (pad_type == 0) esai_mod(1,1,1) = 0._r8

    status = nf_close(ncid)

  end subroutine readCLMdata

  !-----------------------------------------------------------------------
  subroutine readCLMhyd (ncfilename, strt, h2osoi_mod)
    !
    ! !DESCRIPTION:
    ! Read volumetric soil water from CLM netcdf history file
    !
    ! !USES:
    use clm_varpar, only : nlevgrnd
    !
    ! !ARGUMENTS:
    implicit none
    character(len=*), intent(in)  :: ncfilename  ! netcdf filename
    integer, intent(in) :: strt                  ! Time slice of data to retrieve
    real(r8), intent(out) :: h2osoi_mod(1,1,nlevgrnd) ! Volumetric soil moisture for CLM soil layers (m3/m3)
    !
    ! !LOCAL VARIABLES:
    integer  :: ncid                             ! netcdf file ID
    integer  :: status                           ! Function return status
    integer  :: varid                            ! netcdf variable id
    integer  :: start3(3), count3(3)             ! Start and count arrays for reading 3-D data from netcdf files
    !---------------------------------------------------------------------

    status = nf_open(ncfilename, nf_nowrite, ncid)
    if (status /= nf_noerr) call handle_err(status, ncfilename)

    ! Dimensions in FORTRAN are in column major order: the first array index varies the most rapidly.
    ! In NetCDF file the dimensions appear in the opposite order: lat, lon (2-D); time, lat, lon (3-D);
    ! time, levgrnd, lat, lon (4-D)

    start3 = (/ 1,  1, strt /)
    count3 = (/ 1, nlevgrnd, 1 /)

    ! Read h2osoi(nlndgrid, nlevgrnd, ntime): volumetric soil water

    status = nf_inq_varid(ncid, "H2OSOI", varid)
    if (status /= nf_noerr) call handle_err(status, "H2OSOI")

    status = nf_get_vara_double(ncid, varid, start3, count3, h2osoi_mod)
    if (status /= nf_noerr) call handle_err(status, "h2osoi_mod")

    status = nf_close(ncid)

  end subroutine readCLMhyd

end module clmDataMod
