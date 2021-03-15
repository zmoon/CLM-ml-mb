# ###############################################
# Makefile
#
# To run: type 'make' or 'make -f Makefile'
#
# To remove obj files and mod files: type 'make clean'
# ###############################################
 
# LIB_NETCDF=     /usr/local/netcdf/lib
# MOD_NETCDF=     /usr/local/netcdf/lib
LIB_NETCDF = /usr/include
MOD_NETCDF = /usr/include

# Source code directories
 
bgp = ../biogeophys
can = ../canopy
shr = ../csm_share
drv = ../driver
main = ../main
util = ../utils
 
# Directory for object files

O = obj

# Executable is "prgm.exe"

prgm = prgm.exe

# Fortran 95 compiler

# cmplr = pgf90
# cmplr = pgf90 -C -Minform=inform -Ktrap=fp
cmplr = gfortran -fno-range-check -std=f2008 -ffree-line-length-none

# Specify module files (that all other files use). Order matters, because some files
# are used by other files in mods list

mods = ${O}/shr_kind_mod.o ${O}/clm_varcon.o ${O}/clm_varctl.o ${O}/clm_varpar.o \
	${O}/CanopyStateType.o ${O}/EnergyFluxType.o ${O}/FrictionVelocityMod.o ${O}/SoilStateType.o \
	${O}/SolarAbsorbedType.o ${O}/SurfaceAlbedoType.o ${O}/TemperatureType.o ${O}/WaterFluxType.o \
	${O}/WaterStateType.o ${O}/atm2lndType.o ${O}/ColumnType.o ${O}/PatchType.o ${O}/ncdio_pio.o \
	${O}/decompMod.o ${O}/pftconMod.o ${O}/abortutils.o \
	${O}/shr_file_mod.o ${O}/clm_time_manager.o ${O}/WaterVaporMod.o \
	${O}/SoilTexMod.o ${O}/TowerDataMod.o ${O}/clm_varorb.o ${O}/shr_orb_mod.o \
	${O}/histFileMod.o ${O}/restUtilMod.o ${O}/CanopyFluxesMultilayerType.o ${O}/filterMod.o \
	${O}/MathToolsMod.o

# Specify all other files. Order matters, because some files are used by other files in src1 list

src1 = ${O}/SurfaceAlbedoMod.o ${O}/clm_instMod.o ${O}/controlMod.o \
	${O}/BandDiagonalMod.o ${O}/SoilTemperatureMod.o \
	${O}/LeafBoundaryLayerMod.o ${O}/LeafPhotosynthesisMod.o ${O}/LeafTemperatureMod.o \
	${O}/StomataOptimizationMod.o ${O}/LeafFluxesMod.o ${O}/LongwaveRadiationMod.o ${O}/PlantHydraulicsMod.o \
	${O}/SoilFluxesMultilayerMod.o ${O}/MultibandSolarMod.o ${O}/SolarRadiationMod.o ${O}/SurfaceResistanceMod.o \
	${O}/TowerMetMod.o ${O}/clmDataMod.o ${O}/CanopyWaterMod.o ${O}/CanopyNitrogenProfileMod.o \
	${O}/CanopyTurbulenceMod.o

src2 = ${O}/CanopyFluxesMultilayerMod.o ${O}/CLMml_initializeMod.o ${O}/clm_driver.o ${O}/CLMml_driver.o ${O}/CLMml.o

# Object files

objs = $(mods) $(src1) $(src2)

# Compile code into an executable file defined by "prgm"

$(prgm): $(objs)
	$(cmplr) -o $(prgm) $(objs) -L$(LIB_NETCDF) -lnetcdff -lblas -llapack -lm -L/usr/lib64

# mods files first, because they are used by other files

${O}/shr_kind_mod.o: ${shr}/shr_kind_mod.F90
	$(cmplr) -c ${shr}/shr_kind_mod.F90
	mv shr_kind_mod.o ${O}

${O}/shr_file_mod.o: ${shr}/shr_file_mod.F90 ${O}/shr_kind_mod.o ${O}/clm_varctl.o  ${O}/abortutils.o 
	$(cmplr) -c ${shr}/shr_file_mod.F90
	mv shr_file_mod.o ${O}

${O}/shr_orb_mod.o: ${shr}/shr_orb_mod.F90 ${O}/shr_kind_mod.o ${O}/clm_varcon.o
	$(cmplr) -c ${shr}/shr_orb_mod.F90
	mv shr_orb_mod.o ${O}

${O}/clm_varcon.o: ${main}/clm_varcon.F90 ${O}/shr_kind_mod.o
	$(cmplr) -c ${main}/clm_varcon.F90
	mv clm_varcon.o ${O}

${O}/clm_varctl.o: ${main}/clm_varctl.F90 ${O}/shr_kind_mod.o
	$(cmplr) -c ${main}/clm_varctl.F90
	mv clm_varctl.o ${O}

${O}/clm_varpar.o: ${main}/clm_varpar.F90 ${O}/shr_kind_mod.o
	$(cmplr) -c ${main}/clm_varpar.F90
	mv clm_varpar.o ${O}

${O}/decompMod.o: ${main}/decompMod.F90 ${O}/shr_kind_mod.o
	$(cmplr) -c ${main}/decompMod.F90
	mv decompMod.o ${O}

${O}/pftconMod.o: ${main}/pftconMod.F90 ${O}/shr_kind_mod.o ${O}/clm_varpar.o ${O}/clm_varctl.o
	$(cmplr) -c ${main}/pftconMod.F90
	mv pftconMod.o ${O}

${O}/abortutils.o: ${main}/abortutils.F90 ${O}/shr_kind_mod.o ${O}/clm_varctl.o
	$(cmplr) -c ${main}/abortutils.F90 -I$(MOD_NETCDF)
	mv abortutils.o ${O}

${O}/clm_time_manager.o: ${util}/clm_time_manager.F90 ${O}/shr_kind_mod.o ${O}/clm_varctl.o ${O}/abortutils.o
	$(cmplr) -c ${util}/clm_time_manager.F90
	mv clm_time_manager.o ${O}

${O}/clm_varorb.o: ${util}/clm_varorb.F90 ${O}/shr_kind_mod.o
	$(cmplr) -c ${util}/clm_varorb.F90
	mv clm_varorb.o ${O}

${O}/MathToolsMod.o: ${can}/MathToolsMod.F90 ${O}/shr_kind_mod.o ${O}/clm_varctl.o ${O}/abortutils.o
	$(cmplr) -c ${can}/MathToolsMod.F90
	mv MathToolsMod.o ${O}

${O}/WaterVaporMod.o: ${can}/WaterVaporMod.F90 ${O}/shr_kind_mod.o ${O}/clm_varctl.o ${O}/clm_varcon.o
	$(cmplr) -c ${can}/WaterVaporMod.F90
	mv WaterVaporMod.o ${O}

${O}/SoilTexMod.o: ${drv}/SoilTexMod.F90 ${O}/shr_kind_mod.o
	$(cmplr) -c ${drv}/SoilTexMod.F90
	mv SoilTexMod.o ${O}

${O}/TowerDataMod.o: ${drv}/TowerDataMod.F90 ${O}/shr_kind_mod.o
	$(cmplr) -c ${drv}/TowerDataMod.F90
	mv TowerDataMod.o ${O}

${O}/atm2lndType.o: ${main}/atm2lndType.F90 ${O}/shr_kind_mod.o ${O}/clm_varpar.o ${O}/decompMod.o
	$(cmplr) -c ${main}/atm2lndType.F90
	mv atm2lndType.o ${O}

${O}/ColumnType.o: ${main}/ColumnType.F90 ${O}/shr_kind_mod.o ${O}/clm_varpar.o ${O}/clm_varcon.o
	$(cmplr) -c ${main}/ColumnType.F90
	mv ColumnType.o ${O}

${O}/PatchType.o: ${main}/PatchType.F90 ${O}/shr_kind_mod.o ${O}/clm_varcon.o
	$(cmplr) -c ${main}/PatchType.F90
	mv PatchType.o ${O}

${O}/ncdio_pio.o: ${main}/ncdio_pio.F90
	$(cmplr) -c ${main}/ncdio_pio.F90
	mv ncdio_pio.o ${O}

${O}/histFileMod.o: ${main}/histFileMod.F90 ${O}/shr_kind_mod.o ${O}/decompMod.o
	$(cmplr) -c ${main}/histFileMod.F90
	mv histFileMod.o ${O}

${O}/restUtilMod.o: ${util}/restUtilMod.F90 ${O}/shr_kind_mod.o ${O}/clm_varpar.o ${O}/ncdio_pio.o
	$(cmplr) -c ${util}/restUtilMod.F90
	mv restUtilMod.o ${O}

${O}/CanopyStateType.o: ${bgp}/CanopyStateType.F90 ${O}/shr_kind_mod.o ${O}/clm_varcon.o ${O}/clm_varpar.o ${O}/decompMod.o
	$(cmplr) -c ${bgp}/CanopyStateType.F90
	mv CanopyStateType.o ${O}

${O}/EnergyFluxType.o: ${bgp}/EnergyFluxType.F90 ${O}/shr_kind_mod.o ${O}/clm_varcon.o ${O}/clm_varpar.o ${O}/decompMod.o
	$(cmplr) -c ${bgp}/EnergyFluxType.F90
	mv EnergyFluxType.o ${O}

${O}/FrictionVelocityMod.o: ${bgp}/FrictionVelocityMod.F90 ${O}/shr_kind_mod.o ${O}/clm_varcon.o ${O}/decompMod.o
	$(cmplr) -c ${bgp}/FrictionVelocityMod.F90
	mv FrictionVelocityMod.o ${O}

${O}/SoilStateType.o: ${bgp}/SoilStateType.F90 ${O}/shr_kind_mod.o ${O}/clm_varcon.o ${O}/clm_varpar.o ${O}/decompMod.o
	$(cmplr) -c ${bgp}/SoilStateType.F90
	mv SoilStateType.o ${O}

${O}/SolarAbsorbedType.o: ${bgp}/SolarAbsorbedType.F90 ${O}/shr_kind_mod.o ${O}/clm_varcon.o ${O}/clm_varpar.o ${O}/decompMod.o
	$(cmplr) -c ${bgp}/SolarAbsorbedType.F90
	mv SolarAbsorbedType.o ${O}

${O}/SurfaceAlbedoType.o: ${bgp}/SurfaceAlbedoType.F90 ${O}/shr_kind_mod.o ${O}/clm_varcon.o ${O}/clm_varpar.o ${O}/decompMod.o
	$(cmplr) -c ${bgp}/SurfaceAlbedoType.F90
	mv SurfaceAlbedoType.o ${O}

${O}/TemperatureType.o: ${bgp}/TemperatureType.F90 ${O}/shr_kind_mod.o ${O}/clm_varcon.o ${O}/clm_varpar.o ${O}/decompMod.o
	$(cmplr) -c ${bgp}/TemperatureType.F90
	mv TemperatureType.o ${O}

${O}/WaterFluxType.o: ${bgp}/WaterFluxType.F90 ${O}/shr_kind_mod.o ${O}/clm_varcon.o ${O}/clm_varpar.o ${O}/decompMod.o
	$(cmplr) -c ${bgp}/WaterFluxType.F90
	mv WaterFluxType.o ${O}

${O}/WaterStateType.o: ${bgp}/WaterStateType.F90 ${O}/shr_kind_mod.o ${O}/clm_varcon.o ${O}/clm_varpar.o ${O}/decompMod.o
	$(cmplr) -c ${bgp}/WaterStateType.F90
	mv WaterStateType.o ${O}

${O}/CanopyFluxesMultilayerType.o: ${can}/CanopyFluxesMultilayerType.F90 ${O}/shr_kind_mod.o ${O}/clm_varcon.o ${O}/clm_varpar.o ${O}/decompMod.o ${O}/histFileMod.o ${O}/restUtilMod.o ${O}/ncdio_pio.o ${O}/shr_file_mod.o
	$(cmplr) -c ${can}/CanopyFluxesMultilayerType.F90
	mv CanopyFluxesMultilayerType.o ${O}

${O}/filterMod.o: ${main}/filterMod.F90 ${O}/shr_kind_mod.o ${O}/clm_varcon.o
	$(cmplr) -c ${main}/filterMod.F90
	mv filterMod.o ${O}

# src1 files, which use mod files

${O}/clm_instMod.o: ${main}/clm_instMod.F90 $(mods)
	$(cmplr) -c ${main}/clm_instMod.F90
	mv clm_instMod.o ${O}

${O}/controlMod.o: ${drv}/controlMod.F90 $(mods)
	$(cmplr) -c ${drv}/controlMod.F90
	mv controlMod.o ${O}

${O}/SurfaceAlbedoMod.o: ${bgp}/SurfaceAlbedoMod.F90 $(mods)
	$(cmplr) -c ${bgp}/SurfaceAlbedoMod.F90
	mv SurfaceAlbedoMod.o ${O}

${O}/BandDiagonalMod.o: ${bgp}/BandDiagonalMod.F90 $(mods)
	$(cmplr) -c ${bgp}/BandDiagonalMod.F90
	mv BandDiagonalMod.o ${O}

${O}/SoilTemperatureMod.o: ${bgp}/SoilTemperatureMod.F90 $(mods)
	$(cmplr) -c ${bgp}/SoilTemperatureMod.F90
	mv SoilTemperatureMod.o ${O}

${O}/TowerMetMod.o: ${drv}/TowerMetMod.F90 $(mods)
	$(cmplr) -c ${drv}/TowerMetMod.F90 -I$(MOD_NETCDF)
	mv TowerMetMod.o ${O}

${O}/clmDataMod.o: ${drv}/clmDataMod.F90 $(mods)
	$(cmplr) -c ${drv}/clmDataMod.F90 -I$(MOD_NETCDF)
	mv clmDataMod.o ${O}

${O}/SurfaceResistanceMod.o: ${bgp}/SurfaceResistanceMod.F90 $(mods)
	$(cmplr) -c ${bgp}/SurfaceResistanceMod.F90
	mv SurfaceResistanceMod.o ${O}

${O}/MultibandSolarMod.o: ${can}/MultibandSolarMod.F90 $(mods)
	$(cmplr) -c ${can}/MultibandSolarMod.F90
	mv MultibandSolarMod.o ${O}

${O}/SolarRadiationMod.o: ${can}/SolarRadiationMod.F90 $(mods)
	$(cmplr) -c ${can}/SolarRadiationMod.F90
	mv SolarRadiationMod.o ${O}

${O}/CanopyTurbulenceMod.o: ${can}/CanopyTurbulenceMod.F90 $(mods)
	$(cmplr) -c ${can}/CanopyTurbulenceMod.F90
	mv CanopyTurbulenceMod.o ${O}

${O}/LeafBoundaryLayerMod.o: ${can}/LeafBoundaryLayerMod.F90 $(mods)
	$(cmplr) -c ${can}/LeafBoundaryLayerMod.F90
	mv LeafBoundaryLayerMod.o ${O}

${O}/LeafFluxesMod.o: ${can}/LeafFluxesMod.F90 $(mods) ${O}/StomataOptimizationMod.o \
	${O}/LeafPhotosynthesisMod.o ${O}/LeafTemperatureMod.o
	$(cmplr) -c ${can}/LeafFluxesMod.F90
	mv LeafFluxesMod.o ${O}

${O}/LeafTemperatureMod.o: ${can}/LeafTemperatureMod.F90 $(mods)
	$(cmplr) -c ${can}/LeafTemperatureMod.F90
	mv LeafTemperatureMod.o ${O}

${O}/LeafPhotosynthesisMod.o: ${can}/LeafPhotosynthesisMod.F90 $(mods)
	$(cmplr) -c ${can}/LeafPhotosynthesisMod.F90
	mv LeafPhotosynthesisMod.o ${O}

${O}/LongwaveRadiationMod.o: ${can}/LongwaveRadiationMod.F90 $(mods)
	$(cmplr) -c ${can}/LongwaveRadiationMod.F90
	mv LongwaveRadiationMod.o ${O}

${O}/PlantHydraulicsMod.o: ${can}/PlantHydraulicsMod.F90 $(mods)
	$(cmplr) -c ${can}/PlantHydraulicsMod.F90 -I$(MOD_NETCDF)
	mv PlantHydraulicsMod.o ${O}

${O}/SoilFluxesMultilayerMod.o: ${can}/SoilFluxesMultilayerMod.F90 $(mods)
	$(cmplr) -c ${can}/SoilFluxesMultilayerMod.F90
	mv SoilFluxesMultilayerMod.o ${O}

${O}/StomataOptimizationMod.o: ${can}/StomataOptimizationMod.F90 $(mods)
	$(cmplr) -c ${can}/StomataOptimizationMod.F90
	mv StomataOptimizationMod.o ${O}

${O}/CanopyWaterMod.o: ${can}/CanopyWaterMod.F90 $(mods)
	$(cmplr) -c ${can}/CanopyWaterMod.F90
	mv CanopyWaterMod.o ${O}

${O}/CanopyNitrogenProfileMod.o: ${can}/CanopyNitrogenProfileMod.F90 $(mods)
	$(cmplr) -c ${can}/CanopyNitrogenProfileMod.F90
	mv CanopyNitrogenProfileMod.o ${O}

# src2 files, which use mod files and src1 files

${O}/CanopyFluxesMultilayerMod.o: ${can}/CanopyFluxesMultilayerMod.F90 $(mods) $(src1)
	$(cmplr) -c ${can}/CanopyFluxesMultilayerMod.F90
	mv CanopyFluxesMultilayerMod.o ${O}

${O}/clm_driver.o: ${main}/clm_driver.F90 $(mods) $(src1) ${O}/CanopyFluxesMultilayerMod.o
	$(cmplr) -c ${main}/clm_driver.F90
	mv clm_driver.o ${O}

${O}/CLMml_initializeMod.o: ${drv}/CLMml_initializeMod.F90 $(mods)
	$(cmplr) -c ${drv}/CLMml_initializeMod.F90 -I$(MOD_NETCDF)
	mv CLMml_initializeMod.o ${O}

${O}/CLMml_driver.o: ${drv}/CLMml_driver.F90 $(mods) $(src1) ${O}/CLMml_initializeMod.o ${O}/clm_driver.o
	$(cmplr) -c ${drv}/CLMml_driver.F90
	mv CLMml_driver.o ${O}

${O}/CLMml.o: ${drv}/CLMml.F90 $(mods) $(src1) ${O}/CLMml_driver.o
	$(cmplr) -c ${drv}/CLMml.F90
	mv CLMml.o ${O}

# clean up object and mod files

clean:
	rm ${O}/*.o
	rm *.mod
