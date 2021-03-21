project('CLM-ml-zm', 'fortran',
  default_options: [
    # Compiler-agnostic dialect settings supported by Meson
    'fortran_std=f2008',
    'warning_level=0',  # 1 for -Wall
  ]
)

# Compiler-specific dialect options
fc = meson.get_compiler('fortran')
if fc.get_id() == 'gcc'
  dialect = ['-ffree-line-length-none', '-fno-range-check']
else
  error('not supporting other compilers at this time...')
endif

# Detect if running on ACI
hostname = run_command('hostname').stdout().strip()
on_aci = hostname.contains('aci.ics.psu.edu')
if on_aci
  message('PSU ACI detected')
endif

# Find the location of `netcdf.inc` (needed to compile some of the objects)
r = run_command('nf-config', '--includedir')
if r.returncode() != 0
  error(
    'trying to get the netcdf-fortran include dir with `nf-config` failed. ' +
    'Set `f_inc_netcdf` manually in `meson.build`?'
  )
endif
f_inc_netcdf = '-I' + r.stdout().strip()

# Set options for compiling the objects
dialect += f_inc_netcdf
add_project_arguments(dialect, language: 'fortran')

# Define sources
# Note:
# * these can instead be in a `meson.build` of each subdir and loaded with `subdir()`
# * these are in alphabetical order (`ls` order)--Meson figures out the correct build order
src_shr = [
  'csm_share/shr_file_mod.F90',
  'csm_share/shr_kind_mod.F90',
  'csm_share/shr_orb_mod.F90',
]
src_main = [
  'main/ColumnType.F90',
  'main/PatchType.F90',
  'main/abortutils.F90',
  'main/atm2lndType.F90',
  'main/clm_driver.F90',
  'main/clm_varcon.F90',
  'main/clm_varctl.F90',
  'main/clm_varpar.F90',
  'main/clm_instMod.F90',
  'main/decompMod.F90',
  'main/filterMod.F90',
  'main/histFileMod.F90',
  'main/ncdio_pio.F90',
  'main/pftconMod.F90',
]
src_util = [
  'utils/clm_time_manager.F90',
  'utils/clm_varorb.F90',
  'utils/restUtilMod.F90',
]
src_bgp = [
  'biogeophys/BandDiagonalMod.F90',
  'biogeophys/CanopyStateType.F90',
  'biogeophys/EnergyFluxType.F90',
  'biogeophys/FrictionVelocityMod.F90',
  'biogeophys/SoilStateType.F90',
  'biogeophys/SoilTemperatureMod.F90',
  'biogeophys/SolarAbsorbedType.F90',
  'biogeophys/SurfaceAlbedoMod.F90',
  'biogeophys/SurfaceAlbedoType.F90',
  'biogeophys/SurfaceResistanceMod.F90',
  'biogeophys/TemperatureType.F90',
  'biogeophys/WaterFluxType.F90',
  'biogeophys/WaterStateType.F90',
]
src_mb = [
  'canopy/multiband/MultibandSolarDataMod.f90',
  'canopy/multiband/MultibandSolarMod.f90',
]
src_canopy = [
  'canopy/CanopyFluxesMultilayerMod.F90',
  'canopy/CanopyFluxesMultilayerType.F90',
  'canopy/CanopyNitrogenProfileMod.F90',
  'canopy/CanopyTurbulenceMod.F90',
  'canopy/CanopyWaterMod.F90',
  'canopy/LeafBoundaryLayerMod.F90',
  'canopy/LeafFluxesMod.F90',
  'canopy/LeafPhotosynthesisMod.F90',
  'canopy/LeafTemperatureMod.F90',
  'canopy/LongwaveRadiationMod.F90',
  'canopy/MathToolsMod.F90',
  'canopy/PlantHydraulicsMod.F90',
  'canopy/SoilFluxesMultilayerMod.F90',
  'canopy/SolarRadiationMod.F90',
  'canopy/StomataOptimizationMod.F90',
  'canopy/WaterVaporMod.F90',
]
src_drv = [
  'driver/CLMml.F90',
  'driver/CLMml_driver.F90',
  'driver/CLMml_initializeMod.F90',
  'driver/SoilTexMod.F90',
  'driver/TowerDataMod.F90',
  'driver/TowerMetMod.F90',
  'driver/clmDataMod.F90',
  'driver/controlMod.F90',
]
src_all = src_shr + src_main + src_util + src_bgp + src_canopy + src_drv + src_mb

# Create clm-ml executable
deps = [
  # dependency('netcdf', language: 'fortran', required: true),
  fc.find_library('netcdff'),  # ^ above not working
  fc.find_library('blas'),
  fc.find_library('lapack'),  # TODO: look for MKL first?
  fc.find_library('m'),  # doesn't seem to actually be used
]
if on_aci
  # ACI's libnetcdff needs more supporting libs (not all covered by `nf-config --flibs`)
  # Note: presumably could use `find_library` for this too
  link_args = ['-lhdf5', '-lhdf5_hl', '-lsz', '-lnetcdf']
else
  link_args = []
endif
executable('clm-ml',
  sources: files(src_all),
  link_args: link_args,
  dependencies: deps,
)

# Create multiband test program and mark it as a test
test_mb_exe = executable('test_mb', src_mb + ['canopy/multiband/test_MultibandSolarMod.f90'])
test('test_mb', test_mb_exe)
