# CLM-ml-mb
Running Gordon Bonan's CLM-ml, with
* more options
* Python
* multiple sub-bands per waveband in the solar RT (hence `-mb`)

Note
* Bonan has the source code for the multilayer canopy part (here in [`f/src/canopy`](f/src/canopy)) at [gbonan/CLM-ml_v0](https://github.com/gbonan/CLM-ml_v0),
  which goes with a 2018 GMD paper: doi:[10.5194/gmd-11-1467-2018](https://doi.org/10.5194/gmd-11-1467-2018)


## Dependencies

* Cloning the repo
  - [Git LFS](https://git-lfs.github.com/) (to download the `.nc` files)
* Build/compile
  - build system: [Meson](https://mesonbuild.com/) + [Ninja](https://ninja-build.org/)
  - compiler: Fortran compiler (currently only GNU `gfortran` is supported in [the `meson.build`](f/src/meson.build))
  - libraries: NetCDF-Fortran, LAPACK, BLAS
* Running with Python
  - see [the CI Python requirements](.github/workflows/ci_requirements.txt)

### Ubuntu

(or other Linux with `apt`)

```bash
sudo apt install gfortran libnetcdf-dev libnetcdff-dev libblas-dev liblapack-dev meson ninja-build
```

### Windows (MSYS2)

The `meson.build` assumes that `C:/msys64` is the location of the [MSYS2](https://www.msys2.org/) installation.

1. To get `gfortran`(`.exe`) and most of the Fortran dependencies, open MSYS2 and:
   ```bash
   pacman -S mingw64/mingw-w64-x86_64-gcc-fortran mingw64/mingw-w64-x86_64-netcdf mingw64/mingw-w64-x86_64-lapack
   ```
   That LAPACK includes `libblas`.

2. [Download](https://www.unidata.ucar.edu/downloads/netcdf/) NetCDF-Fortran and [build from source](https://www.unidata.ucar.edu/software/netcdf/docs/building_netcdf_fortran.html) with `NCDIR=/mingw64`, `FC=gfortran`, `NFDIR=$NCDIR`.

3. Install Meson and Ninja to Windows. For example, into the Conda env.

For me, the model runs faster on [WSL](https://docs.microsoft.com/en-us/windows/wsl/about), though.


<!-- TODO: Conda, Brew  -->


## Building the model

1. Navigate into `f/src` and invoke `meson setup build`
2. Enter `build` and invoke `meson compile` or `ninja`

Or, without leaving the repo root:
```bash
meson setup f/src/build f/src && meson compile -C f/src/build
```


## Original Bonan notes

*Provided by Bonan in a Word document along with the driver code, modified slightly by me.*

The code is configured to run for one particular tower site (US-UMB) for one year (2006) and one month (July) and has the necessary meteorological forcing and CLM datasets needed to run the model. The directory `output` contains model output for this simulation. These files can be used to verify the model is running correctly.

The model is compiled and run from the directory `exe`.

Executing `make` within directory `exe` creates the executable `prgm.exe`. The code is run using:
```bash
./prgm.exe < namelists/nl.US-UMB.2006
```

Note that the directory `tower-fluxes` is not required. This is just an example of tower fluxes (derived from measurements) used to test the model.

`driver/CLMml_driver.F90` – This is the main time-stepping driver for the model. The model processes one tower site for a specified month and for a specified year.

`driver/controlMod.F90` – This routine sets run-time options used with the model. Model options are controlled by a namelist file (required). **However, the namelist is not up-to-date. This routine overwrites and hardwires model options.** There are multiple options related to solar radiation, stomatal conductance, photosynthesis, and other things. These were used during model development and testing. The default values are set in this routine. Other options should only be used after verifying they work. This routine also sets three directories for input/output:

1. tower meteorological forcing: `diratm = '../tower-forcing/'`
2. CLM4.5 forcing: `dirclm = '../clm4_5/'`
3. Model ouput: `dirout = '../output/'`

Four important variables set in the namelist file that control the tower data are:

1. Tower site: `tower = 'US-UMB'`
2. Tower year: `tower_yrbeg = 2006`
3. Tower year: `tower_yrend = 2006`
4. Tower month: `tower_month = 7`

`driver/TowerDataMod.F90` – This sets parameters associated with each tower site.
