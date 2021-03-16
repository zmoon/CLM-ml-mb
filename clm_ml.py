"""
Run the model (only default setup/case so far)
"""
import contextlib
import os
import subprocess
from pathlib import Path

import f90nml
import numpy as np
import pandas as pd
import xarray as xr


REPO_BASE = Path(__file__).parent
EXE_DIR = REPO_BASE / "exe"
OUT_DIR = REPO_BASE / "output"

# See `driver/CLMml_driver.F90` for `write` calls and variable definitions
# nout1 - *_flux.out - canopy fluxes
FLUX_OUT_VARS = [
    # name, long_name, units (CF format)
    ('rn', 'Net radiation', 'W m-2'),
    ('stc', 'Canopy storage heat flux', 'W m-2'),
    ('sh', 'Sensible heat flux', 'W m-2'),
    ('lh', 'Latent heat flux', 'W m-2'),
    ('gpp', 'Gross primary production by vegetation', 'umol CO2 m-2 s-1'),
    ('ustar', 'Friction velocity', 'm s-1'),
    ('swup', 'Reflected solar radiation', 'W m-2'),
    ('ircan', 'Upward longwave radiation above canopy', 'W m-2'),
    ('taf', 'Air temperature at canopy top', 'K'),
    ('gs', 'Soil heat flux', 'W m-2'),
    ('rns', 'Ground net radiation', 'W m-2'),
    ('shs', 'Ground sensible heat flux', 'W m-2'),
    ('lhs', 'Ground latent heat flux', 'W m-2')
]

# nout2 - *_aux.out
AUX_OUT_VARS = [
    ('btran', 'Ball-Berry soil wetness factor', ''),
    ('lsc', 'Leaf-specific conductance of canopy layer', 'mml H2O (m2 leaf)-1 s-1 MPa-1'),
    ('psis', 'Weighted soil water potential', 'MPa'),
    ('lwptop', 'Leaf water potential at top-of-canopy', 'MPa'),
    ('lwpmid', 'Leaf water potential at mid-canopy', 'MPa'),
    ('fracminlwp', 'Fraction of canopy with lwp < minlwp', '1'),
]

# nout3 - *_profile.out - in-canopy profiles
PROFILE_OUT_VARS = [
    ('calday', 'Currently calendar day', ''),
    ('zs', 'Canopy height (for scalar conc./source)', 'm'),
    ('dpai', 'Layer plant area index', 'm2 leaf m-2'),
    ('rnl', 'Canopy layer net radiation', 'W m-2'),
    ('shl', 'Canopy layer sensible heat flux', 'W m-2'),
    ('lhl', 'Canopy layer latent heat flux', 'W m-2'),
    ('fcl', 'Canopy layer CO2 flux', 'umol CO2 m-2 s-1'),
    ('apar', 'Leaf absorbed PAR', 'umol photon (m2 leaf)-1 s-1'),
    ('gs', 'Leaf stomatal conductance', 'mol H2O (m2 leaf)-1 s-1'),
    ('lwp', 'Leaf water potential', 'MPa'),
    ('tveg', 'Vegetation temperature', 'K'),
    ('wind', 'Wind speed', 'm s-1'),
    ('dta', 'Air temperature difference from reference height', 'K'),
    ('ea', 'Air vapor pressure', 'kPa'),
    ('ras', 'Aerodynamic resistance for scalars', 's m-1'),
]

OUT_VARS = {
    "flux": FLUX_OUT_VARS,
    "aux": AUX_OUT_VARS,
    "profile": PROFILE_OUT_VARS,
}


def _load_out_ds_ts(which):
    # Load data and create series of times
    data = np.loadtxt(OUT_DIR / f"US-UMB_2006-07_{which}.out")
    assert data.shape[1] == len(OUT_VARS[which])
    t = pd.date_range('2006-07-01', '2006-08-01', freq="H", closed="right")
    assert data.shape[0] == t.size

    # Create and return ds
    return xr.Dataset(
        coords={
            "t": ("t", t),
        },
        data_vars={
            info[0]: ('t', col, {'long_name': info[1], 'units': info[2]})
            for info, col in zip(OUT_VARS[which], data.T)
        },
    )


def _load_out_ds_prof(which):
    # Load data and create series of times
    data = np.loadtxt(OUT_DIR / f"US-UMB_2006-07_{which}.out")
    assert data.shape[1] == len(OUT_VARS[which])
    t = pd.date_range('2006-07-01', '2006-08-01', freq="H", closed="right")
    nz = np.unique(data[:,1]).size
    nt = t.size
    assert data.shape[0] == nt * nz

    # Create and return ds
    z = data[:nz,1]
    return xr.Dataset(
        coords={
            "t": ("t", t),
            "z": ("z", z, {"long_name": "Height", "units": "m"}),
        },
        data_vars={
            info[0]: (('t', 'z'), col.reshape((nt, nz)), {'long_name': info[1], 'units': info[2]})
            for info, col in zip(OUT_VARS[which], data.T)
        },
    )


def load_out_ds(which="flux"):
    # Note: currently hard-coded to only work for 2006-07, but that's the only we can run anyway currently...
    
    if which not in OUT_VARS:
        raise ValueError("invalid `which`")

    if which in ("flux", "aux"):
        return _load_out_ds_ts(which)
    else:
        return _load_out_ds_prof(which)


@contextlib.contextmanager
def out_and_back(p):
    """Context manager: change working directory to path `p` but return to original cwd after."""
    cwd = os.getcwd()
    os.chdir(p)
    try:
        yield
    finally:
        os.chdir(cwd)


def build():
    """Build the model with `make`."""
    with out_and_back(EXE_DIR):
        subprocess.run(["make"], check=True)


def run(*, nsb=1):
    # Read default nml
    with open(EXE_DIR / "namelists/nl.US-UMB.2006", "r") as f:
        nml = f90nml.read(f)

    # Update nml with user settings
    nml["clm_inparm"]["nsb"] = nsb

    # Try to run
    s_nml = str(nml) + "\n"  # complains without newline at the end
    with out_and_back(EXE_DIR):
        subprocess.run(["./prgm.exe"], input=s_nml, text=True, check=True)
