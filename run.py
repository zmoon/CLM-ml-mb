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

# See `driver/CLMml_driver.F90` for the definitions
# nout1 - *_flux.out - canopy fluxes
FLUX_OUT_VARS = [
    # name, long_name, units (CF format)
    ('rnet', 'Net radiation', 'W m-2'),
    ('stflx', 'Canopy storage heat flux', 'W m-2'),
    ('shflx', 'Sensible heat flux', 'W m-2'),
    ('lhflx', 'Latent heat flux', 'W m-2'),
    ('gppveg', 'Gross primary production', 'umol CO2 m-2 s-1'),
    ('ustar', 'Friction velocity', 'm s-1'),
    ('swup', 'Reflected solar radiation', 'W m-2'),
    ('ircan', 'Upward longwave radiation above canopy', 'W m-2'),
    ('taf', 'Air temperature at canopy top', 'K'),
    ('gsoi', 'Soil heat flux', 'W m-2'),
    ('rnsoi', 'Ground net radiation', 'W m-2'),
    ('shsoi', 'Ground sensible heat flux', 'W m-2'),
    ('lhsoi', 'Ground latent heat flux', 'W m-2')
]

# nout2 - *_aux.out
AUX_OUT_VARS = [
    ('btran', 'Ball-Berry soil wetness factor', ''),
    ('lsc', 'Leaf-specific conductance of canopy layer', 'mml H2O (m2 leaf)-1 s-1 MPa-1'),
    ('psis', 'Weighted soil water potential', 'MPa'),
    ('lwp', 'Leaf water potential at top-of-canopy', 'MPa'),
    ('lwp', 'Leaf water potential at mid-canopy', 'MPa'),
    ('fracminlwp', 'Fraction of canopy with lwp < minlwp', '1'),
]

# nout3 - *_profile.out - in-canopy profiles

OUT_VARS = {
    "flux": FLUX_OUT_VARS,
    "aux": AUX_OUT_VARS,
}


def load_out_ds(which="flux"):
    # Currently hard-coded to only work for 2006-07, but that's the only we can run anyway currently...
    
    if which not in ("flux", "aux"):
        raise NotImplementedError

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

