"""
Run the model (only default setup/case so far)
"""
import contextlib
import os
import subprocess
from pathlib import Path

import f90nml


REPO_BASE = Path(__file__).parent
EXE_DIR = REPO_BASE / "exe"
OUT_DIR = REPO_BASE / "output"


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

