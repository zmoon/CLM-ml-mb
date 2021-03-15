"""
Run the model (only default setup/case so far)
"""
import contextlib
import os
import subprocess
from pathlib import Path


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
        subprocess.run(["make"])


def run():
    # Read nml to pass to the program
    with open(EXE_DIR / "namelists/nl.US-UMB.2006", "r") as f:
        s_nml = f.read()

    # Try to run
    with out_and_back(EXE_DIR):
        subprocess.run(["./prgm.exe"], input=s_nml, text=True)
