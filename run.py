"""
Run the model (only default setup/case so far)
"""
import os
import subprocess
from pathlib import Path


REPO_BASE = Path(__file__).parent
EXE_DIR = REPO_BASE / "exe"
OUT_DIR = REPO_BASE / "output"


def run():
    # Read nml to pass to the program
    with open(EXE_DIR / "namelists/nl.US-UMB.2006", "r") as f:
        s_nml = f.read()

    # Try to run
    cwd = os.getcwd()
    os.chdir(EXE_DIR)
    try:
        subprocess.run(["./prgm.exe"], input=s_nml, text=True)
    except Exception:
        raise
    finally:
        os.chdir(cwd)

