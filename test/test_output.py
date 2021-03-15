"""
Regression test -- check the output against original provided by Bonan
"""
import datetime
from pathlib import Path

import numpy as np

import run


def test_run_output_reg():
    # Try to run the model (takes ~ 10 s for the month)
    run.build()  # in case of changes
    run.run()

    # Check time stamps of the output files
    now = datetime.datetime.now()
    fouts = list(run.OUT_DIR.glob("*.out"))
    assert len(fouts) == 3
    for fout in fouts:
        mtime = datetime.datetime.fromtimestamp(fout.stat().st_mtime)
        assert (mtime - now).total_seconds() < 2

    # Compare values to original provided by Bonan
    for fout in fouts:
        bn = fout.name
        new = np.loadtxt(fout)
        orig = np.loadtxt(run.REPO_BASE / "test/data/output" / bn)
        np.testing.assert_equal(new, orig)

