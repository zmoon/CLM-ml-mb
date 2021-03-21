import datetime

import numpy as np

import clm_ml


def get_fouts(validate=True):
    """Return list of the output file paths, doing some validation"""
    fouts = list(clm_ml.F_OUTPUT_DIR.glob("*.out"))
    if validate:  # make sure we have the right number and that they are all recent
        now = datetime.datetime.now()
        assert len(fouts) == 3
        for fout in fouts:
            mtime = datetime.datetime.fromtimestamp(fout.stat().st_mtime)
            assert (mtime - now).total_seconds() < 2

    return fouts


def compare_fouts_to_orig(fouts=None):
    """Compare values to original provided by Bonan"""
    if not fouts:
        fouts = get_fouts(validate=True)

    for fout in fouts:
        bn = fout.name
        new = np.loadtxt(fout)
        orig = np.loadtxt(clm_ml.REPO_BASE / "test/data/output" / bn)
        np.testing.assert_equal(new, orig)
