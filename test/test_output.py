"""
Test the output files -- generation (running the model), loading, etc.
"""
import pytest

import run
import utils


@pytest.mark.slow
def test_run_output_reg():
    """We test the outputs with default settings (1 sub-band per waveband)
    against original provided by Bonan as a regression test.
    """
    run.build()  # in case of changes
    run.run()  # takes ~ 10 s for the month
    utils.compare_fouts_to_orig()


@pytest.mark.slow
@pytest.mark.xfail(strict=True)
def test_run_output_sb():
    """With 2 sub-bands, the output will be different."""
    run.build()
    run.run(nsb=2)
    utils.compare_fouts_to_orig()


@pytest.mark.parametrize(
    "which",
    ["flux", "aux", "profile"]
)
def test_load_out_ds(which):
    run.load_out_ds(which)
