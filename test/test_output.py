"""
Test the output files -- generation (running the model), loading, etc.
"""
import pytest

import clm_ml
import utils  # noreorder


@pytest.mark.slow
def test_run_output_reg():
    """We test the outputs with default settings (1 sub-band per waveband)
    against original provided by Bonan as a regression test.
    """
    clm_ml.build()  # in case of changes
    clm_ml.run()  # takes ~ 10 s for the month
    utils.compare_fouts_to_orig()


@pytest.mark.slow
def test_run_output_sb():
    """With 2 sub-bands, the output will be different."""
    clm_ml.build()
    clm_ml.run(nsb=2)
    with pytest.raises(AssertionError):  # "Arrays are not equal"
        utils.compare_fouts_to_orig()


@pytest.mark.parametrize("which", ["flux", "aux", "profile"])
def test_load_out_ds(which):
    clm_ml.load_out_ds(which)
