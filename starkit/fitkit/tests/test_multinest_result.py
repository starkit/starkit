"""Regression tests for reading MultiNest's whitespace-separated output files."""

import numpy as np
import pytest

try:
    from starkit.fitkit.multinest.base import MultiNestResult
except ImportError:
    pytest.skip("needs the multinest extra and the MultiNest library", allow_module_level=True)

PARAMETERS = ["teff_0", "logg_0"]


def _multinest_line(*values):
    # MultiNest writes each value as %28.18E, so every line starts with spaces.
    return "".join("{0:28.18E}".format(v) for v in values) + "\n"


def test_reads_weighted_posterior(tmp_path):
    basename = str(tmp_path / "fit")
    with open(basename + "_.txt", "w") as fh:
        fh.write(_multinest_line(0.25, 1500.0, 4600.0, 2.25))
        fh.write(_multinest_line(0.75, 1502.5, 4650.0, 2.5))

    result = MultiNestResult.from_multinest_basename(basename, PARAMETERS)

    data = result.posterior_data
    assert list(data.columns) == ["weights", "loglikelihood"] + PARAMETERS
    np.testing.assert_array_equal(data["weights"], [0.25, 0.75])
    np.testing.assert_array_equal(data["loglikelihood"], [1500.0, 1502.5])
    np.testing.assert_array_equal(data["teff_0"], [4600.0, 4650.0])
    np.testing.assert_array_equal(data["logg_0"], [2.25, 2.5])
    assert result.parameter_names == PARAMETERS


def test_reads_equally_weighted_posterior(tmp_path):
    basename = str(tmp_path / "fit")
    with open(basename + "_post_equal_weights.dat", "w") as fh:
        for teff, logg, loglike in [(4600.0, 2.25, -750.0), (4650.0, 2.5, -751.25),
                                    (4700.0, 2.75, -752.5), (4750.0, 3.0, -753.75)]:
            fh.write(_multinest_line(teff, logg, loglike))

    result = MultiNestResult.from_multinest_basename(basename, PARAMETERS, equal_weights=True)

    data = result.posterior_data
    np.testing.assert_array_equal(data["teff_0"], [4600.0, 4650.0, 4700.0, 4750.0])
    np.testing.assert_array_equal(data["logg_0"], [2.25, 2.5, 2.75, 3.0])
    np.testing.assert_array_equal(data["weights"], [0.25, 0.25, 0.25, 0.25])
    assert result.parameter_names == PARAMETERS
