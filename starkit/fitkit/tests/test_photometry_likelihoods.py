"""Tests for the likelihoods that compare synthetic and observed magnitudes."""

from types import SimpleNamespace

import numpy as np
import pytest

from starkit.fitkit.likelihoods import (PhotometryColorLikelihood,
                                        RelativePhotometryLikelihood)

MAGNITUDES = [10.0, 9.5, 9.2]
UNCERTAINTIES = [0.1, 0.1, 0.2]
SYNTHETIC = np.array([1.0, 0.4, 0.2])
# Colors do not depend on a common offset of the magnitudes, so their chi^2
# with the full covariance equals the magnitude chi^2 minimized over the
# offset k. Here d = observed - synthetic = (9.0, 9.1, 9.0) with weights
# (100, 100, 25), k = 2035 / 225, and chi^2 = (400 + 625 + 100) / 2025 = 5 / 9.
COLOR_LOGLIKELIHOOD = -5.0 / 18.0


def _minimum_magnitude_chi2(observed, synthetic, uncertainties):
    difference = np.asarray(observed) - np.asarray(synthetic)
    weights = 1.0 / np.asarray(uncertainties) ** 2
    offset = np.sum(weights * difference) / np.sum(weights)
    return np.sum(weights * (difference - offset) ** 2)


def test_color_likelihood_uses_the_covariance_of_adjacent_colors():
    # Adjacent colors share a magnitude, so their errors are correlated.
    observed = [16.2, 15.1, 14.8, 14.75]
    uncertainties = [0.1, 0.05, 0.2, 0.03]
    synthetic = np.array([1.3, 0.25, 0.0, 0.0])
    likelihood = PhotometryColorLikelihood(observed, uncertainties)
    assert likelihood(synthetic) == pytest.approx(
        -0.5 * _minimum_magnitude_chi2(observed, synthetic, uncertainties))


def test_color_likelihood_takes_magnitudes_and_uncertainties():
    likelihood = PhotometryColorLikelihood(MAGNITUDES, UNCERTAINTIES)
    assert likelihood(SYNTHETIC) == pytest.approx(COLOR_LOGLIKELIHOOD)


def test_color_likelihood_takes_a_magnitude_set():
    # wsynphot.MagnitudeSet has these two attributes.
    magnitude_set = SimpleNamespace(magnitudes=tuple(MAGNITUDES),
                                    magnitude_uncertainties=tuple(UNCERTAINTIES))
    likelihood = PhotometryColorLikelihood(magnitude_set)
    assert likelihood(SYNTHETIC) == pytest.approx(COLOR_LOGLIKELIHOOD)


def test_relative_likelihood_takes_magnitudes_and_uncertainties():
    likelihood = RelativePhotometryLikelihood([10.0, 9.5], [0.1, 0.2])
    # -0.5 ((0.1 / 0.1)^2 + (0.2 / 0.2)^2)
    assert likelihood(np.array([10.1, 9.3])) == pytest.approx(-1.0)


def test_likelihood_needs_uncertainties():
    with pytest.raises(TypeError, match='magnitude_uncertainties'):
        RelativePhotometryLikelihood(MAGNITUDES)
