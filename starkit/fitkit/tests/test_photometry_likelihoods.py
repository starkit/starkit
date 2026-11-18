"""Tests for the likelihoods that compare synthetic and observed magnitudes."""

from types import SimpleNamespace

import numpy as np
import pytest

from starkit.fitkit.likelihoods import (PhotometryColorLikelihood,
                                        RelativePhotometryLikelihood)

MAGNITUDES = [10.0, 9.5, 9.2]
UNCERTAINTIES = [0.1, 0.1, 0.2]
SYNTHETIC = np.array([1.0, 0.4, 0.2])
# Observed colors 0.5 and 0.3, synthetic 0.6 and 0.2, with variances 0.02
# and 0.05: -0.5 (0.1^2 / 0.02 + 0.1^2 / 0.05).
COLOR_LOGLIKELIHOOD = -0.35


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
