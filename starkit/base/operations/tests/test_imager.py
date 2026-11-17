"""Tests for the synthetic photometry operation."""

from types import SimpleNamespace

import numpy as np

from starkit.base.operations.imager import Photometry


def test_photometry_integrates_the_flux_through_each_filter():
    # Photometry cannot be constructed without wsynphot, so the integration
    # step is called with filter curves given on the grid wavelengths.
    wavelength = np.linspace(4000.0, 5000.0, 101)
    flux = np.full(wavelength.size, 2.0)
    filters = SimpleNamespace(
        grid_wavelength=wavelength,
        filter_transmission=np.array([np.ones(wavelength.size), np.full(wavelength.size, 0.5)]),
        wavelength_delta=np.array([1000.0, 1000.0]),
        zp_vega=np.array([1.0, 1.0]),
    )
    magnitudes = Photometry.evaluate_grid(filters, wavelength, flux)
    # Mean flux in the filters: 2.0 and 1.0 times the zero point.
    np.testing.assert_allclose(magnitudes, [-2.5 * np.log10(2.0), 0.0], atol=1e-12)
