"""Tests for assemble_model."""

import numpy as np
from astropy import units as u

from starkit.base.assemble_model import assemble_model
from starkit.base.operations.imager import Photometry
from starkit.filters import Filter
from starkit.fix_spectrum1d import SKSpectrum1D

FLAM = u.erg / (u.s * u.cm ** 2 * u.angstrom)
WAVELENGTH = np.linspace(3500.0, 5500.0, 2001)
BOX_WAVELENGTH = np.linspace(4000.0, 5000.0, 101)
BOX = Filter(BOX_WAVELENGTH, np.ones(BOX_WAVELENGTH.size), name='box')


def _grid(make_grid):
    # f_lambda = lambda / 1000 Angstrom.
    return make_grid(WAVELENGTH / 1000.0, WAVELENGTH)


def test_assemble_model_with_only_the_grid(make_grid):
    grid = _grid(make_grid)
    wavelength, flux = assemble_model(grid)()
    expected_wavelength, expected_flux = grid()
    np.testing.assert_array_equal(wavelength, expected_wavelength)
    np.testing.assert_array_equal(flux, expected_flux)


def test_assemble_model_with_only_filters(make_grid):
    grid = _grid(make_grid)
    model = assemble_model(grid, filter_set=[BOX], mag_type='ab')
    expected = grid | Photometry.from_grid(grid, filter_set=[BOX], mag_type='ab')
    assert model.param_names == ('teff_0', 'logg_0')
    np.testing.assert_array_equal(model(), expected())


def test_assemble_model_with_filters_and_extinction(make_grid):
    grid = _grid(make_grid)
    model = assemble_model(grid, filter_set=[BOX], mag_type='ab', a_v=0.5,
                           r_v=4.0)
    assert model.param_names == ('teff_0', 'logg_0', 'a_v_1', 'r_v_1')
    assert (model.a_v_1.value, model.r_v_1.value) == (0.5, 4.0)


def test_assemble_model_with_a_spectrum_and_no_stellar_operation(make_grid):
    grid = _grid(make_grid)
    observed_wavelength = WAVELENGTH[(WAVELENGTH > 4000) & (WAVELENGTH < 5000)]
    observed_flux = observed_wavelength / 1000.0
    spectrum = SKSpectrum1D.from_array(
        observed_wavelength * u.angstrom, observed_flux * FLAM,
        uncertainty=0.01 * observed_flux * FLAM)
    wavelength, flux = assemble_model(grid, spectrum=spectrum,
                                      normalize_npol=1)()
    np.testing.assert_array_equal(wavelength, observed_wavelength)
    # The model is the observed spectrum, so the normalization leaves it.
    np.testing.assert_allclose(flux, observed_flux, rtol=1e-10)
