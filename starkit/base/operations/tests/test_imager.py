"""Tests for the synthetic photometry operation."""

import sys
from types import SimpleNamespace

import numpy as np
import pytest
from astropy import constants as const
from astropy import units as u

from starkit import filters
from starkit.base.operations.imager import Photometry
from starkit.filters import Filter

FLAM = u.erg / (u.s * u.cm ** 2 * u.angstrom)
WAVELENGTH = np.linspace(3500.0, 5500.0, 2001)
BOX_WAVELENGTH = np.linspace(4000.0, 5000.0, 101)
PHOTON_BOX = Filter(BOX_WAVELENGTH, np.ones(BOX_WAVELENGTH.size), name='photon box')
ENERGY_BOX = Filter(BOX_WAVELENGTH, np.ones(BOX_WAVELENGTH.size),
                    detector='energy', name='energy box')
# Mean of f_lambda = lambda / 1000 Angstrom over 4000 to 5000 Angstrom, for
# energy integration and for photon counting:
# ((5000^3 - 4000^3) / 3000) / ((5000^2 - 4000^2) / 2).
ENERGY_MEAN = 4.5
PHOTON_MEAN = ((5000.0 ** 3 - 4000.0 ** 3) / 3000.0
               / ((5000.0 ** 2 - 4000.0 ** 2) / 2.0))


@pytest.fixture
def linear_grid(make_grid):
    # f_lambda = lambda / 1000 Angstrom.
    return make_grid(WAVELENGTH / 1000.0, WAVELENGTH)


def test_photometry_from_grid_gives_vega_magnitudes(flat_vega, linear_grid):
    model = linear_grid | Photometry.from_grid(
        linear_grid, filter_set=[PHOTON_BOX, ENERGY_BOX])
    # Vega has f_lambda = 1, so the zero points are 1.
    np.testing.assert_allclose(
        model(), [-2.5 * np.log10(PHOTON_MEAN), -2.5 * np.log10(ENERGY_MEAN)],
        atol=1e-5)


def test_photometry_without_grid_wavelengths_uses_the_spectrum(flat_vega):
    photometry = Photometry(filter_set=[PHOTON_BOX, ENERGY_BOX])
    np.testing.assert_allclose(
        photometry(WAVELENGTH, WAVELENGTH / 1000.0),
        [-2.5 * np.log10(PHOTON_MEAN), -2.5 * np.log10(ENERGY_MEAN)],
        atol=1e-5)


def test_photometry_needs_no_synphot(monkeypatch, make_grid):
    # A None entry in sys.modules makes "import synphot" fail as if the
    # package were not installed.
    monkeypatch.setitem(sys.modules, 'synphot', None)
    filters.vega_spectrum.cache_clear()
    try:
        f_nu_3631 = (3631 * u.Jy * const.c
                     / (WAVELENGTH * u.angstrom) ** 2).to_value(FLAM)
        ab_grid = make_grid(f_nu_3631, WAVELENGTH, name='ab.h5')
        ab = ab_grid | Photometry.from_grid(
            ab_grid, filter_set=[PHOTON_BOX, ENERGY_BOX], mag_type='ab')
        vega_wavelength, vega_flux = filters.vega_spectrum()
        vega_grid = make_grid(
            np.interp(WAVELENGTH, vega_wavelength, vega_flux), WAVELENGTH,
            name='vega.h5')
        vega = vega_grid | Photometry.from_grid(
            vega_grid, filter_set=[PHOTON_BOX, ENERGY_BOX], mag_type='vega')
        # On the 1 Angstrom grid the vertical edges of the box become
        # 1 Angstrom ramps, which changes the magnitudes by 3e-5.
        np.testing.assert_allclose(ab(), [0.0, 0.0], atol=1e-4)
        # Vega has magnitude 0 in the Vega system.
        np.testing.assert_allclose(vega(), [0.0, 0.0], atol=1e-4)
    finally:
        filters.vega_spectrum.cache_clear()


def test_photometry_accepts_svo_ids_and_wsynphot_filters(svo, flat_vega):
    wsynphot_curve = SimpleNamespace(
        wavelength=BOX_WAVELENGTH * u.angstrom,
        transmission_lambda=np.ones(BOX_WAVELENGTH.size), detector_type=0,
        filter_id='Test/Box.E')
    photometry = Photometry(filter_set=['Generic/Johnson.V', wsynphot_curve,
                                        PHOTON_BOX])
    assert [(item.name, item.detector) for item in photometry.filters] == [
        ('Generic/Johnson.V', 'energy'), ('Test/Box.E', 'energy'),
        ('photon box', 'photon')]


def test_photometry_rejects_unknown_magnitude_types():
    with pytest.raises(ValueError, match="'vega' or 'ab'"):
        Photometry(filter_set=[PHOTON_BOX], mag_type='st')


def test_filters_must_lie_within_the_grid(flat_vega, make_grid):
    wavelength = np.linspace(4500.0, 5500.0, 1001)
    grid = make_grid(np.ones(wavelength.size), wavelength)
    with pytest.raises(ValueError, match='photon box'):
        Photometry.from_grid(grid, filter_set=[PHOTON_BOX])
