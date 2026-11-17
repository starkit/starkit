"""Tests for the synthetic photometry operation."""

import sys
from types import ModuleType, SimpleNamespace

import h5py
import numpy as np
import pandas as pd
import pytest
from astropy import units as u
from scipy.integrate import trapezoid

from starkit.base.operations.imager import Photometry
from starkit.gridkit import load_grid

FLAM = u.erg / (u.s * u.cm ** 2 * u.angstrom)
WAVELENGTH = np.linspace(4000.0, 5000.0, 101)
# The values of wsynphot.io.cache_filters.DetectorType.
ENERGY_COUNTER, PHOTON_COUNTER = 0, 1


class _FilterCurve:
    """The attributes of wsynphot.FilterCurve that Photometry reads."""

    def __init__(self, filter_id, detector_type=ENERGY_COUNTER, zp_vega=1.0,
                 zp_ab=1.0):
        self.filter_id = filter_id
        self.wavelength = WAVELENGTH * u.angstrom
        self.transmission_lambda = np.ones(WAVELENGTH.size)
        self.detector_type = detector_type
        self.zp_vega_f_lambda = zp_vega * FLAM
        self.zp_ab_f_lambda = zp_ab * FLAM

    def calculate_wavelength_delta(self):
        # As in wsynphot, the integral of the transmission is weighted by the
        # wavelength for a photon counter.
        if self.detector_type == PHOTON_COUNTER:
            return trapezoid(self.transmission_lambda * self.wavelength,
                             self.wavelength)
        return trapezoid(self.transmission_lambda, self.wavelength)


FILTERS = {
    'Test/zp1': _FilterCurve('Test/zp1', zp_vega=1.0, zp_ab=4.0),
    'Test/zp2': _FilterCurve('Test/zp2', zp_vega=2.0, zp_ab=1.0),
    'Test/photon': _FilterCurve('Test/photon', detector_type=PHOTON_COUNTER),
}


class _FilterSet:
    """wsynphot.FilterSet, which looks filter names up in FILTERS."""

    def __init__(self, filter_set):
        if hasattr(filter_set[0], 'wavelength'):
            self.filter_set = list(filter_set)
        else:
            self.filter_set = [FILTERS[name] for name in filter_set]

    def __iter__(self):
        return iter(self.filter_set)

    def calculate_vega_magnitudes(self, spectrum):
        raise NotImplementedError('only the grid method is tested here')

    def calculate_ab_magnitudes(self, spectrum):
        raise NotImplementedError('only the grid method is tested here')


@pytest.fixture
def wsynphot(monkeypatch):
    """Photometry imports wsynphot, which the test dependencies lack."""
    module = ModuleType('wsynphot')
    module.FilterSet = _FilterSet
    monkeypatch.setitem(sys.modules, 'wsynphot', module)
    return module


def _load_grid(path, flux):
    """Write and load a grid whose four spectra are all ``flux``."""
    meta = pd.Series({'parameters': ['teff', 'logg'], 'grid_type': 'log',
                      'wavelength_unit': 'Angstrom',
                      'wavelength_type': 'vacuum',
                      'flux_unit': 'erg/s/cm^2/angstrom',
                      'R': 3000.0, 'R_sampling': 4, 'format_version': 'v1.0'})
    meta.to_hdf(path, key='meta')
    pd.DataFrame({'teff': [5000.0, 5000.0, 6000.0, 6000.0],
                  'logg': [4.0, 4.5, 4.0, 4.5]}).to_hdf(path, key='index')
    pd.DataFrame({0: WAVELENGTH}).to_hdf(path, key='wavelength')
    with h5py.File(path, 'a') as handle:
        handle['fluxes'] = np.tile(flux, (4, 1))
    return load_grid(path)


@pytest.fixture
def flat_grid(tmp_path):
    return _load_grid(tmp_path / 'grid.h5', np.full(WAVELENGTH.size, 2.0))


def test_photometry_integrates_the_flux_through_each_filter():
    # The integration step, called with filter curves given on the grid
    # wavelengths.
    wavelength = np.linspace(4000.0, 5000.0, 101)
    flux = np.full(wavelength.size, 2.0)
    filters = SimpleNamespace(
        grid_wavelength=wavelength,
        filter_transmission=np.array([np.ones(wavelength.size), np.full(wavelength.size, 0.5)]),
        wavelength_delta=np.array([1000.0, 1000.0]),
        zero_point=np.array([1.0, 1.0]),
    )
    magnitudes = Photometry.evaluate_grid(filters, wavelength, flux)
    # Mean flux in the filters: 2.0 and 1.0 times the zero point.
    np.testing.assert_allclose(magnitudes, [-2.5 * np.log10(2.0), 0.0], atol=1e-12)


def test_photometry_from_grid_gives_magnitudes(wsynphot, flat_grid):
    filter_set = _FilterSet([FILTERS['Test/zp1'], FILTERS['Test/zp2']])
    model = flat_grid | Photometry.from_grid(flat_grid, filter_set=filter_set)
    # Mean flux of 2 in both filters, with Vega zero points of 1 and 2.
    np.testing.assert_allclose(model(), [-2.5 * np.log10(2.0), 0.0],
                               atol=1e-12)


def test_photometry_from_grid_accepts_filter_names(wsynphot, flat_grid):
    model = flat_grid | Photometry.from_grid(
        flat_grid, filter_set=['Test/zp1', 'Test/zp2'])
    np.testing.assert_allclose(model(), [-2.5 * np.log10(2.0), 0.0],
                               atol=1e-12)


def test_photometry_from_grid_uses_ab_zero_points(wsynphot, flat_grid):
    filter_set = _FilterSet([FILTERS['Test/zp1'], FILTERS['Test/zp2']])
    model = flat_grid | Photometry.from_grid(
        flat_grid, filter_set=filter_set, mag_type='ab')
    # Mean flux of 2 in both filters, with AB zero points of 4 and 1.
    np.testing.assert_allclose(
        model(), [-2.5 * np.log10(0.5), -2.5 * np.log10(2.0)], atol=1e-12)


def test_photometry_from_grid_weights_photon_counts_by_wavelength(
        wsynphot, tmp_path):
    # f_lambda = lambda / 1000 Angstrom.
    grid = _load_grid(tmp_path / 'grid.h5', WAVELENGTH / 1000.0)
    filter_set = _FilterSet([FILTERS['Test/zp1'], FILTERS['Test/photon']])
    model = grid | Photometry.from_grid(grid, filter_set=filter_set)
    # Energy counter: mean of f_lambda over 4000 to 5000 Angstrom, 4.5.
    # Photon counter: int f_lambda lambda dlambda / int lambda dlambda
    # = ((5000^3 - 4000^3) / 3000) / ((5000^2 - 4000^2) / 2) = 4.5185.
    photon_mean = ((5000.0 ** 3 - 4000.0 ** 3) / 3000.0
                   / ((5000.0 ** 2 - 4000.0 ** 2) / 2.0))
    # The trapezoid rule overestimates the integral of lambda^2 by 1e-6.
    np.testing.assert_allclose(
        model(), [-2.5 * np.log10(4.5), -2.5 * np.log10(photon_mean)],
        atol=1e-5)


def test_photometry_rejects_unknown_magnitude_types(wsynphot):
    with pytest.raises(ValueError, match="'vega' or 'ab'"):
        Photometry(filter_set=['Test/zp1'], mag_type='st')

