"""Tests of the comparison of grid fluxes with a blackbody."""

import h5py
import numpy as np
import pandas as pd
import pytest
from scipy.integrate import trapezoid
from astropy import constants as const

from starkit.gridkit.flux_check import (blackbody_surface_flux,
                                        check_flux_scale,
                                        format_flux_scale_check)

WAVELENGTH = np.geomspace(3000.0, 9000.0, 2000)
TEFFS = [4000.0, 4000.0, 6000.0, 6000.0, 9000.0]
LOGGS = [4.0, 4.5, 4.0, 4.5, 4.0]


def _write_grid(path, fluxes, flux_unit):
    meta = pd.Series({'parameters': ['teff', 'logg'], 'grid_type': 'log',
                      'wavelength_unit': 'Angstrom',
                      'wavelength_type': 'vacuum', 'flux_unit': flux_unit,
                      'R': 3000.0, 'R_sampling': 4, 'format_version': 'v1.0'})
    meta.to_hdf(path, key='meta')
    pd.DataFrame({'teff': TEFFS, 'logg': LOGGS}).to_hdf(path, key='index')
    pd.DataFrame({0: WAVELENGTH}).to_hdf(path, key='wavelength')
    with h5py.File(path, 'a') as handle:
        handle['fluxes'] = np.asarray(fluxes)


def _surface_fluxes():
    return [blackbody_surface_flux(WAVELENGTH, teff) for teff in TEFFS]


def test_blackbody_surface_flux_integrates_to_sigma_t4():
    wavelength = np.geomspace(500.0, 1e7, 200000)
    total = trapezoid(blackbody_surface_flux(wavelength, 5772.0), wavelength)
    expected = const.sigma_sb.cgs.value * 5772.0 ** 4
    assert total == pytest.approx(expected, rel=1e-3)


def test_surface_flux_with_surface_flux_label(tmp_path):
    path = tmp_path / 'grid.h5'
    _write_grid(path, _surface_fluxes(), 'erg/s/cm^2/angstrom')
    check = check_flux_scale(path)
    assert [model.parameters['teff'] for model in check.models] == [4000.0, 6000.0, 9000.0]
    assert all(model.ratio == pytest.approx(1.0, rel=1e-9) for model in check.models)
    assert check.interpretation == 'surface flux'
    assert check.consistent_with_label is True


def test_surface_flux_labeled_as_luminosity_density(tmp_path):
    # The labeling of StarKit grids written between December 2017 and
    # January 2025.
    path = tmp_path / 'grid.h5'
    _write_grid(path, _surface_fluxes(), 'erg/s/angstrom')
    check = check_flux_scale(path)
    assert check.label_kind == 'luminosity density'
    assert check.interpretation == 'surface flux'
    assert check.consistent_with_label is False
    assert 'do not agree' in format_flux_scale_check(check)


def test_surface_flux_divided_by_pi(tmp_path):
    path = tmp_path / 'grid.h5'
    _write_grid(path, [flux / np.pi for flux in _surface_fluxes()],
                'erg/s/cm^2/angstrom')
    check = check_flux_scale(path)
    assert check.interpretation == 'surface flux divided by pi'
    assert check.consistent_with_label is False


def test_luminosity_density_gives_the_radius(tmp_path):
    area = 4 * np.pi * const.R_sun.cgs.value ** 2
    path = tmp_path / 'grid.h5'
    _write_grid(path, [flux * area for flux in _surface_fluxes()],
                'erg/s/angstrom')
    check = check_flux_scale(path)
    assert check.interpretation == 'luminosity density'
    assert check.consistent_with_label is True
    assert all(model.implied_radius == pytest.approx(1.0, rel=1e-9)
               for model in check.models)


@pytest.mark.parametrize('label, consistent', [('', True),
                                               ('erg/s/angstrom', False)])
def test_continuum_normalized_spectra(tmp_path, label, consistent):
    path = tmp_path / 'grid.h5'
    fluxes = np.full((len(TEFFS), len(WAVELENGTH)), 0.97)
    fluxes[:, 1000] = 0.4
    _write_grid(path, fluxes, label)
    check = check_flux_scale(path)
    assert check.interpretation == 'continuum-normalized'
    assert check.consistent_with_label is consistent


def test_model_with_bad_fluxes_is_invalid_and_not_used(tmp_path):
    fluxes = _surface_fluxes()
    fluxes[0] = fluxes[0].copy()
    fluxes[0][10] = np.nan
    fluxes[0][20] = -1e30
    path = tmp_path / 'grid.h5'
    _write_grid(path, fluxes, 'erg/s/cm^2/angstrom')
    check = check_flux_scale(path)
    assert check.models[0].interpretation.startswith('invalid')
    assert check.models[0].n_nonfinite == 1
    assert check.interpretation == 'surface flux'
    assert check.consistent_with_label is True
    assert '1 of the 3 models checked are invalid' in format_flux_scale_check(check)


def test_grid_without_teff_is_refused(tmp_path):
    path = tmp_path / 'grid.h5'
    _write_grid(path, _surface_fluxes(), 'erg/s/cm^2/angstrom')
    with pytest.raises(ValueError, match='teff_column'):
        check_flux_scale(path, teff_column='temperature')
