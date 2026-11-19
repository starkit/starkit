"""Checks of the solar spectrum shipped in starkit/data."""

from pathlib import Path

import numpy as np
from astropy import units as u
from astropy.io import fits
from astropy.table import Table

import starkit

SUN_FILE = Path(starkit.__file__).parent / 'data' / 'sun_kurucz2005_r1500.fits'


def test_solar_spectrum_columns_and_units():
    table = Table.read(SUN_FILE, hdu=1)
    assert table['WAVELENGTH'].unit == u.angstrom
    assert table['FLUX'].unit == u.W / u.m**2 / u.nm
    assert len(table) == 13144


def test_solar_spectrum_range_and_sampling():
    table = Table.read(SUN_FILE, hdu=1)
    wavelength = np.asarray(table['WAVELENGTH'])
    flux = np.asarray(table['FLUX'])
    assert np.all(np.diff(wavelength) > 0)
    assert 3009 < wavelength[0] < 3012
    assert 8999 < wavelength[-1] <= 9000
    # The wavelengths of the test grid: R = 3000 with four points per
    # resolution element.
    step = np.median(np.diff(np.log(wavelength)))
    np.testing.assert_allclose(step, 1 / 12000, rtol=1e-3)
    assert np.all(np.isfinite(flux)) and np.all(flux > 0)


def test_solar_spectrum_header():
    header = fits.getheader(SUN_FILE, 1)
    assert header['WAVETYPE'] == 'vacuum'
    assert header['RESOLUTN'] == 1500
    assert 'Kurucz' in header['REFERENC']
