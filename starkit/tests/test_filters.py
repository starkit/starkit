"""Tests for the filter curves used in synthetic photometry."""

import sys

import numpy as np
import pytest
from astropy import constants as const
from astropy import units as u

from starkit import filters
from starkit.filters import Filter

FLAM = u.erg / (u.s * u.cm ** 2 * u.angstrom)
# A box from 4000 to 5000 Angstrom, finely sampled.
BOX_WAVELENGTH = np.linspace(4000.0, 5000.0, 1001)
BOX = np.ones(BOX_WAVELENGTH.size)
SVO_ERROR = """<?xml version="1.0"?>
<VOTABLE version="1.1" xmlns:xsi="http://www.w3.org/2001/XMLSchema-instance">
  <INFO name="QUERY_STATUS" value="ERROR">
    <DESCRIPTION>Filter not found: </DESCRIPTION>
  </INFO>
</VOTABLE>
"""


@pytest.mark.parametrize('detector, expected', [
    # f_lambda = lambda / 1000 Angstrom. Photon counting:
    # int f lambda dlambda / int lambda dlambda
    # = ((5000^3 - 4000^3) / 3000) / ((5000^2 - 4000^2) / 2).
    ('photon', (5000.0 ** 3 - 4000.0 ** 3) / 3000.0
               / ((5000.0 ** 2 - 4000.0 ** 2) / 2.0)),
    # Energy integration: the mean of f_lambda, 4.5.
    ('energy', 4.5),
])
def test_mean_flux_density_weights_by_detector_type(detector, expected):
    box = Filter(BOX_WAVELENGTH * u.angstrom, BOX, detector=detector)
    wavelength = np.linspace(3500.0, 5500.0, 20001)
    assert box.mean_flux_density(wavelength, wavelength / 1000.0) == (
        pytest.approx(expected, rel=1e-6))


def test_filters_count_photons_by_default():
    assert Filter(BOX_WAVELENGTH, BOX).detector == 'photon'


def test_filter_rejects_unknown_detector_types():
    with pytest.raises(ValueError, match="'photon' or 'energy'"):
        Filter(BOX_WAVELENGTH, BOX, detector='bolometer')


def test_filter_reads_a_text_file_in_microns(tmp_path):
    path = tmp_path / 'nb217.txt'
    path.write_text('# wavelength (um)  transmission (%)\n'
                    '2.15 0.0\n2.16 50.0\n2.17 80.0\n2.18 50.0\n2.19 0.0\n')
    nb217 = Filter.read(path, wavelength_unit=u.um)
    np.testing.assert_allclose(nb217.wavelength.to_value(u.angstrom),
                               [21500.0, 21600.0, 21700.0, 21800.0, 21900.0])
    np.testing.assert_allclose(nb217.transmission, [0.0, 50.0, 80.0, 50.0, 0.0])
    assert nb217.detector == 'photon'
    assert nb217.name == 'nb217'


def test_filter_from_svo_reads_the_detector_type(svo):
    johnson_v = Filter.from_svo('Generic/Johnson.V')
    ks = Filter.from_svo('2MASS/2MASS.Ks')
    assert johnson_v.detector == 'energy'
    assert ks.detector == 'photon'
    assert johnson_v.name == 'Generic/Johnson.V'
    assert johnson_v.wavelength.size == 15
    assert johnson_v.wavelength[0] == 4600.0 * u.angstrom
    assert svo['urls'] == [
        'https://svo2.cab.inta-csic.es/theory/fps/fps.php?ID=Generic/Johnson.V',
        'https://svo2.cab.inta-csic.es/theory/fps/fps.php?ID=2MASS/2MASS.Ks']


def test_filter_from_svo_can_override_the_detector_type(svo):
    # SVO labels some curves of the filter alone as energy responses.
    assert Filter.from_svo('Generic/Johnson.V', detector='photon').detector == (
        'photon')


def test_filter_from_svo_accepts_wsynphot_filter_ids(svo):
    # wsynphot named the filters Facility/Instrument/Band.
    assert Filter.from_svo('Generic/Johnson/V').name == 'Generic/Johnson.V'


def test_filter_from_svo_reports_unknown_filters(svo, tmp_path):
    svo['error_file'] = str(tmp_path / 'error.xml')
    (tmp_path / 'error.xml').write_text(SVO_ERROR)
    with pytest.raises(ValueError, match='Nope/Nope.X'):
        Filter.from_svo('Nope/Nope.X')
    # The error response is not kept in the download cache.
    assert svo['cleared'] == [
        'https://svo2.cab.inta-csic.es/theory/fps/fps.php?ID=Nope/Nope.X']


@pytest.mark.parametrize('detector', ['photon', 'energy'])
def test_flat_f_nu_of_3631_jy_has_ab_magnitude_zero(detector):
    box = Filter(BOX_WAVELENGTH, BOX, detector=detector)
    wavelength = np.linspace(3500.0, 5500.0, 20001)
    f_lambda = (3631 * u.Jy * const.c / (wavelength * u.angstrom) ** 2).to_value(FLAM)
    assert box.mean_flux_density(wavelength, f_lambda) == (
        pytest.approx(box.zero_point('ab'), rel=1e-6))


def test_vega_zero_point_is_the_mean_flux_density_of_vega(monkeypatch):
    # A power law in place of the Vega spectrum: f_lambda = 1e-9 (lambda / 5000)^-2.
    wavelength = np.linspace(3000.0, 10000.0, 70001)
    monkeypatch.setattr(filters, 'vega_spectrum',
                        lambda: (wavelength, 1e-9 * (wavelength / 5000.0) ** -2))
    box = Filter(BOX_WAVELENGTH, BOX)
    # Photon counting: 1e-9 5000^2 ln(5000 / 4000) / ((5000^2 - 4000^2) / 2).
    expected = 1e-9 * 5000.0 ** 2 * np.log(1.25) / 4.5e6
    assert box.zero_point('vega') == pytest.approx(expected, rel=1e-6)


def test_vega_zero_points_report_how_to_install_synphot(monkeypatch):
    # A None entry in sys.modules makes "import synphot" fail as if the
    # package were not installed.
    monkeypatch.setitem(sys.modules, 'synphot', None)
    filters.vega_spectrum.cache_clear()
    with pytest.raises(ModuleNotFoundError, match=r"starkit\[photometry\]"):
        Filter(BOX_WAVELENGTH, BOX).zero_point('vega')
    filters.vega_spectrum.cache_clear()


def test_filter_must_lie_within_the_spectrum():
    box = Filter(BOX_WAVELENGTH, BOX)
    wavelength = np.linspace(4200.0, 6000.0, 1801)
    with pytest.raises(ValueError, match='outside'):
        box.mean_flux_density(wavelength, np.ones(wavelength.size))
