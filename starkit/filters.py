"""
Filter curves and zero points for synthetic photometry.

A `Filter` holds a response curve and its detector type. The mean flux
density of a spectrum through a filter is

    <f_lambda> = int f_lambda S(lambda) lambda dlambda / int S(lambda) lambda dlambda

for a photon-counting curve S, and

    <f_lambda> = int f_lambda S'(lambda) dlambda / int S'(lambda) dlambda

for an energy response curve S', as published for photomultiplier
photometry (Bessell & Murphy 2012, PASP 124, 140, Appendix A). Both forms are
computed here with the energy-weighted response, S lambda or S'.
"""

from functools import lru_cache
from pathlib import Path

import numpy as np
from astropy import constants as const
from astropy import units as u
from astropy.io import fits, votable
from astropy.utils.data import clear_download_cache, download_file
from scipy.integrate import trapezoid

__all__ = ['Filter', 'as_filters', 'vega_spectrum']

FLAM = u.erg / (u.s * u.cm ** 2 * u.angstrom)
AB_ZERO_POINT = 3631 * u.Jy
SVO_URL = 'https://svo2.cab.inta-csic.es/theory/fps/fps.php?ID={0}'
DETECTORS = ('photon', 'energy')
# DetectorType in the SVO Filter Profile Service and in wsynphot.
DETECTOR_TYPES = {0: 'energy', 1: 'photon'}
# Largest fraction of a filter's response allowed outside a spectrum.
COVERAGE_LIMIT = 0.01


# CALSPEC spectrum of Vega (Bohlin, Gordon & Tremblay 2014, PASP 126, 711),
# from https://ssb.stsci.edu/trds/calspec/alpha_lyr_stis_011.fits. Another
# CALSPEC FITS table, as a path or a URL, can be set here; call
# vega_spectrum.cache_clear() after changing it.
VEGA_FILE = str(Path(__file__).parent / 'data' / 'alpha_lyr_stis_011.fits')


@lru_cache(maxsize=None)
def vega_spectrum():
    """
    Spectrum of Vega used for the Vega zero points, read from `VEGA_FILE`.

    Returns
    -------
    wavelength : numpy.ndarray
        wavelengths in Angstrom
    flux : numpy.ndarray
        f_lambda in erg/(s cm2 Angstrom)
    """
    with fits.open(VEGA_FILE, cache=True) as hdus:
        table = hdus[1].data
        return (np.array(table['WAVELENGTH'], float),
                np.array(table['FLUX'], float))


class Filter(object):
    """
    Response curve of a photometric filter.

    Parameters
    ----------
    wavelength : astropy.units.Quantity or array
        wavelengths of the curve, increasing; numbers without a unit are
        taken as Angstrom
    transmission : array
        response at these wavelengths, on any scale
    detector : str
        'photon' (the default) for a curve giving the relative number of
        photons detected, as for a filter used with a CCD or infrared array,
        including measured filter transmissions; 'energy' for an energy
        response function, as published for photomultiplier photometry
    name : str, optional
    """

    def __init__(self, wavelength, transmission, detector='photon', name=None):
        if detector not in DETECTORS:
            raise ValueError("detector must be 'photon' or 'energy', not "
                             "{0!r}".format(detector))
        wavelength = np.asarray(
            u.Quantity(wavelength, u.angstrom).to_value(u.angstrom), float)
        transmission = np.asarray(transmission, float)
        if (wavelength.ndim != 1 or wavelength.shape != transmission.shape
                or wavelength.size < 2):
            raise ValueError('wavelength and transmission must be 1-d arrays '
                             'of the same length, with at least two points')
        if not np.all(np.diff(wavelength) > 0):
            raise ValueError('the wavelengths of a filter must increase')
        if not (np.all(np.isfinite(transmission)) and transmission.max() > 0):
            raise ValueError('the transmission must be finite and positive '
                             'somewhere')
        self.wavelength = wavelength * u.angstrom
        self.transmission = transmission
        self.detector = detector
        self.name = name

    def __repr__(self):
        return '<Filter {0} ({1}), {2:g} to {3:g} Angstrom>'.format(
            self.name, self.detector, self.wavelength[0].value,
            self.wavelength[-1].value)

    @classmethod
    def read(cls, path, wavelength_unit=u.angstrom, detector='photon',
             name=None):
        """
        Read a filter curve from a text file.

        The file has two columns, wavelength and transmission, separated by
        white space. Lines starting with ``#`` are skipped.

        Parameters
        ----------
        path : str or pathlib.Path
        wavelength_unit : astropy.units.Unit
            unit of the wavelength column [default Angstrom]
        detector : str
            'photon' (the default) or 'energy'; see `Filter`
        name : str, optional
            the file name without its extension by default
        """
        data = np.loadtxt(path, comments='#', usecols=(0, 1), ndmin=2)
        if name is None:
            name = Path(path).stem
        return cls(data[:, 0] * u.Unit(wavelength_unit), data[:, 1],
                   detector=detector, name=name)

    @classmethod
    def from_svo(cls, filter_id, detector=None, cache=True):
        """
        Download a filter curve from the SVO Filter Profile Service.

        Parameters
        ----------
        filter_id : str
            SVO filter ID, such as '2MASS/2MASS.Ks'; wsynphot-style IDs such
            as '2MASS/2MASS/Ks' are also accepted
        detector : str, optional
            'photon' or 'energy' in place of the detector type given by SVO
        cache : bool
            keep the curve in the astropy download cache [default True]
        """
        if filter_id.count('/') == 2:
            facility, instrument, band = filter_id.split('/')
            filter_id = '{0}/{1}.{2}'.format(facility, instrument, band)
        url = SVO_URL.format(filter_id)
        response = votable.parse(download_file(url, cache=cache),
                                 verify='ignore')
        status = [info.value for info in response.infos
                  if info.name == 'QUERY_STATUS']
        tables = list(response.iter_tables())
        if status != ['OK'] or not tables:
            if cache:
                clear_download_cache(url)
            raise ValueError('The SVO Filter Profile Service has no filter '
                             '{0!r}.'.format(filter_id))
        table = tables[0].to_table()
        wavelength_unit = table['Wavelength'].unit or u.angstrom
        try:
            detector_type = int(
                response.get_field_by_id_or_name('DetectorType').value)
        except KeyError:
            detector_type = 1
        if detector is None:
            detector = DETECTOR_TYPES.get(detector_type, 'photon')
        return cls(np.asarray(table['Wavelength'], float) * wavelength_unit,
                   np.asarray(table['Transmission'], float),
                   detector=detector, name=filter_id)

    def _energy_weighted(self, wavelength, transmission):
        if self.detector == 'photon':
            return transmission * wavelength
        return transmission

    def energy_response(self, wavelength):
        """
        Energy-weighted response on the given wavelengths.

        This is the transmission times the wavelength for a photon-counting
        filter, and the transmission for an energy response, interpolated
        linearly and zero outside the curve.

        Raises
        ------
        ValueError
            if more than 1% of the filter's response lies outside the
            given wavelengths
        """
        wavelength = np.asarray(
            u.Quantity(wavelength, u.angstrom).to_value(u.angstrom), float)
        own_wavelength = self.wavelength.value
        own_response = self._energy_weighted(own_wavelength, self.transmission)
        clipped = np.clip(own_wavelength, wavelength.min(), wavelength.max())
        inside = trapezoid(np.interp(clipped, own_wavelength, own_response),
                           clipped)
        outside = 1.0 - inside / trapezoid(own_response, own_wavelength)
        if outside > COVERAGE_LIMIT:
            raise ValueError(
                '{0:.1%} of the response of filter {1} lies outside the '
                'wavelengths {2:g} to {3:g} Angstrom.'.format(
                    outside, self.name, wavelength.min(), wavelength.max()))
        transmission = np.interp(wavelength, own_wavelength,
                                 self.transmission, left=0.0, right=0.0)
        return self._energy_weighted(wavelength, transmission)

    def mean_flux_density(self, wavelength, flux):
        """
        Mean f_lambda of a spectrum through the filter.

        Parameters
        ----------
        wavelength : astropy.units.Quantity or array
            wavelengths of the spectrum; numbers are taken as Angstrom
        flux : array
            f_lambda of the spectrum

        Returns
        -------
        float
            mean f_lambda, in the unit of ``flux``
        """
        wavelength = np.asarray(
            u.Quantity(wavelength, u.angstrom).to_value(u.angstrom), float)
        response = self.energy_response(wavelength)
        return (trapezoid(response * np.asarray(flux, float), wavelength)
                / trapezoid(response, wavelength))

    @property
    def pivot_wavelength(self):
        """
        Pivot wavelength, which converts the mean f_nu to the mean f_lambda.

        For a photon-counting curve S, lambda_p^2 = int S lambda dlambda /
        int S / lambda dlambda (Bessell & Murphy 2012, equation A16); for an
        energy response S', int S' dlambda / int S' / lambda^2 dlambda.
        """
        # The curve is interpolated linearly, as on a spectrum, onto a fine
        # sampling: SVO curves can have steps of 200 Angstrom.
        own_wavelength = self.wavelength.value
        wavelength = np.union1d(own_wavelength, np.linspace(
            own_wavelength[0], own_wavelength[-1], 10001))
        response = self._energy_weighted(
            wavelength,
            np.interp(wavelength, own_wavelength, self.transmission))
        return np.sqrt(trapezoid(response, wavelength)
                       / trapezoid(response / wavelength ** 2,
                                   wavelength)) * u.angstrom

    def zero_point(self, mag_type):
        """
        Mean f_lambda of a magnitude of 0, in erg/(s cm2 Angstrom).

        Parameters
        ----------
        mag_type : str
            'vega' for the mean flux density of Vega (see `vega_spectrum`),
            'ab' for a source with f_nu = 3631 Jy
        """
        if mag_type == 'ab':
            return (AB_ZERO_POINT * const.c
                    / self.pivot_wavelength ** 2).to_value(FLAM)
        if mag_type == 'vega':
            vega_wavelength, vega_flux = vega_spectrum()
            low, high = self.wavelength[[0, -1]].value
            in_filter = (vega_wavelength > low) & (vega_wavelength < high)
            wavelength = np.union1d(self.wavelength.value,
                                    vega_wavelength[in_filter])
            return self.mean_flux_density(
                wavelength, np.interp(wavelength, vega_wavelength, vega_flux))
        raise ValueError("mag_type must be 'vega' or 'ab', not "
                         "{0!r}".format(mag_type))


def as_filters(filter_set):
    """
    Make a list of `Filter` objects.

    Parameters
    ----------
    filter_set : str, Filter, or a list of them
        filters, with strings taken as SVO filter IDs; wsynphot ``FilterSet``
        and ``FilterCurve`` objects are also accepted

    Returns
    -------
    list of Filter
    """
    if (isinstance(filter_set, (str, Filter))
            or hasattr(filter_set, 'transmission_lambda')):
        filter_set = [filter_set]
    elif hasattr(filter_set, 'filter_set'):
        filter_set = filter_set.filter_set
    return [_as_filter(item) for item in filter_set]


def _as_filter(item):
    if isinstance(item, Filter):
        return item
    if isinstance(item, str):
        return Filter.from_svo(item)
    if hasattr(item, 'transmission_lambda'):
        # wsynphot.FilterCurve. Versions before February 2022 had no
        # detector_type; from October 2017 they counted photons for every
        # curve.
        detector_type = int(getattr(item, 'detector_type', 1))
        return Filter(item.wavelength, item.transmission_lambda,
                      detector=DETECTOR_TYPES.get(detector_type, 'photon'),
                      name=getattr(item, 'filter_id', None))
    raise TypeError('Cannot make a filter from {0!r}; give a Filter or an '
                    'SVO filter ID.'.format(item))
