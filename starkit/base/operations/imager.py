import numpy as np

from astropy import units as u
from starkit.base.operations.base import InstrumentOperationModel
from starkit.filters import as_filters

# numpy 2.5 removed np.trapz, and numpy 1.x has no np.trapezoid.
from scipy.integrate import trapezoid

class ImagerInstrumentOperation(InstrumentOperationModel):
    pass

__all__ = ['Photometry']

class Photometry(ImagerInstrumentOperation):
    """
    Magnitudes of a spectrum in a set of filters.

    The flux is taken to be f_lambda in erg/(s cm2 Angstrom) at the observer.

    Parameters
    ----------
    filter_set : str, starkit.filters.Filter, or a list of them
        filters; strings are SVO filter IDs (see `starkit.filters.as_filters`)
    mag_type : str
        'vega' (the default) or 'ab'
    grid_wavelength : astropy.units.Quantity or array, optional
        wavelengths of every spectrum the operation will receive, such as
        those of a grid; the filter curves are then interpolated onto them
        once instead of at every evaluation. Numbers are taken as Angstrom.
    """
    inputs = ('wavelength', 'flux')
    outputs = ('photometry',)

    @classmethod
    def from_grid(cls, grid, **kwargs):
        return cls(grid_wavelength=grid.wavelength, **kwargs)

    def __init__(self, filter_set, mag_type='vega', grid_wavelength=None):
        super(Photometry, self).__init__()
        if mag_type not in ('vega', 'ab'):
            raise ValueError(
                "mag_type must be 'vega' or 'ab', not {0!r}".format(mag_type))
        self.mag_type = mag_type
        self.filters = as_filters(filter_set)
        self.zero_point = np.array([item.zero_point(mag_type)
                                    for item in self.filters])

        if grid_wavelength is not None:
            self.grid_wavelength = u.Quantity(
                grid_wavelength, u.angstrom).to_value(u.angstrom)
            self.grid_response = self.filter_responses(self.grid_wavelength)
            self.method = 'grid'
        else:
            self.method = 'spectrum'

    def filter_responses(self, wavelength):
        """
        Energy-weighted response of each filter on the given wavelengths,
        as an array of shape (number of filters, number of wavelengths).
        """
        return np.array([item.energy_response(wavelength)
                         for item in self.filters])

    def magnitudes(self, wavelength, flux, response):
        mean_flux = (trapezoid(response * flux, x=wavelength, axis=1)
                     / trapezoid(response, x=wavelength, axis=1))
        return -2.5 * np.log10(mean_flux / self.zero_point)

    def evaluate(self, wavelength, flux):
        if self.method == 'grid':
            return self.magnitudes(self.grid_wavelength, flux,
                                   self.grid_response)
        wavelength = np.asarray(wavelength, float)
        return self.magnitudes(wavelength, flux,
                               self.filter_responses(wavelength))
