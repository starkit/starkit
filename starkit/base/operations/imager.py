import numpy as np

from astropy import units as u
from starkit.fix_spectrum1d import SKSpectrum1D
from starkit.base.operations.base import InstrumentOperationModel

from scipy import interpolate
# numpy 2.5 removed np.trapz, and numpy 1.x has no np.trapezoid.
from scipy.integrate import trapezoid

FLAM = u.erg / (u.s * u.cm ** 2 * u.angstrom)
# wsynphot.io.cache_filters.DetectorType.PHOTON_COUNTER
PHOTON_COUNTER = 1

class ImagerInstrumentOperation(InstrumentOperationModel):
    pass

__all__ = ['Photometry']

class Photometry(ImagerInstrumentOperation):
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
        try:
            from wsynphot import FilterSet
        except ModuleNotFoundError as exc:
            if exc.name != 'wsynphot':
                raise
            raise ModuleNotFoundError(
                'Photometry needs wsynphot, which is not on PyPI. Install it '
                'from https://github.com/starkit/wsynphot; see '
                'https://starkit.github.io/starkit/installation.html.',
                name='wsynphot',
            ) from exc

        if hasattr(filter_set, 'calculate_{0}_magnitudes'.format(mag_type)):
            self.filter_set = filter_set
        else:
            self.filter_set = FilterSet(filter_set)


        if grid_wavelength is not None:
            # Grids give their wavelengths as a Quantity. The filter
            # integral needs plain values, in Angstrom like the zero points.
            self.grid_wavelength = u.Quantity(grid_wavelength,
                                              u.angstrom).value
            (self.zero_point, self.filter_transmission,
             self.wavelength_delta) = (
                self.interpolate_filters_to_grid_wavelength(self.filter_set))
            self.method = 'grid'

        else:
            self.method = 'slow'
            self.calculate_magnitudes = getattr(
                self.filter_set, 'calculate_{0}_magnitudes'.format(mag_type))



    def interpolate_filters_to_grid_wavelength(self, filter_set):
        """
        Interpolate the filter curves onto the grid wavelengths.

        The magnitudes are computed as in wsynphot: the mean flux density in
        a filter is the integral of the transmission times the flux, divided
        by the integral of the transmission. For a photon-counting filter,
        both integrands are also multiplied by the wavelength.

        Returns
        -------
        zero_point : numpy.ndarray
            Vega or AB zero point of each filter, in erg/(s cm2 Angstrom)
        filter_transmission : numpy.ndarray
            transmission of each filter on the grid wavelengths, times the
            wavelength in Angstrom for a photon-counting filter
        wavelength_delta : numpy.ndarray
            integral of the transmission, in Angstrom, or of the transmission
            times the wavelength, in Angstrom^2, for a photon-counting filter
        """
        filter_transmission = np.empty((len(filter_set.filter_set),
                                        self.grid_wavelength.shape[0]))
        zero_point = np.empty(len(filter_set.filter_set))
        wavelength_delta = np.empty(len(filter_set.filter_set))
        for i, item in enumerate(filter_set):
            transmission = interpolate.interp1d(
                item.wavelength.to_value(u.angstrom),
                item.transmission_lambda, bounds_error=False,
                fill_value=0.0)(self.grid_wavelength)

            if getattr(item, 'detector_type', None) == PHOTON_COUNTER:
                filter_transmission[i] = transmission * self.grid_wavelength
                wavelength_delta[i] = item.calculate_wavelength_delta(
                    ).to_value(u.angstrom ** 2)
            else:
                filter_transmission[i] = transmission
                wavelength_delta[i] = item.calculate_wavelength_delta(
                    ).to_value(u.angstrom)

            zero_point[i] = getattr(
                item, 'zp_{0}_f_lambda'.format(self.mag_type)).to_value(FLAM)
        return zero_point, filter_transmission, wavelength_delta


    def evaluate_grid(self, wavelength, flux):
        return -2.5 * np.log10(trapezoid(self.filter_transmission * flux,
                                 x=self.grid_wavelength, axis=1)
                        / self.wavelength_delta / self.zero_point)

    def evaluate_slow(self, wavelength, flux):

        spec = SKSpectrum1D.from_array(wavelength * u.angstrom,
                                       flux * u.erg / u.s / u.cm ** 2 / u.angstrom)
        return np.array(u.Quantity(self.calculate_magnitudes(spec)).value)


    def evaluate(self, wavelength, flux):
        if self.method == 'slow':
            return self.evaluate_slow(wavelength, flux)
        elif self.method == 'grid':
            return self.evaluate_grid(wavelength, flux)
        else:
            raise ValueError('method attribute is only allowed to be grid or slow')



