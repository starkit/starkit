import numpy as np

from astropy import units as u
from starkit.fix_spectrum1d import SKSpectrum1D
from starkit.base.operations.base import InstrumentOperationModel

from scipy import interpolate

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
        try:
            from wsynphot import FilterSet
        except ImportError:
            raise ImportError('The photometry plugin needs wsynphot')

        if hasattr(filter_set, 'calculate_{0}_magnitudes'.format(mag_type)):
            filter_set = filter_set
        else:
            filter_set = FilterSet(filter_set)


        if grid_wavelength is not None:
            self.grid_wavelength = grid_wavelength
            self.zp_vega, self.filter_transmission, self.wavelength_delta = (
                self.interpolate_filters_to_grid_wavelength(filter_set))
            self.method = 'grid'

        else:
            self.method = 'slow'
            self.calculate_magnitudes = getattr(
                self.filter_set, 'calculate_{0}_magnitudes'.format(mag_type))
            self.filter_set = filter_set



    def interpolate_filters_to_grid_wavelength(self, filter_set):
        filter_transmission = np.empty((len(filter_set.filter_set),
                                        self.grid_wavelength.shape[0]))
        zp_vega = np.empty(len(filter_set.filter_set))
        wavelength_delta = np.empty(len(filter_set.filter_set))
        for i, item in enumerate(filter_set):
            filter_transmission[i] = interpolate.interp1d(
                item.wavelength.to('angstrom').value,
                item.transmission_lambda, bounds_error=False,
                fill_value=0.0)(self.grid_wavelength)

            #TODO
            ### be careful the units here are erg / (angstrom cm2 s) ##
            zp_vega[i] = item.zp_vega_f_lambda.value
            wavelength_delta[i] = item.calculate_wavelength_delta().to(
                'angstrom').value
        return zp_vega, filter_transmission, wavelength_delta


    def evaluate_grid(self, wavelength, flux):
        return -2.5 * np.log10(np.trapz(self.filter_transmission * flux,
                                 x=self.grid_wavelength, axis=1)
                        / self.wavelength_delta / self.zp_vega)

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



