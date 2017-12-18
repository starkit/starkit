#from specutils.spectrum1d import Spectrum1D as SUSpectrum1D
from astropy import units as u

class SKSpectrum1D(object):

    @classmethod
    def from_array(cls, wavelength, flux, uncertainty):
        """

        Parameters
        ----------
        wavelength: astropy quantity
        flux: astropy quantity
        uncertainty: astropy quantity

        Returns
        -------
            : SKSpectrum1D object

        """
        return cls(wavelength, flux, uncertainty)

    def __init__(self, wavelength, flux, uncertainty):
        self.wavelength = wavelength
        self.flux = flux
        self.uncertainty = uncertainty


    def uncertainty_getter(self):
        return self._uncertainty

    def uncertainty_setter(self, value):
        if value is None:
            self._uncertainty = None
        else:
            self._uncertainty = u.Quantity(value, self.flux.unit)

    uncertainty = property(uncertainty_getter, uncertainty_setter)

    def slice_index(self, start=None, stop=None, step=None):
        spectral_slice = slice(start, stop, step)
        if self.uncertainty is not None:
            new_uncertainty = self.uncertainty[spectral_slice]
        else:
            new_uncertainty = None

        return SKSpectrum1D(self.wavelength[spectral_slice], self.flux[spectral_slice],
                            new_uncertainty)