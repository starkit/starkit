#from specutils.spectrum1d import Spectrum1D as SUSpectrum1D
from astropy import units as u

class Spectrum1D(object):


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
