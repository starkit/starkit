from starkit.base.model import StarKitModel
from astropy import modeling
from collections import OrderedDict
from itertools import chain


class SpectralOperationModel(StarKitModel):

    inputs = ('wavelength', 'flux')
    outputs = ('wavelength', 'flux')

class InstrumentOperationModel(SpectralOperationModel):
    pass


class DoubleSpectrum(SpectralOperationModel):
    inputs = ('wavelength', 'flux')
    outputs = ('wavelength', 'flux', 'wavelength', 'flux')


    def evaluate(self, wavelength, flux):
        return wavelength, flux, wavelength, flux

