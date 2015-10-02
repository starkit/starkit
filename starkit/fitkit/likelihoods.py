from astropy import units as u, constants as const
from astropy import modeling
from starkit.base.model import StarKitModel
import numpy as np

class SpectralChi2Likelihood(StarKitModel):
    inputs = ('wavelength', 'flux')
    outputs = ('loglikelihood', )

    def __init__(self, observed):
        super(SpectralChi2Likelihood, self).__init__()
        self.observed_wavelength = observed.wavelength.to(u.angstrom).value
        self.observed_flux = observed.flux.value
        self.observed_uncertainty = getattr(observed, 'uncertainty', None)
        if self.observed_uncertainty is not None:
            self.observed_uncertainty = self.observed_uncertainty.value
        else:
            self.observed_uncertainty = np.ones_like(self.observed_wavelength)


    def evaluate(self, wavelength, flux):
        loglikelihood =  -0.5 * np.sum(
            ((self.observed_flux - flux) / self.observed_uncertainty)**2)
        if np.isnan(loglikelihood):
            return -1e300
        return loglikelihood

class PhotometryColorLikelihood(StarKitModel):
    inputs = ('photometry',)
    outputs = ('loglikelihood',)

    def __init__(self, magnitude_set):
        super(PhotometryColorLikelihood, self).__init__()
        self.colors = (magnitude_set.magnitudes[:-1] -
                       magnitude_set.magnitudes[1:])
        self.color_uncertainties = np.sqrt(
            magnitude_set.magnitude_uncertainties[:-1]**2
            + magnitude_set.magnitude_uncertainties[1:]**2)

    def evaluate(self, photometry):
        synth_colors = photometry[:-1] - photometry[1:]
        loglikelihood = -0.5 * np.sum(((self.colors - synth_colors)
                                       / self.color_uncertainties)**2)
        return loglikelihood


class RelativePhotometryLikelihood(StarKitModel):

    inputs = ('photometry', )
    outputs = ('loglikelihood', )

    def __init__(self, magnitude_set):
        super(RelativePhotometryLikelihood, self).__init__()
        self.magnitudes = magnitude_set.magnitudes
        self.magnitude_uncertainties = magnitude_set.magnitude_uncertainties

    def evaluate(self, photometry):
        loglikelihood = -0.5 * np.sum(
            ((self.magnitudes - photometry) /
             self.magnitude_uncertainties)**2)

        return loglikelihood

class SpectroPhotometryColorLikelihood(StarKitModel):
    inputs = ('wavelength', 'flux', 'photometry')
    pass



class Addition(StarKitModel):

    inputs = ('a', 'b')
    outputs = ('x', )

    def __init__(self):
        super(Addition, self).__init__()

    @staticmethod
    def evaluate(a, b):
        return a + b
