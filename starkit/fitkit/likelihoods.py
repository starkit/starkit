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
    
class SpectralChi2LikelihoodAddErr(StarKitModel):
    ## additive error model
    
    inputs = ('wavelength', 'flux')
    outputs = ('loglikelihood', )

    add_err = modeling.Parameter(default=0.0)
    
    def __init__(self, observed, add_err = 0.0):
        super(SpectralChi2LikelihoodAddErr, self).__init__(add_err=add_err)
        self.observed_wavelength = observed.wavelength.to(u.angstrom).value
        self.observed_flux = observed.flux.value
        self.observed_uncertainty = getattr(observed, 'uncertainty', None)
        if self.observed_uncertainty is not None:
            self.observed_uncertainty = self.observed_uncertainty.value
        else:
            self.observed_uncertainty = np.ones_like(self.observed_wavelength)


    def evaluate(self, wavelength, flux, add_err):
        norm = 1.0/np.sqrt(2.0*np.pi*(self.observed_uncertainty**2+add_err**2))
        loglikelihood =  np.sum(np.log(norm) + -0.5 * (
            ((self.observed_flux - flux)**2 / (self.observed_uncertainty**2 + add_err**2))))
        if np.isnan(loglikelihood):
            return -1e300
        return loglikelihood

    
class SpectralL1Likelihood(StarKitModel):
    # this likelihood is for the L1 norm, which is appropriate for
    # Laplacian noise, or to have less sensitvity to outliers. 
    inputs = ('wavelength', 'flux')
    outputs = ('loglikelihood', )

    def __init__(self, observed):
        super(SpectralL1Likelihood, self).__init__()
        self.observed_wavelength = observed.wavelength.to(u.angstrom).value
        self.observed_flux = observed.flux.value
        self.observed_uncertainty = getattr(observed, 'uncertainty', None)
        if self.observed_uncertainty is not None:
            self.observed_uncertainty = self.observed_uncertainty.value
        else:
            self.observed_uncertainty = np.ones_like(self.observed_wavelength)


    def evaluate(self, wavelength, flux):
        loglikelihood =  -1.0*np.sum(
            np.abs(self.observed_flux - flux) / self.observed_uncertainty)
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
