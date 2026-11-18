from astropy import units as u, constants as const
from astropy import modeling
from starkit.base.model import StarKitModel
import numpy as np

class SpectralChi2Likelihood(StarKitModel):
    inputs = ('wavelength', 'flux')
    outputs = ('loglikelihood', )

    def __init__(self, observed, *args, **kwargs):
        super(SpectralChi2Likelihood, self).__init__(*args, **kwargs)
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

    
class SpectralL1Likelihood(SpectralChi2Likelihood):
    # this likelihood is for the L1 norm, which is appropriate for
    # Laplacian noise, or to have less sensitvity to outliers. 
    inputs = ('wavelength', 'flux')
    outputs = ('loglikelihood', )


    def evaluate(self, wavelength, flux):
        loglikelihood =  -1.0*np.sum(
            np.abs(self.observed_flux - flux) / self.observed_uncertainty)
        if np.isnan(loglikelihood):
            return -1e300
        return loglikelihood

class SpectralScaledChi2Likelihood(SpectralChi2Likelihood):
    inputs = ('wavelength', 'flux')
    outputs = ('loglikelihood', )
    lnf = modeling.Parameter()

    def __init__(self, observed, lnf=0.0):
        super(SpectralScaledChi2Likelihood, self).__init__(observed, lnf=lnf)


    def evaluate(self, wavelength, flux, lnf):
        """
        Calculating likelihood and scaling the uncertainties (as shown in http://dfm.io/emcee/current/user/line/)
        This likelihood function is simply a Gaussian where the variance is underestimated by some fractional amount: f.
        One can fit for the natural logarithm of f.

        Parameters
        ----------
        wavelength : numpy.ndarray
            wavelength
        flux : numpy.ndarray
            flux

        Returns
        -------
            : float
            log likelihood
        """
        inv_sigma2 = 1.0/(self.observed_uncertainty**2 + self.observed_flux**2 * np.exp(2 * lnf))
        loglikelihood = -0.5 * (np.sum((flux - self.observed_flux)**2 * inv_sigma2 - np.log(inv_sigma2)))

        if np.isnan(loglikelihood):
            return -1e300
        else:
            return loglikelihood

def _observed_magnitudes(magnitudes, magnitude_uncertainties):
    """
    Observed magnitudes and their uncertainties as arrays, from two arrays or
    from an object with the attributes ``magnitudes`` and
    ``magnitude_uncertainties``, such as a wsynphot ``MagnitudeSet``.
    """
    if magnitude_uncertainties is None:
        if not hasattr(magnitudes, 'magnitude_uncertainties'):
            raise TypeError('Give the magnitude_uncertainties of the '
                            'observed magnitudes.')
        magnitude_uncertainties = magnitudes.magnitude_uncertainties
        magnitudes = magnitudes.magnitudes
    magnitudes = np.asarray(magnitudes, float)
    magnitude_uncertainties = np.asarray(magnitude_uncertainties, float)
    if magnitudes.shape != magnitude_uncertainties.shape:
        raise ValueError('magnitudes and magnitude_uncertainties must have '
                         'the same shape')
    return magnitudes, magnitude_uncertainties


class PhotometryColorLikelihood(StarKitModel):
    """
    Gaussian likelihood of the colors between consecutive filters.

    Parameters
    ----------
    magnitudes : array
        observed magnitudes, in the order of the filters of the Photometry
        operation; an object with the attributes ``magnitudes`` and
        ``magnitude_uncertainties`` is also accepted
    magnitude_uncertainties : array
    """
    inputs = ('photometry',)
    outputs = ('loglikelihood',)

    def __init__(self, magnitudes, magnitude_uncertainties=None):
        super(PhotometryColorLikelihood, self).__init__()
        magnitudes, magnitude_uncertainties = _observed_magnitudes(
            magnitudes, magnitude_uncertainties)
        self.colors = magnitudes[:-1] - magnitudes[1:]
        self.color_uncertainties = np.sqrt(
            magnitude_uncertainties[:-1]**2 + magnitude_uncertainties[1:]**2)

    def evaluate(self, photometry):
        synth_colors = photometry[:-1] - photometry[1:]
        loglikelihood = -0.5 * np.sum(((self.colors - synth_colors)
                                       / self.color_uncertainties)**2)
        return loglikelihood


class RelativePhotometryLikelihood(StarKitModel):
    """
    Gaussian likelihood of the magnitudes.

    Parameters
    ----------
    magnitudes : array
        observed magnitudes, in the order of the filters of the Photometry
        operation; an object with the attributes ``magnitudes`` and
        ``magnitude_uncertainties`` is also accepted
    magnitude_uncertainties : array
    """

    inputs = ('photometry', )
    outputs = ('loglikelihood', )

    def __init__(self, magnitudes, magnitude_uncertainties=None):
        super(RelativePhotometryLikelihood, self).__init__()
        self.magnitudes, self.magnitude_uncertainties = _observed_magnitudes(
            magnitudes, magnitude_uncertainties)

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
