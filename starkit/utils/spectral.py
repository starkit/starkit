import numpy as np
from starkit.fix_spectrum1d import SKSpectrum1D

FWHM2SIGMA_CONST = 2 * np.sqrt(2*np.log(2))

def fwhm2sigma(fwhm):
    """
    Convert full-width half maximum to sigma

    Parameters
    ----------
    fwhm: float

    Returns
    -------
    sigma: float
    """
    return fwhm / FWHM2SIGMA_CONST

def sigma2fwhm(sigma):
    """
    Convert full-width half maximum to sigma

    Parameters
    ----------
    sigma: float

    Returns
    -------
    fwhm: float
    """

    return sigma * FWHM2SIGMA_CONST

def prepare_observed(observed):
    """
    Preparing an observed spectrum

    Parameters
    ----------
    observed: Spectrum1D object

    """

    wavelength, flux = observed.wavelength, observed.flux

    if getattr(observed, 'uncertainty', None) is None:
        uncertainty = np.ones_like(flux)
    else:
        uncertainty = getattr(observed.uncertainty, 'array',
                                       observed.uncertainty).value

    spec = SKSpectrum1D(wavelength, flux, uncertainty)

    return spec