import numpy as np
from starkit.fix_spectrum1d import SKSpectrum1D


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