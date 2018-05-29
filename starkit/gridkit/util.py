## File for Utility functions dealing with grids
import numpy as np
from scipy import ndimage as nd

def convolve_grid_to_R(flux, source_R, source_R_sampling, target_R):
    """
    Requires the grid to be in logarithmic wavelength scaling. Here we define
    the resolution R as lambda/delta_lambda

    Parameters
    ----------
    flux: numpy.ndarray
        flux in numpy ndarray
    source_R: float
        R of the input flux
    source_R_sampling: float
        pixels per resolution element for the input flux
    target_R: float
        R of the output flux
    Returns
    -------
        : numpy.ndarray
        flux convolved to the new resolution target_R with target_R_sampling pixels
        per resolution element
    """

    assert target_R < source_R, ("Requested resolution {target_R} must be "
                                 "smaller than given resolution source "
                                 "resolution {source_R}".format(target_R=target_R,
                                                                source_R=source_R))

    rescale_R = 1 / np.sqrt((1 / target_R) ** 2 - (1 / source_R) ** 2)
    sigma = ((source_R / rescale_R) * source_R_sampling
             / (2 * np.sqrt(2 * np.log(2))))

    convolved_flux = nd.gaussian_filter1d(flux, sigma)

    return convolved_flux
