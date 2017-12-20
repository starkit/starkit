import os

import numpy as np
import pandas as pd
from scipy import ndimage as nd
from scipy.interpolate import interp1d

from astropy import units as u
from astropy.io import fits

from starkit.gridkit.io.process import BaseProcessGrid
from starkit.gridkit.io.bosz.base import convert_bz2_memmap

class BOSZProcessGrid(BaseProcessGrid):
    """

    """

    R_initial = 300000
    R_initial_sampling=2
    def __init__(self, index, input_wavelength, meta, wavelength_start=0*u.angstrom, wavelength_stop=np.inf*u.angstrom,
                 R=5000.0, R_sampling=4):
        """

        Parameters
        ----------
        index
        input_wavelength
        meta
        wavelength_start
        wavelength_stop
        R
        R_sampling
        pre_sampling: Select the sampling that you want to have before convolution (there are strange jumps in the diff)
        """
        super(BOSZProcessGrid, self).__init__(index, input_wavelength, meta, wavelength_start=wavelength_start,
                                                 wavelength_stop=wavelength_stop, R=R, R_sampling=R_sampling)


    def interp_wavelength(self, flux):
        cut_flux = flux[self.start_idx:self.stop_idx]
        rescaled_R = 1 / np.sqrt((1/self.R)**2 - (1/self.R_initial)**2)
        sigma = ((self.R_initial / rescaled_R) * self.R_initial_sampling
                 / (2 * np.sqrt(2 * np.log(2))))

        processed_flux = nd.gaussian_filter1d(cut_flux, sigma)
        output_flux = interp1d(self.cut_wavelength, processed_flux)(
            self.output_wavelength)
        return output_flux


    def load_flux(self, fname):
        """
        load the 1D flux array from file

        Parameters
        ----------
        fname : str

        Returns
        -------
            : numpy.ndarray
        """
        fname_npy = fname.replace('.bz2', '.v1.npy')
        if not os.path.exists(fname_npy):
            convert_bz2_memmap(fname)
        flux = np.load(fname_npy)
        return flux
