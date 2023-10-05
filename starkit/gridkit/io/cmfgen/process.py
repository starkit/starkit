import os

import numpy as np
from scipy.interpolate import interp1d

from astropy import units as u

from starkit.gridkit.io.process import BaseProcessGrid
from starkit.gridkit.io.bosz.base import convert_bz2_memmap
from starkit.gridkit.util import convolve_to_resolution

class CMFGENProcessGrid(BaseProcessGrid):
    """

    """

    R_initial = 300000
    R_initial_sampling=1
    def __init__(self, index, input_wavelength, meta, wavelength_start=0*u.angstrom, wavelength_stop=np.inf*u.angstrom,
                 R=5000.0, R_sampling=4):
        """

        Parameters
        ----------
        index : pandas.DataFrame
        input_wavelength : astropy.units.Quantity
        meta : pandas.Series
        wavelength_start : astropy.units.Quantity
        wavelength_stop : astropy.units.Quantity
        R : float
        R_sampling : integer
        """
        super(CMFGENProcessGrid, self).__init__(index, input_wavelength, meta, wavelength_start=wavelength_start,
                                                 wavelength_stop=wavelength_stop, R=R, R_sampling=R_sampling)




    def interp_wavelength(self, flux):
        cut_flux = flux[self.start_idx:self.stop_idx]
        processed_flux = convolve_to_resolution(cut_flux, self.R_initial, self.R_initial_sampling, self.R)
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
        return flux * np.pi
