import os

import numpy as np
from scipy.interpolate import interp1d
from tqdm import tqdm
from astropy import units as u
from astropy.table import Table

from starkit.gridkit.io.process import BaseProcessGrid
from starkit.gridkit.io.bosz.base import convert_bz2_memmap
from starkit.gridkit.io.gotberg23.base import convert_sed_memmap
from starkit.gridkit.util import convolve_to_resolution

class Gotberg23ProcessGrid(BaseProcessGrid):
    """

    """

    R_initial = 300000
    R_initial_sampling=1
    def __init__(
            self, index, input_wavelength, meta,
            wavelength_start=0*u.angstrom, wavelength_stop=np.inf*u.angstrom,
            R=5000.0, R_sampling=4,
        ):
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
        super(Gotberg23ProcessGrid, self).__init__(
            index, input_wavelength, meta,
            wavelength_start=wavelength_start, wavelength_stop=wavelength_stop,
            R=R, R_sampling=R_sampling,
        )
    
    def get_fluxes(self):
        """
        Load the fluxes from the files given in the raw index, process them, interpolate them and return them.
        This return float64 arrays
        Returns
        -------
            fluxes : numpy.ndarray
        """
        fluxes = np.empty((len(self.index), len(self.output_wavelength)), dtype=np.float64)

        for i, fname in tqdm(enumerate(self.index.filename),
                             total=len(self.index)):
            (cut_wavelength, flux) = self.load_flux(fname)
            fluxes[i] = self.interp_wavelength(flux, cut_wavelength)
        return fluxes
    
    def interp_wavelength(self, flux, cut_wavelength):
        cut_flux = flux[self.start_idx:self.stop_idx]
        processed_flux = convolve_to_resolution(
            cut_flux, self.R_initial, self.R_initial_sampling, self.R,
        )
        
        output_flux = interp1d(cut_wavelength, processed_flux)(
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
        
        sed_tab = Table.read(
            fname,
            format='ascii.commented_header',
        )
        
        # Crop wavelength
        input_wavelength = sed_tab['Wavelength'] * u.angstrom
        self.start_idx, self.stop_idx = self._get_wavelength_bounds(
            input_wavelength, self.wavelength_start, self.wavelength_stop,
        )
        cut_wavelength = input_wavelength[
            self.start_idx:self.stop_idx]
        
        # Get flux
        fname_npy = fname.replace('.txt', '.v1.npy')
        if not os.path.exists(fname_npy):
            convert_sed_memmap(fname)
        flux = np.load(fname_npy)
        
        # Return cut wavelength and loaded flux
        return (cut_wavelength, flux * np.pi)
