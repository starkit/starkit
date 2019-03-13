import os

import h5py
import pandas as pd
import numpy as np
import uuid
from astropy import units as u
from tqdm import tqdm
from scipy import ndimage as nd
from scipy.interpolate import interp1d
import starkit

class BaseProcessGrid(object):

    def __init__(self, index, input_wavelength, meta, wavelength_start=0*u.angstrom, wavelength_stop=np.inf*u.angstrom,
                 R=None, R_sampling=4):

        self.index = index
        self.meta = meta
        self.R = float(R)
        self.R_sampling = int(R_sampling)


        wavelength_start = np.max(u.Quantity([wavelength_start, np.min(input_wavelength)]))
        wavelength_stop = np.min(u.Quantity([wavelength_stop, np.max(input_wavelength)]))
        self.wavelength_start = wavelength_start
        self.wavelength_stop = wavelength_stop

        input_wavelength = input_wavelength.to('angstrom')

        self.start_idx, self.stop_idx = self._get_wavelength_bounds(
            input_wavelength, wavelength_start, wavelength_stop)
        self.cut_wavelength = input_wavelength[
                              self.start_idx:self.stop_idx]

        self.output_wavelength = self._logrange(
            wavelength_start, wavelength_stop, R, R_sampling)






    @staticmethod
    def _get_wavelength_bounds(wavelength, wavelength_start, wavelength_stop):
        start_idx = np.max([wavelength.searchsorted(wavelength_start) - 5, 0])
        stop_idx = np.min([wavelength.searchsorted(wavelength_stop) + 5, len(wavelength)])
        return start_idx, stop_idx

    @staticmethod
    def _logrange(wavelength_start, wavelength_stop, R, sampling):
        return np.exp(
            np.arange(np.log(u.Quantity(wavelength_start, 'angstrom').value),
                      np.log(u.Quantity(wavelength_stop, 'angstrom').value),
                      1 / (sampling * float(R))))

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
            flux = self.load_flux(fname)
            fluxes[i] = self.interp_wavelength(flux)
        return fluxes

    def interp_wavelength(self, flux):
        cut_flux = flux[self.start_idx:self.stop_idx]
        rescaled_R = 1 / np.sqrt((1/self.R)**2 - (1/self.R_initial)**2)
        sigma = ((self.R_initial / rescaled_R) * self.R_initial_sampling
                 / (2 * np.sqrt(2 * np.log(2))))

        processed_flux = nd.gaussian_filter1d(cut_flux, sigma)
        output_flux = interp1d(self.cut_wavelength, processed_flux)(
            self.output_wavelength)
        return output_flux


    def get_meta(self):
        """
        Get the meta data and add processing information on

        Returns
        -------
            : pandas.Series
        """
        meta = self.meta.copy()
        parameters = []

        for param in meta['parameters']:
            if len(self.index[param].unique()) > 1:
                parameters.append(param)
        meta['parameters'] = parameters
        meta['flux_unit'] = 'erg/s/angstrom'
        meta['R'] = self.R
        meta['R_sampling'] = self.R_sampling
        meta['grid_type'] = 'log'
        meta['uuid'] = str(uuid.uuid4())
        meta['format_version'] = starkit.gridkit.FORMAT_VERSION
        return meta

    def get_index(self):
        """
        Remove the filenames from the index
        Returns
        -------
            : pandas DataFrame
        """

        index = self.index.copy()
        return index.drop('filename', axis=1)

    def to_hdf(self, fname, overwrite=False):
        """

        Parameters
        ----------
        fname
        overwrite : bool

        Returns
        -------

        """

        if os.path.exists(fname) and not overwrite:
            raise IOError('File {0} exists and overwrite=False ... aborting'.format(fname))
        if os.path.exists(fname) and overwrite:
            os.system('rm -f {0}'.format(fname))

        fluxes = self.get_fluxes()
        meta = self.get_meta()
        index = self.get_index()

        with h5py.File(fname, mode='w') as fh:
            fh['fluxes'] = fluxes

        meta.to_hdf(fname, 'meta')
        index.to_hdf(fname, 'index')
        pd.DataFrame(self.output_wavelength).to_hdf(fname, 'wavelength')
        
