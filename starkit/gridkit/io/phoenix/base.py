import os

import fnmatch

import logging

logger = logging.getLogger(__name__)

import h5py
from astropy.io import fits
from astropy import units as u

from starkit.gridkit.io.base import BaseSpectralGridIO, BaseSpectralGridDBIO
from starkit.gridkit.io.phoenix.alchemy import (Spectrum, ParameterSet,
                                                PhoenixBase)

from starkit.gridkit.io.phoenix.plugin import PhoenixProcess

class PhoenixSpectralGridIO(BaseSpectralGridDBIO):

    GridBase = PhoenixBase

    spectrum_table = Spectrum
    parameter_set_table = ParameterSet

    wavelength_unit = u.angstrom

    def __init__(self, db_url, base_dir,
                 wavelength_fname='WAVE_PHOENIX-ACES-AGSS-COND-2011.fits'):
        super(PhoenixSpectralGridIO, self).__init__(
            db_url=db_url, base_dir=base_dir, 
            wavelength_fname=os.path.join(base_dir, wavelength_fname))

    def _load_wavelength_solution(self, fname):
        return fits.getdata(fname, ext=0) * self.wavelength_unit

    def ingest(self):
        self.read_spectra()

    def read_spectra(self):
        spectral_files = self.get_spectral_files()
        no_of_spectra = len(spectral_files)
        for i, relative_path in enumerate(spectral_files):
            fname = os.path.basename(relative_path)
            fpath = os.path.relpath(os.path.dirname(relative_path),
                                    self.base_dir)
            logger.info("{0} / {1} importing {2}".format(i,
                                                         no_of_spectra, fname))
            spectrum = Spectrum(fname=fname, fpath=fpath)
            parameter_set = ParameterSet.from_file(spectrum.full_path)
            parameter_set.spectrum = spectrum
            self.session.add_all([spectrum, parameter_set])
        self.session.commit()

    def get_spectral_files(self):
        spectral_files = []
        for root, dirs, files in os.walk(self.base_dir):
            if 'Z' not in root:
                continue
            for filename in fnmatch.filter(files, 'lte*.fits'):
                spectral_files.append(os.path.join(root, filename))

        return spectral_files

    def to_hdf(self, fname, filter_tuple, R, wavelength_range,
               pre_sampling=2, sampling=4, clobber=False):
        """
        Writing the grid to HDF5

        Parameters

        """
        plugin = PhoenixProcess(self.wavelength.value, R, wavelength_range,
                                pre_sampling=pre_sampling, sampling=sampling)
        super(PhoenixSpectralGridIO, self).to_hdf(
            fname, filter_tuple, plugin, 'log', u.erg / u.s / u.angstrom,
            R, sampling, clobber)






