import os
import logging
from abc import ABCMeta, abstractmethod

import numpy as np

from sqlalchemy import create_engine
from sqlalchemy.orm import sessionmaker

import pandas as pd
import h5py
from specutils import Spectrum1D

from astropy import units as u

logger = logging.getLogger(__name__)

class BaseSpectralGridIO(object, metaclass=ABCMeta):
    """
    Base class for reading in the spectral grid information
    """

    GridBase = None




    def _set_grid_base_dir(self, base_dir):
        setattr(self.spectrum_table, 'base_dir', base_dir)
        self.base_dir = base_dir


    def _set_spectrum_wavelength(self, wavelength):
        setattr(self.spectrum_table, 'wavelength', wavelength)


    @staticmethod
    def _load_wavelength_solution():
        raise NotImplementedError('This needs to be implemented in subclasses')


    @staticmethod
    def _create_compound_model(models):
        if len(models) == 0:
            return lambda spectrum: spectrum
        else:
            compound_model = models[0]
            for model in models[1:]:
                compound_model = compound_model | model
            return compound_model

    def get_spectrum_query(self, filter_tuple):
            return self.session.query(self.spectrum_table).join(
            self.parameter_set_table).filter(*filter_tuple)

    
    def to_hdf(self, fname, index, fluxes, wavelength, grid_type,
               R, R_sampling, clobber=False):
        """
        Writing the grid to the Starkit HDF5 format
        """

        assert len(index) == len(fluxes)
        assert fluxes.shape[1] == len(wavelength)

        if os.path.exists(fname):
            if clobber:
                os.remove(fname)
            else:
                raise IOError('File {0} exists - '
                              'if you want overwrite set clobber=True'.format(fname))

        index.to_hdf(fname, 'index')

        with h5py.File(fname, 'a') as fh:
            fh['fluxes'] = fluxes
            fh['fluxes'].attrs['unit'] = str(fluxes.unit)
            fh['wavelength'] = wavelength.value
            fh['wavelength'].attrs['unit'] = str(wavelength.unit)
            fh['wavelength'].attrs['grid'] = 'log'
            fh['wavelength'].attrs['R'] =  R
            fh['wavelength'].attrs['R_sampling'] = R_sampling



class BaseSpectralGridDBIO(BaseSpectralGridIO):
    """
    Base spectral grid IO with using database connections
    """

    def __init__(self, db_url, base_dir,
                 wavelength_fname):
        self.engine = create_engine(db_url)
        self.GridBase.metadata.create_all(self.engine)
        self.GridBase.metadata.bind = self.engine

        self.session = sessionmaker(bind=self.engine)()

        self._set_grid_base_dir(base_dir)

        self.wavelength = self._load_wavelength_solution(wavelength_fname)

        self._set_spectrum_wavelength(self.wavelength)

    def get_query_data(self, filter_tuple, plugin,
                       warning_threshold=1 * u.gigabyte):
        """
        Write spectra to disk
        :param filter_tuple:
        :param models:
        :return:
        """

        query = self.get_spectrum_query(filter_tuple)

        sample_spectrum_row = query.first()
        sample_spectrum_flux = plugin(sample_spectrum_row.get_spectrum1d().flux)

        no_spectra = query.count()

        size_of_spectra = (query.count() *
                           len(sample_spectrum_flux)) * 8 * u.byte

        if size_of_spectra > warning_threshold:
            continue_query = input('The size of the spectra are {0:.2f}. '
                                       'Continue [y/N]'.format(
                size_of_spectra.to(warning_threshold.unit)))
            if continue_query.strip().lower() != 'y':
                raise ValueError('Size of requested grid ({:.2f}) to '
                                 'large for user ... aborting'.format(
                    size_of_spectra.to(warning_threshold.unit)))

        fluxes = np.empty((query.count(),
                          len(sample_spectrum_flux)))
        parameters = []
        param_names = [item.name
                       for item in sample_spectrum_row.parameter_set.parameters]

        for i, spectrum_row in enumerate(query):
            logger.info("{0} {1}/{2}".format(spectrum_row, i, no_spectra))
            spectrum = spectrum_row.get_spectrum1d()
            fluxes[i] = plugin(spectrum.flux)
            parameters.append([getattr(spectrum_row.parameter_set, key)
                               for key in param_names])

        parameters = pd.DataFrame(parameters, columns= param_names)
        output_sample_spectrum = Spectrum1D.from_array(
            plugin.output_wavelength * u.angstrom, sample_spectrum_flux)

        return output_sample_spectrum, parameters, fluxes

    def to_hdf(self, fname, filter_tuple, plugin, grid_type, flux_unit,
               R, R_sampling, clobber=False):

        sample_spectrum, index, fluxes = self.get_query_data(filter_tuple,
                                                                  plugin)

        super(BaseSpectralGridDBIO, self).to_hdf(
            fname, index, u.Quantity(fluxes, flux_unit, copy=False), sample_spectrum.wavelength,
            'log', R, R_sampling, clobber)
