import os

import numpy as np
import pandas as pd
from astropy import units as u

import starkit
from starkit.gridkit.io.base import BaseSpectralGridIO

gaia_data_dir = os.path.join(starkit.__path__[0], 'gridkit', 'io', 'gaia',
                             'data')

speclib_column_description = pd.read_csv(os.path.join(
    gaia_data_dir, 'speclib_column_description.csv'), header=0)

speclib_dict = {name:index for name, index in zip(
        speclib_column_description.name, speclib_column_description.index)}

class GaiaSpecLibIO(BaseSpectralGridIO):

    R=None
    R_sampling=4

    flux_unit = u.W / u.m**2 / u.nm
    wavelength_unit = u.angstrom

    def __init__(self, fname):
        self.fname = fname
        self.raw_data = open(fname).readlines()


    def get_fluxes(self):
        fluxes = [tuple(float(item) for item in self.raw_data[i].split())
                  for i in range(1, len(self.raw_data), 2)
                  if self.raw_data[i].strip() != '']
        return np.array(fluxes) * u.W / u.m**2 / u.nm



    def get_grid_wavelength(self):
        raw_index = self.raw_data[0].split()
        wavelength_start = float(raw_index[speclib_dict['lambda0']]) * 10
        wavelength_end = float(raw_index[speclib_dict['lambdaEnd']]) * 10
        wavelength_step = float(raw_index[speclib_dict['dlambda']]) * 10

        return np.arange(wavelength_start, wavelength_end + wavelength_step,
                         wavelength_step) * self.wavelength_unit


    def get_index(self, index_columns=['teff', 'logg', 'zmetal']):
        raw_indices = [self.raw_data[i]
                       for i in range(0, len(self.raw_data), 2)
                       if self.raw_data[i].strip() != '']
        index = []

        for row in raw_indices:
            index.append([float(row.split()[speclib_dict[meta_param]])
                          for meta_param in index_columns])

        index = pd.DataFrame(data=index, columns=index_columns)

        return index.replace(-999.0, np.nan)



    def to_hdf(self, fname, index_columns, grid_type='log', clobber=False):
        index = self.get_index(index_columns)
        fluxes = self.get_fluxes()

        wavelength = self.get_wavelength()

        super(GaiaSpecLibIO, self).to_hdf(fname, index, fluxes, wavelength,
                                          R=self.R, R_sampling=self.R_sampling,
                                          clobber=clobber, grid_type=grid_type)





