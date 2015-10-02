import os
import re

from glob import glob
import logging
import pandas as pd
import numpy as np


from scipy import interpolate
from astropy import units as u, constants as const
from astropy.modeling import Parameter

from scipy.spatial import cKDTree as KDTree
from starkit.base.model import StarKitModel
from starkit.gridkit.base import load_grid


logger = logging.getLogger(__name__)


def _read_parsec_fname(fname):
    fname_pattern = re.compile(
        '(Z\d+\.\d+)(Y\d+\.\d+).+\_(M\d+\.\d+).+\.DAT')

    Z_str, Y_str, M_str = fname_pattern.match(
        os.path.basename(fname)).groups()
    parsec_data = pd.read_csv(fname, delim_whitespace=True)

    Z_solar = 0.01524
    Y_solar = 0.2485 + 1.78 * 0.0152
    X_solar = 1 - (Y_solar + Z_solar)
    Z = float(Z_str[1:])
    Y = 0.2485 + 1.78 * Z
    X = 1 - (Y + Z)

    R = u.Quantity(10 ** parsec_data['LOG_R'].values, u.cm)
    M = u.Quantity(parsec_data['MASS'].values, u.M_sun)

    log_g = np.log10((const.G * M / R**2).cgs.value)
    parsec_data['X'] = X
    parsec_data['Y'] = Y
    parsec_data['Z'] = Z
    parsec_data['MH'] = np.log10((Z / X) / (Z_solar / X_solar))
    parsec_data['LOG_G'] = log_g
    parsec_data['TEFF'] = 10**parsec_data['LOG_TE']

    return Z_str + Y_str + M_str, parsec_data


def create_parsec_store(fname, parsec_dir, format='raw'):

    if os.path.exists(fname):
        raise IOError('Parsec HDF5 already exists {0}'.format(fname))

    parsec_store = pd.HDFStore(fname, mode='w')
    if format == 'joined':
        full_ev_data = None
    for parsec_dir in glob(os.path.join(parsec_dir, 'Z*Y*')):
        logger.info('Working in directory {0}'.format(parsec_dir))
        if not os.path.isdir(parsec_dir):
            logger.info('{0} not a directory - skipping'.format(parsec_dir))
            continue
        for parsec_fname in glob(os.path.join(parsec_dir, '*.DAT')):
            if 'ADD' in parsec_fname:
                continue
            if format == 'raw':
                key, ev_data = _read_parsec_fname(parsec_fname)
                parsec_store[os.path.join('parsec', key)] = ev_data
                parsec_store.flush()
            elif format == 'joined':

                key, ev_data = _read_parsec_fname(parsec_fname)
                if full_ev_data is None:
                    full_ev_data = ev_data
                else:
                    full_ev_data = pd.concat([full_ev_data, ev_data])
            else:
                raise ValueError('format only allows for raw or joined')
    if format == 'joined':
        parsec_store['parsec'] = full_ev_data
    parsec_store.close()
    logger.info('Wrote Database to {0}'.format(fname))

class Parsec(StarKitModel):

    mh = Parameter()
    mass = Parameter()
    age = Parameter()

    inputs = ()
    outputs = ('teff', 'logg', 'lum')

    def __init__(self, parsec_store, mh=0.0, mass=1.0, age=5e9):
        super(Parsec, self).__init__(mh, mass, age)
        try:
            self.parsec_store = pd.HDFStore(parsec_store)
        except TypeError:
            self.parsec_store = parsec_store

        self.evolution_data = [self.parsec_store[key]
                               for key in self.parsec_store.keys()]
        self.parsec_store.close()
        self.mh_mass = np.empty((len(self.evolution_data), 2))

        for i, ev_data in enumerate(self.evolution_data):
            mh = ev_data['MH'][0]
            mass = ev_data['MASS'][0]
            self.mh_mass[i] = mh, mass

        self.mh_mass_kd_tree = KDTree(self.mh_mass)


    def evaluate(self, mh, mass, age):
        distance, idx = self.mh_mass_kd_tree.query(
            np.array([mh, mass]).squeeze())
        ev_data = self.evolution_data[idx]
        age_ev_data = ev_data['AGE'].values
        out_ev_data = ev_data[['TEFF', 'LOG_G', 'LOG_L']].values
        teff, logg, log_l = interpolate.interp1d(
            age_ev_data, out_ev_data.T, bounds_error=False)(np.squeeze(age))
        return teff, logg, 10**log_l


class ParsecSpectralGrid(Parsec):

    outputs = ('wavelength', 'flux')

    def __init__(self, parsec_store, spectral_grid):
        super(ParsecSpectralGrid, self).__init__(parsec_store)
        self.spectral_grid = load_grid(spectral_grid)




