
import os
import gzip
import logging

from collections import OrderedDict

import numpy as np
import pandas as pd

from astropy import units as u

from starkit.gridkit.io.base import BaseSpectralGridIO
from starkit.gridkit.io.process import BaseProcessGrid



logger = logging.getLogger(__name__)

import re

class Munari1AProcess(BaseProcessGrid):
    pass


class Munari1A(BaseSpectralGridIO):

    R = 10000
    R_sampling = 4
    wavelength_fname = 'LAMBDA_D01.DAT'
    def __init__(self, base_dir):
        self.base_dir = base_dir

        if not os.path.exists(os.path.join(base_dir, self.wavelength_fname)):
            raise IOError('Please make sure that LAMBDA_D01.DAT file exists in '
                          'the root directory (downloadable at: '
                          'http://archives.pd.astro.it/2500-10500/SPECTRA/D01/'
                          'LAMBDA_D01.DAT)')
        else:
            self.wavelength = u.angstrom * np.loadtxt(os.path.join(base_dir,
                                                      self.wavelength_fname))


        self.fname_pattern = re.compile(
            'T(\d+)G(\d+)([PM])(\d+)V(\d+)K(\d+)([SA])(OD|NW)NVD01F\.ASC\.gz')
        self.process = Munari1AProcess(self.wavelength, self.R,
                                       sampling=self.R_sampling)

    def get_wavelength(self):
        np.loadtxt(os.path.join(self.base_dir, ''))

    def get_index(self, odf='prefer_new', alpha=0.0, k=2):
        """

        Parameters
        ----------

        odf: ~str
            "prefer_new": will remove duplicates by dropping the ones with
                old ODFs in favours of the new ones and then dropping the column
                ODF
            None: will not remove any duplicates

        """
        file_list = self._generate_fname_list()

        column_names = ['teff', 'logg', 'mh', 'alpha', 'k', 'odf_new', 'fname']

        index_raw = OrderedDict([(col_name, []) for col_name in column_names])

        for fname in file_list:
            teff, logg, mh, alpha, k, odf_new = self._parse_filename(fname)
            index_raw['teff'].append(teff)
            index_raw['logg'].append(logg)
            index_raw['mh'].append(mh)
            index_raw['alpha'].append(alpha)
            index_raw['k'].append(k)
            index_raw['odf_new'].append(odf_new)
            index_raw['fname'].append(fname)

        index = pd.DataFrame(index_raw)
        index.sort(index.columns.tolist()[:-1], inplace = True)

        if odf == 'prefer_new':
            index = index.drop_duplicates(['teff', 'logg', 'mh', 'alpha', 'k'], take_last=False)
        elif odf is None:
            pass
        else:
            raise NotImplementedError('Method {0} not implemented'.format(odf))

        if alpha is not None:
            index = index[index.alpha == alpha].drop('alpha', 1)

        if k is not None:
            index = index[index.k == k].drop('k', 1)


        return index


    def get_fluxes(self, index):
        file_list = index.fname
        fluxes = []
        len_file_list = len(file_list)
        for i, fname in enumerate(file_list):
            logger.info('[{0}/{1}] Reading spectrum {2}'.format(i,
                                                                len_file_list,
                                                                fname))
            with gzip.GzipFile(fname) as fh:
                fluxes.append(self.process.process(np.loadtxt(fh)))

        return np.array(fluxes)


    def _parse_filename(self, fname):
        teff, logg, mh_sign, mh, vrot, k, alpha_enhanced, odf_str = (
            self.fname_pattern.match(os.path.basename(fname)).groups())

        teff = float(teff)
        logg = float(logg) / 10
        mh = (-1 if mh_sign == 'M' else 1) * float(mh) / 10
        vrot = float(vrot)
        k = float(k)
        if alpha_enhanced == 'A':
            alpha = 0.4
        elif alpha_enhanced == 'S':
            alpha = 0.0
        else:
            raise ValueError('Got {0} for alphaenhanced string - only '
                             'A and S recognized'.format(alpha_enhanced))

        if odf_str == 'NW':
            odf_new = True
        elif odf_str == 'OD':
            odf_new = False
        else:
            raise ValueError('Got {0} for ODF modeltype string - '
                             'only NW and OD recognized')

        return teff, logg, mh, alpha, k, odf_new

    def get_wavelength(self):
        return self.process.output_wavelength * u.angstrom

    def _generate_fname_list(self):
        fname_list = []
        for dirpath, subdirs, files in os.walk(self.base_dir):
            for fname in files:
                if 'V000' in fname:
                    fname_list.append(os.path.join(dirpath, fname))

        return fname_list

    def to_hdf(self, fname, index_odf='prefer_new', clobber=False):
        index = self.get_index(odf=index_odf)
        fluxes = self.get_fluxes(index)

        default_values = {}
        if 'alpha' in index.columns:
            default_values['alpha'] = 0.0

        if 'k' in index.columns:
            default_values['k'] = 2.0

        super(Munari1A, self).to_hdf(fname, index.drop(['fname', 'odf_new'], 1),
                                     fluxes, self.get_wavelength(),
                                     'log', self.R, self.R_sampling,
                                     clobber=clobber,
                                     default_values=default_values)
