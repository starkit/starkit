import logging

from astropy import modeling
from astropy import units as u, constants as const
import pandas as pd
import h5py
from scipy import interpolate
import numpy as np

from astropy.modeling import Parameter
from starkit.fitkit.priors import UniformPrior
from starkit.utils.vacuumair_conversion import (convert_air2vacuum,
                                                convert_vacuum2air)


wavelength_conversions = {'air2vacuum':convert_air2vacuum,
                          'vacuum2air':convert_vacuum2air}
logger = logging.getLogger(__name__)

class BaseSpectralGrid(modeling.Model):
    inputs = tuple()
    outputs = ('wavelength', 'flux')

    def __init__(self, wavelength, index, fluxes, **kwargs):
        self.R = kwargs.pop('R', None)
        self.R_sampling = kwargs.pop('R_sampling', None)
        self.flux_unit = kwargs.pop('flux_unit', None)
        self.index = index
        self.fluxes = fluxes

        super(BaseSpectralGrid, self).__init__(**kwargs)
        self.interpolator = self._generate_interpolator(index, fluxes)
        self.wavelength = wavelength

    def get_grid_extent(self):
        extents = []
        for i, param_name in enumerate(self.param_names):
            extents.append((self.interpolator.points[:, i].min(),
                            self.interpolator.points[:, i].max()))
        return extents

    def get_grid_uniform_priors(self):
        extents = self.get_grid_extent()
        priors = []
        for extent in extents:
            priors.append(UniformPrior(*extent))

        return priors

    def evaluate(self, *args):

        return self.wavelength.value, self.interpolator(np.array(
            args).reshape(len(self.param_names)))[0]

    @staticmethod
    def _generate_interpolator(index, fluxes):
        return interpolate.LinearNDInterpolator(index, fluxes)

    @property
    def velocity_per_pix(self):
        if self.R is None:
            raise ValueError('R and R_sampling for the current grid not known')

        else:
            return const.c / self.R / self.R_sampling

    def _renormalize_grid(self):
        for i in xrange(self.fluxes.shape[0]):
            self.fluxes[i] /= np.trapz(self.fluxes[i], self.wavelength)

class BaseTelluricGrid(BaseSpectralGrid):
    input = ('wavelength', 'flux')

    vrad = modeling.Parameter()

    def __init__(self, wavelength, index, fluxes, vrad=0.0, **kwargs):
        super(BaseTelluricGrid, self).__init__(wavelength, index, fluxes, vrad=vrad, **kwargs)
        self.raw_fluxes = self.fluxes
    def evaluate(self, *args):
        return self.wavelength.value, self.interpolator(np.array(
            args).reshape(len(self.param_names)))[0]

def load_telluric_grid(hdf_fname, wavelength_type=None):
    return load_grid(hdf_fname, wavelength_type=wavelength_type,
                     base_class=BaseTelluricGrid)

def load_grid(hdf_fname, wavelength_type=None, base_class=BaseSpectralGrid):
    """
    Load the grid from an HDF file

    Parameters
    ----------
    hdf_fname: ~str
        filename and path to the HDF file
    wavelength_type: str
        use 'air' or 'vacuum' wavelength and convert if necessary (by inspecting
        what the grid uses in meta)

    Returns
    -------
        : SpectralGrid object

    """
    logger.info('Reading index')
    index = pd.read_hdf(hdf_fname, 'index')
    meta = pd.read_hdf(hdf_fname, 'meta')
    logger.info('Discovered columns {0}'.format(', '.join(meta['parameters'])))
    interpolate_parameters = meta['parameters']

    with h5py.File(hdf_fname) as fh:
        logger.info('Reading Fluxes')
        fluxes = fh['fluxes'].__array__()
    logger.info('Fluxes shape {0}'.format(fluxes.shape))
    flux_unit = u.Unit(meta['flux_unit'])
    wavelength = pd.read_hdf(hdf_fname, 'wavelength').values[:, 0]
    wavelength = u.Quantity(wavelength, meta['wavelength_unit'])

    if wavelength_type is None:
        logger.warn("**** NO WAVELENGTH TYPE SET DEFAULTING TO GRID ({0}) ****\n\n".format(
            meta['wavelength_type']))
        wavelength_type = meta['wavelength_type']
    if wavelength_type not in ['air', 'vacuum']:
        raise ValueError("Wavelength_type can either be 'vacuum' or 'air' not "
                         "{0}".format(wavelength_type))
    if not meta['wavelength_type'] == wavelength_type:
        wave_conv = wavelength_conversions['{0}2{1}'.format(
            meta['wavelength_type'], wavelength_type)]

        wavelength = wave_conv(wavelength)


    if meta['grid_type'] == 'log':
        R = meta['R']
        R_sampling = meta['R_sampling']
    else:
        raise ValueError('No other grid_type than log is supported')

    parameter_defaults = {param: index.loc[index.index[0], param]
                            for param in interpolate_parameters}

    class_dict = {}
    for param in interpolate_parameters:
        if parameter_defaults[param] is None:
            param_descriptor = Parameter()
        else:
            param_descriptor = Parameter(default=parameter_defaults[param])

        class_dict[param] = param_descriptor

    class_dict['__init__'] = base_class.__init__

    SpectralGrid = type('SpectralGrid', (base_class, ), class_dict)

    initial_parameters = {item: index[item].iloc[0]
                          for item in interpolate_parameters}
    logger.info('Initializing spec grid')
    spec_grid = SpectralGrid(wavelength, index[interpolate_parameters], fluxes,
                        R=R, R_sampling=R_sampling, flux_unit=flux_unit,
                        **initial_parameters)

    logger.info('Setting grid extent')
    for param in interpolate_parameters:
        uniform_prior = UniformPrior(index[param].min(), index[param].max())
        parameter = getattr(spec_grid, param)
        parameter.prior = uniform_prior

    return spec_grid
