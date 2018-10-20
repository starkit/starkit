import logging

from astropy import modeling
from astropy import units as u, constants as const
import pandas as pd
from exceptions import ValueError, DeprecationWarning
import h5py
from scipy import interpolate
import numpy as np
import warnings
from astropy.modeling import Parameter
from starkit.fitkit.priors import UniformPrior
from starkit.utils.vacuumair_conversion import (convert_air2vacuum,
                                                convert_vacuum2air)
from starkit.gridkit.util import convolve_to_resolution
import starkit

wavelength_conversions = {'air2vacuum':convert_air2vacuum,
                          'vacuum2air':convert_vacuum2air}
logger = logging.getLogger(__name__)

class BaseSpectralGrid(modeling.Model):
    inputs = tuple()
    outputs = ('wavelength', 'flux')

    def __init__(self, wavelength, index, fluxes, meta, wavelength_type=None, **kwargs):
        super(BaseSpectralGrid, self).__init__(**kwargs)
        if meta['grid_type'] != 'log':
            raise ValueError('No other grid_type than log is supported - meta specifies {0} however'.format(
                meta['grid_type']))
        self.R = meta['R']
        self.R_sampling = meta['R_sampling']

        self.meta_grid = meta
        self.index = index
        self.fluxes = fluxes

        self.setup_wavelength_type(wavelength, wavelength_type)

        self.interpolator = self._generate_interpolator(index, fluxes)

    @property
    def wavelength_type(self):
        return self._wavelength_type

    @wavelength_type.setter
    def wavelength_type(self, wavelength_type):
        if wavelength_type not in ['air', 'vacuum']:
            raise ValueError("Wavelength_type can either be 'vacuum' or 'air' not "
                             "{0}".format(wavelength_type))
        self._wavelength_type = wavelength_type

        self.wavelength = getattr(self, 'wavelength_{0}'.format(self._wavelength_type)).copy()

    def setup_wavelength_type(self, wavelength, wavelength_type):
        """
        Setting the initial wavelength_type and setting the attributes
        wavelength_air and wavelength_vacuum.

        Parameters
        ----------
        wavelength: astropy.units.Quantity
        wavelength_type: str

        Returns
        -------
            : None

        """



        if wavelength_type is None:
            logger.warn("**** NO WAVELENGTH TYPE SET DEFAULTING TO GRID ({0}) ****\n\n".format(
                self.meta_grid['wavelength_type']))
            wavelength_type = self.meta_grid['wavelength_type']

        setattr(self, 'wavelength_{0}'.format(self.meta_grid['wavelength_type']), wavelength)

        wavelength_types = {'air', 'vacuum'}
        wavelength_types.remove(self.meta_grid['wavelength_type'])
        non_grid_wavelength_type = wavelength_types.pop()
        non_grid_wavelength = wavelength_conversions['{0}2{1}'.format(self.meta_grid['wavelength_type'],
                                                                      non_grid_wavelength_type)](wavelength)

        setattr(self, 'wavelength_{0}'.format(non_grid_wavelength_type), non_grid_wavelength)
        self.wavelength_type = wavelength_type

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
        wavelength = self.wavelength.value
        flux = np.squeeze(self.interpolator(np.squeeze(np.array(args))))
        return wavelength, flux

    @staticmethod
    def _generate_interpolator(index, fluxes):
        return interpolate.LinearNDInterpolator(index, fluxes.value)

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

    def __init__(self, wavelength, index, fluxes, meta, wavelength_type=None, target_wavelength=None,
                 target_wavelength_type=None, target_R=None, vrad_telluric=0.0, **kwargs):

        super(BaseTelluricGrid, self).__init__(wavelength, index, fluxes, meta, wavelength_type=wavelength_type,
                                               vrad_telluric=vrad_telluric, **kwargs)
        self.raw_wavelength = wavelength.copy()
        self.raw_fluxes = fluxes.copy()

        if target_R is not None:
            self.reconvolve_fluxes(target_R)

        if target_wavelength is not None:
            if target_wavelength_type is None:
                raise ValueError('target_wavelength_type can not by None if target_wavelength is specified')
            self.wavelength_type = target_wavelength_type
            self.resample_fluxes(target_wavelength)

        self.interpolator = self._generate_interpolator(index, self.fluxes)
        self.c_in_kms = const.c.to(u.km / u.s).value

    def reconvolve_fluxes(self, target_R):
        """
        Reconvolve Fluxes to target

        Parameters
        ----------
        target_R: float

        Returns
        -------

        """

        new_fluxes = np.empty_like(self.raw_fluxes.value)
        for i, flux in enumerate(self.raw_fluxes.value):
            new_fluxes[i] = convolve_to_resolution(flux, self.R, self.R_sampling, target_R)

        self.fluxes = new_fluxes * self.fluxes.unit
        self.R = target_R
        self.R_sampling = None


    def resample_fluxes(self, target_wavelength):
        """
        Resample Telluric absorption on target wavelength_grid.

        Parameters
        ----------
        target_wavelength: astropy.units.Quantity
            target wavelength_grid

        """
        wavelength_interp = interpolate.interp1d(self.wavelength, self.raw_fluxes, bounds_error=False,
                                                 fill_value=np.nan)

        self.fluxes = wavelength_interp(target_wavelength) * self.fluxes.unit
        self.setup_wavelength_type(target_wavelength, self.wavelength_type)


    def evaluate_raw(self, *args):
        return self.wavelength.value, np.squeeze(self.interpolator(np.squeeze(np.array(args))))

    def evaluate_transmission(self, wavelength, flux, vrad_telluric, *args):


        beta = vrad_telluric / self.c_in_kms
        doppler_factor = np.sqrt((1 + beta) / (1 - beta))

        transmission = np.squeeze(self.interpolator(np.squeeze(np.array(args))))
        doppler_transmission_interp = interpolate.interp1d(wavelength * doppler_factor, transmission,
                                                           bounds_error=False, fill_value=np.nan)

        return wavelength, flux * doppler_transmission_interp(wavelength)






def construct_grid_class_dict(meta, index, class_dict=None):
    """
    Construct a SpectralGrid class with the parameter attributes and initial parameters set.

    Parameters
    ----------
    meta: pandas.Series
        normal grid metadata
    index: pandas.Dataframe
        Grid index with parameters
    base_class: type
        Base class to use to construct the grid class

    Returns
    -------
    class_dict : dict
    initial_parameters: dict
        dictionary for initializing the base_class
    """

    interpolation_parameters = meta['parameters']
    initial_parameters = {item: index[item].iloc[0]
                          for item in interpolation_parameters}

    parameter_defaults = {param: index.loc[index.index[0], param]
                          for param in interpolation_parameters}

    parameter_bounds = {param: (index[param].min(), index[param].max())
                        for param in interpolation_parameters}

    if class_dict is None:
        class_dict = {}

    for param in interpolation_parameters:
        cur_param_default = parameter_defaults.get(param, None)
        cur_param_bound = parameter_bounds.get(param, None)
        param_descriptor = Parameter(default=cur_param_default,
                                     bounds=cur_param_bound)
        class_dict[param] = param_descriptor
    if u.Unit(meta['flux_unit']) == u.erg / u.s / u.cm ** 2 / u.angstrom:
        radius_descriptor = Parameter(default=1, bounds=[1e-6, 1e9])
        class_dict['radius'] = radius_descriptor


    return class_dict, initial_parameters


def read_grid(hdf_fname):
    """
    Load the grid information from hdf

    Parameters
    ----------
    hdf_fname: str
        filename and path to the HDF file

    Returns
    -------
        wavelength : astropy.units.Quantity
        meta : pandas.Series
        index : pandas.DataFrame
        fluxes : astropy.units.Quantity

    """

    logger.info('Reading index')
    index = pd.read_hdf(hdf_fname, 'index')
    meta = pd.read_hdf(hdf_fname, 'meta')
    logger.info('Discovered columns {0}'.format(', '.join(meta['parameters'])))

    with h5py.File(hdf_fname) as fh:
        logger.info('Reading Fluxes')
        fluxes = fh['fluxes'].__array__()

    logger.info('Fluxes shape {0}'.format(fluxes.shape))
    flux_unit = u.Unit(meta['flux_unit'])
    wavelength = pd.read_hdf(hdf_fname, 'wavelength').values[:, 0]
    wavelength = u.Quantity(wavelength, meta['wavelength_unit'])

    return wavelength, meta, index, fluxes * flux_unit

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

    wavelength, meta, index, fluxes = read_grid(hdf_fname)
    class_dict = {}
    class_dict['__init__'] = base_class.__init__

    class_dict, initial_parameters = construct_grid_class_dict(meta, index, class_dict=class_dict)


    SpectralGrid = type('SpectralGrid', (base_class,), class_dict)

    logger.info('Initializing spec grid')

    spec_grid = SpectralGrid(wavelength, index[meta['parameters']].values, fluxes, meta,
                             wavelength_type=wavelength_type,  **initial_parameters)

    if 'format_version' not in spec_grid.meta_grid.keys():
        logger.warn('No format_version in meta data for this grid. Please get an'+\
                    ' updated grid. This will fail in the future.')
        warnings.warn('No format_version in meta data for this grid. Please get '+\
                      'an updated grid. This will fail in the future.', DeprecationWarning)
    else:
        current_version = starkit.gridkit.FORMAT_VERSION
        grid_version = spec_grid.meta_grid['format_version']
        if current_version[1] != grid_version[1]:
            raise ValueError('Grid major versions do not match! Curent code '+\
                             'format version: {0} grid format version: {1} '.format(current_version, grid_version))
        
    return spec_grid

def load_telluric_grid(hdf_fname, stellar_grid=None, wavelength_type=None, base_class=BaseTelluricGrid):
    """
    Load the grid from an HDF file

    Parameters
    ----------
    hdf_fname: ~str
        filename and path to the HDF file

    stellar_grid: BaseSpectralGrid
        spectral_grid to adapt to

    wavelength_type: str
        use 'air' or 'vacuum' wavelength and convert if necessary (by inspecting
        what the grid uses in meta)

    Returns
    -------
        : SpectralGrid object

    """

    wavelength, meta, index, fluxes = read_grid(hdf_fname)
    class_dict = {}
    class_dict['vrad_telluric'] = modeling.Parameter(default=0.0)
    class_dict, initial_parameters = construct_grid_class_dict(meta, index, class_dict=class_dict)
    class_dict['__init__'] = base_class.__init__

    if stellar_grid is not None:
        class_dict['inputs'] = ('wavelength', 'flux')
        class_dict['evaluate'] = base_class.evaluate_transmission
        initial_parameters['target_wavelength'] = stellar_grid.wavelength
        initial_parameters['target_R'] = stellar_grid.R
        if wavelength_type is not None and stellar_grid is not None:
            logger.warn('Ignoring requested wavelength_type in favour of spectral_grid wavelength_type')
        initial_parameters['wavelength_type'] = meta['wavelength_type']
        initial_parameters['target_wavelength_type'] = stellar_grid.meta_grid['wavelength_type']

    else:
        class_dict['inputs'] = tuple()
        class_dict['evaluate'] = base_class.evaluate_raw
        initial_parameters['wavelength_type'] = wavelength_type

    TelluricGrid = type('TelluricGrid', (base_class,), class_dict)

    logger.info('Initializing spec grid')

    telluric_grid = TelluricGrid(wavelength, index[meta['parameters']].values, fluxes, meta, **initial_parameters)

    return telluric_grid
