import os
import time
import types
import tempfile
from collections import OrderedDict
from logging import getLogger
import shutil

from scipy import stats
import pandas as pd

from starkit.fitkit.priors import PriorCollection

logger = getLogger(__name__)

import numpy as np


try:
    import pymultinest
except ImportError:
    multinest_available = False
    raise
else:
    multinest_available = True

def multinest_evaluate(self, model_param, ndim, nparam):
    # returns the likelihood of observing the data given the model param_names
    model_param = np.array([model_param[i] for i in xrange(nparam)])
    parameters = self.parameters.copy()
    parameters[~self.fixed_mask()] = model_param
    loglikelihood = self.evaluate(*parameters)
    return float(loglikelihood)

def fixed_mask(self):
    return np.array([getattr(self, param_name).fixed
                      for param_name in self.param_names])



class MultiNestResult(object):


    @classmethod
    def from_multinest_basename(cls, basename, parameter_names,equal_weights=False):
        """
        Reading a MultiNest result from a basename

        Parameters
        ----------

        basename: str
            basename (path + prefix) for a multinest run

        Keywords
        --------
        equal_weights - load the equally weighted chains instead

        Returns
            : ~MultinestResult
        """
        if equal_weights:
            posterior_data = cls.read_equal_posterior_data(basename, parameter_names)
        else:
            posterior_data = cls.read_posterior_data(basename, parameter_names)

        return cls(posterior_data)

    @classmethod
    def from_hdf5(cls, h5_fname, key='multinest'):
        """
        Reading a Multinest result from its generated HDF5 file

        Parameters
        ----------

        h5_fname: ~str
            HDF5 filename

        key: ~str
            group identifier in the store
        """

        posterior_data = pd.read_hdf(h5_fname, key)

        return cls(posterior_data)

    @staticmethod
    def read_posterior_data(basename, parameter_names):
        """
        Reading the posterior data into a pandas dataframe

        Multinest weighted posterior file blah_.txt has the following format
        weights, log likelihood, parameters
        """

        posterior_data = pd.read_csv(
            '{0}_.txt'.format(basename),
            delim_whitespace=True, names=['weights']+['loglikelihood']+parameter_names)
        posterior_data.index = np.arange(len(posterior_data))
        return posterior_data

    @staticmethod
    def read_equal_posterior_data(basename, parameter_names):
        """
        Reading the posterior data into a pandas dataframe

        """

        posterior_data = pd.read_csv(
            '{0}_post_equal_weights.dat'.format(basename),
            delim_whitespace=True, names=parameter_names + ['x'])
        posterior_data.index = np.arange(len(posterior_data))

        # since the chain is equally weighted, we should just put equal weights
        # in
        posterior_data['weights'] = np.zeros(len(posterior_data))+1.0/float(len(posterior_data))

        return posterior_data

    def __init__(self, posterior_data):
        self.posterior_data = posterior_data
        self.parameter_names = [col_name for col_name in posterior_data.columns
                                if col_name not in ['x','weights','loglikelihood']]


    @property
    def mean(self):
        return self.posterior_data[self.parameter_names].mean()

    @property
    def median(self):
        return self.posterior_data[self.parameter_names].median()

    @property
    def maximum(self):
        # returns the maximum in the posterior
        max_ind = self.posterior_data.loglikelihood.argmin() # this should also be the maxium in weight for the non-equal weighted points
        return self.posterior_data[self.parameter_names].iloc[max_ind]

    def __repr__(self):
        return "<MultiNest Result (median)\n{0}>".format(self.median.__repr__())

    def calculate_sigmas(self, sigma_number):
        sigma_dict = []
        for param_name in self.parameter_names:

            # sort the parameter in order to create the CDF
            param_x = np.copy(self.posterior_data[param_name])

            weights = np.copy(self.posterior_data['weights'])
            ind = np.argsort(param_x)
            param_x = np.array(param_x[ind])
            weights = np.array(weights[ind])
            #k = [np.sum(weights[0:i+1]) for i in xrange(len(weights))]

            # make CDF of the weights to determine sigmas later
            k = np.cumsum(weights)

            sigma_lower = np.interp(stats.norm.cdf(-sigma_number), k, param_x)
            sigma_upper = np.interp(stats.norm.cdf(sigma_number), k, param_x)
            sigma_dict.append((param_name, (sigma_lower, sigma_upper)))
        return OrderedDict(sigma_dict)

    def plot_triangle(self, parameters = None, **kwargs):
        '''
        Produce a corner plot of the chains posterior.

        Keywords
        --------
        parameters - a list of paramters to plot. By default, it will plot
                     all fit parameters. This is useful if you run into problems
                     where one of the fit paramters is fixed and corner.py does
                     not work on it
        '''
        try:
            from corner import corner
        except ImportError:
            raise ImportError('Plotting requires corner.py')
        if parameters is None:
            corner(self.posterior_data[self.parameter_names],
                   labels=self.parameter_names,
                   weights=self.posterior_data['weights'], **kwargs)
        else:
            corner(self.posterior_data[parameters],
                   labels=parameters,
                   weights=self.posterior_data['weights'], **kwargs)

    def to_hdf(self, fname_or_buf, key='multinest'):
        """
        Writing the MultiNest result out to HDF5.

        Parameters
        ----------

        fname_or_buf: ~str
            filename or buffer

        key: ~str
            key to save it under default='multinest'
        """

        self.posterior_data.to_hdf(fname_or_buf, key=key)




class MultiNest(object):
    """
    Use multinest to fit a spectrum using a grid of models generated by specgrid.

    Parameters
    ----------


    likelihood: ~Likelihood object, optional
        By default uses the Likelihood object which uses the chi-square for the
        likelihood of observing the data given the model param_names

    run_dir:

    """


    def __init__(self, likelihood, priors, run_dir=None, prefix='specgrid_multinest'):

        self.run_dir = run_dir
        self.prefix = prefix
        self.likelihood = likelihood
        self.likelihood.multinest_evaluate = types.MethodType(
            multinest_evaluate, self.likelihood)

        self.likelihood.fixed_mask = types.MethodType(fixed_mask,
                                                      self.likelihood)
        if not hasattr(priors, 'prior_transform'):
            self.priors = PriorCollection(priors)
        else:
            self.priors = priors




    @property
    def n_params(self):
        return np.sum(~self.likelihood.fixed_mask())

    @property
    def basename_(self):
        return '{0}_'.format(self.basename)

    @property
    def posterior_data(self):
        if self._posterior_data is None:
            self._posterior_data = self.read_posterior_data()

        return self._posterior_data

    def prepare_fit_directory(self, run_dir, prefix):
        if not os.path.exists(run_dir):
            os.mkdir(run_dir)

        # checking if previous chains already exist
        return os.path.join(run_dir, prefix)

    def run(self, clean_up=None, **kwargs):

        if clean_up is None:
            if self.run_dir is None:
                clean_up = True
            else:
                clean_up = False

        if self.run_dir is None:
            run_dir = tempfile.mkdtemp()
        else:
            run_dir = self.run_dir

        basename = self.prepare_fit_directory(run_dir, self.prefix)

        start_time = time.time()

        logger.info('Starting fit in {0} with prefix {1}'.format(run_dir, self.prefix))
        pymultinest.run(self.likelihood.multinest_evaluate, self.priors.prior_transform,
                        self.n_params,
                        outputfiles_basename='{0}_'.format(basename),
                        **kwargs)

        logger.info("Fit finished - took {0:.2f} s"
                    .format(time.time() - start_time))
        fitted_parameter_names = [item for item in self.likelihood.param_names
                                  if not self.likelihood.fixed[item]]

        self.result = MultiNestResult.from_multinest_basename(
            basename, fitted_parameter_names)

        if clean_up:
            logger.info("Cleaning up - deleting {0}".format(run_dir))
            shutil.rmtree(run_dir)
        else:
            logger.info("Multinest files can be found in {0}".format(run_dir))

        self.likelihood.parameters[~self.likelihood.fixed_mask()] = (
            self.result.median.values)
        return self.result





    def __repr__(self):
        return "{0}\n\n{1}".format(
            self.likelihood, self.priors)
