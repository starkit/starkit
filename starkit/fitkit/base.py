import numpy as np
from scipy import stats
import logging

logger = logging.getLogger(__name__)

class ParameterExplorer(object):


    def __init__(self, model, likelihood, priors={}):
        """

        Parameters
        ----------
        model: astropy.modeling.Model
        likelihood: astropy.modeling.Model
        priors: dict
        """

        self.model = model
        self.likelihood = likelihood

        self.full_likelihood = model | likelihood


        self.param_names_raw = self.full_likelihood.param_names

        self.param_names = np.array([item[:item.rfind('_')]
                                     for item in self.param_names_raw])
        self.param_fixed = np.array(
            [getattr(self.full_likelihood, param_name).fixed
             for param_name in self.param_names_raw])

        self.param_values = np.array(
            [getattr(self.full_likelihood, param_name).value
             for param_name in self.param_names_raw])

        self.priors = self._generate_grid_uniform_priors()
        self.priors.update(priors)


        param_no_prior = set(self.param_names) - set(self.priors.keys())

        if len(param_no_prior) > 0:
            logger.warn('Prior not set for parameters {0}'.format(
                ', '.join(param_no_prior)))



    def _generate_grid_uniform_priors(self):
        """
        Get the uniform priors from the grid
        """
        priors = {}
        for pname, pname_raw in zip(self.param_names, self.param_names_raw):
            bound = self.full_likelihood.bounds[pname_raw]
            if bound[0] is None or bound[1] is None:
                continue
            priors[pname] = stats.uniform(loc=bound[0],
                                          scale=bound[1]-bound[0]).ppf

        return priors

    def update_prior_transform(self):
        self._prior_transform_list = [self.priors[pname]
                                     for pname in self.active_param_names]


    def prior_transform(self, parameters):
        transformed_params = np.array(
            [transform(param) for transform, param in
             zip(self._prior_transform_list, parameters)])
        return transformed_params


    @property
    def model_param_values(self):
        return self.param_values[:len(self.model.param_names)]

    @property
    def n_model_param(self):
        return len(self.model_param_values)

    @property
    def model_param_fixed(self):
        return self.param_fixed[:len(self.model.param_names)]

    @property
    def active_param_names(self):
        return [pname for pname, fixed in
                zip(self.param_names, self.param_fixed) if not fixed]

    def model_eval(self, *params):
        model_parameters = self.model_param_values.copy()
        model_parameters[~self.model_param_fixed] = params
        return self.model.evaluate(*model_parameters)

    def full_likelihood_eval(self, *params):
        parameters = self.param_values.copy()
        parameters[~self.param_fixed] = params
        return self.full_likelihood.evaluate(*parameters)


