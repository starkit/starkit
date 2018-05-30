import numpy as np
from scipy import stats

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
        self.priors = priors
        self.priors.update(self._generate_grid_uniform_priors())

        self.full_likelihood = model | likelihood
        self.param_names_raw = self.full_likelihood.param_names
        self.param_names = [item[:item.rfind('_')]
                            for item in self.param_names_raw]
        self.param_fixed = np.array([getattr(self.full_likelihood, param_name).fixed
                            for param_name in self.param_names_raw])
        self.param_values = np.array([getattr(self.full_likelihood, param_name).value
                            for param_name in self.param_names_raw])


    def _generate_grid_uniform_priors(self):
        spectral_grid = self.model[0]
        param_names = spectral_grid.param_names
        bounds = spectral_grid.get_grid_extent()
        priors = {pname:stats.uniform(loc=bound[0], scale=bound[1]-bound[0]).ppf
                  for pname, bound in zip(param_names, bounds)}
        return priors

    def update_prior_transform(self):
        pass

    @property
    def model_param_values(self):
        return self.param_values[:len(self.model.param_names)]

    @property
    def model_param_fixed(self):
        return self.param_fixed[:len(self.model.param_names)]


    def model_eval(self, *params):
        model_parameters = self.model_param_values.copy()
        model_parameters[~self.model_param_fixed] = params
        return self.model.evaluate(*model_parameters)

    def full_likelihood_eval(self, *params):
        parameters = self.param_values.copy()
        parameters[~self.param_fixed] = params
        return self.full_likelihood.evaluate(*parameters)


