import types
import logging
import time

import numpy as np
from jbopt.de import de
from jbopt.classic import classical

from starkit.fitkit.priors import PriorCollection

logger = logging.getLogger(__name__)


def fit_evaluate(self, model_param):
    # returns the likelihood of observing the data given the model param_names
    parameters = self.parameters.copy()
    parameters[~self.fixed_mask()] = model_param


    loglikelihood = self.evaluate(*parameters)
    return float(loglikelihood)

def fixed_mask(self):
    return np.array([getattr(self, param_name).fixed
                      for param_name in self.param_names])


class JBOptPriorCollection(PriorCollection):

    def prior_transform(self, cube):
        cube = np.asarray(cube)
        super(JBOptPriorCollection, self).prior_transform(cube,
                                                          None, len(cube))
        return cube

class JBOpt(object):

    def __init__(self, likelihood, priors, output_basename='test_all'):

        self.likelihood = likelihood
        self.likelihood.fit_evaluate = types.MethodType(
            fit_evaluate, self.likelihood)

        self.likelihood.fixed_mask = types.MethodType(fixed_mask,
                                                      self.likelihood)
        if not hasattr(priors, 'prior_transform'):
            self.priors = JBOptPriorCollection(priors)
        else:
            self.priors = priors

        self.fit_parameter_names = [
            item for i, item in enumerate(self.likelihood.param_names)
            if not self.likelihood.fixed_mask()[i]]

        self.args = dict(loglikelihood=self.likelihood.fit_evaluate,
                         transform=self.priors.prior_transform,
                         prior=lambda x: 0,
                         parameter_names=self.fit_parameter_names,
                         )


    def run(self, output_basename, method='de', start=None, nsteps=2000, verbose=0):
        if start is None:
            start = [0.5] * len(self.fit_parameter_names)

        self.args['start'] = start
        self.args['nsteps'] = nsteps
        self.args['disp'] = verbose
        self.args['output_basename'] = output_basename

        start_time = time.time()
        if method == 'de':
            self.result = self._run_de()
        elif method in ('cobyla', 'ralg', 'mma', 'auglag', 'minuit',
                        'neldermead'):
            self.result = self._run_classical(method)


        logger.info('Fit took {0:.2f}s'.format(time.time() - start_time))
        self.result['best_values'] = self.priors.prior_transform(
            self.result['start'])
        self.likelihood.parameters[~self.likelihood.fixed_mask()] = (
            self.result['best_values'])

        return self.result

    def _run_de(self):

        return de(**self.args)

    def _run_classical(self, method):
        return classical(method=method, **self.args)
