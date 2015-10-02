from astropy.modeling import Parameter

class StarKitParameter(Parameter):

    constraints = tuple(list(Parameter.constraints) + ['prior'])

    @property
    def prior(self):
        if self._model is None:
            return getattr(self, '_prior', None)
        else:
            prior = self._model._constraints.get('prior', {})
            return prior.get(self._name, None)

    @prior.setter
    def prior(self, value):
        if self._model is None:
            raise AttributeError("can't set attribute 'prior' on Parameter "
                                 "definition")
        else:
            if 'prior' not in self._model._constraints:
                self._model._constraints['prior'] = {}
            self._model._constraints['prior'][self._name] = value

        self._prior = value
