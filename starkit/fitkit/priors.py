from scipy import stats


class UniformPrior(object):
    """
    A Uniform distribution prior

    Parameters
    ----------

    lbound: ~float
        lower bound

    ubound: ~float
        upper bound
    """

    def __init__(self, lbound, ubound):
        self.lbound = lbound
        self.ubound = ubound

    def __call__(self, cube):
        return cube * (self.ubound - self.lbound) + self.lbound

    def __repr__(self):
        return "<UniformPrior lbound {0} ubound {1}>".format(self.lbound,
                                                            self.ubound)

class GaussianPrior(object):
    """
    A gaussian prior

    Parameters
    ----------

    m: ~float
        mean of the distribution

    sigma: ~float
        sigma of the distribution

    """

    def __init__(self, m, sigma):
        self.m = m
        self.sigma = sigma

    def __call__(self, cube):
        return stats.norm.ppf(cube,scale=self.sigma,loc=self.m)

    def __repr__(self):
        return "<GaussianPrior mean {0} std {1}>".format(self.m, self.sigma)


class PoissonPrior(object):
    """
    A Poisson prior

    Parameters
    ----------

    m: ~float
        mean of the distribution


    """
    def __init__(self, m):
        self.m = m

    def __call__(self,cube):
        return stats.poisson.ppf(cube,loc=self.m)

    def __repr__(self):
        return "poisson prior: loc {0}".format(self.m)

class FixedPrior(object):
    """
    A fixed value

    Parameters
    ----------

    val: ~float
        fixed value
    """

    def __init__(self, val):
        self.val = val

    def __call__(self, cube):
        return self.val

    def __repr__(self):
        return "fixed prior: {0}".format(self.val)





class PriorCollection(object):
    """
    A collection of prior objects that will be evaluated
    """
    def __init__(self, priors_list):

        self.priors = priors_list

        for prior in self.priors:
            if not hasattr(prior, '__call__'):
                raise TypeError('Given prior {0} is not callable'.format(prior))

    def prior_transform(self, cube, ndim, nparam):
        # will be given an array of values from 0 to 1 and transforms it
        # according to the prior distribution

        for i in range(nparam):
            cube[i] = self.priors[i](cube[i])

    def _generate_prior_str(self):
        return [repr(item) for item in self.priors]

    def __repr__(self):
        return "Priors in collection:\n---\n{0}\n---".format(
            '\n'.join(self._generate_prior_str()))

