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
        return poisson.ppf(cube,loc=self.m)

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





class PriorCollections(object):


    def __init__(self, prior_dict, parameter_order=[]):

        self.priors = prior_dict

        for prior in self.priors.values():
            if not hasattr(prior, '__call__'):
                raise TypeError('Given prior {0} is not callable'.format(prior))

    def prior_transform(self, cube, ndim, nparam):
        # will be given an array of values from 0 to 1 and transforms it
        # according to the prior distribution

        for i in xrange(nparam):
            cube[i] = self.priors.values()[i](cube[i])

