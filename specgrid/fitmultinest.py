from scipy.stats import norm, poisson
import pymultinest

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
        return norm.ppf(cube,scale=self.sigma,loc=self.m)
        

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

class SpectralData(object):
    pass

class Likelihood(object):
    
    def __init__(model_star):
        pass
    
    def __call__(spectral_data, *pars):
        pass
        #return float

class FitMultiNest(object):
    
    def __init__(spectrum, model_star, priors, likelihood=None):
        
        self.priors = PriorCollection(priors)
        if likelihood is None:
            self.likelihood = Likelihood(model_star)
        self.data = spectrum
        
    def run(self, **kwargs):
        # run pymultinest on the data object
        
    

class PriorCollection(object):
    def __init__(prior_dict):
        #something
        pass
        
        
    def prior(cube, ndim, x):
        for current_prior, c in zip(self.priors, cube):
            c = current_prior(c)
            
        #change the cube
        