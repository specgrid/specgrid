from scipy.stats import norm, poisson

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
    
#def multinest_spectrum(grid, spectrum, 
#priors={'teff':UniformPrior(4000, 9000), 'logg':GaussianPrior(3, 0.5), 'feh':FixedPrior(-0.5))

class MultiNest(object):
    def __init__(spectrum, model_star, priors):
        pass

class PriorCollections(object):
    def __init__(prior_dict):
        #something
        pass
        
        
    def prior(cube, ndim, x):
        for current_prior, c in zip(self.priors, cube):
            c = current_prior(c)
            
        #change the cube
        