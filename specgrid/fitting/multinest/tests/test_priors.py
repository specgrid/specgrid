import numpy.testing as nptest
import numpy as np

from collections import OrderedDict

from specgrid.fitting.multinest import priors

def test_uniform_prior():
    uniform_prior = priors.UniformPrior(lbound=3000, ubound=9000)
    
    nptest.assert_almost_equal(uniform_prior(0.0), uniform_prior.lbound)
    nptest.assert_almost_equal(uniform_prior(0.0), uniform_prior.lbound)
    

def test_gaussian_prior():
    gaussian_prior = priors.GaussianPrior(m=3, sigma=2)
    
    nptest.assert_almost_equal(gaussian_prior(0.0), -np.inf)
    nptest.assert_almost_equal(gaussian_prior(1.0), np.inf)
    nptest.assert_almost_equal(gaussian_prior(0.5), 3)
    nptest.assert_almost_equal(gaussian_prior(0.5 + 0.682689492 / 2), 5)
    
def test_fixed_prior():
    fixed_prior = priors.FixedPrior(2.0)
    
    nptest.assert_almost_equal(fixed_prior(.5), 2.0)
    
def test_prior_collections():
    prior_dict=OrderedDict([('teff', priors.UniformPrior(4000, 9000)),
            ('logg', priors.GaussianPrior(3, 0.5)),
            ('feh', priors.FixedPrior(-0.5))])
    
    col = priors.PriorCollection(prior_dict.values())
    c = np.array([0.5, 0.5, 1.0])
    col.prior_transform(c, 3, 3)
    
    nptest.assert_almost_equal(c, [(9000.0 + 4000.0) * 0.5, 3, -0.5])


    
