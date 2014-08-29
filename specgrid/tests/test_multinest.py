from specgrid import fitmultinest
import numpy.testing as nptest
import numpy as np


def test_uniform_prior():
    uniform_prior = fitmultinest.UniformPrior(lbound=3000, ubound=9000)
    
    nptest.assert_almost_equal(uniform_prior(0.0), uniform_prior.lbound)
    nptest.assert_almost_equal(uniform_prior(0.0), uniform_prior.lbound)
    

def test_gaussian_prior():
    gaussian_prior = fitmultinest.GaussianPrior(m=3, sigma=2)
    
    nptest.assert_almost_equal(gaussian_prior(0.0), -np.inf)
    nptest.assert_almost_equal(gaussian_prior(1.0), np.inf)
    nptest.assert_almost_equal(gaussian_prior(0.5), 3)
    nptest.assert_almost_equal(gaussian_prior(0.5 + 0.682689492 / 2), 5)