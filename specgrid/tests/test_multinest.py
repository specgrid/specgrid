import numpy.testing as nptest
import numpy as np
import os
from collections import OrderedDict

from specgrid import fitmultinest
from specgrid.fix_spectrum1d import Spectrum1D
from specgrid.specgrid import BaseSpectralGrid
from specgrid.composite import ModelStar

import specgrid

def data_path(filename):
    return os.path.join(specgrid.__path__[0], 'data', filename)


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
    
def test_fixed_prior():
    fixed_prior = fitmultinest.FixedPrior(2.0)
    
    nptest.assert_almost_equal(fixed_prior(.5), 2.0)
    
def test_prior_collections():
    priors=OrderedDict([('teff', fitmultinest.UniformPrior(4000, 9000)), 
            ('logg', fitmultinest.GaussianPrior(3, 0.5)), 
            ('feh', fitmultinest.FixedPrior(-0.5))])
    
    col = fitmultinest.PriorCollections(priors)
    c = np.array([0.5, 0.5, 1.0])
    col.prior_transform(c, 3, 3)
    
    nptest.assert_almost_equal(c, [(9000.0+4000.0) * 0.5, 3, -0.5])

class TestSimpleMultinest(object):
    def setup(self):
        self.spec_grid = BaseSpectralGrid(data_path('munari_small.h5'))
        self.model_star = ModelStar([self.spec_grid])
        self.priors=OrderedDict([('teff', fitmultinest.UniformPrior(4000, 9000)), 
            ('logg', fitmultinest.GaussianPrior(3, 0.5)), 
            ('feh', fitmultinest.FixedPrior(-0.5))])
        self.fit_multinest = fitmultinest.FitMultinest(None, self.priors, self.model_star)
    

        
        
class TestLikelihood(object):
    def setup(self):
        self.spec_grid = BaseSpectralGrid(data_path('munari_small.h5'))
        self.model_star = ModelStar([self.spec_grid])
        self.priors=OrderedDict([('teff', fitmultinest.UniformPrior(4000, 9000)), 
            ('logg', fitmultinest.GaussianPrior(3, 0.5)), 
            ('feh', fitmultinest.FixedPrior(-0.5))])
        self.likelihood = fitmultinest.Likelihood(self.model_star, self.priors.keys())
        self.data = self.model_star.evaluate(teff=5780,logg=4.14,feh=0.0)
        self.data.uncertainty = np.ones(self.data.flux.shape) * self.data.flux.unit
        
    def test_likelihood(self):
        likelihood = self.likelihood(self.data,[5780,4.14,0.0])
        nptest.assert_almost_equal(likelihood, 0.0)
        
        l2 = self.likelihood(self.data,[5781,4.14,0.0])
        nptest.assert_almost_equal(l2, -14008598406.946743)

    
