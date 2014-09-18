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


class TestLikelihood(object):
    def setup(self):
        self.spec_grid = BaseSpectralGrid(data_path('munari_small.h5'))
        self.model_star = ModelStar([self.spec_grid])
        self.priors=OrderedDict([('teff', fitmultinest.UniformPrior(4000, 9000)), 
            ('logg', fitmultinest.GaussianPrior(3, 0.5)), 
            ('feh', fitmultinest.FixedPrior(-0.5))])
        spectrum = self.model_star.evaluate(teff=5780,logg=4.14,feh=0.0)
        spectrum.uncertainty = np.ones(spectrum.flux.shape) * spectrum.flux.unit
        self.likelihood = fitmultinest.Likelihood(spectrum, self.model_star, self.priors.keys())
        
    def test_likelihood(self):
        likelihood = self.likelihood([5780,4.14,0.0], 4,4 )
        nptest.assert_almost_equal(likelihood, 0.0)
        
        l2 = self.likelihood([5781,4.14,0.0], 4, 4)
        nptest.assert_almost_equal(l2, -14008598406.946743)

    
class TestSimpleMultinest(object):
    def setup(self):
        self.spec_grid = BaseSpectralGrid(data_path('munari_small.h5'))
        self.model_star = ModelStar([self.spec_grid])
        spectrum = self.model_star.evaluate(teff=5780,logg=4.14,feh=0.0)
        spectrum.uncertainty = (np.ones(spectrum.flux.shape)+np.sqrt(spectrum.flux.value)) * spectrum.flux.unit
        self.priors=OrderedDict([('teff', fitmultinest.UniformPrior(5000, 6000)), 
            ('logg', fitmultinest.GaussianPrior(4.3,0.3)), 
            ('feh', fitmultinest.FixedPrior(0.05))])
        self.fit_multinest = fitmultinest.FitMultinest(spectrum, self.priors, self.model_star)
    
    def test_multinest(self):
         # test a run of multinest
        try:
            os.mkdir('chains')
        except IOError:
            pass
        except OSError:
            pass
        
        self.fit_multinest.run(seed=741761)
         
        nptest.assert_almost_equal(self.fit_multinest.mean,[5779.616533128825, 4.135490388533586, 0.04999999999999998])
         
        nptest.assert_almost_equal(self.fit_multinest.sigma1,
        [[5779.607572222019, 5779.625214566039], [4.135304645124507, 4.135682928510821], [0.05, 0.05]])
         
        nptest.assert_almost_equal(self.fit_multinest.sigma3, 
        [[5779.590363406517, 5779.6427247009115], [4.134914605525979, 4.13604688052353], [0.05, 0.05]])
         
