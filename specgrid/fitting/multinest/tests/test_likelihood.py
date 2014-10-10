from collections import OrderedDict
from specgrid import ModelStar
import numpy as np
from specgrid.fitting.multinest import priors, SimpleSpectrumLikelihood
import pytest

@pytest.fixture()
def test_multinest_model_star(test_specgrid):
    return ModelStar(test_specgrid)

@pytest.fixture()
def test_multinest_spectrum(test_multinest_model_star):
    spectrum = test_multinest_model_star.evaluate(teff=5780., logg=4.14, feh=0.0)
    spectrum.uncertainty = np.ones(spectrum.flux.shape)

    return spectrum


def test_simple_spec_likelihood1(test_multinest_model_star,
                                 test_multinest_spectrum):
    likelihood = SimpleSpectrumLikelihood(test_multinest_spectrum,
                                          test_multinest_model_star,
                                          parameter_names=['lkj', 'dasdas'])

"""
class TestSimpleSpectrumLikelihood(object):
    def setup(self, test_spec_grid):
        self.spec_grid = test_spec_grid
        self.model_star = ModelStar(self.spec_grid)
        self.priors = OrderedDict([('teff', priors.UniformPrior(4000, 9000)),
            ('logg', priors.GaussianPrior(3, 0.5)),
            ('feh', priors.FixedPrior(-0.5))])
        spectrum = self.model_star.evaluate(teff=5780,logg=4.14,feh=0.0)
        spectrum.uncertainty = np.ones(spectrum.flux.shape) * spectrum.flux.unit
        self.likelihood = SimpleSpectrumLikelihood(spectrum, self.model_star, self.priors.keys())

    def test_likelihood(self):
        likelihood = self.likelihood([5780,4.14,0.0], 4,4 )
        nptest.assert_almost_equal(likelihood, 0.0)

        l2 = self.likelihood([5781,4.14,0.0], 4, 4)
        nptest.assert_almost_equal(l2, -14008598406.946743)
"""