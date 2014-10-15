from collections import OrderedDict
from specgrid import ModelStar
import numpy as np
from specgrid.fitting.multinest import priors, SimpleSpectrumLikelihood
import pytest

import numpy.testing as nptesting



def test_simple_spec_likelihood1(test_multinest_model_star,
                                 test_multinest_spectrum):
    with pytest.raises(ValueError):
        likelihood = SimpleSpectrumLikelihood(test_multinest_spectrum,
                                              test_multinest_model_star,
                                              parameter_names=['lkj', 'dasdas'])

def test_simple_spec_likelihood2(test_multinest_model_star,
                                 test_multinest_spectrum):

    likelihood = SimpleSpectrumLikelihood(test_multinest_spectrum,
                                              test_multinest_model_star,
                                              parameter_names=['teff', 'logg',
                                                               'feh'])

    likelihood_eval = likelihood([5780, 4.14, 0.0], 4, 4)
    nptesting.assert_almost_equal(likelihood_eval, 0.0)

    likelihood_eval2 = likelihood([5781, 4.14, 0.0], 4, 4)
    nptesting.assert_almost_equal(likelihood_eval2, -6464.8036889299492)

