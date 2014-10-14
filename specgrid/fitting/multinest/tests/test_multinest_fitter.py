import pytest
import numpy as np
from collections import OrderedDict
from specgrid.fitting.multinest import priors, fit_simple_spectrum_multinest
import numpy as np

import numpy.testing as nptesting




try:
    import pymultinest
except ImportError:
    pymultinest_available = False
else:
    pymultinest_available = True
    from specgrid import fitmultinest

pytestmark = pytest.mark.skipif(not pymultinest_available,
                                reason='pymultinest not available')



def test_simple_multinest(test_multinest_model_star, test_multinest_spectrum):
    prior_dict = OrderedDict([('teff', priors.UniformPrior(5000, 6000)),
            ('logg', priors.GaussianPrior(4.3, 0.3)),
            ('feh', priors.FixedPrior(0.05))])
    multinest_fitter = fit_simple_spectrum_multinest(test_multinest_spectrum,
                                  test_multinest_model_star, prior_dict)

    multinest_fitter.run(clean_up=False, seed=741761)

    nptesting.assert_allclose(multinest_fitter.result.mean.values(),
                                  [5779.616533128825, 4.135490388533586,
                                   0.04999999999999998],
                                  rtol=1e-5)

#    nptesting.assert_almost_equal(self.fit_multinest.sigma1,
#        [[5779.607572222019, 5779.625214566039], [4.135304645124507, 4.135682928510821], [0.05, 0.05]])

#    nptesting.assert_almost_equal(self.fit_multinest.sigma3,
#        [[5779.590363406517, 5779.6427247009115], [4.134914605525979, 4.13604688052353], [0.05, 0.05]])



