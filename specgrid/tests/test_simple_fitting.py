import pytest
from specgrid.fitting import fit_spectrum
from specgrid import assemble_observation
from specgrid import Spectrum1D
from astropy import units as u
import numpy as np

import numpy.testing as nptesting

@pytest.fixture
def test_model_star(test_specgrid, test_spectrum):
    model_star = assemble_observation(
        test_specgrid, plugin_names=['doppler', 'rotation', 'resolution'],
        spectrum = test_spectrum,
        normalize_npol=4)
    return model_star


def test_simple_fit_crash1(test_spectrum, test_model_star):
    with pytest.raises(ValueError):
        fit_spec = fit_spectrum(test_spectrum, test_model_star)


def test_simple_fit_crash2(test_spectrum, test_model_star):
    with pytest.raises(ValueError):
        fit_spec = fit_spectrum(test_spectrum, test_model_star, wrong_param=4)

def test_simple_fit1(test_spectrum, test_model_star):
    fit_result = fit_spectrum(test_spectrum, test_model_star, teff=5200, logg=4.1,
                            feh=0.1)

    nptesting.assert_allclose(fit_result.best_fit_values[0], 5779.88207292)
    nptesting.assert_allclose(fit_result.best_fit_values[1], 4.40831046502)
    nptesting.assert_allclose(fit_result.best_fit_values[2], 0.328212333713)

@pytest.mark.xfail
def test_simple_fit1(test_spectrum, test_model_star):
    fit_result = fit_spectrum(test_spectrum, test_model_star, teff=5600, logg=4.3,
                            feh=0.2, fitter='Nelder-Mead')
    assert False, "finish nelder mead test"




