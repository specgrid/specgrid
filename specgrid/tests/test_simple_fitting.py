import pytest
from specgrid.fitting import fit_spectrum
from specgrid import assemble_model_star
from specgrid import Spectrum1D
from astropy import units as u
import numpy as np
import numpy.testing as nptesting

@pytest.fixture
def test_model_star(test_specgrid, test_spectrum):
    model_star = assemble_model_star(
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
    fit_spec = fit_spectrum(test_spectrum, test_model_star, teff=5000, logg=4.1,
                            feh=0.3)

    nptesting.assert_allclose(fit_spec[0]['teff'], 5098.85879801)
    nptesting.assert_allclose(fit_spec[0]['logg'], 3.77366711329)
    nptesting.assert_allclose(fit_spec[0]['feh'], 2.47711074314e-07)





