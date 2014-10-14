# this contains imports plugins that configure py.test for astropy tests.
# by importing them here in conftest.py they are discoverable by py.test
# no matter how it is invoked within the source tree.

from astropy.tests.pytest_plugins import *
import pytest
## Uncomment the following line to treat all DeprecationWarnings as
## exceptions
# enable_deprecations_as_exceptions()

import h5py
import numpy as np
from astropy import units as u

import specgrid
from specgrid import SpectralGrid
from specgrid import Spectrum1D, ModelStar

def data_path(filename):
    return os.path.join(specgrid.__path__[0], 'data', filename)


@pytest.fixture()
def test_specgrid():
    return SpectralGrid(data_path('munari_small.h5'))

@pytest.fixture(scope='session')
def h5_test_data():
    return h5py.File(data_path('test_data.h5'), mode='r')

@pytest.fixture(scope='session')
def test_spectrum():
    wave, flux = np.loadtxt(data_path('test_spec.txt'), unpack=True)
    return Spectrum1D.from_array(wave * u.angstrom,
                                 flux * u.erg / u.s / u.cm**2 / u.angstrom)

#### MULTINEST Fixtures ###

@pytest.fixture()
def test_multinest_model_star(test_specgrid):
    return ModelStar(test_specgrid)

@pytest.fixture()
def test_multinest_spectrum(test_multinest_model_star):
    spectrum = test_multinest_model_star.evaluate(teff=5780., logg=4.14, feh=0.0)
    spectrum.uncertainty = (np.ones(spectrum.flux.shape) +
                            np.sqrt(spectrum.flux.value))

    return spectrum
