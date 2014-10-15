import os
import specgrid
from specgrid import SpectralGrid
import numpy.testing as nptesting
import numpy as np
import h5py
import pytest


def data_path(filename):
    return os.path.join(specgrid.__path__[0], 'data', filename)

@pytest.fixture
def h5_test_data():
    return h5py.File(data_path('test_data.h5'), mode='r')



def test_default_call(test_specgrid):
    default_interp_spec = test_specgrid()
    nptesting.assert_allclose(default_interp_spec.wavelength,
                              test_specgrid.wavelength)
    nptesting.assert_allclose(default_interp_spec.flux.value,
                              test_specgrid.fluxes[0])

def test_evalulate_method(test_specgrid):
    spec = test_specgrid.evaluate(teff=5780., logg=4.4)
    assert hasattr(spec, 'flux')
    assert hasattr(spec, 'wavelength')
    nptesting.assert_allclose(test_specgrid.teff, 5780.)
    nptesting.assert_allclose(test_specgrid.logg, 4.4)

def test_simple_interpolation(test_specgrid, h5_test_data):
    spec = test_specgrid.evaluate(teff=5780., logg=4.4, feh=0.0)
    comp_flux = h5_test_data['test_data']['simple_interpolation_sun'].__array__()
    nptesting.assert_allclose(spec.flux.value, comp_flux)

