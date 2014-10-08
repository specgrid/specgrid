# this contains imports plugins that configure py.test for astropy tests.
# by importing them here in conftest.py they are discoverable by py.test
# no matter how it is invoked within the source tree.

from astropy.tests.pytest_plugins import *
import pytest
## Uncomment the following line to treat all DeprecationWarnings as
## exceptions
# enable_deprecations_as_exceptions()

import specgrid
from specgrid import SpectralGrid
import h5py

def data_path(filename):
    return os.path.join(specgrid.__path__[0], 'data', filename)


@pytest.fixture()
def test_specgrid():
    return SpectralGrid(data_path('munari_small.h5'))

@pytest.fixture
def h5_test_data():
    return h5py.File(data_path('test_data.h5'), mode='r')

