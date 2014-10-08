import os
import specgrid
from specgrid import BaseSpectralGrid
import numpy.testing as nptesting
import numpy as np
import pytest

def data_path(filename):
    return os.path.join(specgrid.__path__[0], 'data', filename)

class TestSimpleSpectralGrid():
    def setup(self):
        self.spec_grid = BaseSpectralGrid(data_path('munari_small.h5'))

    def test_default_call(self):
        spec_grid = BaseSpectralGrid(data_path('munari_small.h5'))
        default_interp_spec = spec_grid()
        nptesting.assert_allclose(default_interp_spec.wavelength,
                                  spec_grid.wavelength)
        nptesting.assert_allclose(default_interp_spec.flux.value,
                                  spec_grid.fluxes[0])

    def test_evalulate_method(self):
        spec = self.spec_grid.evaluate(teff=5780., logg=4.4)
        assert hasattr(spec, 'flux')
        assert hasattr(spec, 'wavelength')
        nptesting.assert_allclose(self.spec_grid.teff, 5780.)
        nptesting.assert_allclose(self.spec_grid.logg, 4.4)
        nptesting.assert_allclose(self.spec_grid.feh,
                                  self.spec_grid.index['feh'][0])
