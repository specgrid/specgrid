from specgrid.fix_spectrum1d import Spectrum1D
import numpy as np
from astropy import units as u



def test_flux_fix_spectrum1d():
    test_spec = Spectrum1D.from_array(np.arange(3000, 9000), np.random.random(6000),
                          dispersion_unit='angstrom',
                          unit='erg/s')

    assert test_spec.flux.unit == u.Unit('erg/s')

