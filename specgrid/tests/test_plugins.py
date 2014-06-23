import numpy as np
from specutils import Spectrum1D
from astropy import units as u
from specgrid.plugins import Convolve
from scipy.integrate import simps

my_spec = Spectrum1D.from_array(np.linspace(6000, 9000, 3000), np.ones(3000), 
                                dispersion_unit='angstrom')
                                
                                
my_spec.flux[1000] = 0.0

def test_convolution():
    convolve_plugin = Convolve(R, central_wavelength)
    conv_spectrum = convolve_plugin(my_spec)
    assert conv_spectrum.flux[999] < 1.
    assert conv_spectrum.flux[1000] > 0

    integral_pre_convolved = simps(my_spec.flux, my_spec.wavelength.value)
    intergral_convolved = simps(conv_spectrum.flux, conv_spectrum.wavelength.value)


    np.testing.assert_allclose(integral_pre_convolved,integral_convolved,rtol=10,atol=0.0002)
    
