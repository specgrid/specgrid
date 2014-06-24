import numpy as np
import numpy.testing as npt
from specutils import Spectrum1D, rvmeasure
from astropy import units as u
from specgrid.plugins import Convolve, Normalize
from scipy.integrate import simps

"""
The two lines of code below are for the optional plotting feature (see further for more details). Remove the ''' above and below the two lines of code to turn on this feature of the test
"""

'''
import matplotlib.pyplot as plt
'''

def test_convolution():
    my_spec = Spectrum1D.from_array(np.linspace(6000, 8000, 2000), np.ones(2000), dispersion_unit = 'Angstrom', unit = ' erg / (cm^2 s Angstrom)')
    my_spec.flux[1000] = 0.0

    R = 5000.
    central_wavelength = 7000. * u.Angstrom
    
    convolve_plugin = Convolve(R, central_wavelength)
    conv_spectrum = convolve_plugin(my_spec)
    assert conv_spectrum.flux[999] < 1.
    assert conv_spectrum.flux[1000] > 0

    integral_pre_convolved = simps(my_spec.flux, my_spec.wavelength.value)
    integral_convolved = simps(conv_spectrum.flux, conv_spectrum.wavelength.value)
    #This allows you to manually check the values of the integrals. They should be very close in value if the gaussian convolution was sucessful. 
    print "integral_pre_convolved = {0}".format(integral_pre_convolved)
    print "integral_convolved = {0}".format(integral_convolved)
    #When I tested this, I got the values: integral_pre_convolved = 1998.99949975 and integral_convolved = 1998.99949975.

    npt.assert_allclose(integral_pre_convolved, integral_convolved, rtol = 1.0, atol=0.0002)


    """
    Optional plotting, to get a feel for why this test works conceptually. To use this feature, remove the ''' above and below the code, and zoom into the plot at x = 7000 to see the desired effect.
    """
    '''
    plt.clf()
    plt.figure(1)
    plt.xlim(5500, 8500)
    plt.ylim(-0.25, 1.25)
    plt.plot(my_spec.wavelength.value, my_spec.flux, color='blue', label = 'my_spec')
    plt.plot(conv_spectrum.wavelength.value, conv_spectrum.flux, color = 'red', label = 'conv_spectrum')
    plt.legend(loc = 'upper right')
    plt.show()
    '''


def test_normalization():
    test_spec = Spectrum1D.from_array(np.arange(1,1001), np.arange(1,1001), dispersion_unit = 'Angstrom', unit = ' erg / (cm^2 s Angstrom)')
    ON = orderofnormalization = 1
    normalization_plugin = Normalize(ON)
    norm_spec = normalization_plugin(test_spec)
    norm_flux = norm_spec.flux
    norm_flux_mean = np.mean(norm_flux)
    
    #This allows you to manually check the mean value of the notmalized flux. This should be very close to 1 if the normalization was sucessful.
    print "norm_flux_mean = {0}".format(norm_flux_mean)
    #When I tested this, I got: norm_flux_mean =  1.0
    
    npt.assert_almost_equal(norm_flux_mean, 1, decimal = 3) 

    """
    Optional plotting, to get a feel for why this test works conceptually. To use this feature, remove the ''' below, and zoom near y = 1 to see the normalization.
    """
    '''
    plt.clf()
    plt.figure(1)
    plt.xlim(-100, 1100)
    plt.ylim(-100, 1100)
    plt.plot(test_spec.wavelength.value, test_spec.flux, color = 'blue', label = 'test_spec')
    plt.plot(norm_spec.wavelength.value, norm_spec.flux, color = 'red', label = 'norm_spec')
    plt.legend(loc = 'upper left')
    plt.show()
    '''
