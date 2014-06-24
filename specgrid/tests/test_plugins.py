import numpy as np
import numpy.testing as npt
from specutils import Spectrum1D
from astropy import units as u
from specgrid.plugins import Convolve, Interpolate
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
    plt.plot(my_spec.wavelength.value, my_spec.flux, color= 'blue', label = 'my_spec')
    plt.plot(conv_spectrum.wavelength.value, conv_spectrum.flux, color = 'red', label = 'conv_spectrum')
    plt.legend(loc = 'upper right')
    plt.show()
    '''


def test_interpolate():
    obs_spectrum = Spectrum1D.from_array(np.arange(1,26), np.zeros(25), dispersion_unit = 'micron', unit = 'erg / (cm^2 s Angstrom)')
    test_spectrum = Spectrum1D.from_array(np.arange(0.5,26.5,1), np.arange(0.5,26.5,1), dispersion_unit = 'micron', unit = 'erg / (cm^2 s Angstrom)')
    
    interpolate_plugin = Interpolate(obs_spectrum)
    interpolated_test_spectrum = interpolate_plugin(test_spectrum)

    expected_interpolated_flux = np.arange(1,26)
    interpolated_test_flux = interpolated_test_spectrum.flux
    
    #This allows you to manually check the mean value of the interpolated test flux. This should be an array of 1 to 25 (spacing of 1).
    print "interpolated_test_flux = {0}".format(interpolated_test_flux)
    #When I tested this, I got: interpolated_test_flux = [  1.   2.   3.   4.   5.   6.   7.   8.   9.  10.  11.  12.  13.  14.  15.  16.  17.  18.  19.  20.  21.  22.  23.  24.  25.] ... As expected.

    npt.assert_array_almost_equal(interpolated_test_flux, expected_interpolated_flux, decimal = 6)

    """
    Optional plotting, to get a feel for why this test works conceptually. To use this feature, remove the ''' below, and take notice that the green and red markers lie on the same x-values, while the green and blue markers DO NOT lie on the same x-values (given that the interpolate worked correctly).
    """
    '''
    plt.clf()
    plt.figure(1)
    plt.xlim(0, 26)
    plt.ylim(-1, 26)
    plt.plot(obs_spectrum.wavelength.value, obs_spectrum.flux, marker = 'o', linestyle = 'None', color = 'green', label = 'obs_spectrum')
    plt.plot(test_spectrum.wavelength.value, test_spectrum.flux, marker = 'o', linestyle = 'None', color = 'blue', label = 'test_spectrum')
    plt.plot(interpolated_test_spectrum.wavelength.value, interpolated_test_spectrum.flux, marker = 'o', linestyle = 'None', color = 'red', label = 'interpolated_test_spectrum')
    plt.legend(loc = 'upper left', numpoints = 1)
    plt.grid()
    plt.xticks([i for i in range(0,27)])
    plt.yticks([i for i in range(0,27)])
    plt.show()
    '''
