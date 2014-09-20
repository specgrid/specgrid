import numpy as np
import numpy.testing as npt
from numpy.polynomial import Polynomial
from specgrid.fix_spectrum1d import Spectrum1D
#from specutils import Spectrum1D
from astropy import units as u
from specgrid.plugins import Convolve, Interpolate, Normalize, NormalizeParts
from scipy.integrate import simps


def test_convolution():
    my_spec = Spectrum1D.from_array(np.linspace(6000, 8000, 2000),
                                    np.ones(2000), dispersion_unit='Angstrom',
                                    unit='erg/(cm^2 s Angstrom)')
    my_spec.flux[1000] = 0.0

    R = 5000.
    central_wavelength = 7000. * u.Angstrom

    convolve_plugin = Convolve(R, central_wavelength)
    conv_spectrum = convolve_plugin(my_spec)
    assert conv_spectrum.flux[999].value < 1.
    assert conv_spectrum.flux[1000].value > 0

    integral_pre_convolved = simps(my_spec.flux, my_spec.wavelength.value)
    integral_convolved = simps(conv_spectrum.flux,
                               conv_spectrum.wavelength.value)

    npt.assert_allclose(integral_pre_convolved, integral_convolved,
                        rtol=1.0, atol=0.0002)


def test_interpolate():
    obs_spectrum = Spectrum1D.from_array(np.arange(1,26), np.zeros(25),
                                         dispersion_unit='micron',
                                         unit='erg/(cm^2 s Angstrom)')
    test_spectrum = Spectrum1D.from_array(np.arange(0.5, 26.5, 1),
                                          np.arange(0.5, 26.5, 1),
                                          dispersion_unit='micron',
                                          unit='erg/(cm^2 s Angstrom)')

    interpolate_plugin = Interpolate(obs_spectrum)
    interpolated_test_spectrum = interpolate_plugin(test_spectrum)

    expected_interpolated_flux = np.arange(1,26)
    interpolated_test_flux = interpolated_test_spectrum.flux

    npt.assert_array_almost_equal(interpolated_test_flux.value,
                                  expected_interpolated_flux, decimal=6)


def test_normalize():
    w = np.arange(4000., 5000., 10.)
    coeff = [1., 0.1, 0.1]
    pol = Polynomial(coeff)
    obs_spectrum = Spectrum1D.from_array(w, pol(w), dispersion_unit=u.AA,unit='erg/(cm^2 s Angstrom)')
    norm = Normalize(obs_spectrum, npol=3)
    model = Spectrum1D.from_array(w, np.ones_like(w), dispersion_unit=u.AA,unit='erg/(cm^2 s Angstrom)')
    fit = norm(model)
    npt.assert_allclose(fit.flux.value, obs_spectrum.flux.value)
    npt.assert_allclose(norm.polynomial(obs_spectrum.wavelength.value),
                        obs_spectrum.flux.value)
    npt.assert_allclose(norm.polynomial.convert().coef,
                        np.array(coeff + [0.]), rtol=1e-5, atol=1.e-10)


def test_normalize_parts():
    w = np.arange(4000., 5000., 10.)
    coeff = [1., 0.1, 0.1]
    pol = Polynomial(coeff)
    obs_spectrum = Spectrum1D.from_array(w, pol(w), dispersion_unit=u.AA)
    parts = [slice(0, 30), slice(30, None)]
    norm_parts = NormalizeParts(obs_spectrum, parts, npol=[3, 3])
    model = Spectrum1D.from_array(w, np.ones_like(w), dispersion_unit=u.AA)
    fit = norm_parts(model)
    npt.assert_allclose(fit.flux.value, obs_spectrum.flux.value)
    for normalizer in norm_parts.normalizers:
        npt.assert_allclose(normalizer.polynomial.convert().coef,
                            np.array(coeff + [0.]), rtol=1e-3, atol=1.e-5)
    # try also with boolean arrays
    parts = [np.where(w < 4400.), np.where(w >= 4400.)]
    norm_parts = NormalizeParts(obs_spectrum, parts, npol=[3, 3])
    model = Spectrum1D.from_array(w, np.ones_like(w), dispersion_unit=u.AA)
    fit = norm_parts(model)
    npt.assert_allclose(fit.flux, obs_spectrum.flux)
    for normalizer in norm_parts.normalizers:
        npt.assert_allclose(normalizer.polynomial.convert().coef,
                            np.array(coeff + [0.]), rtol=1e-3, atol=1.e-5)
