import scipy.ndimage as nd
import numpy as np

import astropy.units as u
import astropy.constants as const

from specutils import Spectrum1D


class RotationalBroadening(object):
    vrot = 1. * u.km / u.s
    resolution = (20 * u.km / u.s / const.c).to(1)
    limb_darkening = 0.6
    parameters = ['vrot']

    def rotational_profile(self):
        vrot_by_c = (self.vrot / const.c).to(1)
        half_width = np.round(vrot_by_c / self.resolution).astype(int)
        profile_velocity = np.linspace(-half_width, half_width,
                                       2 * half_width + 1) * self.resolution
        profile = np.maximum(0.,
                             1. - (profile_velocity / vrot_by_c) ** 2)
        profile = ((2 * (1-self.limb_darkening) * np.sqrt(profile) +
                    0.5 * np.pi * self.limb_darkening * profile) /
                   (np.pi * vrot_by_c * (1.-self.limb_darkening/3.)))
        return profile/profile.sum()

    def __call__(self, spectrum):
        wavelength, flux = spectrum.wavelength.value, spectrum.flux.value
        log_grid_log_wavelength = np.arange(np.log(wavelength.min()),
                                            np.log(wavelength.max()),
                                            self.resolution.to(1).value)
        log_grid_wavelength = np.exp(log_grid_log_wavelength)
        log_grid_flux = np.interp(log_grid_wavelength, wavelength, flux)
        profile = self.rotational_profile()
        log_grid_convolved = nd.convolve1d(log_grid_flux, profile)
        convolved_flux = np.interp(wavelength, log_grid_wavelength,
                                   log_grid_convolved)
        return Spectrum1D.from_array(spectrum.wavelength,
                                     convolved_flux * spectrum.flux.unit,
                                     dispersion_unit=spectrum.wavelength.unit)


class DopplerShift(object):

    vrad = 0. * u.km / u.s
    parameters = ['vrad']

    def __call__(self, spectrum):
        doppler_factor = 1. + self.vrad / const.c
        return Spectrum1D.from_array(spectrum.wavelength * doppler_factor,
                                     spectrum.flux,
                                     dispersion_unit=spectrum.wavelength.unit)
