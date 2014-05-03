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
        vrot_by_c = (np.maximum(0.1 * u.m / u.s, np.abs(self.vrot)) /
                     const.c).to(1)
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
        wavelength, flux = spectrum.wavelength.value, spectrum.flux
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
                                     convolved_flux,
                                     dispersion_unit=spectrum.wavelength.unit,
                                     unit=spectrum.unit)


class DopplerShift(object):

    vrad = 0. * u.km / u.s
    parameters = ['vrad']

    def __call__(self, spectrum):
        doppler_factor = 1. + self.vrad / const.c
        return Spectrum1D.from_array(spectrum.wavelength * doppler_factor,
                                     spectrum.flux,
                                     dispersion_unit=spectrum.wavelength.unit)


class Interpolate(object):

    parameters = []

    def __init__(self, observed):
        self.observed = observed

    def __call__(self, spectrum):
        wavelength, flux = spectrum.wavelength.value, spectrum.flux.value

        interpolated_flux = np.interp(self.observed.wavelength.value,
                                      wavelength, flux)
        return Spectrum1D.from_array(
            self.observed.wavelength,
            interpolated_flux * self.observed.flux,
            dispersion_unit=self.observed.wavelength.unit,
            unit=self.spectrum.unit)

class CCM89Extinction(object):
    parameters = ['a_v', 'r_v']

    def __init__(self, a_v=0.0, r_v=3.1):
        self.a_v = a_v
        self.r_v = r_v

    def __call__(self, spectrum):

        from specutils import extinction

        extinction_factor = 10**(-0.4*extinction.extinction_ccm89(
            spectrum.wavelength, a_v=self.a_v, r_v=self.r_v))


        return Spectrum1D.from_array(
            spectrum.wavelength.value,
            extinction_factor * spectrum.flux,
            dispersion_unit=spectrum.wavelength.unit, unit=spectrum.unit)



def observe(model, wgrid, slit, seeing, overresolve, offset=0.):
    """Convolve a model with a seeing profile, truncated by a slit, & pixelate

    Parameters
    ----------
    model: Table (or dict-like)
       Holding wavelengths and fluxes in columns 'w', 'flux'
    wgrid: array
       Wavelength grid to interpolate model on
    slit: float
       Size of the slit in wavelength units
    seeing: float
       FWHM of the seeing disk in wavelength units
    overresolve: int
       Factor by which detector pixels are overresolved in the wavelength grid
    offset: float, optional
       Offset of the star in the slit in wavelength units (default 0.)

    Returns
    -------
    Convolved model: Table
       Holding wavelength grid and interpolated, convolved fluxes
       in columns 'w', 'flux'
    """
    # make filter
    wgridres = np.min(np.abs(np.diff(wgrid)))
    filthalfsize = np.round(slit/2./wgridres)
    filtgrid = np.arange(-filthalfsize,filthalfsize+1)*wgridres
    # sigma ~ seeing-fwhm/sqrt(8*ln(2.))
    filtsig = seeing/np.sqrt(8.*np.log(2.))
    filt = np.exp(-0.5*((filtgrid-offset)/filtsig)**2)
    filt /= filt.sum()
    # convolve with pixel width
    filtextra = int((overresolve-1)/2+0.5)
    filt = np.hstack((np.zeros(filtextra), filt, np.zeros(filtextra)))
    filt = nd.convolve1d(filt, np.ones(overresolve)/overresolve)
    mint = np.interp(wgrid, model['w'], model['flux'])
    mconv = nd.convolve1d(mint, filt)
    return Table([wgrid, mconv], names=('w','flux'), meta={'filt': filt})
