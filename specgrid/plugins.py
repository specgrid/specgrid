import warnings
from collections import OrderedDict


import scipy.ndimage as nd
from scipy.ndimage.filters import gaussian_filter1d
import numpy as np
from numpy.polynomial import Polynomial

import astropy.units as u
import astropy.constants as const
from fix_spectrum1d import Spectrum1D

class RotationalBroadening(object):

    @property
    def vrot(self):
        return getattr(self, '_vrot', 1. * u.km / u.s)

    @vrot.setter
    def vrot(self, value):
        self._vrot = u.Quantity(value, u.km / u.s)


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

    @property
    def vrad(self):
        return getattr(self, '_vrad', 0. * u.km / u.s)

    @vrad.setter
    def vrad(self, value):
        self._vrad = u.Quantity(value, u.km / u.s)

    parameters = ['vrad']

    def __call__(self, spectrum):
        doppler_factor = 1. + self.vrad / const.c
        return Spectrum1D.from_array(spectrum.wavelength * doppler_factor,
                                     spectrum.flux,
                                     dispersion_unit=spectrum.wavelength.unit)


class InstrumentConvolve(object):
    """
    Convolve with a gaussian with given resolution to mimick an instrument

    Parameters
    ----------

    R : float or astropy.units.Quantity (unitless)
        resolution of the spectrum R = lambda/delta_lambda

    sampling: float
        number of pixels per resolution element (default=2.)

    """

    @property
    def R(self):
        return getattr(self, '_R', 0. * u.km / u.s)

    @R.setter
    def R(self, value):
        self._R = u.Quantity(value, u.Unit(1))


    parameters = ['R']

    def __init__(self, R, sampling=2.):
        self.R = u.Quantity(R, u.Unit(1))

        self.sampling = float(sampling)

    def __call__(self, spectrum):
        wavelength, flux = spectrum.wavelength.value, spectrum.flux
        log_grid_log_wavelength = np.arange(np.log(wavelength.min()),
                                            np.log(wavelength.max()),
                                            1 / (self.sampling *
                                                 self.R.to(1).value))
        log_grid_wavelength = np.exp(log_grid_log_wavelength)
        log_grid_flux = np.interp(log_grid_wavelength, wavelength, flux)
        sigma = self.sampling / (2 * np.sqrt(2 * np.log(2)))
        log_grid_convolved = nd.gaussian_filter1d(log_grid_flux, sigma)
        convolved_flux = np.interp(wavelength, log_grid_wavelength,
                                   log_grid_convolved)

        return Spectrum1D.from_array(spectrum.wavelength,
                                     convolved_flux,
                                     dispersion_unit=spectrum.wavelength.unit,
                                     unit=spectrum.unit)



class Interpolate(object):

    """
    This class can be called to do a interpolation on a given spectrum.
    You must initialize it with the observed spectrum. The output will be a
    Spectrum1D object.

    Parameters
    ----------
    observed: Spectrum1D object
        This is the observed spectrum which you want to interpolate your
        (model) spectrum to.
    """

    parameters = []

    def __init__(self, observed):
        self.observed = observed

    def __call__(self, spectrum):
        wavelength, flux = spectrum.wavelength.value, spectrum.flux
        interpolated_flux = np.interp(self.observed.wavelength.value,
                                      wavelength, flux)
        return Spectrum1D.from_array(
            self.observed.wavelength,
            interpolated_flux,
            dispersion_unit=self.observed.wavelength.unit,
            unit=self.observed.unit)


class Normalize(object):
    """Normalize a model spectrum to an observed one using a polynomial

    Parameters
    ----------
    observed : Spectrum1D object
        The observed spectrum to which the model should be matched
    npol : int
        The degree of the polynomial
    """

    parameters = []

    def __init__(self, observed, npol):
        if getattr(observed, 'uncertainty', None) is None:
            self.uncertainty = 1. * observed.flux.unit
        else:
            self.uncertainty = getattr(observed.uncertainty, 'array',
                                       observed.uncertainty)
        self.signal_to_noise = observed.flux / self.uncertainty
        self.flux_unit = observed.unit
        self._rcond = (len(observed.flux) *
                       np.finfo(observed.flux.dtype).eps)
        self._Vp = np.polynomial.polynomial.polyvander(
            observed.wavelength/observed.wavelength.mean() - 1., npol)
        self.domain = u.Quantity([observed.wavelength.min(),
                                  observed.wavelength.max()])
        self.window = self.domain/observed.wavelength.mean() - 1.

    def __call__(self, model):
        # V[:,0]=mfi/e, Vp[:,1]=mfi/e*w, .., Vp[:,npol]=mfi/e*w**npol
        
        V = self._Vp * (model.flux / self.uncertainty)[:, np.newaxis]
        # normalizes different powers
        scl = np.sqrt((V*V).sum(0))
        if np.isfinite(scl[0].value):  # check for validity before evaluating
            sol, resids, rank, s = np.linalg.lstsq(V/scl, self.signal_to_noise,
                                                   self._rcond)
            sol = (sol.T / scl).T
            if rank != self._Vp.shape[-1] - 1:
                msg = "The fit may be poorly conditioned"
                warnings.warn(msg)

            fit = np.dot(V, sol) * self.uncertainty
            # keep coefficients in case the outside wants to look at it
            self.polynomial = Polynomial(sol, domain=self.domain.value,
                                         window=self.window.value)
            return Spectrum1D.from_array(
                model.wavelength.value,
                fit)
        else:
            return model


class NormalizeParts(object):
    """Normalize a model spectrum to an observed one in multiple parts

    Here, different parts could, e.g., be different echelle orders or
    different chips for GMOS spectra.

    Parameters
    ----------
    observed : Spectrum1D object
        The observed spectrum to which the model should be matched
    parts : list of slices, index arrays, or boolean arrays
        Different parts of the observed spectrum which should be normalized
        separately.
    npol : list of int
        Polynomial degrees for the different parts
    """
    parameters = []

    def __init__(self, observed, parts, npol):
        self.parts = parts
        self.normalizers = []
        try:
            if len(parts) != len(npol):
                raise ValueError("List of parts should match in length to "
                                 "list of degrees")
        except TypeError:  # npol is single number
            npol = [npol]*len(parts)

        for part, _npol in zip(parts, npol):
            self.normalizers.append(
                Normalize(self.spectrum_1d_getitem(observed, part), _npol))

    @staticmethod
    def spectrum_1d_getitem(observed, part):
        observed_part = Spectrum1D.from_array(
            observed.wavelength[part],
            observed.flux[part])
        if getattr(observed, 'uncertainty', None) is not None:
            observed_part.uncertainty = getattr(observed.uncertainty, 'array',
                                                observed.uncertainty)[part]
        return observed_part

    def __call__(self, model):
        fit = np.zeros_like(model.wavelength.value)
        for part, normalizer in zip(self.parts, self.normalizers):
            fit[part] = normalizer(self.spectrum_1d_getitem(model, part)).flux

        return Spectrum1D.from_array(
            model.wavelength.value,
            fit, unit=self.normalizers[0].flux_unit,
            dispersion_unit=model.wavelength.unit)


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


stellar_physics_plugins = OrderedDict([('rotation', RotationalBroadening),
                                       ('doppler', DopplerShift),
                                       ('ccm89', CCM89Extinction),])

instrument_physics_plugins = OrderedDict([
    ('spec_resolution', InstrumentConvolve)])