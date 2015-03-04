from collections import OrderedDict

from astropy import units as u
import numpy as np

from scipy import optimize

from specgrid.fitting.base import (SimpleSpectrumFitnessFunction,
                                   SimpleLeastsqFitResult, SimpleFitResult)


class SpectroPhotometryFitnessFunction(SimpleSpectrumFitnessFunction):
    def __init__(self, spectrum, model_observation, magnitude_set, fill_value=1e99):
        self.spectrum = spectrum
        self.model_observation = model_observation
        self.magnitude_set = magnitude_set
        self.fill_value = fill_value
        self.spectrum_uncertainty = self._get_spectrum_uncertainty()


    def fitness_function(self, param_values, param_names,
                         return_square_sum=True):
        """

        :param param_values:
        :param param_names:
        :param return_square_sum:
        :return:
        """
        parameter_dict = OrderedDict(zip(param_names, param_values))
        self.model_observation._set_parameters(**parameter_dict)

        model_star_spectrum = self.model_observation.model_star()

        model_observation_spectrum = self.model_observation.model_instrument(
            model_star_spectrum)

        model_observation_flux = self._replace_nan_flux(
            model_observation_spectrum)

        spectrum_chi2 = np.empty(len(model_observation_flux) +
                                 len(self.magnitude_set.magnitudes) - 1)

        spectrum_chi2[:len(model_observation_flux)] = (
            (self.spectrum.flux - model_observation_flux)
            / self.spectrum_uncertainty)

        spectrum_chi2[len(model_observation_flux):] = (
            self._calculate_magnitude_chi2(model_star_spectrum))

        if return_square_sum:
            return (spectrum_chi2**2).sum()
        else:
            return spectrum_chi2


    def _calculate_magnitude_chi2(self, model_star_spectrum):
        synth_mag = u.Quantity(self.magnitude_set.calculate_vega_magnitudes(
            model_star_spectrum))

        synth_colors = synth_mag[0:-1] - synth_mag[1:]

        colors = (self.magnitude_set.magnitudes[0:-1] -
                  self.magnitude_set.magnitudes[1:])
        color_uncertainties = np.sqrt(self.magnitude_set.magnitudes[0:-1]**2 +
                                      self.magnitude_set.magnitudes[1:]**2)

        chi2 = ((colors - synth_colors) / color_uncertainties)

        return chi2


def fit_spectro_photometry(spectrum, magnitude_set, model_observation,
                           method='leastsq', fill_value=1e99, **guesses):
    """
    Parameters
    ----------

    spectrum: ~specutils.Spectrum1D
        spectrum to be fit
    magnitude_set: wsynphot.MagnitudeSet
        set of magnitudes and filters
    model_observation: ~specgrid.ModelObservation
        model of observation, that returns a spectrum
    method: str
        method name one of the scipy minimize options or leastsq
    fill_value: ~float
    valid_slice: ~slice
        slice for valid data
    :return:
    """

    if len(guesses) == 0:
        raise ValueError('fit_spectrum requires initial guessses '
                         'for the param_names that should be fit. '
                         'e.g. fit_spectrum(spec, mstar, teff=5000, logg=4)')

    if not set(guesses.keys()).issubset(set(model_observation.param_names)):
        raise ValueError('Some parameter guesses are not param_names'
                         ' in the current model star ({0})'.format(
            ', '.join(set(guesses.keys()) - set(model_observation.param_names))))


    parameter_guesses = OrderedDict((key, guesses[key])
                                  for key in model_observation.param_names
                                  if key in guesses)

    if np.isnan(model_observation.evaluate(**parameter_guesses).flux.value[0]):
        raise ValueError('Initial guess ({0}) is outside the confines '
                         'of the grid -- aborting'.format(parameter_guesses))

    fitness_func = SpectroPhotometryFitnessFunction(spectrum, model_observation,
                                                    magnitude_set, fill_value)

    return_square_sum = not (method == 'leastsq')
    if method == 'leastsq':

        fit = optimize.leastsq(fitness_func.fitness_function,
                               np.array(parameter_guesses.values()),
                               args=(parameter_guesses.keys(),
                                     return_square_sum),
                               full_output=True)

        bestfit_spectrum = model_observation()

        result = SimpleLeastsqFitResult(fit[0], fit[1], parameter_guesses, bestfit_spectrum,
                                 spectrum, model_observation, fit)
        return result
    else:
        fit = optimize.minimize(fitness_func.fitness_function,
                               np.array(parameter_guesses.values()),
                               args=(parameter_guesses.keys(),
                                     return_square_sum))

        best_fit_spectrum = model_observation()

        result = SimpleFitResult(fit.x, parameter_guesses, best_fit_spectrum,
                                 spectrum, model_observation, fit, fit.success)

    return result

