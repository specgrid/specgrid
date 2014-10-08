from scipy import optimize
import numpy as np
from collections import OrderedDict


def fit_spectrum(spectrum, model_star, fitter='leastsq',
                 fill_value=1e99, valid_slice=slice(None), **guesses):
    """

    :param spectrum:
    :param model_star:
    :param fitter:
    :param fill_value:
    :param valid_slice:
    :return:
    """

    if len(guesses) == 0:
        raise ValueError('fit_spectrum requires initial guessses '
                         'for the parameters that should be fit. '
                         'e.g. fit_spectrum(spec, mstar, teff=5000, logg=4)')

    if not set(guesses.keys()).issubset(set(model_star.parameters)):
        raise ValueError('Some parameter guesses are not parameters'
                         ' in the current model star ({0})'.format(
            ', '.join(set(guesses.keys()) - set(model_star.parameters))))


    parameter_guesses = OrderedDict((key, guesses[key])
                                  for key in model_star.parameters
                                  if key in guesses)


    def spectral_model_fit(parameter_values):

        parameter_dict = {key:value for key, value in zip(
            parameter_guesses.keys(), parameter_values)}
        model = model_star.evaluate(**parameter_dict)

        if getattr(spectrum, 'uncertainty', None) is not None:
            uncertainty = spectrum.uncertainty
        else:
            uncertainty = np.ones_like(spectrum.flux)

        if np.isnan(model.flux[0]):
            model_flux = np.ones_like(model.flux) * fill_value
        else:
            model_flux = model.flux

        quality = ((spectrum.flux - model_flux) / uncertainty)[valid_slice]
        return quality if fitter == 'leastsq' else quality.sum()

    if fitter == 'leastsq':

        fit = optimize.leastsq(spectral_model_fit,
                               np.array(parameter_guesses.values()),
                               full_output=True)

        stellar_params = OrderedDict((key, par) for key, par in zip(
            parameter_guesses.keys(), fit[0]))

        if fit[1] is not None:
            stellar_params_uncertainty = OrderedDict(
                (key, np.sqrt(par)) for key, par in
                zip(parameter_guesses.keys(), np.diag(fit[1])))
        else:
            stellar_params_uncertainty = OrderedDict((key, None)
                                                     for key in parameter_guesses.keys())
    else:
        fit = optimize.minimize(spectral_model_fit,
                                np.array(parameter_guesses.values()),
                                method=fitter)

        stellar_params = OrderedDict(
            (key, par) for key, par in zip(parameter_guesses.keys(), fit['x']))
        stellar_params_uncertainty = OrderedDict(
            (key, None) for key, par in zip(parameter_guesses.keys(), fit['x']))

    return stellar_params, stellar_params_uncertainty, fit
