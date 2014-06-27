from scipy import optimize
import numpy as np
from collections import OrderedDict


def fit_spectrum(spectrum, guess, model_star, fitter='leastsq',
                 fill_value=1e99, valid_slice=slice(None)):

    def spectral_model_fit(pars):
        pardict = OrderedDict()
        for key, par in zip(guess.keys(), pars):
            pardict[key] = par

        model = model_star.eval(**pardict)

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
        fit = optimize.leastsq(spectral_model_fit, np.array(guess.values()),
                               full_output=True)

        stellar_params = OrderedDict((key, par)
                                     for key, par in zip(guess.keys(), fit[0]))
        if fit[1] is not None:
            stellar_params_uncertainty = OrderedDict(
                (key, np.sqrt(par)) for key, par in
                zip(guess.keys(), np.diag(fit[1])))
        else:
            stellar_params_uncertainty = OrderedDict((key, None)
                                                     for key in guess.keys())
    else:
        fit = optimize.minimize(spectral_model_fit, np.array(guess.values()),
                                method=fitter)
        stellar_params = OrderedDict(
            (key, par) for key, par in zip(guess.keys(), fit['x']))
        stellar_params_uncertainty = OrderedDict(
            (key, None) for key, par in zip(guess.keys(), fit['x']))

    return stellar_params, stellar_params_uncertainty, fit
