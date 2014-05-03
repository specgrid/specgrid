from scipy import optimize
import numpy as np
from collections import OrderedDict

def fit_spectrum(spectrum, guess, model_star, fitter='leastsq'):
    def spectral_model_fit(pars):
        pardict = OrderedDict()
        for key, par in zip(guess.keys(), pars):
            pardict[key] = par

        model = model_star.eval(**pardict)

        if hasattr(spectrum, 'uncertainty') and spectrum.uncertainty is not None:
            uncertainty = spectrum.uncertainty
        else:
            uncertainty = np.ones_like(spectrum.flux)

        if np.isnan(model.flux[0]):
            return np.inf
        else:
            if fitter == 'leastsq':
                return ((spectrum.flux - model.flux) / uncertainty)
            else:
                return np.sum(((spectrum.flux - model.flux) / uncertainty)**2)


    if fitter == 'leastsq':
        fit = optimize.leastsq(spectral_model_fit, np.array(guess.values()),
                            full_output=True)

        stellar_params = OrderedDict((key, par) for key, par in zip(guess.keys(), fit[0]))
        if fit[1] is not None:
            stellar_params_uncertainty = OrderedDict(
                (key, np.sqrt(par)) for key, par in
                 zip(guess.keys(), np.diag(fit[1])))
        else:
            stellar_params_uncertainty = OrderedDict((key, None) for key in guess.keys())
    else:
        fit =  optimize.minimize(spectral_model_fit, np.array(guess.values()), method=fitter)
        stellar_params = OrderedDict((key, par) for key, par in zip(guess.keys(), fit['x']))
        stellar_params_uncertainty = OrderedDict((key, None) for key, par in zip(guess.keys(), fit['x']))

    return stellar_params, stellar_params_uncertainty, fit
