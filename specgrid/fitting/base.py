from scipy import optimize
import numpy as np
from collections import OrderedDict\

from specgrid import Spectrum1D

class BaseFitResult():
    pass

class SimpleFitResult(BaseFitResult):
    """
    Simple fit result object

    :param fit_parameter_names:
    :param best_fit_values:
    :param best_fit_spectrum:
    :param spectrum:
    :param model:
    :param full_output:
    :return:
    """

    def __init__(self, best_fit_values, parameter_guesses, best_fit_spectrum,
                 spectrum, model, full_output, success):
        self.parameter_guesses = parameter_guesses
        self.best_fit_values = best_fit_values

        self.best_fit_spectrum = best_fit_spectrum
        self.spectrum = spectrum
        self.model = model
        self.full_output = full_output
        self.success = success and not np.all(np.isnan(best_fit_spectrum.flux))

    @property
    def residuals(self):
        residuals = getattr(self, '_residuals', None)
        if residuals is None:
            residuals = Spectrum1D.from_array(self.spectrum.wavelength,
                                              self.best_fit_spectrum.flux -
                                              self.spectrum.flux)
            self._residuals = residuals

        return residuals


    @property
    def chi2(self):
        """
        :math:`\\chi2`
        """
        chi2 = getattr(self, '_chi2', None)
        if chi2 is None:
            if self.spectrum.uncertainty is None:
                uncertainty = np.ones(self.residuals.flux.shape)
            else:
                uncertainty = self.spectrum.uncertainty
            chi2 = np.sum((self.residuals.flux.value / uncertainty) **2)
            self._chi2 = chi2
        return chi2

    @property
    def reduced_chi2(self):
        """
        Compute reduced :math:`\\chi^2_\\textrm{red}`
        """
        return self.chi2 / (len(self.residuals.flux)
                            - len(self.parameter_guesses))

    def _generate_fit_repr(self):
        repr_str = ""

        for i, parameter_name in enumerate(self.parameter_guesses.keys()):
            repr_str += '{0} {1:.3f} ({2:.3f})\n'.format(
                parameter_name, self.best_fit_values[i],
                self.parameter_guesses[parameter_name])
        return repr_str



    def __repr__(self):
        representation = 'Fit Result:\n'
        if self.success:
            representation += 'Fit successful\n'
        else:
            representation += 'Fit failed\n\n'

        representation += self._generate_fit_repr()

        return representation


class SimpleLeastsqFitResult(SimpleFitResult):
    def __init__(self, best_fit_values, covariance_matrix, parameter_guesses, best_fit_spectrum,
                 spectrum, model, full_leastsq):
        self.parameter_guesses = parameter_guesses
        self.best_fit_values = best_fit_values

        self.best_fit_spectrum = best_fit_spectrum
        self.spectrum = spectrum
        self.model = model
        self.covariance_matrix = covariance_matrix
        self.full_output = full_leastsq[2]
        self.full_output['mesg'] = full_leastsq[3]
        if full_leastsq[4] in [1, 2, 3, 4]:
            self.success = not np.all(np.isnan(best_fit_spectrum.flux))
        else:
            self.success = False

    @property
    def best_fit_uncertainties(self):
        if self.covariance_matrix is not None:
            return np.sqrt(np.diag(self.covariance_matrix))

    def _generate_fit_repr(self):
        repr_str = ""

        for i, parameter_name in enumerate(self.parameter_guesses.keys()):
            repr_str += '{0} {1:.3f} +/- {2:.3f} ({3:.3f})\n'.format(
                parameter_name, self.best_fit_values[i],
                self.best_fit_uncertainties[i],
                self.parameter_guesses[parameter_name])
        return repr_str

class SimpleSpectrumFitnessFunction(object):
    def __init__(self, spectrum, model_observation, fill_value = 1e99):
        self.spectrum = spectrum
        self.model_observation = model_observation
        self.fill_value = fill_value
        

    def _get_spectrum_uncertainty(self):
        if getattr(self.spectrum, 'uncertainty', None) is not None:
            return self.spectrum.uncertainty
        else:
            return np.ones_like(self.spectrum.flux)

    def _get_model_flux(self, param_values, param_names):
        parameter_dict = OrderedDict(zip(param_names, param_values))

        model_spectrum = self.model_observation.evaluate(**parameter_dict)

        return self._replace_nan_flux(model_spectrum)

    def _replace_nan_flux(self, model_spectrum):
        if np.isnan(model_spectrum.flux[0]):
            model_flux = np.ones_like(model_spectrum.flux) * self.fill_value
        else:
            model_flux = model_spectrum.flux

        return model_flux

    def fitness_function(self, param_values, param_names, return_square_sum):

        model_flux = self._get_model_flux(param_values, param_names)

        uncertainty = self._get_spectrum_uncertainty()



        quality = ((self.spectrum.flux - model_flux) / uncertainty)


        return quality if return_square_sum else (quality**2).sum()





def fit_spectrum(spectrum, model_observation, method='leastsq',
                 fill_value=1e99, valid_slice=slice(None), **guesses):
    """
    Parameters
    ----------

    spectrum: ~specutils.Spectrum1D
        spectrum to be fit
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

    fitness_func = SimpleSpectrumFitnessFunction(spectrum, model_observation,
                                                 fill_value=fill_value)

    return_square_sum = not (method == 'leastsq')
    if method == 'leastsq':

        fit = optimize.leastsq(fitness_func,
                               np.array(parameter_guesses.values()),
                               args=(parameter_guesses.keys(),
                                     return_square_sum),
                               full_output=True)

        bestfit_spectrum = model_observation()

        result = SimpleLeastsqFitResult(fit[0], fit[1], parameter_guesses, bestfit_spectrum,
                                 spectrum, model_observation, fit)
        return result
    else:
        fit = optimize.minimize(fitness_func,
                               np.array(parameter_guesses.values()),
                               args=(parameter_guesses.keys(),
                                     return_square_sum))

        best_fit_spectrum = model_observation()

        result = SimpleFitResult(fit.x, parameter_guesses, best_fit_spectrum,
                                 spectrum, model_observation, fit, fit.success)

    return result

