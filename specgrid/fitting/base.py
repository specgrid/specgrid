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


def fit_spectrum(spectrum, model_spectrum, fitter='leastsq',
                 fill_value=1e99, valid_slice=slice(None), **guesses):
    """
    Parameters
    ----------

    spectrum: ~specutils.Spectrum1D
        spectrum to be fit
    model_spectrum: ~specgrid.ModelObservation
        model of observation, that returns a spectrum
    fitter: str
        fitter name one of the scipy minimize options or leastsq
    fill_value: ~float
    valid_slice: ~slice
        slice for valid data
    :return:
    """

    if len(guesses) == 0:
        raise ValueError('fit_spectrum requires initial guessses '
                         'for the param_names that should be fit. '
                         'e.g. fit_spectrum(spec, mstar, teff=5000, logg=4)')

    if not set(guesses.keys()).issubset(set(model_spectrum.param_names)):
        raise ValueError('Some parameter guesses are not param_names'
                         ' in the current model star ({0})'.format(
            ', '.join(set(guesses.keys()) - set(model_spectrum.param_names))))


    parameter_guesses = OrderedDict((key, guesses[key])
                                  for key in model_spectrum.param_names
                                  if key in guesses)

    if np.isnan(model_spectrum.evaluate(**parameter_guesses).flux.value[0]):
        raise ValueError('Initial guess ({0}) is outside the confines '
                         'of the grid -- aborting'.format(parameter_guesses))
    def spectral_model_fit(parameter_values):

        parameter_dict = {key:value for key, value in zip(
            parameter_guesses.keys(), parameter_values)}
        model = model_spectrum.evaluate(**parameter_dict)

        if getattr(spectrum, 'uncertainty', None) is not None:
            uncertainty = spectrum.uncertainty
        else:
            uncertainty = np.ones_like(spectrum.flux)

        if np.isnan(model.flux[0]):
            model_flux = np.ones_like(model.flux) * fill_value
        else:
            model_flux = model.flux

        quality = ((spectrum.flux - model_flux) / uncertainty)[valid_slice]
        return quality if fitter == 'leastsq' else (quality**2).sum()

    if fitter == 'leastsq':

        fit = optimize.leastsq(spectral_model_fit,
                               np.array(parameter_guesses.values()),
                               full_output=True)
        bestfit_spectrum = model_spectrum()

        result = SimpleLeastsqFitResult(fit[0], fit[1], parameter_guesses, bestfit_spectrum,
                                 spectrum, model_spectrum, fit)
        return result
    else:
        fit = optimize.minimize(spectral_model_fit,
                                np.array(parameter_guesses.values()),
                                method=fitter)
        best_fit_spectrum = model_spectrum()

        result = SimpleFitResult(fit.x, parameter_guesses, best_fit_spectrum,
                                 spectrum, model_spectrum, fit, fit.success)

    return result
