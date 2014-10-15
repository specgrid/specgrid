from collections import OrderedDict
from logging import getLogger
import numpy as np

logger = getLogger(__name__)

class SimpleSpectrumLikelihood(object):
    """
    A spectrum likelihood object. Initialize with the spectral grid model.
    Calls will returns the log-likelihood for observing a spectrum given the model
    param_names and the uncertainty.

    Parameters
    ----------

    observation: ~model from specpgrid
        model object containing the spectral grid to evaluate


    """
    def __init__(self, spectrum, observation, parameter_names):

        if not set(parameter_names).issubset(set(observation.param_names)):
            raise ValueError('Some of the parameter names ({0}) are not valid '
                             'parameters for the model_observation'.format(
                ', '.join(set(parameter_names) - set(observation.param_names))))


        if set(parameter_names) != set(observation.param_names):
            fixed_parameters = ['{0} {1:.2f}'.format(param_name,
                                                     getattr(observation,
                                                             param_name))
                                 for param_name in (set(observation.param_names)
                                                    -  set(parameter_names))]

            logger.warning('Not all parameters of the observation model '
                           'have an associated prior and are thus fixed:\n {0}'
                           '\n-------'.format('\n'.join(
                fixed_parameters)))
        self.observation = observation


        self.parameter_names = sorted(
            parameter_names,
            key=lambda item: observation.param_names.index(item))

        self.spectrum = spectrum

    @property
    def param_names(self):
        return self.parameter_names


    @property
    def fixed_parameters(self):
        if getattr(self, '_fixed_parameters', None) is None:
            self._fixed_parameters = OrderedDict(
                [(param_name,getattr(self.observation, param_name))
                for param_name in (set(self.observation.param_names)
                                   - set(self.parameter_names))])
        return self._fixed_parameters

    def __call__(self, model_param, ndim, nparam):
        # returns the likelihood of observing the data given the model param_names
        param_dict = OrderedDict([(key, value) for key, value in
                                  zip(self.parameter_names, model_param)])

        m = self.observation.evaluate(**param_dict)

        # log-likelhood for chi-square
        return (-0.5 * ((self.spectrum.flux.value - m.flux.value) /
        self.spectrum.uncertainty.value)**2).sum()


    def __repr__(self):
        repr_str = ("Likelihood Parameters\n----------\n{0}\n\n"
                    "Fixed Parameters\n----------------\n{1}")

        fixed_parameter_str = ['{0} {1}'.format(p_name, p_value)
                               for p_name, p_value
                               in self.fixed_parameters.items()]

        return repr_str.format('\n'.join(self.parameter_names),
                               '\n'.join(fixed_parameter_str))


class IsochroneSpectrumLikelihood(SimpleSpectrumLikelihood):

    def __init__(self, spectrum, observation_model, v, bv, isochrone_interpolator,
                 spectrum_parameter_names,
                 isochrone_parameter_names=['age', 'feh', 'mass', 'distance',
                                            'ebv']):
        super(IsochroneSpectrumLikelihood, self).__init__(
            spectrum, observation_model, spectrum_parameter_names)

        self.v = v
        self.bv = bv
        self.isochrone_interpolator = isochrone_interpolator
        self.parameter_names = isochrone_parameter_names


    def __call__(self, model_param, ndim, nparam):
        isochrone_values = self.isochrone_interpolator.interpolate_point(
            model_param[0], model_param[1], model_param[2])

        v_likelihood = self._calculate_v_log_likelihood(isochrone_values.m_v,
                                                        model_param[3])
        likelihood = v_likelihood

        bv_likelihood = self._calculate_bv_log_likelihood(isochrone_values.bv,
                                                          model_param[4])

        likelihood += bv_likelihood

        spectrum_param_dict = OrderedDict([('teff', isochrone_values.teff),
                                           ('logg', isochrone_values.logg),
                                           ('feh', model_param[1])])

        m = self.observation.evaluate(**spectrum_param_dict)
        spec_likelihood = (-0.5 * ((self.spectrum.flux.value - m.flux.value) /
        self.spectrum.uncertainty.value)**2).sum() / self.spectrum.flux.shape[0]

        likelihood += spec_likelihood

        return likelihood


    def _calculate_v_log_likelihood(self, isochrone_m_v, distance):
        model_v = isochrone_m_v + (5*np.log10(distance) - 5)
        return -0.5 * ((model_v - self.v[0]) / self.v[1])**2

    def _calculate_bv_log_likelihood(self, isochrone_bv, ebv):
        model_bv = isochrone_bv + ebv
        return -0.5 * ((model_bv - self.bv[0]) / self.bv[1])**2

