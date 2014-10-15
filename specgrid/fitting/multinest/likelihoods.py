from collections import OrderedDict
from logging import getLogger
#likelihood = SimpleLikelihood(self.spectrum, observation, self.parameter_names)

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