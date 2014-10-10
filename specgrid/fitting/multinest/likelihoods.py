from collections import OrderedDict

#likelihood = SimpleLikelihood(self.spectrum, observation, self.parameter_names)


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
                             'parameters for the model_observation')

        self.observation = observation


        self.parameter_names = parameter_names
        self.spectrum = spectrum

    @property
    def param_names(self):
        return self.parameter_names

    def __call__(self, model_param, ndim, nparam):
        # returns the likelihood of observing the data given the model param_names
        param_dict = OrderedDict([(key, value) for key, value in
                                  zip(self.parameter_names, model_param)])

        m = self.observation.evaluate(**param_dict)

        # log-likelhood for chi-square
        return (-0.5 * ((self.spectrum.flux.value - m.flux.value) /
        self.spectrum.uncertainty.value)**2).sum()


