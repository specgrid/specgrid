from collections import OrderedDict

#likelihood = SimpleLikelihood(self.spectrum, model_star, self.parameter_names)


class SimpleSpectrumLikelihood(object):
    """
    A spectrum likelihood object. Initialize with the spectral grid model.
    Calls will returns the log-likelihood for observing a spectrum given the model
    param_names and the uncertainty.

    Parameters
    ----------

    model_star: ~model from specpgrid
        model object containing the spectral grid to evaluate


    """
    def __init__(self, spectrum, model_star, parameter_names):
        self.model_star = model_star
        self.parameter_names = parameter_names
        self.data = spectrum


    def __call__(self, model_param, ndim, nparam):
        # returns the likelihood of observing the data given the model param_names
        param_dict = OrderedDict([(key, value) for key, value in
                                  zip(self.parameter_names, model_param)])

        m = self.model.evaluate(**param_dict)

        # log-likelhood for chi-square
        return (-0.5 * ((self.data.flux.value - m.flux.value) /
        self.data.uncertainty.value)**2).sum()


