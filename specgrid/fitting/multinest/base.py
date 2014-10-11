from specgrid.fitting.multinest.priors import PriorCollection
from specgrid.fitting.multinest.likelihoods import SimpleSpectrumLikelihood

from specgrid.fitting.multinest import BaseMultinestFitter


def fit_simple_spectrum_multinest(spectrum, observation, priors):
    """

    Parameters
    ----------

    spectrum: ~Spectrum1D
        Spectrum1D object with the observed spectrum

    observation: ~specgrid.Observation
        an observation model

    priors: ~dict
        A dictionary with the param_names to fit as well as their priors as
        implemented in the prior classes available (UniformPrior, GaussianPrior,
        PoissonPrior,FixedPrior)
    """

    if not set(priors.keys()).issubset(set(observation.param_names)):
        raise ValueError('Given priors ({0}) are not valid parameters in '
                         'the observation model'.format(
            set(priors.keys()) - set(observation.param_names)))

    likelihood = SimpleSpectrumLikelihood(spectrum, observation, priors.keys())
    priors_list = [priors[param_name]
                   for param_name in observation.param_names
                   if param_name in priors]

    prior_collection = PriorCollection(priors_list)

    return BaseMultinestFitter(likelihood, prior_collection)






