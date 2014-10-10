from specgrid.multinest.priors import PriorCollections
from specgrid.fitting.multinest.likelihoods import SimpleSpectrumLikelihood


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
                         'the observation model')

    likelihood = SimpleSpectrumLikelihood(spectrum, observation)




