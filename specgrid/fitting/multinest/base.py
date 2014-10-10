from specgrid.multinest.priors import PriorCollections
from specgrid.multinest.likelihoods import SimpleSpectrumLikelihood


def fit_simple_spectrum_multinest(spectrum, specgrid, prior_grid):
    """

    Parameters
    ----------

    spectrum: ~Spectrum1D from specpgrid
        Spectrum1D object with the observed spectrum

    priors: ~dict
        A dictionary with the param_names to fit as well as their priors as
        implemented in the prior classes available (UniformPrior, GaussianPrior,
        PoissonPrior,FixedPrior)


    """

    model_star

