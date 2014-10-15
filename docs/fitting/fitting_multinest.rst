***********************************************
Simple fit of a stellar spectrum with MultiNest
***********************************************

Multinest is a very useful tool when fitting spectra. This first important step
to use specgrid with MultiNest is to install `PyMultiNest <http://johannesbuchner.github.io/PyMultiNest/>`_

Once PyMultiNest is running, fitting a spectrum is relatively easy to fit a spectrum with it::

    >>> from specgrid import Spectrum1D, SpectralGrid, assemble_observation, fitting
    >>> from specgrid.fitting.multinest import priors, fit_simple_spectrum_multinest

    >>> spec_grid = SpectralGrid('munari.h5')
    >>> model_observation = assemble_observation(spec_grid, plugin_names=['doppler', 'rotation', 'resolution'])

    >>> my_spectrum = model_observation.evaluate(teff=4580., logg=3.0, feh=0.0)

    >>> my_spectrum.data *= random.normal(1, 0.3, my_spectrum.flux.shape) #making my own spectrum
    >>> my_spectrum.uncertainty = np.sqrt(my_spectrum.flux.value)

    >>> prior_dict = dict(teff=priors.UniformPrior(5000, 6000), logg=priors.GaussianPrior(4.3, 0.3), feh=priors.FixedPrior(0.05))
    >>> multinest_fitter = fit_simple_spectrum_multinest(my_spectrum,
                                  model_observation, prior_dict)
    specgrid.fitting.multinest.likelihoods - WARNING - Not all parameters of the observation model have an associated prior and are thus fixed:
     R inf
    vrad 0.00 km / s
    vrot 0.00 km / s
    -------
    specgrid.fitting.multinest.multinest_fitter - INFO - Starting fit in /var/folders/v2/z0nm04hj495cj432_4xwvzmh0000gn/T/tmpiH4uVM with prefix specgrid_multinest
     MultiNest Warning: no resume file found, starting from scratch
     *****************************************************
     MultiNest v3.7
     Copyright Farhan Feroz & Mike Hobson
     Release June 2014

     no. of live points =  400
     dimensionality =    3
     *****************************************************




