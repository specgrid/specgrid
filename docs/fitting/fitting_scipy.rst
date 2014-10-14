****************************************************
Simple fit of a stellar spectrum with scipy.optimize
****************************************************

The simplest way to fit a spectrum is to use one of the scipy minimzers.
It uses the minimizers given in the `scipy.optimize <http://docs.scipy.org/doc/scipy-0.14.0/reference/optimize.html>`_ package.

 ::

    >>> from specgrid import Spectrum1D, SpectralGrid, assemble_observation, fitting
    >>> spec_grid = SpectralGrid('munari.h5')
    >>> model_observation = assemble_observation(spec_grid, plugin_names=['doppler', 'rotation', 'resolution'])

    >>> my_spectrum = model_observation.evaluate(teff=4580., logg=3.0, feh=0.0)

    >>> my_spectrum.data *= random.normal(1, 0.3, my_spectrum.flux.shape) #making my own spectrum

    >>> result = fitting.fit_spectrum(my_spectrum, model_observation, teff=5780., logg=4.4, feh=-1.)
    >>> result
    Fit Result:
    Fit successful
    teff 4585.617 +/- 0.000 (5780.000)
    logg 3.020 +/- 0.000 (4.400)
    feh 0.038 +/- 0.000 (-1.000)

    >>> result = fitting.fit_spectrum(my_spectrum, model_observation, teff=5780., logg=4.4, feh=-1., method='Nelder-Mead')

.. automodapi:: specgrid.fitting.base
    :no-inheritance-diagram: