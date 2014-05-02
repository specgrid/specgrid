Composite stellar model
=======================


A composite stellar model allows us to create a model that has in addition to
 the grid parameters the ability for other parameters like extinction, rotation,
 doppler shifts, etc::

    >>> from specgrid import specgrid
    >>> from specgrid import plugins
    >>> from specgrid import composite
    >>> spectral_grid = specgrid.MunariGrid('munari_new.h5')

    >>> spectral_grid.parameters = ['teff', 'logg', 'feh']
    >>> rb = plugins.RotationalBroadening() # module for broadening
    >>> ds = plugins.DopplerShift() # module for doppler shift

    # making a model_star out of the grid module followed by rotational
    # broadening and dopplershift
    >>> model_star = composite.ModelStar([spectral_grid, rb, ds])
    >>> model_star.parameters
    ['logg', 'teff', 'feh', 'vrot', 'vrad']

Using `model_star.eval` with different arguments one can get a spectrum for
different parameters. All other parameters are then taken from the
attributes (e.g. `model_star.logg`::

    >>> model_star.eval(teff=5000)
    Spectrum1D(




