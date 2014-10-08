Model star
==========


A composite stellar model allows us to create a model that has in addition to
the grid parameters (often Teff, log(g) and [Fe/H]) the ability to add other
astrophysical transformations like doppler shift, extinction, etc. as well as
instrumental properties like resolution, normalization, etc. Here is a simple
way to create such a model star::

    >>> from specgrid import SpectralGrid, assemble_model_star
    >>> spec_grid = SpectralGrid('munari.h5')
    >>> model_star = assemble_model_star(spec_grid, plugin_names=['doppler', 'rotation', 'resolution'])
    >>> model_star.parameters
    ['teff', 'logg', 'feh', 'vrot', 'vrad', 'R']

Similar to the basic spectral grid, calling the model_star will evaluate it with
the current parameters

    >>> model_star()
    Spectrum1D([  7.78958400e-01,   9.68851000e-01,   1.73870800e+00, ...,
                  1.70572100e+05,   1.72696300e+05,   1.67735400e+05])





