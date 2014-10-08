**************************
Munari Specgrid Operations
**************************

The first step is too obtain an appropriate specgrid from the author. The
munari specgrid is publicly available
`here <http://moria.astro.utoronto.ca/~wkerzend/files/munari.h5>`_.

We will be dealing with the munari spectral grid here.
It is very easy to open the spectral grid

    >>> from specgrid import SpectralGrid
    >>> munari_grid = SpectralGrid('munari.h5')

The default values are Teff=5780, log(g)=4.4, [Fe/H] = -1., so if we just
call the grid object we will get a sun-like star back (as a `Spectrum1D`-object)::

    >>> my_spectrum = munari_grid
