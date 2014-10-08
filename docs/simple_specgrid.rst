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
    Opening munari.h5 in read-only mode
    >>> munari.parameters
    ['teff', 'logg', 'feh']
    >>> munari.teff
    3500.

The easiest way do an interpolation is to use the `SpectralGrid.evaluate` method::

    >>> munari.evaluate(5780, 4.4, 0.0)
    Spectrum1D([ 587668.18 ,  578633.564,  678022.4  , ...,  920261.52 ,
             912980.964,  885117.936])
    >>> munari.teff
    5780.0

As can be seen this also sets the parameters. The most generic way to
interpolate is to just call the object

The default values are Teff=5780, log(g)=4.4, [Fe/H] = -1., so if we just
call the grid object we will get a sun-like star back (as a `Spectrum1D`-object)::

    >>> my_spectrum = munari_grid
