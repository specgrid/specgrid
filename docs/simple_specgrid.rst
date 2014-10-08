*********************************
Simple spectral grid interpolation
*********************************

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

    >>> spec = munari.evaluate(5780, 4.4, 0.0)
    >>> spec
    Spectrum1D([ 587668.18 ,  578633.564,  678022.4  , ...,  920261.52 ,
                 912980.964,  885117.936])
    >>> spec.wavelength
    <Quantity [  2500.5,  2501.5,  2502.5,...,  10497.5, 10498.5, 10499.5] Angstrom>
    >>> spec.flux
    <Quantity [ 587668.18 , 578633.564, 678022.4  ,...,  920261.52 ,
            912980.964, 885117.936] erg / (Angstrom cm2 s)>
    >>> munari.teff
    5780.0


As can be seen this also sets the parameters. The most generic way to
interpolate is to just call the object, which will use the currently set
attributes as the point of interpolation (currently set to solar values from the
last example)::

    >>> munari()
    Spectrum1D([ 587668.18 ,  578633.564,  678022.4  , ...,  920261.52 ,
                 912980.964,  885117.936])

The final way to interpolate is to use the evaluate with keyword arguments. This
will only update the selected keyword values and leave the rest the same::

    >>> munari.evaluate(teff=7000)
    Spectrum1D([ 4185839.8,  4061805.8,  4112774.6, ...,  1457224.2,
                 1444824. ,  1400608.6])
    >>> munari.teff
    7000.0
    >>> munari.logg
    4.4
    >>> munari.feh
    0.0

