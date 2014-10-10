********************************
Simple fit of a stellar spectrum
********************************

There are a variety of ways to simply fit spectra with specgrid. The simplest
way is to use one of the scipy minimzers::

    >>> from specgrid import Spectrum1D,
    >>> from specgrid import Spectrum1D, SpectralGrid, assemble_observation
    >>> spec_grid = SpectralGrid('munari.h5')
    >>> model_star = assemble_observation(spec_grid, plugin_names=['doppler', 'rotation', 'resolution'])

    >>>
    >>> from specgrid import plugins
    >>> from specgrid import composite
    >>> from specgrid import fitting
    >>> spectral_grid = specgrid.MunariGrid('munari_new.h5')

    >>> spectral_grid.parameters = ['teff', 'logg', 'feh']

    >>> model_star = composite.ModelStar([spectral_grid])


    >>> my_spectrum = model_star.eval(teff=4580., logg=3.0, feh=0.0)

    >>> my_spectrum.data *= random.normal(1, 0.3, my_spectrum.flux.shape)

    >>> fitting.fit_spectrum(my_spectrum, dict(teff=5780., logg=4.4, feh=-1.), model_star)
    (OrderedDict([('logg', 3.0000126710686539), ('teff', 4576.0798460409605), ('feh', 0.010206474079974394)]),
     OrderedDict([('logg', 5.5321554013186524e-07), ('teff', 1.6274642182531751e-05), ('feh', 1.6156094042917857e-07)]),
     (array([  3.00001267e+00,   4.57607985e+03,   1.02064741e-02]),
      array([[  3.06047434e-13,   2.53529237e-12,   1.26205508e-14],
             [  2.53529237e-12,   2.64863978e-10,  -1.12665953e-12],
             [  1.26205508e-14,  -1.12665953e-12,   2.61019375e-14]]),
      {'fjac': array([[ -7.18599597e+06,  -2.25474628e+04,   1.37052309e+03, ...,
                -6.56633417e-03,  -8.69272172e-03,  -8.85865277e-03],
              [  5.96878859e+05,   1.88384566e+06,  -5.51902911e+02, ...,
                -2.19135566e-03,  -2.49293194e-03,  -2.55636110e-03],
              [ -3.62806347e+04,  -1.80322729e+04,   6.14452833e+04, ...,
                 1.81251095e-03,   4.53151590e-04,   1.57638325e-04]]),
       'fvec': array([   -72.65181911,   -225.683385  ,    117.16965947, ...,
               31788.9843088 , -70767.51640272,  -7729.5426484 ]),
       'ipvt': array([3, 1, 2], dtype=int32),
       'nfev': 63,
       'qtf': array([-24452.54863629,  40649.33237999,  15681.20523754])},
      'Both actual and predicted relative reductions in the sum of squares\n  are at most 0.000000',
      1))