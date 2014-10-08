from collections import OrderedDict
from astropy.units import Quantity

from specgrid import plugins


class ModelStar(object):
    """
    A model star combines a normal spectral grid (~specgrid.SpectralGrid) with a
    number of plugins that allow to add further physical and instrumental effects
    (e.g. normalizing, dopplershift, rotation, etc.)

    Parameters
    ----------

    spectral_grid: ~specgrid.SpectralGrid

    """

    param2model = OrderedDict()

    def __init__(self, spectral_grid, astrophysics_plugins=[],
                 instrument_plugins=[]):
        self.spectral_grid = spectral_grid
        self.astrophysics_plugins = astrophysics_plugins
        self.instrument_plugins = instrument_plugins


        self.param2model = OrderedDict()

        for model in ([spectral_grid] + self.all_plugins):
            self.param2model.update(OrderedDict([(param, model)
                                          for param in model.parameters]))
        self.parameters = self.param2model.keys()

    @property
    def all_plugins(self):
        return self.astrophysics_plugins + self.instrument_plugins

    def __getattr__(self, item):
        if item in self.param2model:
            return getattr(self.param2model[item], item)
        else:
            return super(ModelStar, self).__getattribute__(item)

    def __setattr__(self, item, value):
        if item in self.param2model:
            return setattr(self.param2model[item], item, value)
        else:
            super(ModelStar, self).__setattr__(item, value)

    def _call_astrophysics_plugins(self, spectrum):
        """
        Evaluate the composite model of all astrophysics plugins

        Parameters
        ----------

        spectrum: ~specutils.Spectrum1D
        """

        for model in self.astrophysics_plugins:
            spectrum = model(spectrum)

        return spectrum

    def __call_instrument_plugins(self, spectrum):
        """
        Evaluate the composite model of all instrument plugins

        Parameters
        ----------

        spectrum: ~specutils.Spectrum1D
        """

        for model in self.instrument_plugins:
            spectrum = model(spectrum)

        return spectrum

    def __call__(self):
        """
        Evaluate the spectrum with the current parameters
        """
        spectrum = self.spectral_grid()

        spectrum = self._call_astrophysics_plugins(spectrum)
        spectrum = self.__call_instrument_plugins(spectrum)

        return spectrum
    
    def evaluate(self, *args, **kwargs):
        """
        Interpolating on the grid to the necessary parameters

        Examples
        --------

        This can either be called with arguments ``specgrid.evaluate(5780, 4.4, -1)`` or
        using keyword way of calling (then not all parameters have to be given)
        ``specgrid.evaluate(logg=4.4)``
        """


        if len(args) > 0:
            if len(kwargs) > 0:
                raise ValueError('One can either use arguments or '
                                 'keyword arguments not both')
            if len(args) != len(self.parameters):
                raise ValueError(
                    'evaluate() takes {0} arguments '
                    'for each parameter ({1}) - {2} given'.format(
                        len(self.parameters),
                        ', '.join(self.parameters),
                        len(args)))
            for param_name, value in zip(self.parameters, args):
                setattr(self, param_name, value)

        for key in kwargs:
            if key not in self.parameters:
                raise ValueError('{0} not a parameter of the current '
                                 'model_star (parameters are {1})'.format(
                    key, ','.join(self.parameters)))
            setattr(self, key, kwargs[key])

        return self.__call__()



def assemble_model_star(spectral_grid, spectrum=None, normalize_npol=None, plugin_names=[]):
    """

    Parameters
    ----------

    spectral_grid: ~specgrid.SpectralGrid
        spectral grid to be used in model_star
    spectrum: ~specutils.Spectrum1D
        spectrum to be used for interpolation, if None neither interpolation nor
            will be performed [default None]

    normalize_npol: int
        degree of polynomial to be used for interpolation, only if not None and
        spectrum is not None will the normalization plugin be used [default None]

    plugin_names: ~list of ~str
        select between the following available plugin choices:
        {plugin_names}

    Returns
    -------
        : ~ModelStar


    """

    astrophysics_plugins = []
    instrument_plugins = []

    for plugin_name in plugin_names:
        if plugin_name in plugins.astrophysics_plugins:
            current_plugin = plugins.astrophysics_plugins[plugin_name]()
            astrophysics_plugins.append(current_plugin)

        elif plugin_name in plugins.instrument_plugins:
            current_plugin = plugins.instrument_plugins[plugin_name]()

            instrument_plugins.append(current_plugin)

    astrophysics_plugins = sorted(
        astrophysics_plugins,
        key=lambda item: plugins.astrophysics_plugins.values().index(item.__class__))

    instrument_plugins = sorted(
        instrument_plugins,
        key=lambda item: plugins.instrument_plugins.values().index(item.__class__))

    if spectrum is not None:
        instrument_plugins += [plugins.Interpolate(spectrum)]

    if not (spectrum is None or normalize_npol is None):
        instrument_plugins += [plugins.Normalize(spectrum, npol=normalize_npol)]

    model_star = ModelStar(spectral_grid, astrophysics_plugins,
                           instrument_plugins)

    return model_star

assemble_model_star.__doc__ = assemble_model_star.__doc__.format(
    plugin_names=', '.join(plugins.astrophysics_plugins.keys() +
                           plugins.instrument_plugins.keys()))