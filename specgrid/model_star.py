from collections import OrderedDict
from astropy.units import Quantity

from specgrid import plugins

class SpecGridCompositeModel(object):
    param2model = OrderedDict()

    def __init__(self, models=[]):
        self.models = models

        self.param2model = OrderedDict()

        for model in models:
            self.param2model.update(OrderedDict(
                [(param, model) for param in model.param_names]))

    @property
    def param_names(self):
        return self.param2model.keys()

    @property
    def param_dict(self):
        return OrderedDict([(param_name, getattr(self, param_name))
                            for param_name in self.param_names])

    def __getattr__(self, item):
        if item in self.param2model:
            return getattr(self.param2model[item], item)
        else:
            return super(SpecGridCompositeModel, self).__getattribute__(item)

    def __setattr__(self, item, value):
        if item in self.param2model:
            return setattr(self.param2model[item], item, value)
        else:
            super(SpecGridCompositeModel, self).__setattr__(item, value)


    def _set_parameters(self, *args, **kwargs):
        """
        Setting the parameters

        Examples
        --------

        This can either be called with arguments ``specgrid._set_parameters(5780, 4.4, -1)`` or
        using keyword way of calling (then not all param_names have to be given)
        ``specgrid._set_parameters(logg=4.4)``
        """



        if len(args) > 0:
            if len(kwargs) > 0:
                raise ValueError('One can either use arguments or '
                                 'keyword arguments not both')
            if len(args) != len(self.param_names):
                raise ValueError(
                    'evaluate() takes {0} arguments '
                    'for each parameter ({1}) - {2} given'.format(
                        len(self.param_names),
                        ', '.join(self.param_names),
                        len(args)))
            for param_name, value in zip(self.param_names, args):
                setattr(self, param_name, value)

        for key in kwargs:
            if key not in self.param_names:
                raise ValueError('{0} not a parameter of the current '
                                 'observation (param_names are {1})'.format(
                    key, ','.join(self.param_names)))
            setattr(self, key, kwargs[key])

    def __call__(self, spectrum):

        for model in self.models:
            spectrum = model(spectrum)

        return spectrum

class ModelStar(SpecGridCompositeModel):
    """
    A model star combines a normal spectral grid (~specgrid.SpectralGrid) with a
    number of plugins that allow to add further physical and instrumental effects
    (e.g. normalizing, dopplershift, rotation, etc.)

    Parameters
    ----------

    spectral_grid: ~specgrid.SpectralGrid

    """

    def __init__(self, spectral_grid, plugins=[]):
        self.spectral_grid = spectral_grid
        super(ModelStar, self).__init__([spectral_grid] + plugins)
        self.models = self.models[1:]

    def __call__(self):
        """
        Evaluate the spectrum with the current param_names
        """
        spectrum = self.spectral_grid()

        for model in self.models:
            spectrum = model(spectrum)

        return spectrum

    def evaluate(self, *args, **kwargs):
        """
        Interpolating on the grid to the necessary param_names

        Examples
        --------

        This can either be called with arguments ``specgrid.evaluate(5780, 4.4, -1)`` or
        using keyword way of calling (then not all param_names have to be given)
        ``specgrid.evaluate(logg=4.4)``
        """

        self._set_parameters(*args, **kwargs)

        return self.__call__()

    def __repr__(self):

        param_str = '\n'.join(['{0} {1}'.format(key, value)
                     for key, value in self.param_dict.items()])
        return "Model Star Parameters:\n\n{0}".format(param_str)


class ModelInstrument(SpecGridCompositeModel):
    def evaluate(self, spectrum, *args, **kwargs):
        """
        Interpolating on the grid to the necessary param_names

        Examples
        --------

        This can either be called with arguments ``specgrid.evaluate(5780, 4.4, -1)`` or
        using keyword way of calling (then not all param_names have to be given)
        ``specgrid.evaluate(logg=4.4)``
        """

        self._set_parameters(*args, **kwargs)

        return self.__call__(spectrum)

    def __repr__(self):

        param_str = '\n'.join(['{0} {1}'.format(key, value)
                     for key, value in self.param_dict.items()])
        return "Model Instrument Parameters:\n\n{0}".format(param_str)


class Observation(SpecGridCompositeModel):
    def __init__(self, model_star, model_instrument):
        self.model_star = model_star
        self.model_instrument = model_instrument

        self.param2model = self.model_star.param2model.copy()
        self.param2model.update(self.model_instrument.param2model.copy())

    def __call__(self):
        spectrum = self.model_star()
        return self.model_instrument(spectrum)

    @property
    def all_plugins(self):
        return self.model_star.models + self.model_instrument.models

    def __repr__(self):
        return "Model Observation:\n\n{0}\n\n{1}".format(self.model_star, self.model_instrument)

    def evaluate(self, *args, **kwargs):
        self._set_parameters(*args, **kwargs)
        return self.model_instrument(self.model_star())

def assemble_observation(spectral_grid, plugin_names=[], spectrum=None, normalize_npol=None, ):
    """

    Parameters
    ----------

    spectral_grid: ~specgrid.SpectralGrid
        spectral grid to be used in observation
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
        : ~Observation


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

    model_star = ModelStar(spectral_grid, astrophysics_plugins)

    model_instrument = ModelInstrument(instrument_plugins)

    return Observation(model_star, model_instrument)

assemble_observation.__doc__ = assemble_observation.__doc__.format(
    plugin_names=', '.join(plugins.astrophysics_plugins.keys() +
                           plugins.instrument_plugins.keys()))