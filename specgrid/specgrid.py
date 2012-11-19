#specgrid class
from scipy import interpolate
import numpy as np
try:
    from pyspec import oned
    pyspec_available = True
except ImportError:
    pyspec_available = False
    
class specgrid(object):
    def __init__(self, params, fluxes, wave, param_names,
                 interpolator = interpolate.LinearNDInterpolator,
                 normalizer = None, convolver=None, metric={}, inverse_metric={}, fill_value=-1):
        
        self.param_names = param_names
        self.wave = wave


        #setting up limits
        self.param_mins = params.min(axis=0)
        self.param_maxs = params.max(axis=0)

        #setting up metrics
        self.metric = dict([(item, lambda x: x) for item in self.param_names])
        self.inverse_metric = dict([(item, lambda x: x) for item in self.param_names])
        
        
        #inserting metrics
        self.metric.update(metric)
        self.inverse_metric.update(inverse_metric)

        self.params = params
        
        #applying metric
        for i, param_name in enumerate(self.param_names):
            self.params[:,i] = self.inverse_metric[param_name](self.params[:,i])
            
        self.fluxes = fluxes
        self.normalizer = normalizer
        self.convolver = convolver
        self.interpolated_grid = interpolator(params, fluxes, fill_value=fill_value)
        
        #creating sets
        self.param_name_set = set(self.param_names)
        self.all_param_name_set = set(self.param_names)
        
        self.plugins = {}
        
        
    def add_plugin(self, plugin_class, **kwargs):
        plugin_object = plugin_class(self.wave, normalizer=self.normalizer, **kwargs)
        self.plugins[plugin_object.param_name] = plugin_object
        self.all_param_name_set.add(plugin_object.param_name)
    
    
    def interpolate_spectrum(self, **kwargs):
        if pyspec_available:
            flux = self.interpolate(**kwargs)
            wave = self.wave
            return oned.onedspec(wave, flux, mode='waveflux')
        else:
            raise NotImplementedError('PySpec is not available')
        
    def interpolate(self, **kwargs):
        #checking if the kwargs are in the param_names
        req_param_names = set(kwargs.keys())
        
        if req_param_names < self.param_name_set:
            raise TypeError('Missing interpolation arguments need to supply keywords: ' + ','.join(self.param_names))
        
        
        if not req_param_names.issubset(self.all_param_name_set):
            raise TypeError('Requesting parameters that are not grid')
        
        
        req_param_values = [self.inverse_metric[name](kwargs[name]) for name in self.param_names]
        
        interp_flux = self.interpolated_grid(*req_param_values)
        
        for key, value in kwargs.items():
            if key in self.plugins.keys():
                plugin = self.plugins[key]
                #plugin.inverse_metric
                interp_flux = plugin(interp_flux, value)
        return interp_flux