#specgrid class
from scipy import interpolate
import numpy as np
try:
    from pyspec import oned
    pyspec_available = True
except ImportError:
    pyspec_available = False

import pandas as pd
import os
import h5py
from astropy import units as u

from astropy import modeling


class SpectralGrid(modeling.ParametricModel):
    """
    A SpectralGrid interpolation class. Can serve as a model maybe

    Parameters
    ----------

    grid_hdf5_fname: filename for HDF5 File
    """

    def __init__(self, grid_hdf5_fname, interpolator = interpolate.LinearNDInterpolator):
        super(SpectralGrid, self).__init__()
        if not os.path.exists(grid_hdf5_fname):
            raise ValueError('{0} does not exists'.format(grid_hdf5_fname))
        self.grid_hdf5_fname = grid_hdf5_fname
        self.index = pd.read_hdf(self.grid_hdf5_fname, 'index')
        self.params = None
        self.flux = None

        self.interpolator = interpolator
        self.interpolate_grid = None

        self.grid_param_names = self.index.columns

        #cleaning param_names of
        self.grid_param_names = self.grid_param_names.drop(['id', 'fname'])

        with h5py.File(self.grid_hdf5_fname, 'r') as h5file:
            self.hidden_parameters = h5file['spectra'].attrs['hidden_parameters']
            self.wave = h5file['spectra'].attrs['wave'] * u.Unit(h5file['spectra'].attrs['wave_unit'])

        self.grid_param_names = self.grid_param_names.drop(self.hidden_parameters)

        for param_name in self.grid_param_names:
            self.__setattr__(param_name, modeling.Parameter(param_name, default=0))


    @property
    def current_parameters(self):
        return [self.__getattribute__(param_name) for param_name in self.grid_param_names]

    def load_dataset(self, query_string):
        """
        Load the dataset

        Parameters
        ----------

        query_string: str
            similar as pandas ???????

        """
        selected_index = self.index.query(query_string).index

        with h5py.File(self.grid_hdf5_fname, 'r') as h5file:
            self.flux = h5file['spectra/block0_values'][list(selected_index)]

        self.params = self.index[self.grid_param_names].ix[selected_index].values

        self.interpolate_grid = self.interpolator(self.params, self.flux)

    def __call__(self, wavelength):
        wavelength = u.Quantity(wavelength, u.angstrom)

        return self.interpolate_grid(*self.current_parameters)


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