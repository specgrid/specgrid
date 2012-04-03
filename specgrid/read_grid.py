#Reading Grid into memory and grid information

import os
import ConfigParser
import numpy as np
import sqlite3
import numpy as np
import specgrid
from glob import glob
from pyspec import oned
try:
    import sqlparse
    sqlparse_available = True
except ImportError:
    sqlparse_available = False
    
    
metric_dict = {'lin':lambda x:x,
          'log10': np.log10}

inverse_metric_dict = {'lin':lambda x:x,
                  'log10':lambda x:10**x}

def read_general_config(conf_path='~/.specgrid/specgrid.ini', general_dict={}):
    conf_path = os.path.expanduser(conf_path)
    if not os.path.exists(conf_path):
        raise IOError('Config directory does not exist at %s' % conf_path)

    conf = ConfigParser.ConfigParser()
    conf.read(conf_path)
    
    for item, value in conf.items('general'):
        if item == 'warn_thresh':
            value = int(value) * 1024**2
        general_dict[item] = value
    return general_dict

def read_grid_config(fname, grid_dict={}):
    fname = os.path.expanduser(fname)
    if not os.path.exists(fname):
        raise IOError('No file found at %s' % fname)
    
    conf = ConfigParser.ConfigParser()
    conf.read(fname)
    for grid in conf.sections():
        grid_dict[grid] = {}
        for item, value in conf.items(grid):
            if item == 'datatype':
                if value == 'float64': value = np.float64
                else: raise TypeError('Type %s not supported' % value)
            
            elif item == 'specsize':
                value = int(value)
            
            elif item == 'indexdb' or item == 'datadir':
                value = os.path.expanduser(value)
            
            elif item == 'r':
                value = np.float64(value)
            
            elif item.startswith('default_'):
                #item = item.replace('default_', '')
                try:
                    value = np.float64(value)
                except ValueError:
                    value = value
            grid_dict[grid][item] = value
    return grid_dict

def read_standard_grid_configs(conf_path='~/.specgrid/', grid_dict={}):
    conf_path = os.path.expanduser(conf_path)
    if not os.path.exists(conf_path):
        raise IOError('Config directory does not exist at %s' % conf_path)
    conf = ConfigParser.ConfigParser()
    for fname in glob(os.path.join(conf_path,'*.ini')):
        if 'specgrid.ini' in fname: continue
        print fname
        grid_dict = read_grid_config(fname, grid_dict)
    
    return grid_dict


def read_grid(grid_name, **kwargs):
    """Reading in the spectral grid and creating a specgrid dir
    Parameters:
    -----------
    grid_name: str
        Name of the grid to be read
        
    grid_dict: dict
        Specific grid configuration
    general_dict: dict
        General configuration for all grids,
    
    warn_size: bool
        Warn if requested dataset is above warning threshold. Interactive keypress required
        Threshold can be set in specgrid.ini
        
    ignore: tuple
        ignore the parameter in the select and grid.
        
    normalizer: normalize object
        normalizer object"""
    
    #Readint grid configuration
    
    warn_size = kwargs.pop('warn_size', True)
    grid_dict = kwargs.pop('grid_dict', None)
    general_dict = kwargs.pop('general_dict', None)
    normalizer = kwargs.get('normalizer', None)
    convolver = kwargs.get('convolver', None)
    ignore = kwargs.pop('ignore', ())
    
    if general_dict == None:
        general_dict = read_general_config()
    if grid_dict == None:
        grid_dict = read_standard_grid_configs()
    
    
    #Reading configuration dictionary for specific grid
    config_dict = grid_dict[grid_name]
    table_name = config_dict['table']
    conn = sqlite3.connect(config_dict['indexdb'])
    
    #getting parameter names
    param_names = map(str, zip(*conn.execute('pragma table_info(%s)' % table_name).fetchall())[1])
    param_names.remove('id')
    param_names.remove('fname')
    param_values = []
    requested_param_names = []
    #extracting the parameters for grid (using defaults if none are given)
    for param_name in param_names:
        
        default_param_value = config_dict.get('default_' + param_name, None)
        
        if param_name in ignore:
            continue
        
        
        
        
        if param_name in kwargs.keys():
            requested_param_names.append(param_name)
            param_value = kwargs[param_name]

        elif default_param_value == None:
            requested_param_names.append(param_name)
            param_value = None
        else:
            param_value = default_param_value
            
        param_values.append(param_value)
        
    #building sqlite query
    
    grid_query_conditions = []

    #now add the limitations
    
    for param_name, param_value in zip(param_names, param_values):
        #specifically requesting a None selection on a parameter
                
        if param_value == None:
            continue
        
        if param_name in ignore:
            continue

        
        if isinstance(param_value, basestring):
            grid_query_conditions.append("%s='%s'" % (param_name, param_value))
            try:
                requested_param_names.remove(param_name)
            except ValueError:
                pass
            
            continue
            
        param_len = None
        try:
            param_len = len(param_value)
        except TypeError:
            pass
        #detected value selection
        if param_len == None:
            grid_query_conditions.append('%s=%s' % (param_name, param_value))
            if param_name in requested_param_names:
                requested_param_names.remove(param_name)
        #detected range selection
        elif param_len == 2:
            
            if param_value.count(None) == 0:
                grid_query_conditions.append('%s between %s and %s' % (param_name, param_value[0], param_value[1]))
                
            elif param_value.count(None) == 1:
                if param_value[0] == None:
                    grid_query_conditions.append('%s <= %s' % (param_name, param_value[1]))
                elif param_value[1] == None:
                    grid_query_conditions.append('%s >= %s' % (param_name, param_value[0]))
            
            elif param_value.count(None) > 1:
                raise ValueError('Range tuple can only have one \'None\' value')
            
            
        else: raise ValueError('params only support tuple or single numbers')


    grid_query = 'select %s, fname from %s' % (','.join(requested_param_names), table_name)
    
    if len(grid_query_conditions) > 0:
        grid_query += ' where %s' % ' and '.join(grid_query_conditions)
    
    grid_query += ' order by %s' % ','.join(requested_param_names)
    
    
    print 'Constructed database query out of input:\n' + '-'*40 
    if sqlparse_available:
        print sqlparse.format(grid_query, reindent=True, keyword_case='upper')
    else:
        print grid_query
    print '-'*40 + '\n'
    
    raw_data = conn.execute(grid_query).fetchall()
    requested_params = np.array(zip(*zip(*raw_data)[:-1]))

    fnames = zip(*raw_data)[-1]
    
    grid_points = len(fnames)
    
    if grid_points == 0:
        raise ValueError('Your query did not return any gridpoints. Data requested lies outside grid.')
    
    print "Found %d gridpoints" % grid_points
    data_size = grid_points * config_dict['specsize']
    
    
    
    if data_size < 1024:
        data_size_string = '%d bytes' % data_size
    elif data_size < 1024**2:
        data_size_string = '%.2f kilobytes' % (data_size / 1024.)
    elif data_size < 1024**3:
        data_size_string = '%.2f megabytes' % (data_size / 1024.**2)
    else:
        data_size_string = '%.2f gigabytes' % (data_size / 1024.**3)
    
    print "Loading %s into memory" % data_size_string
    
    
    if data_size > general_dict['warn_thresh'] and warn_size:
        print "WARNING " * 3
        print "Requested data is above specified warning threshold of %.2f MB" % (general_dict['warn_thresh'] / 1024**2)
        if not raw_input('Continue [Y/N]').lower().startswith('y'):
            return
    
    
    #reading metrics
    metric = {}
    inverse_metric = {}
    for key, value in config_dict.items():
        if key.startswith('metric_'):
            metric[key[7:]] = metric_dict[value]
            inverse_metric[key[7:]] = inverse_metric_dict[value]
        
    
    wave = np.fromfile(config_dict['datadir'] + 'wave.npmap')
    fluxes = load_spectra(config_dict, fnames, wave, **kwargs)
    
    
    return specgrid.specgrid(requested_params, fluxes, wave,
                             requested_param_names, normalizer=normalizer, convolver=convolver,
                             metric=metric, inverse_metric=inverse_metric)

def load_spectra(config_dict, fnames, wave, **kwargs):
    
    normalizer = kwargs.get('normalizer', None)
    convolver = kwargs.get('convolver', None)
    
    if normalizer != None:
        normalizer.wave = wave
        normalizer.calculate_idx()
    
    if convolver != None:
        convolver.wave = wave

    if config_dict['datatype'] == np.float64:
        pixel_per_spectrum = config_dict['specsize'] / 8
    else:
        raise ValueError("Datatype %s not implemented for spectra" % config_dict['datatype'])

    #initializing spectral grid
    fluxes = np.empty((len(fnames), pixel_per_spectrum))
    
    
    for i, fname in enumerate(fnames):
        flux = np.fromfile(config_dict['datadir'] + fname, dtype = config_dict['datatype'])
        
        if convolver != None:
            flux = convolver.convolve_grid(flux)
        if normalizer != None:
            flux = normalizer.normalize_grid(flux)
        
        
            
        fluxes[i] = flux
    return fluxes
    
    

class NormRange(object):
    def __init__(self, norm_range, wave=None):
        self.norm_range = norm_range
        self.wave = wave
        if self.wave != None:
            self.calculate_idx()
        
        
        
    def calculate_idx(self):
        self.min_idx = self.wave.searchsorted(self.norm_range[0])
        self.max_idx = self.wave.searchsorted(self.norm_range[1])
    
    
    def normalize_grid(self, flux):
        norm_factor = np.mean(flux[self.min_idx:self.max_idx])
        return flux / norm_factor
    
    def normalize_spectrum(self, spectrum):
        norm_factor = spectrum[slice(*self.norm_range)].flux.mean()
        return spectrum / norm_factor
    
class ConvolveResolution(object):
    def __init__(self, requested_resolution, initial_resolution=np.inf, wave=None):
        self.requested_resolution = requested_resolution
        self.initial_resolution = initial_resolution
    
    def set_wave(self, wave):
        self.wave = wave
    
    def convolve_grid(self, flux):
        tmp_spec = oned.onedspec(self.wave, flux, mode='waveflux')
        convolved_spec = tmp_spec.convolve_profile(self.requested_resolution, self.initial_resolution)
        return convolved_spec.flux