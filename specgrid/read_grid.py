#Reading Grid into memory and grid information

import os
import ConfigParser
import numpy as np
import sqlite3
import numpy as np
import specgrid

try:
    import sqlparse
    sqlparse_available = True
except ImportError:
    sqlparse_available = False
    
def read_grid_config(fname='~/.specgrid/specgrid.ini'):
    fname = os.path.expanduser(fname)
    if not os.path.exists(fname):
        raise IOError('Config file does not exist at %s' % fname)
    conf = ConfigParser.ConfigParser()
    conf.read(fname)
    grid_names = conf.sections()
    grid_names.remove('general')
    
    #print "Found the following grids:\n %s" % ('\n'.join(grid_names))
    
    grid_dict = {}
    
    for grid in grid_names:
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
                value = np.float64(value)    
            grid_dict[grid][item] = value
    
    general_dict = {}
    
    for item, value in conf.items('general'):
        if item == 'warn_thresh':
            value = int(value) * 1024**2
        general_dict[item] = value
    return grid_dict, general_dict


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
        Threshold can be set in specgrid.ini"""
    
    #Readint grid configuration
    
    warn_size = kwargs.get('warn_size', True)
    grid_dict = kwargs.get('grid_dict', None)
    general_dict = kwargs.get('general_dict', None)
    
    if grid_dict == None or general_dict == None:
        grid_dict, general_dict = read_grid_config()
    
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
        
        if param_name in kwargs.keys():
            requested_param_names.append(param_name)
            param_value = kwargs[param_name]
        else:
            param_value = config_dict.get('default_' + param_name, None)
        print param_name, param_value
        param_values.append(param_value)
        
    #building sqlite query
    
    grid_query_conditions = []

    #now add the limitations
    
    for param_name, param_value in zip(param_names, param_values):
        #specifically requesting a None selection on a parameter
        
        if param_value == None:
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
                    
            if param_value.count(None) > 1:
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
    
    
    fluxes = load_spectra(config_dict, fnames, **kwargs)
    
    return specgrid.specgrid(requested_params, fluxes, None, None)

def load_spectra(config_dict, fnames, **kwargs):
    
    if config_dict['datatype'] == np.float64:
        pixel_per_spectrum = config_dict['specsize'] / 8
    else:
        raise ValueError("Datatype %s not implemented for spectra" % config_dict['datatype'])

    #initializing spectral grid
    fluxes = np.zeros((len(fnames), pixel_per_spectrum))
    
    
    for i, fname in enumerate(fnames):
        fluxes[i] = np.fromfile(config_dict['datadir'] + fname, dtype = config_dict['datatype'])
    return fluxes
    
    
    
    ###### DELETE FROM HEREON DOWN
    #reading config
    waveFName = config.get('wave', 'wave')
    specDType = config.get('structure', 'datatype')
    specSize = config.getint('structure', 'specsize')
    specsize = config.getint('wave', 'islog')
    R = config.getfloat('params', 'r')
    
    if kwargs.has_key('smoothres'):
        if kwargs['smoothres']>R:
            raise ValueError('requested resolution (R=%f) higher than'
                             'models intrinsic resolution (R=%f)' % (kwargs['smoothres'], R))
    #Reading wave solution
    wave = np.fromfile(os.path.join(specDataDir, waveFName))
    
    if kwargs.has_key('wave'):
        gridSize = len(fnames) * len(kwargs['wave'].tostring()) / 1024**2
    else:
        gridSize = len(fnames) * specSize / 1024**2
    print "Processing %d spectra" % len(fnames)
    
    print "Processing %.3f MB in to Memory." % gridSize
    
    specs = []

    startTime = time.time()
    
    for i, specFName in enumerate(fnames):
        flux = np.fromfile(os.path.join(specDataDir, 'data', specFName))
        spec = oned.onedspec(wave, flux, mode='waveflux')
        if i%100 == 0:
            print "@%d took %.2f s" % (i, time.time() - startTime)
            startTime = time.time()
        
            
        if kwargs.has_key('normrange'):
            normrange = kwargs['normrange']
            normFac = np.mean(spec[slice(*normrange)].flux)
            spec.flux /= normFac
            
            
            
        if kwargs.has_key('smoothres') or kwargs.has_key('smoothrot'):
            if kwargs.has_key('wave'):
                tmpSpec = spec[float(kwargs['wave'].min()):float(kwargs['wave'].max())]
                logDelta, logSpec = tmpSpec.interpolate_log()
            else:
                logDelta, logSpec = spec.interpolate_log()
            if kwargs.has_key('smoothres'):
                logSpec = logSpec.convolve_profile(kwargs['smoothres'], smallDelta=logDelta)
            if kwargs.has_key('smoothrot'):
                logSpec = logSpec.convolve_rotation(kwargs['smoothrot'], smallDelta=logDelta)
            spec = logSpec
            
        if kwargs.has_key('wave'):
            spec = spec.interpolate(kwargs['wave'])
        else:
            spec = spec.interpolate(wave)

            
        specs.append(spec.flux)
    if kwargs.has_key('wave'):
        return kwargs['wave'], np.array(specs)
    else:
        return wave, np.array(specs)