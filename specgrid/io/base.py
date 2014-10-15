import pandas as pd
import os
import sqlite3
import ConfigParser
import numpy as np
from astropy import units as u, constants as const

try:
    import h5py
except ImportError:
    h5py_available = False
else:
    h5py_available = True




def read_indexdb_to_dataframe(path_to_grid, sql_stmt, indexdb_fname='index.db3'):
    """
    Reading the index db to a pandas DataFrame

    :param path_to_grid:
    :param sql_stmt:
    :param indexdb_fname:
    :return:
    """
    indexdb_path = os.path.join(path_to_grid, indexdb_fname)

    print "Reading database index from {0}".format(indexdb_path)
    if not os.path.exists(indexdb_path):
        raise IOError('The index database {0} does not exist'.format(
            indexdb_path))

    conn = sqlite3.connect(indexdb_path)
    indexdf = pd.read_sql(sql_stmt, conn)
    if 'fname' not in indexdf.columns:
        raise ValueError('Query must contain the fname of the spectra for '
                         'further processing')

    for column in indexdf:
        if indexdf[column].dtype == 'O':
            indexdf[column] = indexdf[column].astype(str)
    conn.close()

    return indexdf


def read_spectra(gridname, indexdf, process_plugins=[]):
    grid_config = ConfigParser.ConfigParser()
    grid_config.read(os.path.join(gridname, 'config.ini'))

    dtype = grid_config.get('structure', 'datatype')
    wave_unit = u.Unit(grid_config.get('wave', 'wave_unit'))
    flux_unit = u.Unit(grid_config.get('flux', 'flux_unit'))

    spec0 = np.fromfile(os.path.join(gridname, 'data', indexdf['fname'].loc[0]),
                        dtype=dtype)

    no_of_spectra = len(indexdf)
    fluxes = np.empty((len(indexdf), len(spec0)))

    for i, fname in enumerate(indexdf['fname']):
        fname = os.path.join(gridname, 'data', fname)
        print "[{}/{}] Reading {}".format(i, no_of_spectra, fname)
        fluxes[i] = np.fromfile(fname, dtype)


    wave = np.fromfile(os.path.join(gridname, 'wave.npmap'), dtype) * wave_unit


    return wave, fluxes, flux_unit



def make_hdf5(gridname, sql_stmt, h5_fname, ignore_columns=[]):
    """
    Making an HDF5 File from a databased grid

    Parameters
    ----------

    gridname: ~str
        name of grid
    sql_stmt: ~str
        sql statement to query the index with
    h5_fname: ~str
        path to save HDF5 grid to
    """

    index_df = read_indexdb_to_dataframe(gridname, sql_stmt)
    wavelength, fluxes, flux_unit = read_spectra(gridname, index_df)

    data_columns = []
    for column in index_df.columns:
        if column == 'fname' or column in ignore_columns:
            continue
        if len(index_df[column].unique()) > 1:
            data_columns.append(column)

    index_df[data_columns].to_hdf(h5_fname, 'index', mode='w')

    with h5py.File(h5_fname, 'a') as fh:
        fh['fluxes'] = fluxes
        fh['fluxes'].attrs['wavelength'] = wavelength.value
        fh['fluxes'].attrs['wavelength.unit'] = str(wavelength.unit)
        fh['fluxes'].attrs['flux.unit'] = 'erg / (cm^2 s Angstrom)'
