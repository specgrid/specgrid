import pandas as pd
import os
import sqlite3
import ConfigParser


def read_spectra(gridname, index):
    grid_config = ConfigParser.ConfigParser()
    grid_config.read(os.path.join(gridname, 'config.ini'))

    dtype = grid_config.get('structure', 'datatype')
    spec0 = np.fromfile(os.path.join(gridname, 'data', index['fname'].loc[0]), dtype=dtype)
    no_of_spectra = len(index)
    spectra = np.empty((len(index), len(spec0)))

    for i, fname in enumerate(index['fname']):
        fname = os.path.join(gridname, 'data', fname)
        print "[{}/{}] Reading {}".format(i, no_of_spectra, fname)
        spectra[i] = np.fromfile(fname, dtype)

    spectra_df = pd.DataFrame(spectra)
    spectra_df.index = index.index

    wave = np.fromfile(os.path.join(gridname, 'wave.npmap'), dtype)

    return wave, spectra_df


def convert_to_hdf5(grid_fname, gridname):
    index = read_index2df(gridname)
    wave, spectra = read_spectra(gridname, index)

    store = pd.HDFStore(grid_fname)
    store.put('index', index, format='table')
    store.put('spectra', spectra, format='fixed')
    store.get_storer('spectra').attrs.wave = wave

    store.close()

def read_indexdb_to_dataframe(gridname, sql_stmt, index_fname='index.db3'):
    """
    Reading the index db to a pandas DataFrame

    :param gridname:
    :param sql_stmt:
    :param index_fname:
    :return:
    """
    indexdb_path = os.path.join(gridname, index_fname)
    if not os.path.exists(indexdb_path):
        raise IOError('The index database {0} does not exist'.format(
            indexdb_path))

    conn = sqlite3.connect(index_fname)
    index = pd.read_sql(sql_read_stmt, conn)
    index['fname'] = index['fname'].astype(str)
    index['odftype'] = index['odftype'].astype(str)
    conn.close()
    return index


