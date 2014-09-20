import pandas as pd
import os
import sqlite3
import numpy as np
import h5py
from astropy import units as u



sql_stmt = """
    select
        teff, logg, feh, odf_select(odftype, fname) as fname
    from
        munari
    where
        vrot=0 and alpha=0.0 and k=2.0
        and teff between 5000 and 6500
        and logg between 3 and 4.5
        and feh between -0.3 and 1
    group by
        teff, logg, feh
"""

sample_sql_template = """
    select
        teff, logg, feh, odf_select(odftype, fname) as fname
    from
        munari
    where
        vrot=0 and alpha=0.0 and k=2.0
        and teff between {teff_lower} and {teff_upper}
        and logg between {logg_lower} and {logg_upper}
        and feh between {feh_lower} and {feh_upper}
    group by
        teff, logg, feh
"""

class ODFSelect(object):
    def __init__(self):
        self.odf_type = []
        self.fname = []

    def step(self, odf_type, fname):
        self.odf_type.append(odf_type)
        self.fname.append(fname)

    def finalize(self):
        #print self.odf_type, self.fname

        if 'new' in self.odf_type:
            return self.fname[self.odf_type.index('new')]
        else:
            return self.fname[0]

def read_index2df(gridname, conn):

    sql_read_stmt = 'select id, teff, logg, feh, vrot, k, alpha, odftype, fname from munari where vrot=0'
    index = pd.read_sql(sql_read_stmt, conn)
    index['fname'] = index['fname'].astype(str)
    index['odftype'] = index['odftype'].astype(str)
    conn.close()
    return index



if __name__ == '__main__':
   gridname = 'munari'
   gridhdf5_fname = 'munari_small.h5'
   wavelength = np.fromfile(os.path.join(gridname, 'wave.npmap'), dtype=np.float64) * u.angstrom
   with sqlite3.connect(os.path.join(gridname, 'index.db3')) as conn:
        conn.create_aggregate('odf_select', 2, ODFSelect)
        index = pd.read_sql(sql_stmt, conn)
        index[['teff', 'logg', 'feh']].to_hdf(gridhdf5_fname, 'index', mode='w')
        flux0 = np.fromfile(os.path.join(gridname, 'data', index.fname[0]), dtype=np.float64)

        fluxes = np.empty((len(index.fname), len(flux0)))
        for i, fname in enumerate(index['fname']):
            print '{0} of {1}'.format(i, len(index.fname))
            fluxes[i] = np.fromfile(os.path.join(gridname, 'data', fname), dtype=np.float64)

        with h5py.File(gridhdf5_fname, 'a') as fh:
            fh['fluxes'] = fluxes
            fh['fluxes'].attrs['wavelength'] = wavelength.value
            fh['fluxes'].attrs['wavelength.unit'] = str(wavelength.unit)
            fh['fluxes'].attrs['flux.unit'] = 'erg / (cm^2 s Angstrom)'