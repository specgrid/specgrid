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

from specutils import Spectrum1D


class SpectralGrid(object):
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

        #for param_name in self.grid_param_names:
        #    self.__setattr__(param_name, modeling.Parameter(param_name, default=0))


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

class MunariGrid(SpectralGrid):

    teff = 5780.
    logg = 4.4
    feh = 0.0

    def __init__(self, grid_hdf5_fname,
                 interpolator = interpolate.LinearNDInterpolator,
                 query_string='teff > 0'):

        super(MunariGrid, self).__init__(grid_hdf5_fname, interpolator)
        self.load_dataset(query_string)


    def __call__(self):
        flux = self.interpolate_grid(self.teff, self.logg, self.feh)
        return Spectrum1D.from_array(self.wave, flux * u.Unit(1),
                                     dispersion_unit=u.angstrom)

