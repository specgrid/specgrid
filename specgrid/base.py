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

from fix_spectrum1d import Spectrum1D


class SpectralGrid(object):
    """
    A SpectralGrid interpolation class. Can serve as a model maybe

    Parameters
    ----------

    grid_hdf5_fname: filename for HDF5 File
    """

    param_names = None

    def __init__(self, grid_hdf5_fname,
                 interpolator=interpolate.LinearNDInterpolator):

        super(SpectralGrid, self).__init__()

        if not os.path.exists(grid_hdf5_fname):
            raise IOError('{0} does not exists'.format(grid_hdf5_fname))

        self._load_index(grid_hdf5_fname)
        self._load_fluxes(grid_hdf5_fname)


        self.interpolate_grid = interpolator(self.index.values, self.fluxes)
        self.interpolator = interpolator

    def __call__(self):
        return Spectrum1D.from_array(self.wavelength, self._interpolate_flux())

    def _interpolate_flux(self):
        parameter_values = [getattr(self, item) for item in self.param_names]
        return self.interpolate_grid(parameter_values)[0] * self.flux_unit

    def _load_index(self, grid_hdf5_fname):
        """
        Loading the index from the hdf5 file

        Parameters
        ----------

        grid_hdf5_fname: ~str
            path to HDF5 file
        """

        self.index = pd.read_hdf(grid_hdf5_fname, 'index')

        self.param_names = self.index.columns.tolist()

        for parameter_name in self.param_names:
            setattr(self, parameter_name, self.index[parameter_name].iloc[0])

    def _load_fluxes(self, grid_hdf5_fname):
        """
        Loading the fluxes from the HDF5 file

        Parameters
        ----------

        grid_hdf5_fname: ~str
            path to HDF5 file

        """

        with h5py.File(grid_hdf5_fname, 'r') as h5file:
            wavelength_unit = u.Unit(h5file['fluxes'].attrs['wavelength.unit'])
            self.wavelength = h5file['fluxes'].attrs['wavelength'] * \
                              wavelength_unit
            self.flux_unit = u.Unit(h5file['fluxes'].attrs['flux.unit'])
            self.fluxes = np.array(h5file['fluxes'])

    def evaluate(self, *args, **kwargs):
        """
        Interpolating on the grid to the necessary param_names

        Examples
        --------

        This can either be called with arguments ``specgrid.evaluate(5780, 4.4, -1)`` or
        using keyword way of calling (then not all param_names have to be given)
        ``specgrid.evaluate(logg=4.4)``
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
                                 'spectral grid (param_names are {1})'.format(
                    key, ','.join(self.param_names)))
            setattr(self, key, kwargs[key])

        return self.__call__()

