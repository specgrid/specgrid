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

from . import fix_spectrum1d
from fix_spectrum1d import Spectrum1D


class BaseSpectralGrid(object):
    """
    A BaseSpectralGrid interpolation class. Can serve as a model maybe

    Parameters
    ----------

    grid_hdf5_fname: filename for HDF5 File
    """

    parameters = None

    def __init__(self, grid_hdf5_fname,
                 interpolator=interpolate.LinearNDInterpolator):
        super(BaseSpectralGrid, self).__init__()
        if not os.path.exists(grid_hdf5_fname):
            raise IOError('{0} does not exists'.format(grid_hdf5_fname))

        self.grid_hdf5_fname = grid_hdf5_fname
        self.index = pd.read_hdf(self.grid_hdf5_fname, 'index')

        self.parameters = self.index.columns
        
        for parameter_name in self.parameters:
            setattr(self, parameter_name, None)

        with h5py.File(self.grid_hdf5_fname, 'r') as h5file:
            wavelength_unit = u.Unit(h5file['fluxes'].attrs['wavelength.unit'])
            self.wavelength = h5file['fluxes'].attrs['wavelength'] * \
                              wavelength_unit
            self.flux_unit = u.Unit(h5file['fluxes'].attrs['flux.unit'])
            self.fluxes = np.array(h5file['fluxes'])

        self.interpolate_grid = interpolator(self.index.values, self.fluxes)

    def __call__(self):
        return self.eval(*[getattr(self, item) for item in self.parameters])


    def eval(self, *args):
        """
        Interpolating on the grid to the necessary parameters

        Parameters
        ----------

        teff: float
            effective temperature
        logg: float
            base ten logarithm of surface gravity in cgs
        feh: float
            [Fe/H]

        """
        flux = self.interpolate_grid(*args)
        return Spectrum1D.from_array(self.wavelength.value,
                                     flux,
                                     dispersion_unit=self.wavelength.unit,
                                     unit=self.flux_unit)

class MunariGrid(BaseSpectralGrid):

    teff = 5780.
    logg = 4.4
    feh = 0.0

    parameters = ['teff', 'logg', 'feh']

    def __init__(self, grid_hdf5_fname,
                 interpolator = interpolate.LinearNDInterpolator,
                 query_string='teff > 0'):

        super(MunariGrid, self).__init__(grid_hdf5_fname, interpolator)
        self.load_dataset(query_string)


    def __call__(self):
        return self.eval(self.teff, self.logg, self.feh)


    def eval(self, teff, logg, feh):
        """
        Interpolating on the grid to the necessary parameters

        Parameters
        ----------

        teff: float
            effective temperature
        logg: float
            base ten logarithm of surface gravity in cgs
        feh: float
            [Fe/H]

        """
        flux = self.interpolate_grid(teff, logg, feh)
        return Spectrum1D.from_array(self.wavelength,
                                     flux,
                                     dispersion_unit=u.angstrom,
                                     unit=u.Unit('erg/ (cm2 s Angstrom)'))
