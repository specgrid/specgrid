from specgrid.base import SpectralGrid
import numpy as np
from scipy.optimize import minimize
class TestInterpolation(SpectralGrid):
    pass

    def montecarlo_interpolation_test(self, points_idx=None):
        if points_idx is None:
            points_idx = np.arange(len(self.fluxes))

        all_points_idx = np.arange(len(self.fluxes))

        miss_interpolate_flux = []
        actual_flux = []

        fitted_values = []


        for i, point_id in enumerate(points_idx):
            print "({0}/{1})".format(i, len(points_idx))
            missing_point = self.index.values[point_id]
            exclusion_mask = all_points_idx != point_id
            new_index = self.index.values[exclusion_mask]
            new_index[:,1] = 10** new_index[:,1]
            new_index[:,2] = 10** new_index[:,2]

            new_fluxes = self.fluxes[exclusion_mask]
            self.interpolate_grid = self.interpolator(new_index, new_fluxes)
            miss_interpolate_flux.append(
                self.interpolate_grid(missing_point)[0])
            actual_flux.append(self.fluxes[point_id])

            fitting_func = TestingFittingFunction(self.interpolate_grid, self.fluxes[point_id])
            fit = minimize(fitting_func, tuple(missing_point.tolist()), method='Nelder-Mead')
            fitted_values.append(fit['x'])

        return self.index.values[points_idx], np.array(fitted_values)



class TestingFittingFunction(object):

    def __init__(self, interpolator, real_flux):
        self.interpolator = interpolator
        self.real_flux = real_flux

    def __call__(self, *args):
        interp_point = args[0].copy()
        interp_point[1] = 10**interp_point[1]
        interp_point[2] = 10**interp_point[2]
        return ((self.interpolator(interp_point)[0] - self.real_flux)**2).sum()