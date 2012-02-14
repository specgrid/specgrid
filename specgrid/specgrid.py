#specgrid class
from scipy import interpolate


class specgrid(object):
    def __init__(self, params, fluxes, wave, param_names,
                 interpolator = interpolate.LinearNDInterpolator):
        self.wave = wave
        self.interpolated_grid = interpolator(params, fluxes)
        