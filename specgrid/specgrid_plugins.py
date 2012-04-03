#module implementing the specgrid plugin system. classes that

from scipy import ndimage
import numpy as np
try:
    import pysynphot
    pysynphot_available = True
except ImportError:
    pysynphot_available = False

cardelli_data = '/Users/wkerzend/scripts/python/specgrid/specgrid/data/reddening/milkyway_rv3.1.dat'
class GaussConvolvePlugin(object):
    def __init__(self, wave, **kwargs):
        self.wave = wave
        self.param_name = 'gauss'
        self.kwargs = kwargs
    def __call__(self, flux, sigma):
        return ndimage.gaussian_filter1d(flux, sigma, **self.kwargs)
        



class CardelliReddenPlugin(object):
    def __init__(self, wave, **kwargs):
        self.wave = wave
        self.red_wave, self.red_avscaled = np.loadtxt(cardelli_data, unpack=1)
        self.param_name = 'ebv'
        self.kwargs = kwargs
        self.normalizer = self.kwargs.get('normalizer', None)
        self.interp_avscaled = np.interp(self.wave, self.red_wave, self.red_avscaled)
    def __call__(self, flux, ebv):
        extinct_flux = 10**(-0.4 * ebv * self.interp_avscaled) * flux
        if self.normalizer == None:
            return extinct_flux
        else:
            return self.normalizer.normalize_grid(extinct_flux)
    
        
class reddenGrid(object):
    def __init__(self, wave, enableDIB=False, enableFlux=True, extinctionLaw = 'gal3', normRange=None):
        self.points = np.array(ebvRange)
        values = []
        for ebv in ebvRange:
            extinctFlux = np.ones(wave.shape)
            if enableFlux:
                extinct = pysynphot.Extinction(ebv, extinctionLaw)    
                extinctThroughPut = extinct.GetThroughput()[::-1]
                f = interpolate.interp1d(extinct.wave[::-1], extinctThroughPut)
                extinctFlux *= f(wave)
            
            if enableDIB:
                extinctFlux *= pydib.makeSpectrum(wave, ebv).flux
            if normRange != None:
                extinctFluxSpec = oned.onedspec(wave, extinctFlux, mode='waveflux')
                extinctFlux /= np.mean(extinctFluxSpec[slice(*normRange)].flux)
            values.append(extinctFlux)
        self.values = np.array(values)
        self.interpGrid = interpolate.interp1d(self.points, self.values.transpose(), fill_value=1, bounds_error=False, kind='cubic')
    
    def __call__(self, *args):
        return self.interpGrid(args)[:,0]        