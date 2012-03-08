#module containing code to fit grid based
import numpy as np

def interp_specgrid(specgrid, wave):
    new_fluxes = np.empty((specgrid.fluxes.shape[0], wave.shape[0]))
    for i in xrange(specgrid.fluxes.shape[0]):
        new_fluxes[i] = np.interp(wave, specgrid.wave, specgrid.fluxes[i])
    
    return new_fluxes
    
def fit_grid_simple(specgrid, spectrum, error=None, interp_wave=False):
    if interp_wave:
        fluxes = interp_specgrid(specgrid, spectrum.wave)
    else:
        fluxes = specgrid.fluxes
    
    if specgrid.normalizer != None:
        spectrum = specgrid.normalizer.normalize_spectrum(spectrum)

    return spectrum

    chi2 = calc_chi2(fluxes, spectrum)
    min_param = specgrid.params[np.argmin(chi2)]
    
    return min_param
    
def calc_chi2(model_flux, spectrum):
    error = spectrum.error
    mask = spectrum.mask

    if mask != None:
        spectrum_flux = np.ma.MaskedArray(spectrum.flux, mask)
    else:
        print "no mask"
        spectrum_flux = spectrum.flux

    if error == None:
        chi2 = np.sum((model_flux - spectrum_flux)**2, axis=spectrum.flux.ndim - 1)
    else:
        chi2 = np.sum(((model_flux - spectrum_flux) / error)**2, axis=spectrum.flux.ndim - 1)
    
    return chi2


class minuit_function(object):
    def __init__(self, specgrid, samplespec, priors=None):
        self.specgrid = specgrid

        if not all(samplespec.wave == self.specgrid.wave):
            self.samplespec = samplespec.interpolate(specgrid.wave)
        else:
            self.samplespec = samplespec
            
        if specgrid.normalizer != None:
            self.samplespec = specgrid.normalizer.normalize_spectrum(self.samplespec)
        
        self.priors = priors
        self. param_names = self.specgrid.param_names + self.specgrid.plugins.keys()
        self.varnames(*self.param_names)
        
    class func_code:
        co_varnames = []
        co_argcount = 0

    def varnames(self, *args):
        self.func_code.co_varnames = args
        self.func_code.co_argcount = len(self.func_code.co_varnames)

    def __call__(self, *args):
        if len(args) != self.func_code.co_argcount:
            raise TypeError, "Function takes %d arguments (%d given)" % (self.func_code.co_argcount, len(args))
        
        varnames = self.func_code.co_varnames
        samplespec = self.samplespec
        priors = self.priors
        
        
        
        model_flux = self.specgrid.interpolate(**dict(zip(varnames, args)))
        
        chi2 = calc_chi2(model_flux, self.samplespec)
        return chi2