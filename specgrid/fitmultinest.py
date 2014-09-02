from scipy.stats import norm, poisson
from collections import OrderedDict
from pylab import plt

class UniformPrior(object):
    """
    A Uniform distribution prior
    
    Parameters
    ----------
    
    lbound: ~float
        lower bound
    
    ubound: ~float
        upper bound
    """
    
    def __init__(self, lbound, ubound):
        self.lbound = lbound
        self.ubound = ubound
    
    def __call__(self, cube):
        return cube * (self.ubound - self.lbound) + self.lbound

class GaussianPrior(object):
    """
    A gaussian prior
    
    Parameters
    ----------
    
    m: ~float
        mean of the distribution
        
    sigma: ~float
        sigma of the distribution
    
    """
    
    def __init__(self, m, sigma):
        self.m = m
        self.sigma = sigma
        
    def __call__(self, cube):
        return norm.ppf(cube,scale=self.sigma,loc=self.m)
        

class PoissonPrior(object):
    """
    A Poisson prior
    
    Parameters
    ----------
    
    m: ~float
        mean of the distribution
        
    
    """
    def __init__(self, m):
        self.m = m
        
    def __call__(self,cube):
        return poisson.ppf(cube,loc=self.m)

class FixedPrior(object):
    """
    A fixed value
    
    Parameters
    ----------
    
    val: ~float
        fixed value
    """
    
    def __init__(self, val):
        self.val = val
    
    def __call__(self, cube):
        return self.val
    
#def multinest_spectrum(grid, spectrum, 
#priors={'teff':UniformPrior(4000, 9000), 'logg':GaussianPrior(3, 0.5), 'feh':FixedPrior(-0.5))

class SpectralData(object):
    pass

class Likelihood(object):
    """
    A spectrum likelihood object. Initialize with the spectral grid model.
    Calls will returns the log-likelihood for observing a spectrum given the model
    parameters and the uncertainty.
    
    Parameters
    ----------
    
    model_star: ~model from specpgrid
        model object containing the spectral grid to evaluate
        
    
    """    
    def __init__(self, model_star, parameter_names):
        self.model = model_star
        self.parameter_names = parameter_names
        
    def __call__(self, data, model_param):
        # returns the likelihood of observing the data given the model parameters
        param_dict = OrderedDict([(key, value) for key, value in zip(self.parameter_names, model_param)])
        
        m = self.model.evaluate(**param_dict)
        
        # log-likelhood for chi-square
        return (-0.5 * ((data.flux.value - m.flux.value) / 
        data.uncertainty.value)**2).sum()

class FitMultinest(object):
    """
    Fitting Multinest
    """
    

    def __init__(self,spectrum, priors, model_star, likelihood=None, 
    basename='chains/spectrumFit'):
        self.spectrum = spectrum    # Spectrum1D object with the data
        self.model = model_star     # grid of spectra from specgrid
        parameter_names = sorted(priors.keys(), key=lambda key: model_star.parameters.index(key))
        self.priors = PriorCollections(priors, 
                                parameter_order=model_star.parameters)  
                                # PriorCollection from the input priors
        self.n_params = len(priors)             # number of parameters to fit
        
        if likelihood is None:
            # use the default likelihood for a Spectrum1D object
            self.likelihood = Likelihood(model_star)
            
        # variables that will be filled in after the fit has been run
        self.fit_mean = None    # mean of fitted parameters
        self.sigma = None      # 1 sigma (68% credible intervales) for fitted parameters
        self.evidence = None   # the global evidence value for the best fit
        
    
    def run(self, no_plots=False,**kwargs):
        # runs pymultinest
        l = self.likelihood
        priors = self.priors
        pymultinest.run(l, priors.priorTransform, self.n_params, outputfiles_basename=self.basename,**kwargs)
        json.dump(parameters, open(basename+'params.json', 'w')) # save parameter names
        
        # analyze the output data
        a = pymultinest.Analyzer(outputfiles_basename=self.basename, n_params = self.n_params)
        s = a.get_stats()
        modes = s['modes'][0]   
        self.mean = modes['mean']
        self.sigma = modes['sigma']
        self.evidence = s['global evidence']
        
        if not(no_plots):
            self.mkplots()
        
        
    def mkplots(self):
        # run to make plots of the resulting posteriors. Modified from marginal_plots.py
        # from pymultinest. Produces basename+marg.pdf and basename+marge.png files
        prefix = self.basename
        
        parameters = json.load(file(prefix + 'params.json'))
        n_params = len(parameters)
        
        a = pymultinest.Analyzer(n_params = n_params, outputfiles_basename = prefix)
        s = a.get_stats()
        
        p = pymultinest.PlotMarginal(a)
        
        
        values = a.get_equal_weighted_posterior()
        assert n_params == len(s['marginals'])
        modes = s['modes']
        
        dim2 = os.environ.get('D', '1' if n_params > 20 else '2') == '2'
        nbins = 100 if n_params < 3 else 20
        if dim2:
        	plt.figure(figsize=(5*n_params, 5*n_params))
        	for i in range(n_params):
        		plt.subplot(n_params, n_params, i + 1)
        		plt.xlabel(parameters[i])
        	
        		m = s['marginals'][i]
        		plt.xlim(m['5sigma'])
        	
        		oldax = plt.gca()
        		x,w,patches = oldax.hist(values[:,i], bins=nbins, edgecolor='grey', color='grey', histtype='stepfilled', alpha=0.2)
        		oldax.set_ylim(0, x.max())
        	
        		newax = plt.gcf().add_axes(oldax.get_position(), sharex=oldax, frameon=False)
        		p.plot_marginal(i, ls='-', color='blue', linewidth=3)
        		newax.set_ylim(0, 1)
        	
        		ylim = newax.get_ylim()
        		y = ylim[0] + 0.05*(ylim[1] - ylim[0])
        		center = m['median']
        		low1, high1 = m['1sigma']
        		print center, low1, high1
        		newax.errorbar(x=center, y=y,
        			xerr=numpy.transpose([[center - low1, high1 - center]]), 
        			color='blue', linewidth=2, marker='s')
        		oldax.set_yticks([])
        		#newax.set_yticks([])
        		newax.set_ylabel("Probability")
        		ylim = oldax.get_ylim()
        		newax.set_xlim(m['5sigma'])
        		oldax.set_xlim(m['5sigma'])
        		#plt.close()
        	
        		for j in range(i):
        			plt.subplot(n_params, n_params, n_params * (j + 1) + i + 1)
        			p.plot_conditional(i, j, bins=20, cmap = plt.cm.gray_r)
        			for m in modes:
        				plt.errorbar(x=m['mean'][i], y=m['mean'][j], xerr=m['sigma'][i], yerr=m['sigma'][j])
        			plt.xlabel(parameters[i])
        			plt.ylabel(parameters[j])
        			plt.xlim([m['mean'][i]-5*m['sigma'][i],m['mean'][i]+5*m['sigma'][i]])
        			plt.ylim([m['mean'][j]-5*m['sigma'][j],m['mean'][j]+5*m['sigma'][j]])
        			#plt.savefig('cond_%s_%s.pdf' % (params[i], params[j]), bbox_tight=True)
        			#plt.close()
        
        	plt.savefig(prefix + 'marg.pdf')
        	plt.savefig(prefix + 'marg.png')
        	plt.close()
        else:
        	from matplotlib.backends.backend_pdf import PdfPages
        	print '1dimensional only. Set the D environment variable D=2 to force'
        	print '2d marginal plots.'
        	pp = PdfPages(prefix + 'marg1d.pdf')
        	
        	for i in range(n_params):
        		plt.figure(figsize=(5, 5))
        		plt.xlabel(parameters[i])
        		
        		m = s['marginals'][i]
        		plt.xlim(m['5sigma'])
        	
        		oldax = plt.gca()
        		x,w,patches = oldax.hist(values[:,i], bins=20, edgecolor='grey', color='grey', histtype='stepfilled', alpha=0.2)
        		oldax.set_ylim(0, x.max())
        	
        		newax = plt.gcf().add_axes(oldax.get_position(), sharex=oldax, frameon=False)
        		p.plot_marginal(i, ls='-', color='blue', linewidth=3)
        		newax.set_ylim(0, 1)
        	
        		ylim = newax.get_ylim()
        		y = ylim[0] + 0.05*(ylim[1] - ylim[0])
        		center = m['median']
        		low1, high1 = m['1sigma']
        		print center, low1, high1
        		newax.errorbar(x=center, y=y,
        			xerr=numpy.transpose([[center - low1, high1 - center]]), 
        			color='blue', linewidth=2, marker='s')
        		oldax.set_yticks([])
        		newax.set_ylabel("Probability")
        		ylim = oldax.get_ylim()
        		newax.set_xlim(m['5sigma'])
        		oldax.set_xlim(m['5sigma'])
        		plt.savefig(pp, format='pdf', bbox_inches='tight')
        		plt.close()
        	pp.close()
    
    
class PriorCollections(object):
    
    
    def __init__(self, prior_dict, parameter_order=[]):
        self.priors = prior_dict
        
    def prior_transform(self, cube, ndim, nparam):
        # will be given an array of values from 0 to 1 and transforms it 
        # according to the prior distribution
        
        for i in xrange(len(cube)):
            cube[i] = self.priors.values()[i](cube[i])

        
        