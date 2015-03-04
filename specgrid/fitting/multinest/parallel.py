import os

from logging import getLogger
from copy import deepcopy


from specgrid.model_star import Observation
from specgrid.fitting.multinest.base import fit_simple_spectrum_multinest

from IPython.parallel import interactive

logger = getLogger(__name__)

def set_engines_cpu_affinity():
    import sys
    if sys.platform.startswith('linux'):
        try:
            import psutil
        except ImportError:
            logger.warning('psutil not available - can not set CPU affinity')
        else:
            from multiprocessing import cpu_count
            p = psutil.Process(os.getpid())
            p.set_cpu_affinity(range(cpu_count()))




@interactive
def parallel_simple_spectrum_multinest(spectrum, spectrum_label,
                                       model_instrument, priors,
                                       data_path=None):
    import os
    import pymultinest
    from specgrid.model_star import Observation
    from specgrid.fitting.multinest.base import fit_simple_spectrum_multinest
    fname = '{0}_simple.h5'.format(spectrum_label)
    if data_path is not None:
        fname = os.path.join(data_path, fname)
    if os.path.exists(fname):
        return True
    model_observation = Observation(model_star, model_instrument)
    mn_fit = fit_simple_spectrum_multinest(spectrum, model_observation,
                                           priors)

    mn_fit.result.posterior_data.to_hdf(fname, 'multinest', mode='w')
    return True

@interactive
def parallel_spectro_photometry_multinest(spectrum, spectrum_label, model_instrument, priors,  magnitude_set_dill, data_path=None):
    import os
    import pymultinest
    import dill
    from specgrid.model_star import Observation
    from specgrid.fitting.multinest.base import fit_spectro_photometry_multinest
    fname = '{0}_spectro_photom.h5'.format(spectrum_label)
    magnitude_set = dill.loads(magnitude_set_dill)
    if data_path is not None:
        fname = os.path.join(data_path, fname)
    if os.path.exists(fname):
        return True
    model_observation = Observation(model_star, model_instrument)

    mn_fit = fit_spectro_photometry_multinest(spectrum, magnitude_set,
                                              model_observation, priors)

    mn_fit.result.posterior_data.to_hdf(fname, 'multinest', mode='w')
    return True


class ParallelMultiNestFitter(object):

    def __init__(self, remote_clients, model_star, model_instrument, fitness_function):
        self.model_star = model_star
        self.model_instrument = model_instrument

        self.remote_clients = remote_clients
        self.prepare_remote_clients(remote_clients, model_star)
        self.fitness_function = fitness_function

        self.lbv = remote_clients.load_balanced_view()


    @staticmethod
    def prepare_remote_clients(clients, model_star):
        """
        Preparing the remote clients for computation: Uploading the atomic
        data if available and making sure that the clients can run on different
        CPUs on each Node

        Parameters
        ----------

        clients: IPython.parallel.Client
            remote clients from ipython

        model_star: tardis.atomic.AtomData or None
            remote atomic data, if None each queue needs to bring their own one
        """

        logger.info('Sending model star to remote '
                    'clients and importing specgrid')
        clients.block = True
        for i, client in enumerate(clients):
            logger.info('At client {0} - '.format(i))
            client['model_star'] = model_star
            client.execute('')

        clients.block = False


        for client in clients:
            client.apply(set_engines_cpu_affinity)

    def queue_spectrum(self, spectrum, spectrum_label, priors, data_path=None, **kwargs):
            model_instrument = deepcopy(self.model_instrument)
            for model in model_instrument.models:
                if model.__class__.__name__ == 'NormalizeParts':
                    model._update_observed_spectrum(spectrum, kwargs.pop('parts'))
                elif hasattr(model, '_update_observed_spectrum'):
                    model._update_observed_spectrum(spectrum)
                else:
                    continue
            return self.lbv.apply(self.fitness_function, spectrum, spectrum_label, model_instrument, priors, data_path=data_path, **kwargs)

