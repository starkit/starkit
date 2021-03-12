import logging
import os
from IPython.parallel.util import interactive

logger = logging.getLogger(__name__)

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
            p.set_cpu_affinity(list(range(cpu_count())))


from starkit.fitkit.multinest.base import MultiNest

@interactive
def simple_worker(fitter_factory):
    """
    This is a starkit worker

    Parameters
    ----------

    fitter_factory:

    """
    return fitter_factory.run_fit(spectral_grid)



class BaseFitterFactory(object):
    def __init__(self, model):
         self.model = model

class MultiNestFactory(BaseFitterFactory):
    def __init__(self, model, priors, multinest_kwargs={}):
        super(MultiNestFactory, self).__init__(model)
        self.priors = priors
        self.multinest_kwargs = multinest_kwargs

    def run_fit(self, spectral_grid):
        likelihood = spectral_grid | self.model
        mn = MultiNest(likelihood, self.priors)
        return mn.run(**self.multinest_kwargs)


class BaseGridLauncher(object):
    def __init__(self, remote_clients, spectral_grid_path):
        self.remote_clients = remote_clients
        self.prepare_remote_clients(remote_clients, spectral_grid_path)
        self.lbv = remote_clients.load_balanced_view()



    @staticmethod
    def prepare_remote_clients(clients, spectral_grid_path):
        """
        Preparing the remote clients for computation: Uploading the spectral grid
         if available and making sure that the clients can run on different
        CPUs on each Node

        Parameters
        ----------

        clients: IPython.parallel.Client
            remote clients from ipython

        atom_data: tardis.atomic.AtomData or None
            remote atomic data, if None each queue needs to bring their own one
        """

        logger.info('Loading the spectral grid on each engine')
        clients.block = True
        for client in clients:
            client['spectral_grid_path'] = spectral_grid_path
            client.execute('from starkit.gridkit import load_grid')
            client.execute('spectral_grid = load_grid(spectral_grid_path)')

        clients.block = False

        for client in clients:
            client.apply(set_engines_cpu_affinity)

    def queue_fit(self, fitter_factory):
        """
        Add single parameter set to the queue

        Parameters
        ----------

        parameter_set_dict: ~dict
            a valid configuration dictionary for TARDIS
        """

        return self.lbv.apply(simple_worker, fitter_factory)