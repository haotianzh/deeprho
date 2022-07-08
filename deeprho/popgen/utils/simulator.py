from ..base import Replicate
import msprime as msp
import numpy as np
import warnings


class ExpLogGenerator(object):
    """ A random float generator for uncertain rates in simulation. """

    def __init__(self, start, end):
        self.start = start
        self.end = end

    def __call__(self):
        return np.exp(np.random.uniform(np.log(self.start), np.log(self.end)))

    def __str__(self):
        return "ExpLogGenerator(%.4e, %.4e)" % (self.start, self.end)

    def __repr__(self):
        return "ExpLogGenerator(%.4e, %.4e)" % (self.start, self.end)


class Simulator(object):
    """
    Using for creating simulated data.
    calling method:
        nsam: number of samples
        nreps: number of replicates
    configs:
        rate: mutation rate (per base)
        recombination_rate: recombination rate (per base)
        sequence_length: length of simulated genome (bp)
        population_size(or demography): either one of them should be specified
        ploidy: haploid or diploid
    examples:
        sim_configs = {'rate': 1e-8,
                    'recombination_rate': 1e-8,
                    'sequence_length': 1e5,
                    'ploidy': 1,
                    'population_size': 1e5}
        simulator = Simulator(sim_configs)
        # simulate a genome with 10 individuals
        for genome in simulator(10, 1):
            pass
    """

    def __init__(self, configs=None):
        self.__simulator__ = '{name}/{version}'.format(name='msprime', version=msp.__version__)
        self.configs = None
        self._mutation_configs = None
        self._ancestry_configs = None
        if configs is not None:
            self.set_configs(configs)

    def __call__(self, nsam, nreps):
        """ Run simulation for nsam samples and nreps replicates """
        configs = self.configs
        if configs is None:
            raise Exception("no configuration for simulation found.")
        replicates = msp.sim_ancestry(samples=nsam, num_replicates=nreps, **self._ancestry_configs)
        # simulate mutations once for each of genealogical trees.
        for ts in replicates:
            configs = {'rate': self._mutation_configs['rate']()}
            mts = msp.sim_mutations(ts, model=msp.InfiniteSites(), **configs)
            configs.update(self._ancestry_configs)
            mts.__setattr__('rr', configs['recombination_rate'])
            mts.__setattr__('mr', configs['rate'])
            rep = Replicate(mts, configs)
            yield rep

    def update(self, u):
        assert isinstance(u, dict), Exception('key-value pairs must be provided.')
        configs = self.configs.copy()
        configs.update(u)
        self.set_configs(configs)

    def set_configs(self, configs):
        requires = ['sequence_length']
        if configs is None:
            raise Exception("configuration cannot be empty.")
        for key in requires:
            if key not in configs:
                raise Exception("%s must be set." % key)
        if 'rate' in configs:
            assert isinstance(configs['rate'], float) or isinstance(configs['rate'], list) or isinstance(configs['rate'], ExpLogGenerator), \
                Exception("mutation rate should be either a float number or an interval.")
            if isinstance(configs['rate'], ExpLogGenerator):
                pass
            elif isinstance(configs['rate'], float):
                self._mutation_configs = {'rate': ExpLogGenerator(configs['rate'], configs['rate'])}
            else:
                assert len(configs['rate']) == 2, Exception("length of list must be exactly 2.")
                self._mutation_configs = {'rate': ExpLogGenerator(configs['rate'][0], configs['rate'][1])}
        else:
            self._mutation_configs = {'rate': ExpLogGenerator(1e-8, 1e-8)}
        if 'population_size' in configs and 'demography' in configs:
            warnings.warn("both population size and demography detected, population size will be removed.")
            del(configs['population_size'])
        self._ancestry_configs = {'recombination_rate': 1e-8, 'ploidy': 1}
        self._ancestry_configs.update(configs)
        if 'rate' in self._ancestry_configs:
            del(self._ancestry_configs['rate'])
        if self.configs is None:
            self.configs = dict()
        self.configs.update(self._mutation_configs)
        self.configs.update(self._ancestry_configs)

    def __str__(self):
        return self.__simulator__
