import numpy as np
import tskit


class Haplotype(object):
    """
    Class for haplotype
    """

    def __init__(self, ts=None, positions=None, matrix=None, ancestral_state=0):
        if ts is not None:
            assert isinstance(ts, tskit.TreeSequence), Exception("a tskit.TreeSequence object must be provided.")
            self.matrix = ts.genotype_matrix().T
            self.positions = [int(v.position) for v in ts.variants()]
        elif matrix is not None and positions is not None:
            self.matrix = matrix
            self.positions = positions
        self.ancestral_state = ancestral_state
        length = self.positions[-1] - self.positions[0]
        self.scaled_positions = ((np.array(self.positions) - self.positions[0]) / length * 10000 + 1)

    def cut(self):
        # using tskit.TreeSequence.keep_intervals() to cut genome into fixed-length fragments.
        pass

    @property
    def nsites(self):
        if self.matrix is not None:
            return self.matrix.shape[1]
        else:
            return 0

    @property
    def nsamples(self):
        if self.matrix is not None:
            return self.matrix.shape[0]
        else:
            return 0

    def converse(self):
        # converts ancestral state from 0 to 1.
        pass

    def from_tskit(self, ts):
        # converts from tskit format.
        pass
