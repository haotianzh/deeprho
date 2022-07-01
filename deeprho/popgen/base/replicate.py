from . import Node, BaseTree, Haplotype


class Replicate(object):

    def __init__(self, ts, configs):
        self.configs = configs
        self.ts = ts
        self.haplotype = Haplotype(ts)

    # def genealogies(self, branch_length=False):
    #     genealogies = []
    #     for tree in self.ts.trees():
    #         genealogies.append(tree.newick(include_branch_lengths=branch_length))
    #     return genealogies

    def genealogies(self, branch_length=False):
        breakpoints = list(self.ts.breakpoints())
        trees = []
        for i, tree in enumerate(self.ts.aslist()):
            if tree.num_roots == 1:
                trees.append(tree.newick(include_branch_lengths=branch_length))
            else:
                trees.append(None)
        alls = []
        for variant in self.ts.variants():
            pos = variant.position
            for i in range(len(breakpoints)):
                if breakpoints[i] < pos <= breakpoints[i + 1]:
                    alls.append(trees[i])
        return alls

    def __str__(self):
        return f'{self.haplotype.nsamples} samples, {self.haplotype.nsites} sites, {self.ts.num_trees} topologies'

