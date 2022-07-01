from ..base import Haplotype, Replicate
import os
import tskit
import numpy as np


# for class Haplotype
def filter_none_mutation(hap):
    filtered_matrix = []
    filtered_pos = []
    for col, count in enumerate(hap.matrix.sum(axis=0)):
        if not count == 0 and not count == hap.matrix.shape[0]:
            filtered_matrix.append(hap.matrix[:, col])
            filtered_pos.append(hap.positions[col])
    filtered_matrix = np.array(filtered_matrix).T
    new_hap = Haplotype(matrix=filtered_matrix, positions=filtered_pos)
    return new_hap


# for class Haplotype
def sliding_windows(hap, window_size, step_size=None, drop_last=True):
    if step_size is None:
        step_size = window_size
    windows = []
    positions, matrix = hap.positions, hap.matrix
    for i in range(0, hap.nsites, step_size):
        if drop_last:
            if i + window_size > hap.nsites:
                break
        window_mat = matrix[:, i:i + window_size]
        window_pos = positions[i:i + window_size]
        # length = window_pos[-1] - window_pos[0]
        # scaled_positions = ((window_pos - window_pos[0]) / length * 1e5 + 1).astype(np.int)
        window_hap = Haplotype(matrix=window_mat, positions=window_pos)
        windows.append(window_hap)
    return windows


# for class Replicate
def filter_replicate(replicate):
    assert isinstance(replicate, Replicate), Exception("replicate should be tskit.TreeSequence")
    ts = replicate.ts
    sites_for_delete = []
    for site in ts.variants():
        if 0 not in site.genotypes or 1 not in site.genotypes:
            sites_for_delete.append(site.site.id)
    treeseq = ts.delete_sites(sites_for_delete)
    treeseq.__setattr__('rr', replicate.configs['recombination_rate'])
    treeseq.__setattr__('mr', replicate.configs['rate'])
    filtered_replicate = Replicate(treeseq, replicate.configs)
    return filtered_replicate


# for class Replicate
def cut_replicate(replicate, window_size=50, drop_last=True):
    ts = replicate.ts
    assert ts is not None, Exception("ts shouldn't be none.")
    indices = [0]
    count = 0
    for variant in ts.variants():
        count += 1
        if count == window_size:
            indices.append(variant.position+1)
            count = 0
    replicates = []
    for i in range(1, len(indices)):
        ts_fragment = ts.keep_intervals([[indices[i - 1], indices[i]]])
        ts_fragment.__setattr__('rr', replicate.configs['recombination_rate'])
        ts_fragment.__setattr__('mr', replicate.configs['rate'])
        rep_fragment = Replicate(ts_fragment, replicate.configs)
        replicates.append(rep_fragment)
    if not drop_last:
        if count == 0:
            ts_last = None
        else:
            ts_last = ts.keep_intervals([[indices[-1], replicate.haplotype.positions[-1]+1]])
            ts_last.__setattr__('rr', replicate.configs['recombination_rate'])
            ts_last.__setattr__('mr', replicate.configs['rate'])
            rep_last = Replicate(ts_last, replicate.configs)
        return replicates, rep_last
    return replicates


def write_vcf(trees, name):
    out = open(name, 'w')
    trees.write_vcf(out, ploidy=2)
    out.close()


def write_fasta(trees, name):
    out = open(name, 'w')
    trees.write_fasta(out)
    out.close()


def write_ms(trees, name):
    out = open(name, 'w')
    genotypes = trees.genotype_matrix()
    positions = []
    for v in trees.variants():
        positions.append(v.position)
    positions = np.array(positions).astype(int)
    out.write('ms {size} 1 -t {mu} -r {r} {length}\n'.format(size=100, mu=1e-8, r=1e-8, length=1e5))
    out.write('24263 40612 14324\n\n')
    out.write('//\n')
    out.write('segsites: {0}\n'.format(len(positions)))
    out.write('positions: ')
    out.write(' '.join(['{}'.format(val) for val in positions]))
    out.write('\n')
    for row in genotypes.T:
        out.write(''.join([str(val) for val in row]))
        out.write('\n')
    out.close()
