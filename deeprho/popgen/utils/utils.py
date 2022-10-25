from ..base import Haplotype, Replicate
from . import *
import os
import tskit
import pickle
import numpy as np
import multiprocessing as mp


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


def get_num_bkpts(configs, samples=50, repeats=50):
    points = []
    simulator = popgen.Simulator(configs)
    for rep in simulator(samples, repeats):
        points.append(rep.ts.num_trees + 1)
    return np.mean(points)


SEQUENCE_LENGTH = 2e3
def get_lookup_table(population_size=5e4, demography=None, samples=50, repeats=50,
                          ploidy=2, r_min=0, r_max=1e-6, draws=200, num_thread=8):
    configs = {
        'rate': 2.5e-8,  # mutation rate
        'ploidy': ploidy,  # number of ploidy
        'sequence_length': SEQUENCE_LENGTH,  # genome length
    }
    rates =  np.linspace(r_min, r_max, draws)
    ps = []
    with mp.get_context('spawn').Pool() as pool:
        params = []
        for rate in rates:
            bc = configs.copy()
            bc['recombination_rate'] = rate
            if demography is not None:
                bc['demography'] = demography
            else:
                bc['population_size'] = population_size
            params.append((bc, samples, repeats))
        for res in pool.starmap(get_num_bkpts, params):
            ps.append(res)
    return {'rates': rates, 'bkpts': ps}


def harmonic_sum(n):
    """
        return \sum_1^{n-1} 1/n
    """
    return 2 * sum([1 / i for i in range(2, n)])


def get_r_given_bkpts(n, rates, bkpts):
    '''
        getting the true r given number of bkpts under a certain population model
    '''
    def linear_interpolation(n, i, j):
        a, b, c, d = bkpts[i], bkpts[j], rates[i], rates[j]
        return (n-a)*(d-c)/(b-a)+c
    i, j = 0, len(rates)-1
    while i<j:
        mid = (i+j) // 2
        if bkpts[mid] < n <= bkpts[mid+1]:
            return linear_interpolation(n, mid, mid+1)
        elif bkpts[mid+1] < n:
            i = mid + 1
        else:
            j = mid - 1
    return rates[i]


def get_bkpts(table, paras):
    '''
        Input:
            table: precalculated lookup table under constant model (rhos x lengths)
            paras[0]: rho
            paras[1]: length
        Return:
            number of bkpts
    '''
    def bsearch(li, val):
        list_len = len(li)
        i, j = 0, list_len-1
        while i <= j:
            mid = (i+j) // 2
            a, b = li[min(mid+1, list_len-1)], li[mid]
            if a > val >= b:
                return b, a
            elif val < b:
                j = mid - 1
            else:
                i = mid + 1
        if j < 0:
            return li[i], li[i]
        else:
            return li[j], li[j]
    rhos = table['rho'].unique()
    lens = table['len'].unique()
    table = table.set_index(['rho', 'len']).sort_index()
    bkpts = []
    for rho, l in paras:
        rs= bsearch(rhos, rho)
        ls = bsearch(lens, l)
        vals = []
        for r in rs:
            for ll in ls:
                vals.append(table.loc[(r, ll)].bkpt)
        if rs[0] == rs[1]:
            v1, v2 = vals[0], vals[1]
        else:
            ratio = (rho-rs[0])/(rs[1]-rs[0])
            v1, v2 = ratio*vals[2]+(1-ratio)*vals[0], ratio*vals[3]+(1-ratio)*vals[1]
        if ls[0] == ls[1]:
            res = v1
        else:
            ratio = (l-ls[0])/(ls[1]-ls[0])
            res = ratio*v2 + (1-ratio)*v1
        bkpts.append(res)
    return bkpts


def convert2r(table_constant, table, rhos, positions, window_size=50, step_size=50):
    '''
        convert estimated rho to r based on precomputed lookup table
    '''
    rhos = rhos.reshape(-1)
    positions = positions.reshape(-1)
    lens = []
    rights = []
    lefts = []
    for i, r in enumerate(rhos):
        right = positions[i*step_size + window_size - 1]
        left = positions[i*step_size]
        lens.append(right - left)
        rights.append(right)
        lefts.append(left)
    paras = zip(rhos, lens)
    bkpts = get_bkpts(table_constant, paras)
    rs = []
    for bkpt, l, left, right in zip(bkpts, lens, lefts, rights):
        r_prime = get_r_given_bkpts(bkpt, table['rates'], table['bkpts']) * SEQUENCE_LENGTH / l
        rs.append(r_prime)
    return np.array(rs), np.array(lefts), np.array(rights)


# # --------------------- OTHERS ---------------------------------
# def get_pyrho_map(file, length=100000):
#     pyrho = pd.read_table(file, header=None)
#     pyrate = np.zeros(length, dtype=np.float64)
#     for start, end, r in zip(pyrho[0], pyrho[1], pyrho[2]):
#         pyrate[int(start): int(end)] = r
#     return pyrate
#
#
# def get_deeprho_map(rs, positions, window_size=50, step_size=50):
#     rs = rs.reshape(-1)
#     ls = []
#     rights = []
#     lefts = []
#     for i, r in enumerate(rs):
#         right = positions[i*step_size + window_size - 1]
#         left = positions[i*step_size]
#         l = right - left
#         ls.append(l)
#         rights.append(right)
#         lefts.append(left)
#     res = np.zeros(positions[-1])
#     for r, left, right in zip(rs, lefts, rights):
#         res[left:right] = r
#     return res
#
#
# def get_deeprho_map_1(rates, bounds):
#     starts = [val[0] for val in bounds]
#     ends = [val[1] for val in bounds]
#     begin = starts[0]
#     end = ends[-1]
#     length = end - begin
#     per_site_rates = [[] for i in range(length+1)]
#     for start, end, r in zip(starts, ends, rates):
#         for i in range(start, end):
#             per_site_rates[i-begin].append(r)
#     crates = np.zeros(length+1, np.float64)
#     for i, r in enumerate(per_site_rates):
#         if r:
#             crates[i] = np.mean(r)
#             # if len(r) == 2:
#             #     crates[i] = np.dot(np.array(r), np.array([0.3,0.7]))
#     return crates
#
#
# def get_true_map(file, length=100000):
#     with open(file, 'rb') as stream:
#         mslice = pickle.load(stream)
#     rate = np.zeros(length, dtype=np.float64)
#     for start, end, r in zip(mslice.position[:-1], mslice.position[1:], mslice.rate):
#         rate[int(start): int(end)]= r
#     return rate


# ---------------------- IO PARTS ------------------------------
def write_vcf(trees, name):
    '''
        write to .vcf using tskit.write_vcf
    '''
    out = open(name, 'w')
    trees.write_vcf(out, ploidy=2)
    out.close()


def write_fasta(trees, name):
    '''
        write to .fa using tskit.write_fasta
    '''
    out = open(name, 'w')
    trees.write_fasta(out)
    out.close()


def write_ms(trees, name):
    '''
        write to .ms (random seed doesn't matter)
    '''
    out = open(name, 'w')
    genotypes = trees.genotype_matrix()
    positions = []
    for v in trees.variants():
        positions.append(v.position)
    positions = np.array(positions).astype(int)
    out.write('ms {size} 1 -t {mu} -r {r} {length}\n'.format(size=genotypes.shape[1], mu=1e-8, r=1e-8, length=positions[-1]))
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
