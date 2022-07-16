"""
    Module for data simulation under scenarios users specify.
    Author: ZHT
"""


import os
import logging
from deeprho import LazyLoader
from deeprho.config import CONFIG
import argparse
import coloredlogs
import numpy as np
import pickle
from sklearn.model_selection import train_test_split
from multiprocessing.dummy import Pool
popgen = LazyLoader('deeprho.popgen')
logger = logging.getLogger(__name__)
global_window_size = 1000
window_size = 50


def build_configuration(args):
    configuration = {}
    configuration['sequence_length'] = 5e5
    configuration['rate'] = args.mutation_rate
    configuration['ploidy'] = args.ploidy
    if args.demography is not None:
        configuration['demography'] = popgen.utils.load_demography_from_file(args.demography)
    else:
        configuration['population_size'] = args.ne
    return configuration


def recombination_rate_quadratic_interpolation(nsam, i, rmin, rmax):
    rate = np.square((np.sqrt(rmax)-np.sqrt(rmin))/(nsam-1)*i + np.sqrt(rmin))
    return rate


def recombination_rate_const_interpolation(nsam, i, rmin, rmax):
    rate = (rmax-rmin) / (nsam-1)*i + rmin
    return rate


def save_training_data(path, data):
    assert path is not None, f'no file provided.'
    assert not os.path.exists(path), f'file has already existed.'
    x, y = data
    x_train, x_test, y_train, y_test = train_test_split(x, y, test_size=0.2)
    with open(path, mode='wb') as out:
        pickle.dump((x_train, x_test, y_train, y_test), out)
    logger.info(f'train size: {x_train.shape[0]}. test size: {x_test.shape[0]}')
    print(f'training data has been stored in {path}')


def simulate(configs, args, r):
    haplotypes = []
    genealogies = []
    rhos = []
    configs['recombination_rate'] = r
    simulator = popgen.Simulator(configs)
    for data in simulator(args.npop, 1):
        data = popgen.utils.filter_replicate(data)
        reps = popgen.utils.cut_replicate(data, window_size=global_window_size)
        haps = [rep.haplotype for rep in reps]
        inferred_genealogies = popgen.utils.rentplus(haps, num_thread=args.num_thread)
        for hap, gen in zip(haps, inferred_genealogies):
            for j in range(args.ndraw):
                start = np.random.choice(range(hap.nsites - window_size))
                haplotypes.append(hap.matrix[:, start:start + window_size])
                genealogies.append(gen[start: start + window_size])
                length = hap.positions[start + window_size] - hap.positions[start]
                # !have to redesign when demography included
                if args.demography is not None:
                    pop_size = configs['demography'].events[0].initial_size
                else:
                    pop_size = configs['population_size']
                scaled_rho = 2 * pop_size * r * length
                rhos.append(scaled_rho)
    logger.info(f'complete sampling in {r}, {data}')
    return haplotypes, genealogies, rhos


def run(args):
    assert args.out is not None, f'no output name.'
    assert args.rmax >= args.rmin, f'r_max should be greater than r_min.'
    if args.verbose:
        coloredlogs.install(logger=logger, level='INFO', field_styles=dict(
            asctime={"color": 2},
            message={"color": 6},
            levelname={"color": 3},
            programname={"color": 1}
        ), fmt='%(asctime)s [deeprho_v2] %(programname)s %(levelname)s - %(message)s')
    logger.info(f'----------- simulation -------------')
    logger.info(f'nsam:{args.nsam}, ndraw:{args.ndraw}')
    pool = Pool(args.num_thread // 2)
    haplotypes = []
    genealogies = []
    rhos = []
    # init simulator
    configs = build_configuration(args)
    # generate random haplotypes and infer their genealogies.
    paras = [(configs.copy(), args, recombination_rate_const_interpolation(args.nsam, i, args.rmin, args.rmax)) for i in range(args.nsam)]
    with pool:
        results = pool.starmap(simulate, paras)
    for result in results:
        haplotypes += result[0]
        genealogies += result[1]
        rhos += result[2]

    # build the whole set, train set and test set. compute Linkage Disequilibrium, ld_cluster, Robinson-Foulds distance, and also triplet distance.
    rf_distance = popgen.utils.rf_dist([list(val) for val in genealogies], num_thread=args.num_thread)
    tri_distance = popgen.utils.triplet_dist([list(val) for val in genealogies], num_thread=args.num_thread)
    lds = popgen.utils.linkage_disequilibrium(haplotypes)
    lds = np.expand_dims(np.array(lds), axis=-1)
    rfs = np.expand_dims(np.array(rf_distance), axis=-1)
    tris = np.expand_dims(np.array(tri_distance), axis=-1)
    data = np.concatenate([lds,rfs,tris], axis=-1).astype(np.float64)
    rhos = np.array(rhos).reshape(-1,1)
    save_training_data(args.out, (data, rhos))
    
def gt_args(parser):
    parser.add_argument('--nsam', type=int, help='number of sampling for rhos', default=CONFIG.N_SAMPLE)
    parser.add_argument('--ndraw', type=int, help='number of draws per sample', default=CONFIG.N_DRAW)
    parser.add_argument('--npop', type=int, help='number of individual', default=CONFIG.N_POP)
    parser.add_argument('--ne', type=float, help='effective population size', default=CONFIG.EFFECTIVE_POPULATION_SIZE)
    parser.add_argument('--ploidy', type=int, help='ploidy', default=CONFIG.PLOIDY)
    parser.add_argument('--mutation-rate', type=float, help='mutation rate', default=CONFIG.MUTATION_RATE)
    parser.add_argument('--demography', type=str, help='demography file path', default=CONFIG.DEMOGRAPHY)
    parser.add_argument('--rmin', type=float, help='minimum recombination rate', default=CONFIG.R_MIN)
    parser.add_argument('--rmax', type=float, help='maximum recombination rate', default=CONFIG.R_MAX)
    parser.add_argument('--num-thread', type=int, help='number of threads', default=CONFIG.NUM_THREAD)
    parser.add_argument('--out', type=str, help='output path')
    parser.add_argument('--verbose', help='show loggings', action='store_true')


if __name__ == '__main__':
    parser = argparse.ArgumentParser('data simulator')
    gt_args(parser)
    args = parser.parse_args(['--nsam', '10',
                              '--npop', '100',
                              '--ne', '1e4',
                              '--ploidy', '2',
                              '--rmin', '1e-9',
                              '--rmax', '1e-7',
                              '--out', '../garbo/train.data',
                              '--verbose'
                              ])
    run(args)


    
