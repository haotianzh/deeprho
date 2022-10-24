"""
    Author: Haotian Z
    Generate lookup table for a given demography or an effective population size
"""
import pickle
import os
import argparse
import logging
from deeprho import LazyLoader
from deeprho.config import CONFIG
utils = LazyLoader('deeprho.popgen.utils')
logger = logging.getLogger(__name__)


def run(args):
    assert args.out, f'an output name should be specified.'
    demography = None
    if args.demography:
        assert os.path.exists(args.demography), f'demography {args.demography} not found.'
        demography = utils.load_demography_from_file(args.demography)
    paras = {
        'population_size': args.ne,
        'demography': demography,
        'ploidy': args.ploidy,
        'num_thread': args.num_thread,
        'r_min': args.rmin,
        'r_max': args.rmax,
        'repeats': args.repeat,
        'draws': args.draw,
        'samples': args.npop
    }
    table = utils.get_lookup_table(**paras)
    with open(args.out, 'wb') as out:
        pickle.dump(table, out)



def gt_args(parser):
    parser.add_argument('--num-thread', type=int, help='number of threads', default=CONFIG.NUM_THREAD)
    parser.add_argument('--ne', type=float, help='effective population size', default=None)
    parser.add_argument('--demography', help='demography history', default=None)
    parser.add_argument('--ploidy', type=int, help='ploidy (default 2)', default=CONFIG.PLOIDY)
    parser.add_argument('--npop', type=int, help='number of individuals', default=CONFIG.MT_SAMPLE)
    parser.add_argument('--rmin', type=float, help='min of recombination rate per base per generation', default=CONFIG.MT_R_MIN)
    parser.add_argument('--rmax', type=float, help='max of recombination rate per base per generation', default=CONFIG.MT_R_MAX)
    parser.add_argument('--repeat', type=int, help='number of repeats in simulation', default=CONFIG.MT_REPEAT)
    parser.add_argument('--draw', type=int, help='number of draws', default=CONFIG.MT_DRAW)
    parser.add_argument('--out', type=str, help='table name', default=None)
    parser.add_argument('--verbose', help='show loggings', action='store_true')


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='deeprho estimator')
    gt_args(parser)
    args = parser.parse_args([
                              '--out', 'yri_table',
                              '--ploidy', '2',
                              '--demography', '../examples/YRI_pop_sizes.csv'])
    run(args)