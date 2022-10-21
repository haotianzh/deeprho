"""
    Author: Haotian Z
    Simulate a test case
"""
import argparse
import os
import coloredlogs
import logging
import pickle
import pathlib
import multiprocessing as mp
from deeprho import LazyLoader
popgen = LazyLoader('deeprho.popgen')
logger = logging.getLogger(__name__)


def build_configuration(args):
    configuration = {}
    configuration['rate'] = args.mutation_rate
    configuration['ploidy'] = args.ploidy
    if args.rate_map:
        path = pathlib.Path(args.rate_map)
        if path.suffix == '.pkl':
            map = pickle.load(open(args.rate_map, 'rb'))
        else:
            map = popgen.utils.load_recombination_map_from_file(args.rate_map)
        configuration['recombination_rate'] = map
        configuration['sequence_length'] = map.sequence_length
    else:
        configuration['recombination_rate'] = args.recombination_rate
        configuration['sequence_length'] = args.sequence_length
    if args.demography:
        configuration['demography'] = popgen.utils.load_demography_from_file(args.demography, mode='generation', generation=29)
    else:
        configuration['population_size'] = args.ne
    return configuration

def simulate_single_genome(configs, args):
    simulator = popgen.Simulator(configs)
    logger.info(str(configs))
    genome = next(simulator(args.npop, 1))
    return genome


def run(args):
    if args.verbose:
        coloredlogs.install(logger=logger, level='INFO', field_styles=dict(
            asctime={"color": 2},
            message={"color": 6},
            levelname={"color": 3},
            programname={"color": 1}
        ),  fmt='%(asctime)s [deeprho_v2] %(programname)s %(levelname)s - %(message)s')
    assert args.out, f'no output name.'
    if args.demography:
        assert os.path.exists(args.demography)
    if args.rate_map:
        assert os.path.exists(args.rate_map)
    configs = build_configuration(args)
    genome = simulate_single_genome(configs, args)
    with open(args.out, 'w') as vcf_file:
        genome.ts.write_vcf(output=vcf_file)


def gt_args(parser):
    parser.add_argument('--npop', type=int, help='number of individuals', default=50)
    parser.add_argument('--ne', type=float, help='effective population size', default=5e4)
    parser.add_argument('--ploidy', type=int, help='ploidy', default=2)
    parser.add_argument('--mutation-rate', type=float, help='mutation rate', default=2.5e-8)
    parser.add_argument('--demography', type=str, help='demography file path', default=None)
    parser.add_argument('--recombination-rate', type=float, help='recombination rate')
    parser.add_argument('--sequence-length', type=float, help='sequence length of genome', default=5e5)
    parser.add_argument('--rate-map', type=str, help='recombination rate map', default=None)
    parser.add_argument('--num-thread', type=int, help='number of threads', default=mp.cpu_count()-2)
    parser.add_argument('--out', type=str, help='output path')
    parser.add_argument('--verbose', help='show loggings', action='store_true')


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='simulate whole genome')
    gt_args(parser)
    args = parser.parse_args(['--npop','50',
                              '--ploidy', '2',
                              # '--mutation-rate', '2.5e-8',
                              '--rate-map', '../examples/test_recombination_map.txt',
                              # '--rate-map', '../garbo/map.pkl',
                              '--demography', '../examples/YRI_pop_sizes.csv',
                              # '--demography', 'ms.txt.demo.csv',
                              # '--ne', '50000',
                              '--out', '../garbo/test7.vcf',
                              '--verbose'])
    run(args)