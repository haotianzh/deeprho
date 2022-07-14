"""
    Estimator of recombination rate. Be capable to scale properly even under demography.
    Author: ZHT
"""


import os
import argparse
import pathlib
import coloredlogs
import logging
import numpy as np
import matplotlib.pyplot as plt
from deeprho import LazyLoader
from deeprho.config import CONFIG
os.environ['TF_CPP_MIN_LOG_LEVEL'] = '2'
tf = LazyLoader('tensorflow')
utils = LazyLoader('deeprho.popgen.utils')
logger = logging.getLogger(__name__)


def get_deeprho_map(rates, bounds, length):
    length = int(length)
    starts = [val[0] for val in bounds]
    ends = [val[1] for val in bounds]
    per_site_rates = [[] for i in range(length)]
    for start, end, r in zip(starts, ends, rates):
        for i in range(start, end):
            per_site_rates[i].append(r)
    crates = np.zeros(length, np.float64)
    for i, r in enumerate(per_site_rates):
        if r:
            crates[i] = np.mean(r)
            # if len(r) == 2:
            #     crates[i] = np.dot(np.array(r), np.array([0.3,0.7]))
    return crates


def load_data(file):
    _extensions = {'.ms': utils.load_ms_from_file, '.vcf': utils.load_vcf_from_file}
    ext = pathlib.Path(file).suffix
    assert ext in _extensions, f'only {_extensions} formats are supported.'
    return _extensions[ext](file)


def output(rates, out_name):
    with open(out_name, 'w') as out:
        out.write('Start\tEnd\tRate')
        start = 0
        for i, rate in enumerate(rates):
            if not rate == rates[start]:
                out.write(f'\n{start}\t{i}\t{rates[start]}')
                start = i
    print(f"result is saved as '{os.path.abspath(out_name)}'")


def plot(rates, threshold, out_name):
    plt.figure(figsize=(12,6))
    plt.title('deeprho estimates')
    plt.xlabel('bp')
    plt.ylabel('recombination rate')
    plt.plot(rates, label='deeprho')
    plt.axhline(y=threshold, color='r', linestyle='--', label='hotspot threshold')
    plt.legend()
    plt.savefig(out_name, dpi=100)
    print(f"figure is saved as '{os.path.abspath(out_name)}'")


def estimate(haplotype, model_fine_path, model_large_path, window_size=50, step_size=50,
             sequence_length=None, global_window_size=1000, n_pop=100, ploidy=1, ne=1e5,
             resolution=1e4, threshold=5e-8, num_thread=4):
    haplotype = utils.filter_none_mutation(haplotype)
    haplotypes = utils.sliding_windows(haplotype,
                                 window_size=window_size,
                                 step_size=step_size,
                                 drop_last=True)
    # infer genealogy with global window size as 1e3
    logger.info('inferring local genealogies')
    global_genealogies = []
    global_slices = utils.sliding_windows(haplotype,
                                    window_size=global_window_size,
                                    step_size=global_window_size,
                                    drop_last=False)
    _pre_slices = utils.rentplus(global_slices[:-1], num_thread=num_thread)
    _last_slice = utils.rentplus(global_slices[-1])
    for sli in _pre_slices:
        global_genealogies += sli
    global_genealogies += _last_slice
    _slices = []
    for i in range(0, haplotype.nsites, step_size):
        if i + window_size <= haplotype.nsites:
            _slices.append(global_genealogies[i: i+window_size])
    # calculate distance
    logger.info('calculating distances')
    lds = utils.linkage_disequilibrium(haplotypes)
    rfs = utils.rf_dist(_slices, num_thread=num_thread)
    tris = utils.triplet_dist(_slices, num_thread=num_thread)
    lds = np.expand_dims(np.array(lds, np.float32), axis=-1)
    rfs = np.expand_dims(np.array(rfs, np.float32), axis=-1) / (2*(n_pop-3))
    tris = np.expand_dims(np.array(tris, np.float32), axis=-1) / (n_pop*(n_pop-1)*(n_pop-2)/6)
    x = np.concatenate([lds, rfs, tris], axis=-1)
    # two-stages model
    if sequence_length is None:
        sequence_length = haplotype.positions[-1]
    logger.info('loading two-stages models')
    model_fine = tf.keras.models.load_model(model_fine_path)
    model_large = tf.keras.models.load_model(model_large_path)
    logger.info('predicting')
    rates = model_fine.predict(x, verbose=0)
    rates_large = model_large.predict(x, verbose=0) * CONFIG.SCALE_FACTOR
    scaled_rates, _, __ = utils.stat(rates, haplotype.positions,
                               sequence_length=sequence_length,
                               bin_width=resolution,
                               window_size=window_size,
                               step_size=step_size,
                               ploidy=ploidy,
                               ne=ne)
    scaled_rates_large, _, __ = utils.stat(rates_large, haplotype.positions,
                                     sequence_length=sequence_length,
                                     bin_width=resolution,
                                     window_size=window_size,
                                     step_size=step_size,
                                     ploidy=ploidy,
                                     ne=ne)
    r_fine = get_deeprho_map(scaled_rates, _, length=sequence_length)
    r_large = get_deeprho_map(scaled_rates_large, _, length=sequence_length)
    r_fine[r_fine>threshold] = r_large[r_fine>threshold]
    return r_fine


def run(args):
    assert args.file is not None, f'no input provided.'
    assert args.m1 is not None and args.m2 is not None, f'no specified models, --m1 and --m2 should be specified.'
    assert os.path.exists(args.m1) and os.path.exists(args.m2), f'model file not found.'

    if args.verbose:
        coloredlogs.install(logger=logger, level='INFO', field_styles=dict(
            asctime={"color": 2},
            message={"color": 6},
            levelname={"color": 3},
            programname={"color": 1}
        ),  fmt='%(asctime)s [deeprho_v2] %(programname)s %(levelname)s - %(message)s')
    logger.info(f'loading data from {args.file}')
    haplotype = load_data(args.file)
    logger.debug(f'data: {haplotype.nsites} SNPs, {haplotype.nsamples} individuals')
    length = args.length
    if length is None:
        length = haplotype.positions[-1] - haplotype.positions[0]
    step_size = args.ss
    if step_size is None:
        step_size = args.ws
    ne = args.ne
    if args.demography is not None:
        ne = utils.calculate_average_ne(args.demography)
    paras = {
        'haplotype': haplotype,
        'num_thread': args.num_thread,
        'sequence_length': length,
        'ne': ne,
        'global_window_size': args.gws,
        'window_size': args.ws,
        'step_size': step_size,
        'threshold': args.threshold,
        'model_fine_path': args.m1,
        'model_large_path': args.m2,
        'ploidy': args.ploidy,
        'n_pop': haplotype.nsamples,
        'resolution': args.res
    }
    rates = estimate(**paras)
    output(rates, args.file + '.rate.txt')
    if args.savenp:
        np.save(args.file + '.rate.npy', rates)
        print(f"numpy object is saved as '{os.path.abspath(args.file + '.rate.npy')}'")
    if args.plot:
        plot(rates, args.threshold, args.file + '.rate.png')


def gt_args(parser):
    parser.add_argument('--file', type=str, help='filename')
    parser.add_argument('--num-thread', type=int, help='number of threads', default=CONFIG.NUM_THREAD)
    parser.add_argument('--ne', type=float, help='effective population size', default=CONFIG.EFFECTIVE_POPULATION_SIZE)
    parser.add_argument('--demography', help='demography history', default=CONFIG.DEMOGRAPHY)
    parser.add_argument('--length', type=float, help='genome length', default=CONFIG.LENGTH)
    parser.add_argument('--ploidy', type=int, help='ploidy (default 2)', default=CONFIG.PLOIDY)
    parser.add_argument('--threshold', type=float, help='hotspot threshold', default=CONFIG.THRESHOLD)
    parser.add_argument('--gws', type=int, help='global window size', default=CONFIG.GLOBAL_WINDOW_SIZE)
    parser.add_argument('--ws', type=int, help='window size', default=CONFIG.WINDOW_SIZE)
    parser.add_argument('--ss', type=int, help='step size', default=CONFIG.STEP_SIZE)
    parser.add_argument('--res', type=float, help='resolution(bp)', default=CONFIG.RESOLUTION)
    parser.add_argument('--m1', type=str, help='fine-model path', default=CONFIG.MODEL_FINE)
    parser.add_argument('--m2', type=str, help='large-model path', default=CONFIG.MODEL_LARGE)
    parser.add_argument('--plot', help='plot or not', action='store_true')
    parser.add_argument('--savenp', help='save as numpy object', action='store_true')
    parser.add_argument('--verbose', help='show loggings', action='store_true')


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='deeprho estimator')
    gt_args(parser)
    args = parser.parse_args(['--file', '../garbo/test.vcf',
                              '--demography', '../examples/ACB_pop_sizes.csv',
                              '--ploidy', '2',
                              '--m1', '../models/model_epoch_99.hdf5',
                              '--m2', '../models/model_epoch_123.hdf5',
                              '--res', '1e3',
                              '--plot',
                              '--verbose'])
    run(args)
