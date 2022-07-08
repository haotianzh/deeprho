import os
import pathlib
import argparse
import pathlib
import logging
import numpy as np
import matplotlib.pyplot as plt
os.environ['TF_CPP_MIN_LOG_LEVEL'] = '1'
import tensorflow as tf
# import popgen related lib
from deeprho.popgen.utils import load_vcf_from_file, load_ms_from_file
from deeprho.popgen.utils import filter_none_mutation, sliding_windows, stat
from deeprho.popgen.utils import rentplus, rfdist, triplet_dist, linkage_disequilibrium

# how much the scaling for rate during training. 10 means 10x
scaling_factor = 10


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
    _extensions = {'.ms': load_ms_from_file, '.vcf': load_vcf_from_file}
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
    print(f"result is saved as '{out_name}'")


def plot(rates, threshold, out_name):
    plt.figure(figsize=(12,6))
    plt.title('deeprho estimates')
    plt.xlabel('bp')
    plt.ylabel('recombination rate')
    plt.plot(rates, label='deeprho')
    plt.axhline(y=threshold, color='r', linestyle='--', label='hotspot threshold')
    plt.legend()
    plt.savefig(out_name, dpi=100)
    print(f"figure is saved as '{out_name}'")

def estimate(haplotype, model_fine_path, model_large_path, window_size=50, step_size=50,
             sequence_length=None, global_window_size=1000, n_pop=100, ploidy=1, ne=1e5,
             resolution=1e4, threshold=5e-8, num_thread=4):
    haplotype = filter_none_mutation(haplotype)
    haplotypes = sliding_windows(haplotype,
                                 window_size=window_size,
                                 step_size=step_size,
                                 drop_last=True)
    # infer genealogy with global window size as 1e3
    logging.info('inferring local genealogies')
    global_genealogies = []
    global_slices = sliding_windows(haplotype,
                                    window_size=global_window_size,
                                    step_size=global_window_size,
                                    drop_last=False)
    _pre_slices = rentplus(global_slices[:-1], num_thread=num_thread)
    _last_slice = rentplus(global_slices[-1])
    for sli in _pre_slices:
        global_genealogies += sli
    global_genealogies += _last_slice
    _slices = []
    for i in range(0, haplotype.nsites, step_size):
        if i + window_size <= haplotype.nsites:
            _slices.append(global_genealogies[i: i+window_size])
    # calculate distance
    logging.info('calculating distances')
    lds = linkage_disequilibrium(haplotypes)
    rfs = rfdist(_slices, num_thread=num_thread)
    tris = triplet_dist(_slices, num_thread=num_thread)
    lds = np.expand_dims(np.array(lds, np.float32), axis=-1)
    rfs = np.expand_dims(np.array(rfs, np.float32), axis=-1) / (2*(n_pop-3))
    tris = np.expand_dims(np.array(tris, np.float32), axis=-1) / (n_pop*(n_pop-1)*(n_pop-2)/6)
    x = np.concatenate([lds, rfs, tris], axis=-1)
    # two-stages model
    if sequence_length is None:
        sequence_length = haplotype.positions[-1]
    logging.info('loading two-stages models')
    model_fine = tf.keras.models.load_model(model_fine_path)
    model_large = tf.keras.models.load_model(model_large_path)
    logging.info('predicting')
    rates = model_fine.predict(x)
    rates_large = model_large.predict(x) * scaling_factor
    scaled_rates, _, __ = stat(rates, haplotype.positions,
                               sequence_length=sequence_length,
                               bin_width=resolution,
                               window_size=window_size,
                               step_size=step_size,
                               ploidy=ploidy,
                               ne=ne)
    scaled_rates_large, _, __ = stat(rates_large, haplotype.positions,
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
        logging.basicConfig(format=f'[deeprho_v2] {os.path.basename(__file__)} %(levelname)s %(asctime)s - %(message)s',
                            level=logging.INFO,
                            datefmt='%m/%d %I:%M:%S')
    logging.info(f'loading data from {args.file}')
    haplotype = load_data(args.file)
    logging.debug(f'data: {haplotype.nsites} SNPs, {haplotype.nsamples} individuals')
    length = args.length
    if length is None:
        length = haplotype.positions[-1] - haplotype.positions[0]
    step_size = args.ss
    if step_size is None:
        step_size = args.ws

    paras = {
        'haplotype': haplotype,
        'num_thread': args.num_thread,
        'sequence_length': length,
        'ne': args.ne,
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
    out_prefix = pathlib.Path(args.file).name.split('.')[0]
    output(rates, out_prefix+'.out.txt')
    if args.savenp:
        np.save(out_prefix + '.out.npy', rates)
        print(f"numpy object is saved as '{out_prefix + '.out.npy'}'")
    if args.plot:
        plot(rates, args.threshold, out_prefix + '.out.png')


def gt_args(parser):
    # model_default_dir = pathlib.Path(popgen.__file__).parent.parent.joinpath('models')
    # model_fine_default_path = model_default_dir.joinpath('model_fine.hdf5')
    # model_large_default_path = model_default_dir.joinpath('model_large.hdf5')
    parser.add_argument('--file', type=str, help='filename')
    parser.add_argument('--num-thread', type=int, help='number of threads', default=4)
    parser.add_argument('--ne', type=float, help='effective population size', default=1e5)
    parser.add_argument('--length', type=float, help='genome length', default=None)
    parser.add_argument('--ploidy', type=int, help='ploidy (default 1)', default=1)
    parser.add_argument('--threshold', type=float, help='hotspot threshold', default=5e-8)
    parser.add_argument('--gws', type=int, help='global window size', default=1000)
    parser.add_argument('--ws', type=int, help='window size', default=50)
    parser.add_argument('--ss', type=int, help='step size', default=None)
    parser.add_argument('--res', type=float, help='resolution(bp)', default=1e4)
    parser.add_argument('--m1', type=str, help='fine-model path', default=None)
    parser.add_argument('--m2', type=str, help='large-model path', default=None)
    parser.add_argument('--plot', help='plot or not', action='store_true')
    parser.add_argument('--savenp', help='save as numpy object', action='store_true')
    parser.add_argument('--verbose', help='show loggings', action='store_true')


if __name__ == '__main__':
    logging.basicConfig(format=f'[deeprho_v2] {os.path.basename(__file__)} %(levelname)s %(asctime)s - %(message)s',
                        level=logging.INFO,
                        datefmt='%m/%d %I:%M:%S')
    parser = argparse.ArgumentParser(description='deeprho estimator')
    gt_args(parser)
    args = parser.parse_args(['--file', 'test.vcf',
                              '--ploidy', '2',
                              '--m1', '../models/model_fine.hdf5',
                              '--m2', '../models/model_large.hdf5',
                              '--ne', '138482',
                              '--res', '1e3',
                              '--plot'])
    run(args)
