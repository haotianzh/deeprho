from deeprho import popgen
from deeprho.popgen import utils
import numpy as np
import pandas as pd
import scipy
import matplotlib.pyplot as plt
import multiprocessing as mp
from time import time
import pickle


# create lookup table for rho in range [0, 1000]
rho_min = 0
rho_max = 1000  # 896 is the max in train set.


def gd(n, l, rho):
    samples = 50
    r = rho / 4 / n / l
    simulator = popgen.Simulator()
    configs = {'rate': 2.5e-8, 'ploidy': 2}
    configs['population_size'] = n
    configs['sequence_length'] = l
    configs['recombination_rate'] = r
    simulator.set_configs(configs)
    tmp = []
    for i in range(50):
        try:
            d = next(simulator(samples, 1))
            tmp.append(d.ts.num_trees + 1)
        except:
            pass
    print('_'.join([str(int(rho)), str(int(n)), str(int(l))]))
    return rho, l, n, np.mean(tmp)


def generate_for_human(samples=24):
    lookup_table = {}
    #     configs = {'rate': 2.5e-8, 'ploidy':2}
    N = np.linspace(10000, 100000, 5)
    Rho = np.linspace(rho_min, rho_max, 250)
    L = np.arange(17) * 500 + 2000
    #     simulator = popgen.Simulator()
    paras = []
    for rho in Rho:
        for l in L:
            for n in N:
                paras.append([n, l, rho])

    with mp.Pool(10) as pool:
        res = pool.starmap(gd, paras)

    return res


if __name__ ==  '__main__':
    start = time()
    table = generate_for_human()
    end = time()
    with open('../table/table_500_1000.pkl', 'wb') as out:
        pickle.dump(table, out)
    print(end - start)