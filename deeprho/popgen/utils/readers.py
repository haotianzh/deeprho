"""
    author: Haotian
    created_at: 12/2/2021
    description: loading data from different file types and systems, especially for "ms-format" and "vcf-format"
"""
import vcf
import pandas as pd
import numpy as np
import os
import msprime as msp
from ..base import Haplotype

# define hyperparameter in demography
GENERATION = 25

def check_file_existence(file):
    assert os.path.exists(file), f'file {file} does not exist.'

def load_ms_from_file(file, true_genealogy=False):
    """
        load haplotype data from ms-formatted file.
        Input: file path, true_genealogy = True if "-T" option is specified in simulation. (currently not implemented)
        Output: Haplotype
    """
    check_file_existence(file)
    f = open(file, 'r')
    ms_command = f.readline().strip().split()
    pop_size = int(ms_command[1])
    sample_size = int(ms_command[2])
    random_seeds = f.readline().strip().split()
    data = {'ms_command': ms_command, 'n_pop':pop_size, 'n_hap':sample_size, 'ms_data': []}
    # read for ms outputs
    labels = []
    sample_data = []  # initial list for each sample.
    line = 'a'
    count = 0
    while line:
        line = f.readline()
        if line.strip() == '':
            count += 1
            if line:
                pass
            if sample_data:
                data['ms_data'].append(np.array(sample_data))  # append a new sample into memory
            sample_data = []
            continue
        # tbs option in ms program
        if line.startswith('//'):
            if len(line.strip().split()) == 1:
                pass
            else:
                label = float(line.strip().split()[1])
                labels.append(label)
            continue
        if not line[0].isdigit():
            key, values = line.strip().split(':')
            if key in data:
                data[key].append(np.array([float(val) for val in values.strip().split()]))
            else:
                data[key] = [np.array([float(val) for val in values.strip().split()])]
            continue
        else:
            sample_data.append([int(val) for val in list(line.strip())])
            continue
    haplotype = Haplotype(matrix=data['ms_data'][0], positions=data['positions'][0].astype(np.int))
    return haplotype

def load_vcf_from_file(file):
    """
        load haplotype data from VCF-formatted file.
        Input: file path
        Output: Haplotype
    """
    check_file_existence(file)
    reader = vcf.Reader(open(file, 'r'))
    sites = []
    positions =[]
    for record in reader:
        site = []
        missing = False
        for i, sample in enumerate(record.samples):
            alleles = sample.gt_alleles
            if sample.ploidity > 1:
                assert '|' in sample.data.GT , f'data should be phased when ploidy is greater than 1.'
            if None in alleles: # skip missing calls.
                missing = True
                break
            for allele in alleles:
                site.append(int(allele))
        if not missing:
            sites.append(site)
            positions.append(record.POS)
    matrix = np.array(sites).T
    positions = np.array(positions)
    haplotype = Haplotype(matrix=matrix, positions=positions)
    return haplotype



def load_demography_from_file(file, generation=GENERATION):
    """
        loading demography from file using msprime
        Input: smc++ output  (5 columns)
        ======================================================================
        |label	|x                  |y	                |plot_type	|plot_num|
        |--------------------------------------------------------------------|
        |ACB	|0.0	            |138482.84333082315	|path	    |0       |
        |ACB	|50.0	            |138482.84333082315	|path	    |0       |
        |ACB	|53.97505585700569  |139331.82583178935	|path	    |0       |
        ======================================================================
    """
    check_file_existence(file)
    demography_file = pd.read_csv(file)
    assert demography_file['label'].unique().size == 1, 'there are more than one population.'
    demography = msp.Demography()
    population_name = demography_file['label'][0]
    demography.add_population(name=population_name, initial_size=1e5)
    for i, row in demography_file.iterrows():
        demography.add_population_parameters_change(time=row['x']/generation, initial_size=row['y'])
    return demography


if __name__ == '__main__':
    filename = 'C:/Users/zht/Desktop/deeprho_ results/examples/data.vcf'
    load_vcf_from_file(filename)