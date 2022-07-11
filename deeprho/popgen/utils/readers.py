"""
    author: Haotian
    created_at: 12/2/2021
    description: loading data from different file types and systems, especially for "ms-format" and "vcf-format"
"""
import os
import logging
import pandas as pd
import numpy as np
import re
from collections import deque, namedtuple
import msprime as msp
from ..base import Haplotype

use_pyvcf = False
logger = logging.getLogger(__name__)
try:
    import vcf
    use_pyvcf = True
except:
    logger.warning('PyVCF module not installed, use builtin VCF reader.')


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
    if use_pyvcf:
        return load_vcf_from_file_1(file)
    return load_vcf_from_file_2(file)


# use PyVCF, but issues sometimes happen during pip install.
def load_vcf_from_file_1(file):
    """
        Load haplotype data from VCF file.
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
            if None in alleles:# skip missing calls.
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


def load_vcf_from_file_2(file):
    check_file_existence(file)
    vcf = read_vcf(file)
    df = vcf.data[vcf.data['FILTER'].apply(lambda x: x.lower())=='pass']
    sites = []
    positions = []
    samples = df.columns[9:].to_list()
    for i, row in df.iterrows():
        UNPHASED = False
        MISSING = False
        ref, alt, pos, format = row['REF'], row['ALT'], row['POS'], row['FORMAT']
        if 'GT' not in format:
            raise Exception('no genotype data in VCF file.')
        sep = ' '
        if not format == 'GT':
            sep = format[2]
        site = []
        for sample in samples:
            gt = row[sample].split(sep)[0]
            if '/' in gt:
                UNPHASED = True
                break
            if '.' in gt:
                MISSING = True
                break
            alleles = [0 if val == str(ref) else 1 for val in gt.split('|')]
            site.extend(alleles)
        if MISSING:
            logger.warning('missing alleles are detected, sites with missing alleles are removed.')
            continue
        if UNPHASED:
            logger.warning('unphased genotypes are detected, sites with unphased genotypes are removed.')
            continue
        sites.append(site)
        positions.append(pos)
    haplotype = Haplotype(positions=np.array(positions), matrix=np.array(sites).T.astype(np.int))
    return haplotype


def read_vcf(file):
    """
        Read VCF file.
        VCF v4.2 docs: https://samtools.github.io/hts-specs/VCFv4.2.pdf
        Input: VCF file path
        Output: namedtuplle(metadata, data)
    """
    # check_file_existence(file)
    # read meta
    VCF = namedtuple('VCF', ['metadata', 'data', 'file'])
    count_comments = 0
    meta = {}
    with open(file, 'r') as f:
        line = f.readline().strip()
        while line.startswith('##'):
            count_comments += 1
            info = parse_vcf_metadata(line)
            for key, value in info.items():
                if key not in meta:
                    meta[key] = value
                elif isinstance(meta[key], list):
                    meta[key].append(value)
                else:
                    meta[key] = [meta[key], value]
            line = f.readline().strip()
    vcf_df = pd.read_table(file, skiprows=count_comments)
    return VCF(meta, vcf_df, file)


def parse_vcf_metadata(st):
    """
        Parse string like:
        "##INFO=<ID=ID,Number=number,Type=type,Description="description",Source="source",Version="version">"

        Output: a json string
    """
    entities = deque()
    st = st[2:]
    entity = ''
    pre_symbol = ''
    for i, ch in enumerate(st):
        if ch == '=' or ch == ',' or ch == '>':
            if pre_symbol is not '>':
                entities.append(entity)
                entity = ''
        elif ch == '<':
            entities.append(ch)
        else:
            entity += ch
            if i == len(st)-1:
                entities.append(entity)
        if ch is ',' or ch is '>':
            value = entities.pop()
            key = entities.pop()
            entities.append({key:value})
        if ch is '>':
            values = {}
            value = entities.pop()
            while value is not '<':
                values.update(value)
                value = entities.pop()
            entities.append(values)
        if ch is not ' ':
            pre_symbol = ch
    return {entities[0]: entities[1]}


def load_demography_from_file(file, generation=GENERATION):
    """
        Load demography from file using msprime
        Input: smc++ output  (5 columns)
        ======================================================================
        |label	|x                  |y	                |plot_type	|plot_num|
        |--------------------------------------------------------------------|
        |ACB	|0.0	            |138482.84333082315	|path	    |0       |
        |ACB	|50.0	            |138482.84333082315	|path	    |0       |
        |ACB	|53.97505585700569  |139331.82583178935	|path	    |0       |
        ======================================================================
        Output: msprime.Demography object
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


def load_recombination_map_from_file(file, background_rate=1e-10):
    """
        Load recombination map from file.
        Input: a deeprho output formatted file, details as follows.
        If the intervals are not continuous, the uncovered region will be padded as background_rate.
        --------------------
        Start	End	Rate
        0	4000	1e-9
        4000	9000	1e-9
        9000	11000	1e-8
        11000	20000	1e-7
        20000	30000	1e-9
        30000	100000	1e-8
        ----------------------
    """
    check_file_existence(file)
    rate_map_file = pd.read_csv(file, sep='\t')
    positions = []
    rates = []
    for index, row in rate_map_file.iterrows():
        if index == 0:
            positions.append(row['Start'])
            positions.append(row['End'])
            rates.append(row['Rate'])
        else:
            start, end, rate= row['Start'], row['End'], row['Rate']
            if start == positions[-1]:
                positions.append(end)
                rates.append(rate)
            else:
                positions.append(start)
                rates.append(background_rate)
                positions.append(end)
                rates.append(rate)
    map = msp.RateMap(position=positions, rate=rates)
    return map












