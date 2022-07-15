"""
    author: Haotian
    created_at: 12/2/2021
    description: loading data from different file types and systems, especially for "ms-format" and "vcf-format"
"""
import os
import logging
from deeprho import LazyLoader
from deeprho.config import CONFIG
import numpy as np
from collections import deque, namedtuple
from ..base import Haplotype
from .statistics import _calculate_average_ne
pd = LazyLoader('pandas')
msp = LazyLoader('msprime')


use_pyvcf = False
logger = logging.getLogger(__name__)
try:
    import vcf
    use_pyvcf = True
except:
    logger.warning('PyVCF module not installed, use builtin VCF reader.')


def check_file_existence(file):
    assert os.path.exists(file), f'file {file} does not exist.'


def load_ms_from_file(file, true_genealogy=False):
    """
        Load haplotype data from ms-formatted file.
        Input: file path, true_genealogy = True if "-T" option is specified in simulation. (currently not implemented)
        Return: Haplotype
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
        Read haplotype data from VCF file in two ways, 1. use PyVCF module 2. or use built-in method
        Input: filename
        Return: Haplotype
    """
    if use_pyvcf:
        return _load_vcf_from_file_1(file)
    return _load_vcf_from_file_2(file)


def _load_vcf_from_file_1(file):
    """
        Use PyVCF module for reading VCF, Issues sometimes happen during pip install.
        Input: filename
        Return: Haplotype
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


def _load_vcf_from_file_2(file):
    """
        Built-in module for reading VCF. Missing sites and unphased sites are filtered.
        Input: filename
        Return: Haplotype
    """
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
        Return: namedtuplle(metadata, data)
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
            info = _parse_vcf_metadata(line)
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


def _parse_vcf_metadata(st):
    """
        Parse string like:
        "##INFO=<ID=ID,Number=number,Type=type,Description="description",Source="source",Version="version">"
        Input: a string
        Return: a json string
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


def load_demography_from_file(file, mode='year', generation=CONFIG.GENERATION):
    """
        Load demography from file using msprime
        # xxx_pop_sizes.csv
        --------------------------------------------
        label   x       y       plot_type   plot_num
        ACB     0.0     138482  path        0
        ACB     50.0    138482  path        0
        ACB     53.9    139331  path        0
        -------------------/------------------------
        Input: smc++ output, mode='year' or 'generation', years per generation
        Return: msprime.Demography object
    """
    check_file_existence(file)
    modes = ['year', 'generation']
    assert mode in modes, "mode should be either 'year' or 'generation'."
    demography_file = pd.read_csv(file)
    assert demography_file['label'].unique().size == 1, 'there are more than one population.'
    name = demography_file['label'][0]
    times = []
    sizes = []
    for i, row in demography_file.iterrows():
        times.append(row['x'] / (generation if mode == 'year' else 1))
        sizes.append(row['y'])
    return _load_demography(name, times, sizes)


def _load_demography(name, times, sizes):
    demography = msp.Demography()
    demography.add_population(name=name, initial_size=0)
    for time, size in zip(times, sizes):
        demography.add_population_parameters_change(time=time, initial_size=size)
    return demography


def parse_ms_command(ms_command, initial_size=CONFIG.EFFECTIVE_POPULATION_SIZE):
    """
        Parse ms command, especially, demographic events (In ms, time is scaled by 4N_0)
        Example:
            ms 2 100 -t 81960 -r 13560 30000000 -eN 0.01 0.05 -eN 0.0375 0.5 -eN 1.25 1
        Input: a ms command string
        Return: a dictionary of parameters
    """
    parameters = {'-t': 1, '-r': 2, '-eN': 2}
    ms_to_msprime = {'-t': 'rate', '-r': 'recombination_rate', '-eN': 'demography'}
    commands = ms_command.strip().split()
    program_name = commands[0]
    n_sam = commands[1]
    n_rep = commands[2]
    i = 3
    parse_result = {}
    while i < len(commands):
        command = commands[i]
        args = []
        try:
            for i in range(i+1, i+parameters[command]+1):
                args.append(command[i])
            parse_result[ms_to_msprime[command]] = args
            i += 1
        except KeyError:
            raise Exception('no such arg.')
        except IndexError:
            raise Exception('parse error.')
    sequence_length = parse_result['recombination_rate'][1]
    recombination_rate = parse_result['recombination_rate'][0] / sequence_length
    times = []
    sizes = []
    for i in range(0, len(parse_result['demography']), 2):
        times.append(parse_result['demography'][i] * 4 * initial_size)
        sizes.append(parse_result['demography'][i+1] * initial_size)
    res = {}
    res['program']






def load_recombination_map_from_file(file, background_rate=1e-10):
    """
        Load recombination map from file.
        If the intervals are not consecutive, the uncovered region will be padded as *background_rate*.
        # xxx_recombination_map.txt
        --------------------
        Start	End	    Rate
        0	4000	1e-9
        4000	9000	1e-9
        9000	11000	1e-8
        11000	20000	1e-7
        20000	30000	1e-9
        30000	100000	1e-8
        ---------------------
        Input: a deeprho output formatted file, details shown as above
        Return: a msprime RateMap
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
