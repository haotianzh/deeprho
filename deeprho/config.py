"""
    Global configurations in deeprho.
"""
import os


class CONFIG:
    """
        Default settings in deeprho.
        See https://github.com/haotianzh/deeprho_v2/blob/main/README.md
    """
    PKG_DIR = os.path.dirname(__file__)
    # global default settings
    NUM_THREAD = os.cpu_count() // 2    # Take a half cores
    GENERATION = 20 # Years per generation
    EFFECTIVE_POPULATION_SIZE = 5e4   # Diploid effective population size
    PLOIDY = 2  # Diploid as default
    SCALE_FACTOR = 10   # Normalizing factor for scaled rho to speed up training,\
                        # in training set and test set, rhos are divided by this factor

    # estimate.py
    THRESHOLD = 5e-8    # Hotspot threshold above which is regarded as recombination hotspot
    GLOBAL_WINDOW_SIZE = 1000   # Window size used for inferring local genealogies
    WINDOW_SIZE = 50    # Window size used for estimating rates (sliding window algorithm)
    STEP_SIZE = 25 # Step size (sliding window algorithm)
    RESOLUTION = 1e4    # Rate map resolution, deeprho will estimate rates at 10kb level
    LENGTH = None   # Genome length, None means it will be inferred from data
    MODEL_FINE = os.path.join(PKG_DIR, 'static', 'model_fine.h5')   # Fine model file path (should be sth like '.../models/your_model_fine_name.h5')
    MODEL_LARGE = os.path.join(PKG_DIR, 'static', 'model_large.h5')  # Large model file path (should be sth like '.../models/your_model_large_name.h5')
    CONSTANT_TABLE = os.path.join(PKG_DIR, 'static', 'table_0_1000.csv')
    TABLE = None

    # dpp.py
    DP_N_SAMPLE = 200  # Sampling points uniformly distributed between R_MIN and R_MAX
    DP_N_DRAW = 5  # Sampling times at a particular rate
    DP_N_POP = 50  # Number of individuals in simulation
    DP_MUTATION_RATE = 2.5e-8  # Mutation rate per site per generation
    DEMOGRAPHY = None   # Demography settings, see Docs
    DP_R_MIN = 1e-9    # Minimum recombination rate to simulate
    DP_R_MAX = 5e-7    # Maximum recombination rate to simulate



    # TODO test_provider.py
    # make_table.py
    MT_R_MIN = 0
    MT_R_MAX = 1e-6
    MT_REPEAT = 50
    MT_DRAW = 200
    MT_SAMPLE = 50






