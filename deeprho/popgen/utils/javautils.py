from ..base import Haplotype
import numpy as np
import jpype
import jpype.imports
from jpype.types import *
from labwu.rentplus.main import MPRentPlus
from labwu.py import RFDistance
from labwu.py import TripletDistance

def rentplus(haps, num_thread=1, infer_branch=False):
    matrices = []
    positions = []
    if isinstance(haps, Haplotype):
        matrices.append(haps.matrix)
        positions.append(haps.scaled_positions)
    elif isinstance(haps, list):
        for hap in haps:
            matrices.append(hap.matrix)
            positions.append(hap.scaled_positions)
    matrices = np.array(matrices, np.int)
    positions = np.array(positions, np.int)
    matrices_java = JInt[:][:][:](matrices)
    positions_java = JInt[:][:](positions)
    rent_res = MPRentPlus.pythonMultiProcess(matrices_java, positions_java, num_thread, infer_branch)
    if isinstance(haps, Haplotype):
        return rent_res[0]
    return rent_res


def rfdist(newicks, num_thread=10):
    dist = RFDistance.rfDistance(newicks, num_thread)
    return dist

def triplet_dist(newicks, num_thread=10):
    dist = TripletDistance.tripletDistance(newicks, num_thread)
    return dist
