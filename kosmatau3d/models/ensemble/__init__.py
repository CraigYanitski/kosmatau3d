"""
A subpackage containing information dealing with the probabilities for 
different combinations of clumps contained within a voxel.
This subpackage uses the ISM mass and discrete clump masses to compute the 
combinations of clumps that are likely in a line-of-sight.
"""


# import numpy as np
# from numba import jit_module
import os


if "READTHEDOCS" not in os.environ:
    from .ensemble import *
from kosmatau3d.models import constants


clumpMass = 0

clumpVelocities = [[] for _ in range(len(constants.clump_mass_number))]

clumpNj = [[] for _ in range(len(constants.clump_mass_number))]
clumpDeltaNji = [[] for _ in range(len(constants.clump_mass_number))]
clumpNormalisedNj = [[] for _ in range(len(constants.clump_mass_number))]
clumpNormalisedDeltaNji = [[] for _ in range(len(constants.clump_mass_number))]
clumpSurfaceProbability = [[] for _ in range(len(constants.clump_mass_number))]
clumpProbableNumber = [[] for _ in range(len(constants.clump_mass_number))]
clumpStandardDeviation = [[] for _ in range(len(constants.clump_mass_number))]
CLmaxSurfaceProbability = [[] for _ in range(len(constants.clump_mass_number))]
CLmaxProbableNumber = [[] for _ in range(len(constants.clump_mass_number))]
CLmaxStandardDeviation = [[] for _ in range(len(constants.clump_mass_number))]

clumpNumberRange = [[] for _ in range(len(constants.clump_mass_number))]
CLmaxNumberRange = [[] for _ in range(len(constants.clump_mass_number))]

clumpCombinations = [[] for _ in range(len(constants.clump_mass_number))]
CLmaxCombinations = [[] for _ in range(len(constants.clump_mass_number))]

clumpLargestCombination = [0 for _ in range(len(constants.clump_mass_number))]
clumpLargestIndex = [0 for _ in range(len(constants.clump_mass_number))]

clumpProbability = [[] for _ in range(len(constants.clump_mass_number))]
CLmaxProbability = [[] for _ in range(len(constants.clump_mass_number))]

clumpIndeces = [[] for _ in range(len(constants.clump_mass_number))]

# jit_module(nopython=False)
