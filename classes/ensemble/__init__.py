import numpy as np
from numba import jit_module

import constants
from .ensemble import *

clumpMass = 0

clumpNj = [[] for _ in range(len(constants.clumpMassNumber))]
clumpDeltaNji = [[] for _ in range(len(constants.clumpMassNumber))]
clumpNormalisedNj = [[] for _ in range(len(constants.clumpMassNumber))]
clumpNormalisedDeltaNji = [[] for _ in range(len(constants.clumpMassNumber))]
clumpSurfaceProbability = [[] for _ in range(len(constants.clumpMassNumber))]
clumpProbableNumber = [[] for _ in range(len(constants.clumpMassNumber))]
clumpStandardDeviation = [[] for _ in range(len(constants.clumpMassNumber))]

clumpNumberRange = [[] for _ in range(len(constants.clumpMassNumber))]

clumpCombinations = [[] for _ in range(len(constants.clumpMassNumber))]

clumpLargestCombination = [0 for _ in range(len(constants.clumpMassNumber))]
clumpLargestIndex = [0 for _ in range(len(constants.clumpMassNumber))]

clumpProbability = [[] for _ in range(len(constants.clumpMassNumber))]
CLmaxProbability = [[] for _ in range(len(constants.clumpMassNumber))]

clumpIndeces = [[] for _ in range(len(constants.clumpMassNumber))]

jit_module(nopython=False)
