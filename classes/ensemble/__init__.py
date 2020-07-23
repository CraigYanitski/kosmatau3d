import numpy as np
from numba import jit_module

from .ensemble import *

clumpMass = 0

clumpNj = []
clumpDeltaNji = []

clumpNumberRange = []

clumpCombinations = []

clumpLargestIndex = 0

clumpProbability = []
CLmaxProbability = []

clumpIndeces = []



interclumpMass = 0

interclumpNj = []
interclumpDeltaNji = []

interclumpNumberRange = []

interclumpCombinations = []

interclumpLargestIndex = 0

interclumpProbability = []
ICmaxProbability = []

interclumpIndeces = []

jit_module(nopython=False)
