import numpy as np
from numba import jit_module

from .combination import *
from .. import constants

clumpCombination = [[] for _ in range(len(constants.clumpMassNumber))]

clumpIntensity = [[] for _ in range(len(constants.clumpMassNumber))]
clumpOpticalDepth = [[] for _ in range(len(constants.clumpMassNumber))]

jit_module(nopython=False)
