from numba import jit_module
import numpy as np

from .masspoint import *
from .. import constants

# Properties

logFUV = [0 for _ in range(len(constants.clumpMassNumber))]

clumpLogDensity = [[] for _ in range(len(constants.clumpMassNumber))]
clumpRadius = [[] for _ in range(len(constants.clumpMassNumber))]

# KOSMA-tau outputs

clumpIntensity = [[] for _ in range(len(constants.clumpMassNumber))]
clumpOpticalDepth = [[] for _ in range(len(constants.clumpMassNumber))]

jit_module(nopython=False)
