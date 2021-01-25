from numba import jit_module
import numpy as np

from .masspoint import *
from .. import species
from .. import constants

# Properties

logFUV = [0 for _ in range(len(constants.clumpMassNumber))]

clumpLogDensity = [np.zeros(constants.clumpMassNumber[_]) for _ in range(len(constants.clumpMassNumber))]
clumpRadius = [np.zeros(constants.clumpMassNumber[_]) for _ in range(len(constants.clumpMassNumber))]

# KOSMA-tau outputs

clumpIntensity = [np.zeros((constants.clumpMassNumber[_], len(species.molecules)+constants.wavelengths[constants.nDust].size)) for _ in range(len(constants.clumpMassNumber))]
clumpOpticalDepth = [np.zeros((constants.clumpMassNumber[_], len(species.molecules)+constants.wavelengths[constants.nDust].size)) for _ in range(len(constants.clumpMassNumber))]

# jit_module(nopython=False)
