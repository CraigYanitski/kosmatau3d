import numpy as np
from numba import jit_module

from .masspoint import *
from kosmatau3d.models import species
from kosmatau3d.models import constants


# Properties


logFUV = [0 for _ in range(len(constants.clumpMassNumber))]

clumpLogDensity = [np.zeros(constants.clumpMassNumber[_]) for _ in range(len(constants.clumpMassNumber))]
clumpRadius = [np.zeros(constants.clumpMassNumber[_]) for _ in range(len(constants.clumpMassNumber))]

# KOSMA-tau outputs

clumpIntensity = [np.zeros((constants.clumpMassNumber[_],
                            len(species.molecules)+constants.wavelengths[constants.nDust].size))
                  for _ in range(len(constants.clumpMassNumber))]
clumpOpticalDepth = [np.zeros((constants.clumpMassNumber[_],
                               len(species.molecules)+constants.wavelengths[constants.nDust].size))
                     for _ in range(len(constants.clumpMassNumber))]

# jit_module(nopython=False)
