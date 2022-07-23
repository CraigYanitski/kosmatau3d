import numpy as np
from numba import jit_module

from .combination import *
from kosmatau3d.models import constants

clump_combination = [[] for _ in range(len(constants.clumpMassNumber))]
clump_max_combination = [[] for _ in range(len(constants.clumpMassNumber))]

clump_species_intensity = [[] for _ in range(len(constants.clumpMassNumber))]
clump_species_optical_depth = [[] for _ in range(len(constants.clumpMassNumber))]
clump_dust_intensity = [[] for _ in range(len(constants.clumpMassNumber))]
clump_dust_optical_depth = [[] for _ in range(len(constants.clumpMassNumber))]

# jit_module(nopython=False)
