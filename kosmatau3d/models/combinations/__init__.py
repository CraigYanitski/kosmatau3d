import numpy as np
from numba import jit_module

from .combination import *
from kosmatau3d.models import constants

clump_combination = [[] for _ in range(len(constants.clump_mass_number))]
clump_max_combination = [[] for _ in range(len(constants.clump_mass_number))]

clump_species_intensity = [[] for _ in range(len(constants.clump_mass_number))]
clump_species_optical_depth = [[] for _ in range(len(constants.clump_mass_number))]
clump_dust_intensity = [[] for _ in range(len(constants.clump_mass_number))]
clump_dust_optical_depth = [[] for _ in range(len(constants.clump_mass_number))]

clump_hi_tb = [[] for _ in range(len(constants.clump_mass_number))]
clump_hi_tau = [[] for _ in range(len(constants.clump_mass_number))]

# jit_module(nopython=False)
