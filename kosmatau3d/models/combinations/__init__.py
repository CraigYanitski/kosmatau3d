"""
A subpackage to work with the different combinations of clumps contained within
a voxel.
This subpackage is to keep track of the combinations of clumps used
in the calculation of the ensemble-averaged emission, as well as the
associated emissivity and absorption.
"""


# import numpy as np
# from numba import jit_module
import os


if "READTHEDOCS" not in os.environ:
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
