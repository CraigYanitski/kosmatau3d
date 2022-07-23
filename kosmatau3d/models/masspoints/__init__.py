import numpy as np
from numba import jit_module

from .masspoint import *
from kosmatau3d.models import species
from kosmatau3d.models import constants


# Properties

log_crir = 0
log_fuv = [0 for _ in range(len(constants.clump_mass_number))]

clump_log_density = [np.zeros(constants.clump_mass_number[_]) for _ in range(len(constants.clump_mass_number))]
clump_radius = [np.zeros(constants.clump_mass_number[_]) for _ in range(len(constants.clump_mass_number))]

# KOSMA-tau outputs

clump_intensity = [np.zeros((constants.clump_mass_number[_],
                             len(species.molecules)+constants.wavelengths[constants.n_dust].size))
                   for _ in range(len(constants.clumpMassNumber))]
clump_opticalDepth = [np.zeros((constants.clumpMassNumber[_],
                                len(species.molecules)+constants.wavelengths[constants.n_dust].size))
                      for _ in range(len(constants.clump_mass_number))]

# jit_module(nopython=False)
