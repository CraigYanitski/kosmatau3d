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
clump_t_gas = [np.zeros(constants.clump_mass_number[_]) for _ in range(len(constants.clump_mass_number))]
clump_t_dust = [np.zeros(constants.clump_mass_number[_]) for _ in range(len(constants.clump_mass_number))]
clump_hi_col_dens = [np.zeros(constants.clump_mass_number[_]) for _ in range(len(constants.clump_mass_number))]
clump_h2_col_dens = [np.zeros(constants.clump_mass_number[_]) for _ in range(len(constants.clump_mass_number))]
clump_hi_mass = [np.zeros(constants.clump_mass_number[_]) for _ in range(len(constants.clump_mass_number))]
clump_h2_mass = [np.zeros(constants.clump_mass_number[_]) for _ in range(len(constants.clump_mass_number))]

# KOSMA-tau outputs

clump_intensity = [np.zeros((constants.clump_mass_number[_],
                             len(species.molecules)+constants.wavelengths[constants.n_dust].size))
                   for _ in range(len(constants.clump_mass_number))]
clump_optical_depth = [np.zeros((constants.clump_mass_number[_],
                                 len(species.molecules)+constants.wavelengths[constants.n_dust].size))
                       for _ in range(len(constants.clump_mass_number))]
clump_species_intensity = [np.zeros((constants.clump_mass_number[_], len(species.molecules)))
                           for _ in range(len(constants.clump_mass_number))]
clump_species_optical_depth = [np.zeros((constants.clump_mass_number[_], len(species.molecules)))
                               for _ in range(len(constants.clump_mass_number))]
clump_dust_intensity = [np.zeros((constants.clump_mass_number[_],
                                  constants.wavelengths[constants.n_dust].size))
                        for _ in range(len(constants.clump_mass_number))]
clump_dust_optical_depth = [np.zeros((constants.clump_mass_number[_],
                                     constants.wavelengths[constants.n_dust].size))
                            for _ in range(len(constants.clump_mass_number))]
clump_hi_tb = [np.zeros(constants.clump_mass_number[_]) for _ in range(len(constants.clump_mass_number))]
clump_hi_tau = [np.zeros(constants.clump_mass_number[_]) for _ in range(len(constants.clump_mass_number))]

# jit_module(nopython=False)
