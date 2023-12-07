"""
A subpackage containing information for the discrete clumps within a voxel.
This is a subpackage to appropriately handle all of the information concerning
a KOSMA-:math:`\\tau` clump.
The interpolated values are loaded into some of the subpackage variables.
"""


import numpy as np
import os

# from numba import jit_module


if "READTHEDOCS" not in os.environ:
    from .masspoint import *
from kosmatau3d.models import species
from kosmatau3d.models import constants


# Properties

log_crir = 0
log_fuv = [0 for _ in range(len(constants.clump_mass_number))]

clump_log_density = [
    np.zeros(constants.clump_mass_number[_])
    for _ in range(len(constants.clump_mass_number))
]
clump_log_density_orig = [
    np.zeros(constants.clump_mass_number[_])
    for _ in range(len(constants.clump_mass_number))
]
clump_radius = [
    np.zeros(constants.clump_mass_number[_])
    for _ in range(len(constants.clump_mass_number))
]
clump_t_gas = [
    np.zeros(constants.clump_mass_number[_])
    for _ in range(len(constants.clump_mass_number))
]
clump_t_dust = [
    np.zeros(constants.clump_mass_number[_])
    for _ in range(len(constants.clump_mass_number))
]
clump_hi_col_dens = [
    np.zeros(constants.clump_mass_number[_])
    for _ in range(len(constants.clump_mass_number))
]
clump_h2_col_dens = [
    np.zeros(constants.clump_mass_number[_])
    for _ in range(len(constants.clump_mass_number))
]
clump_hi_mass = [
    np.zeros(constants.clump_mass_number[_])
    for _ in range(len(constants.clump_mass_number))
]
clump_h2_mass = [
    np.zeros(constants.clump_mass_number[_])
    for _ in range(len(constants.clump_mass_number))
]
clump_N_species = [
    np.zeros((constants.clump_mass_number[_], len(constants.abundances)))
    for _ in range(len(constants.clump_mass_number))
]

# KOSMA-tau outputs

clump_intensity = [
    np.zeros(
        (
            constants.clump_mass_number[_],
            len(species.clump_transitions)
            + constants.wavelengths[constants.n_dust].size,
        )
    )
    for _ in range(len(constants.clump_mass_number))
]
clump_optical_depth = [
    np.zeros(
        (
            constants.clump_mass_number[_],
            len(species.clump_transitions)
            + constants.wavelengths[constants.n_dust].size,
        )
    )
    for _ in range(len(constants.clump_mass_number))
]
clump_species_intensity = [
    np.zeros((constants.clump_mass_number[_], len(species.clump_transitions)))
    for _ in range(len(constants.clump_mass_number))
]
clump_species_optical_depth = [
    np.zeros((constants.clump_mass_number[_], len(species.clump_transitions)))
    for _ in range(len(constants.clump_mass_number))
]
clump_dust_intensity = [
    np.zeros(
        (constants.clump_mass_number[_], constants.wavelengths[constants.n_dust].size)
    )
    for _ in range(len(constants.clump_mass_number))
]
clump_dust_optical_depth = [
    np.zeros(
        (constants.clump_mass_number[_], constants.wavelengths[constants.n_dust].size)
    )
    for _ in range(len(constants.clump_mass_number))
]
clump_hi_tb = [
    np.zeros(constants.clump_mass_number[_])
    for _ in range(len(constants.clump_mass_number))
]
clump_hi_tau = [
    np.zeros(constants.clump_mass_number[_])
    for _ in range(len(constants.clump_mass_number))
]

# jit_module(nopython=False)
