import os
import numpy as np
from copy import copy

from numba import jit_module
from numba.core.errors import NumbaDeprecationWarning, NumbaPendingDeprecationWarning
import warnings

from kosmatau3d.models import constants
from kosmatau3d.models import species

warnings.simplefilter('ignore', category=NumbaDeprecationWarning)
warnings.simplefilter('ignore', category=NumbaPendingDeprecationWarning)

'''
This is a script to contain all of the methods needed to change the model parameters.
'''


def change_velocity_number(num):
    num = np.max((num, 2))
    constants.velocity_number = copy(num)
    constants.velocity_range = np.linspace(constants.velocity_bin[0],
                                           constants.velocity_bin[-1],
                                           num=constants.velocity_number)
    constants.velocity_step = constants.velocity_range[1] - constants.velocity_range[0]
    return


def change_velocity_range(range):
    constants.velocity_bin = copy(range)
    constants.velocity_range = np.linspace(constants.velocity_bin[0],
                                          constants.velocity_bin[-1],
                                          num=constants.velocity_number)
    constants.velocity_step = constants.velocity_range[1] - constants.velocity_range[0]
    return


def change_clump_mass_number(num):
    # This will affect all clump sets, so the `num` needs to be a list with length of the number of clump sets.
    if isinstance(num, list) or isinstance(num, np.ndarray):
        constants.clump_mass_number = copy(num)
        constants.clump_log_mass = [[] for _ in range(len(num))]
        for i in range(len(num)):
            constants.clump_log_mass[i] = np.linspace(constants.clump_log_mass_range[i][0],
                                                      constants.clump_log_mass_range[i][-1],
                                                      num=constants.clump_mass_number[i])
            constants.clump_log_mass[i].resize(1, constants.clump_mass_number[i])
    return


def change_clump_mass_range(mass_range):
    # This will affect all clump sets, so the `massRange` needs to be a list with length of the number of clump sets.
    if isinstance(mass_range[0], list) or isinstance(mass_range[0], np.ndarray):
        constants.clump_log_mass_range = copy(mass_range)
        constants.clump_log_mass = [[] for _ in range(len(mass_range))]
        for i in range(len(mass_range)):
            constants.clump_log_mass[i] = np.linspace(min(constants.clump_log_mass_range[i]),
                                                      max(constants.clump_log_mass_range[i]),
                                                      num=constants.clump_mass_number[i])
            constants.clump_log_mass[i].resize(1, constants.clump_mass_number[i])
    return


def add_clumps(mass_range=[], num=0, n_max=1, reset=False):
    # Add another set of clumps to evaluate. Set the density kwarg if you do not want to use the voxel density.
    if reset:
        constants.clump_log_mass_range = []
        constants.clump_mass_number = []
        constants.clump_max_indeces = []
        constants.clump_n_max = []
        constants.clump_log_mass = []
    if isinstance(mass_range[0], int) or isinstance(mass_range[0], float):
        mass_range = [mass_range]
        num = [num]
        n_max = [n_max]
    for i in range(len(num)):
        constants.clump_log_mass_range.append(mass_range[i])
        constants.clump_mass_number.append(num[i])
        constants.clump_max_indeces.append(0)
        constants.clump_n_max.append(n_max[i])
        constants.clump_log_mass.append(np.resize(np.linspace(min(mass_range[i]),
                                                              max(mass_range[i]),
                                                              num=num[i]),
                                                (1,num[i])))
        constants.ensembles = len(num)
    return


def set_interclump_ensemble(idx):
    if isinstance(idx, int):
        i_cl = [False for _ in range(constants.ensembles)]
        i_cl[idx] = True
    elif isinstance(idx, (list, tuple, np.ndarray)):
        i_cl = list(idx)
    else:
        raise TypeError('interclump_idx must be type int or type list.')
    constants.interclump_idx = i_cl
    return


def set_interclump_wnm(idx):
    if isinstance(idx, int):
        i_cl = [False for _ in range(constants.ensembles)]
        i_cl[idx] = True
    elif isinstance(idx, (list, tuple, np.ndarray)):
        i_cl = list(idx)
    else:
        raise TypeError('interclump_idx must be type int or type list.')
    constants.interclump_wnm_idx = i_cl
    return


def reset_clumps():
    # This will restore the clump list to its default.
    constants.clump_mass_number = [3, 1]
    constants.clump_log_mass_range = [[0,2], [-2]]
    constants.clump_log_mass = [[], []]
    # constants.clumpDensity = [None, 1911]
    # constants.clumpFUV = [None, 10]
    for i in range(2):
        constants.clump_log_mass[i] = np.linspace(constants.clump_log_mass_range[i][0],
                                                  constants.clump_log_mass_range[i][-1],
                                                  num=constants.clump_mass_number[i])
        constants.clump_log_mass[i].resize(1,constants.clump_mass_number[i])
    return


def change_mass_function_parameters(alpha=1.84, gamma=2.31):
    # Use this to change the parameters of power-law distribution used
    #  to calculate the clump properties. The default are those from
    #  Heithausen et al. (1998).
    constants.alpha = alpha
    constants.gamma = gamma
    return


def change_directory(direc):

    constants.directory = direc
    if constants.directory[-1] != '/':
        constants.directory = constants.directory + '/'

    directory = constants.HISTORYPATH + constants.directory + constants.history

    # os.chmod(constants.HISTORYPATH, 0o777)
    if not os.path.exists(directory):
        os.makedirs(directory)

    return


def change_dust_wavelengths(limit='all'):
  
    # constants.dustWavelengths = limit
    constants.dust = limit
    
    # Only use the longest wavelength with roughly the same intensity for all the species transitions for a
    #  reduced dust model.
    if isinstance(limit, str):
        if limit == 'reduced':
            limit = ['3.1mm']
        elif limit in constants.dust_names:
            limit = [limit]
      
    # Check if individual wavelengths of the dust continuum are specified
    if isinstance(limit, list):
        constants.dust_wavelengths = copy(limit)
        n_dust = [line == np.asarray(constants.dust_names) for line in limit]
        constants.n_dust = np.any(n_dust, 0)
    # Use PAH to include the PAH features of the dust continuum
    elif limit == 'PAH':
        constants.n_dust = constants.wavelengths>constants.limit_pah
        constants.dust_wavelengths = constants.dust_names[constants.n_dust]
    # Use molecular to use just the section of the dust continuum relevant for the species transitions
    elif limit == 'molecular':
        constants.n_dust = constants.wavelengths>constants.limit_molecular
        constants.dust_wavelengths = constants.dust_names[constants.n_dust]
    # Otherwise include the entire dust continuum.
    else:
        constants.n_dust = constants.wavelengths>0
        constants.dust_wavelengths = constants.dust_names[constants.n_dust]
      
    return


def setup_species(species):
    constants.transitions = species
    constants.transition_number = len(species)
    return


def setup_interclump_species(species):
    constants.interclump_transitions = species
    constants.interclump_transition_number = len(species)
    return


def resort_wavelengths(interclump=False):
    if interclump:
        all_wavelengths = np.append(constants.wavelengths[constants.n_dust], species.interclump_transition_wavelengths)
    else:
        all_wavelengths = np.append(constants.wavelengths[constants.n_dust], species.clump_transition_wavelengths)
    constants.sorted_indeces = all_wavelengths.argsort()
    constants.sorted_wavelengths = all_wavelengths[constants.sorted_indeces]
    return


def dust_wavelengths():
    return constants.wavelengths[constants.n_dust]

# jit_module(nopython=False)
