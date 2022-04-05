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


def changeVelocityNumber(num):
    num = np.max((num, 2))
    constants.velocityNumber = copy(num)
    constants.velocityRange = np.linspace(constants.velocityBin[0],
                                          constants.velocityBin[-1],
                                          num=constants.velocityNumber)
    constants.velocityStep = constants.velocityRange[1] - constants.velocityRange[0]
    return


def changeVelocityRange(range):
    constants.velocityBin = copy(range)
    constants.velocityRange = np.linspace(constants.velocityBin[0],
                                          constants.velocityBin[-1],
                                          num=constants.velocityNumber)
    constants.velocityStep = constants.velocityRange[1] - constants.velocityRange[0]
    return


def changeClumpMassNumber(num):
    # This will affect all clump sets, so the `num` needs to be a list with length of the number of clump sets.
    if isinstance(num, list) or isinstance(num, np.ndarray):
        constants.clumpMassNumber = copy(num)
        constants.clumpLogMass = [[] for _ in range(len(num))]
        for i in range(len(num)):
            constants.clumpLogMass[i] = np.linspace(constants.clumpLogMassRange[i][0],
                                                    constants.clumpLogMassRange[i][-1],
                                                    num=constants.clumpMassNumber[i])
            constants.clumpLogMass[i].resize(1, constants.clumpMassNumber[i])
    return


def changeClumpMassRange(massRange):
    # This will affect all clump sets, so the `massRange` needs to be a list with length of the number of clump sets.
    if isinstance(massRange[0], list) or isinstance(massRange[0], np.ndarray):
        constants.clumpLogMassRange = copy(massRange)
        constants.clumpLogMass = [[] for _ in range(len(massRange))]
        for i in range(len(massRange)):
            constants.clumpLogMass[i] = np.linspace(min(constants.clumpLogMassRange[i]),
                                                    max(constants.clumpLogMassRange[i]),
                                                    num=constants.clumpMassNumber[i])
            constants.clumpLogMass[i].resize(1, constants.clumpMassNumber[i])
    return


def addClumps(massRange=[], num=0, Nmax=1, reset=False):
    # Add another set of clumps to evaluate. Set the density kwarg if you do not want to use the voxel density.
    if reset:
        constants.clumpLogMassRange = []
        constants.clumpMassNumber = []
        constants.clumpMaxIndeces = []
        constants.clumpNmax = []
        constants.clumpLogMass = []
    if isinstance(massRange[0], int) or isinstance(massRange[0], float):
        massRange = [massRange]
        num = [num]
        Nmax = [Nmax]
    for i in range(len(num)):
        constants.clumpLogMassRange.append(massRange[i])
        constants.clumpMassNumber.append(num[i])
        constants.clumpMaxIndeces.append(0)
        constants.clumpNmax.append(Nmax[i])
        constants.clumpLogMass.append(np.resize(np.linspace(min(massRange[i]),
                                                            max(massRange[i]),
                                                            num=num[i]),
                                                (1,num[i])))
        constants.ensembles = len(num)
    return


def resetClumps():
    # This will restore the clump list to its default.
    constants.clumpMassNumber = [3, 1]
    constants.clumpLogMassRange = [[0,2], [-2]]
    constants.clumpLogMass = [[], []]
    # constants.clumpDensity = [None, 1911]
    # constants.clumpFUV = [None, 10]
    for i in range(2):
        constants.clumpLogMass[i] = np.linspace(constants.clumpLogMassRange[i][0],
                                                constants.clumpLogMassRange[i][-1],
                                                num=constants.clumpMassNumber[i])
        constants.clumpLogMass[i].resize(1,constants.clumpMassNumber[i])
    return


def changeMassFunctionParameters(alpha=1.84, gamma=2.31):
    # Use this to change the parameters of power-law distribution used
    #  to calculate the clump properties. The default are those from
    #  Heithausen et al. (1998).
    constants.alpha = alpha
    constants.gamma = gamma
    return


def changeDirectory(direc):

    constants.directory = direc
    if constants.directory[-1] != '/':
        constants.directory = constants.directory + '/'

    directory = constants.HISTORYPATH + constants.directory + constants.history

    # os.chmod(constants.HISTORYPATH, 0o777)
    if not os.path.exists(directory):
        os.makedirs(directory)

    return


def changeDustWavelengths(limit=''):
  
    constants.dustWavelengths = limit
    
    # Only use the longest wavelength with roughly the same intensity for all the species transitions for a
    #  reduced dust model.
    if isinstance(limit, str):
        if limit == 'reduced':
            limit = ['3.1mm']
        elif limit in constants.dustNames:
            limit = [limit]
      
    # Check if individual wavelengths of the dust continuum are specified
    if isinstance(limit, list):
        nDust = [line == np.asarray(constants.dustNames) for line in limit]
        constants.nDust = np.any(nDust, 0)
    # Use PAH to include the PAH features of the dust continuum
    elif limit == 'PAH':
        constants.nDust = constants.wavelengths>constants.limitPAH
    # Use molecular to use just the section of the dust continuum relevant for the species transitions
    elif limit == 'molecular':
        constants.nDust = constants.wavelengths>constants.limitMolecular
    # Otherwise include the entire dust continuum.
    else:
        constants.nDust = constants.wavelengths>0
      
    return


def setupMolecules(species):
    constants.molecules = species
    constants.moleculeNumber = len(species)
    return


def resortWavelengths():
    allWavelengths = np.append(constants.wavelengths[constants.nDust], species.moleculeWavelengths)
    constants.sortedIndeces = allWavelengths.argsort()
    constants.sortedWavelengths = allWavelengths[constants.sortedIndeces]
    return


def dustWavelengths():
    return constants.wavelengths[constants.nDust]

# jit_module(nopython=False)
