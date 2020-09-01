from numba import jit_module
import numpy as np

from numba.errors import NumbaDeprecationWarning, NumbaPendingDeprecationWarning
import warnings
warnings.simplefilter('ignore', category=NumbaDeprecationWarning)
warnings.simplefilter('ignore', category=NumbaPendingDeprecationWarning)

import constants

import species
'''
This is a script to contain all of the methods needed to change the model parameters.
'''
def changeVelocityNumber(num):
  num = np.max((num, 2))
  constants.velocityNumber = num
  constants.velocityRange = np.linspace(constants.velocityBin[0], constants.velocityBin[-1], num=constants.velocityNumber)
  constants.velocityStep = constants.velocityRange[1] - constants.velocityRange[0]
  return

def changeVelocityRange(range):
  constants.velocityBin = range
  constants.velocityRange = np.linspace(constants.velocityBin[0], constants.velocityBin[-1], num=constants.velocityNumber)
  constants.velocityStep = constants.velocityRange[1] - constants.velocityRange[0]
  return

def changeClumpMassNumber(num):
  # This will affect all clump sets, so the `num` needs to be a list with length of the number of clump sets.
  if isinstance(num, list) or isinstance(num, np.ndarray):
    constants.clumpMassNumber = num
    constants.clumpLogMass = [[] for _ in range(len(num))]
    for i in range(len(num)):
      constants.clumpLogMass[i] = np.linspace(constants.clumpLogMassRange[i][0], constants.clumpLogMassRange[i][-1], num=constants.clumpMassNumber[i])
      constants.clumpLogMass[i].resize(1,constants.clumpMassNumber[i])
  return

def changeClumpMassRange(massRange):
  # This will affect all clump sets, so the `massRange` needs to be a list with length of the number of clump sets.
  if isinstance(massRange[0], list) or isinstance(massRange[0], np.ndarray):
    constants.clumpLogMassRange = massRange
    constants.clumpLogMass = [[] for _ in range(len(massRange))]
    for i in range(len(massRange)):
      constants.clumpLogMass[i] = np.linspace(constants.clumpLogMassRange[i][0], constants.clumpLogMassRange[i][-1], num=constants.clumpMassNumber[i])
      constants.clumpLogMass[i].resize(1,constants.clumpMassNumber[i])
  return

def addClumps(massRange=[], num=0, density=None, fuv=None, Nmax=1, reset=False):
  # Add another set of clumps to evaluate. Set the density kwarg if you do not want to use the voxel density.
  if reset:
    constants.clumpLogMassRange = []
    constants.clumpMassNumber = []
    constants.clumpDensity = []
    constants.clumpFUV = []
    constants.clumpMaxIndeces = []
    constants.clumpNmax = []
    constants.clumpLogMass = []
  if isinstance(massRange[0], int) or isinstance(massRange[0], float):
    massRange = [massRange]
    num = [num]
    density = [density]
    fuv = [fuv]
    Nmax = [Nmax]
  for i in range(len(num)):
    constants.clumpLogMassRange.append(massRange[i])
    constants.clumpMassNumber.append(num[i])
    constants.clumpDensity.append(density[i])
    constants.clumpFUV.append(fuv[i])
    constants.clumpMaxIndeces.append(0)
    constants.clumpNmax.append(Nmax[i])
    constants.clumpLogMass.append(np.resize(np.linspace(massRange[i][0], massRange[i][-1], num=num[i]), (1,num[i])))
  return

def resetClumps():
  # This will restore the clump list to its default.
  constants.clumpMassNumber = [3, 1]
  constants.clumpLogMassRange = [[0,2], [-2]]
  constants.clumpLogMass = [[], []]
  constants.clumpDensity = [None, 1911]
  constants.clumpFUV = [None, 10]
  for i in range(2):
    constants.clumpLogMass[i] = np.linspace(constants.clumpLogMassRange[i][0], constants.clumpLogMassRange[i][-1], num=constants.clumpMassNumber[i])
    constants.clumpLogMass[i].resize(1,constants.clumpMassNumber[i])
  return

def changeDirectory(direc):
  constants.directory = direc
  if constants.directory[-1]!='/': constants.directory = constants.directory + '/'
  return

def changeDustWavelengths(limit=''):
  constants.dustWavelengths = limit
  if limit=='PAH':
    constants.nDust = constants.wavelengths>constants.limitPAH
  elif limit=='molecular':
    constants.nDust = constants.wavelengths>constants.limitMolecular
  elif limit=='':
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

jit_module(nopython=False)
