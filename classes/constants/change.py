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
  constants.clumpMassNumber = num
  constants.clumpLogMass = np.linspace(constants.clumpLogMassRange[0], constants.clumpLogMassRange[-1], num=constants.clumpMassNumber)
  constants.clumpLogMass.resize(1,constants.clumpMassNumber)
  return

def changeClumpMassRange(massRange):
  constants.clumpLogMassRange = massRange
  constants.clumpLogMass = np.linspace(constants.clumpLogMassRange[0], constants.clumpLogMassRange[-1], num=constants.clumpMassNumber)
  constants.clumpLogMass.resize(1,constants.clumpMassNumber)
  return

def changeInterclumpMassNumber(num):
  constants.interclumpMassNumber = num
  constants.interclumpLogMass = np.linspace(constants.interclumpLogMassRange[0], constants.interclumpLogMassRange[-1], num=constants.interclumpMassNumber)
  constants.interclumpLogMass.resize(1,constants.interclumpMassNumber)
  return

def changeInterclumpMassRange(massRange):
  constants.interclumpLogMassRange = massRange
  constants.interclumpLogMass = np.linspace(constants.interclumpLogMassRange[0], constants.interclumpLogMassRange[-1], num=constants.interclumpMassNumber)
  constants.interclumpLogMass.resize(1,constants.interclumpMassNumber)
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
