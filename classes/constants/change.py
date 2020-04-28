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
  constants.velocityRange = np.linspace(constants.velocityBin[0], constants.velocityBin[1], num=constants.velocityNumber)
  velocityStep = velocityRange[1] - velocityRange[0]
  return

def changeVelocityRange(range):
  constants.velocityBins = range
  constants.velocityRange = np.linspace(constants.velocityBin[0], constants.velocityBin[1], num=constants.velocityNumber)
  velocityStep = velocityRange[1] - velocityRange[0]
  return

def changeClumpMassNumber(num):
  constants.clumpmassNumber = num
  constants.clumpMass = np.linspace(constants.clumpMassRange[0], constants.clumpMassRange[1], num=constants.clumpMassNumber)
  return

def changeClumpMassRange(range):
  constants.clumpMassRange = range
  constants.clumpMass = np.linspace(constants.clumpMassRange[0], constants.clumpMassRange[1], num=constants.clumpMassNumber)
  return

def changeInterclumpMassNumber(num):
  constants.interclumpMassNumber = num
  constants.interclumpMassRange = np.linspace(constants.interclumpMassRange[0], constants.interclumpMassRange[1], num=constants.interclumpMassNumber)
  return

def changeInterclumpMassRange(range):
  constants.interclumpMassRange = range
  constants.interclumpMass = np.linspace(constants.interclumpMassRange[0], constants.interclumpMassRange[1], num=constants.interclumpMassNumber)
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

def resortWavelengths():
  allWavelengths = np.append(constants.wavelengths, species.moleculeWavelengths)
  constants.sortedIndeces = allWavelengths.argsort()
  constants.sortedWavelengths = allWavelengths[constants.sortedIndeces]

jit_module(nopython=False)
