import numpy as np
from numba import jit_module
import importlib as il

import combinations

import masspoints

'''
This is a class to handle a combination of fractal masses in an ensemble.
It will have its associated probability, which will scale its intrinsic
intensity and optical depth. It returns a tuple of the combination's
probability, intensity, optical depth, and FUV field.
'''


def initialise(clumpCombination=[], interclumpCombination=[]):
  combinations.clumpCombination = clumpCombination
  combinations.interclumpCombination = interclumpCombination
  return

def setClumpCombination(clumpCombination=[]):
  combinations.clumpCombination = clumpCombination
  return

def setInterlumpCombination(interclumpCombination=[]):
  combinations.interclumpCombination = interclumpCombination
  return

def getAfuv(verbose=False):
  # The input is the density in the voxel. The probability should be included outside of this module.
  #  (velocity-independent)
  clAfuv,icAfuv = masspoints.getAfuv()
  clAfuv = clAfuv * combinations.clumpCombination
  icAfuv = icAfuv * combinations.interclumpCombination
  if verbose:
    print('clump: {}'.format(clAfuv))
    print('interclump: {}'.format(icAfuv))
  return (np.exp(-clAfuv.sum(1)), np.exp(-icAfuv.sum(1)))#(self.__probability[self.__probability.nonzero()].prod(),self.__probability[self.__probability.nonzero()].prod()*np.exp(-Afuv))

def calculateEmission(probability=1, debug=False, test=False):
  '''
  This retrieves the emission from the masspoints, with dimensions (masspoints,species,velocity,velocity). It
  sums over the masspoint dimension. The probability will remain dormant until it is further-developed.
  '''
  if debug:
    print('combinations')
    print(combinations.clumpCombination)
    print(combinations.interclumpCombination)
    input()
  CLintensityList = []
  CLopticalDepthList = []
  ICintensityList = []
  ICopticalDepthList = []

  for c in combinations.clumpCombination:
    CLintensityList.append((c*masspoints.clumpIntensity.T).T)
    CLopticalDepthList.append((c*masspoints.clumpOpticalDepth.T).T)
  for c in combinations.interclumpCombination:
    ICintensityList.append((c*masspoints.interclumpIntensity.T).T)
    ICopticalDepthList.append((c*masspoints.interclumpOpticalDepth.T).T)
  if debug:
    print('Combination emissions')
    print('\n', CLintensityList[0], '\n\n', CLopticalDepthList[1])
    print('\n', ICintensityList[0], '\n\n', ICopticalDepthList[1])
    input()

  combinations.clumpIntensity = np.array(CLintensityList).sum(1)
  combinations.clumpOpticalDepth = np.array(CLopticalDepthList).sum(1)
  combinations.interclumpIntensity = np.array(ICintensityList).sum(1)
  combinations.interclumpOpticalDepth = np.array(ICopticalDepthList).sum(1)

  combinations.clumpIntensity[combinations.clumpIntensity<=0] = 1e-100
  combinations.clumpOpticalDepth[combinations.clumpOpticalDepth<=0] = 1e-100
  combinations.interclumpIntensity[combinations.interclumpIntensity<=0] = 1e-100
  combinations.interclumpOpticalDepth[combinations.interclumpOpticalDepth<=0] = 1e-100

  return

jit_module(nopython=False)
