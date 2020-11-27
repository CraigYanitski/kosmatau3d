import numpy as np
from numba import jit_module
import importlib as il

from .. import combinations

from .. import constants
from .. import masspoints

'''
This is a class to handle a combination of fractal masses in an ensemble.
It will have its associated probability, which will scale its intrinsic
intensity and optical depth. It returns a tuple of the combination's
probability, intensity, optical depth, and FUV field.
'''


def initialise(clumpCombination=[], interclumpCombination=[]):
  combinations.clumpCombination = clumpCombination
  return

def setClumpCombination(clumpCombination=[]):
  combinations.clumpCombination = clumpCombination
  return

def getAfuv(verbose=False):
  # The input is the density in the voxel. The probability should be included outside of this module.
  #  (velocity-independent)
  clAfuv = masspoints.getAfuv()
  for ens in range(len(clAfuv)):
    clAfuv[ens] = clAfuv[ens] * combinations.clumpCombination[ens]
  return [np.exp(-clAfuv[ens].sum(1)) for ens in range(len(clAfuv))]

def calculateEmission(probability=1, debug=False, test=False):
  '''
  This retrieves the emission from the masspoints, with dimensions (masspoints,species,velocity,velocity). It
  sums over the masspoint dimension. The probability will remain dormant until it is further-developed.
  '''
  CLintensityList = [[] for _ in range(len(constants.clumpMassNumber))]
  CLopticalDepthList = [[] for _ in range(len(constants.clumpMassNumber))]

  for ens in range(len(constants.clumpMassNumber)):
    
    for c in combinations.clumpCombination[ens]:
      CLintensityList[ens].append((c*masspoints.clumpIntensity[ens].T).T)
      CLopticalDepthList[ens].append((c*masspoints.clumpOpticalDepth[ens].T).T)
    
    combinations.clumpIntensity[ens] = np.array(CLintensityList[ens]).sum(1)
    combinations.clumpOpticalDepth[ens] = np.array(CLopticalDepthList[ens]).sum(1)

  return

def reinitialise():
  # Reinitialise all temporary variables to the correct number of clump sets.

  combinations.clumpCombination = [[] for _ in range(len(constants.clumpMassNumber))]

  combinations.clumpIntensity = [[] for _ in range(len(constants.clumpMassNumber))]
  combinations.clumpOpticalDepth = [[] for _ in range(len(constants.clumpMassNumber))]

  return

# jit_module(nopython=False)
