import numpy as np
from numba import jit_module
import importlib as il

from kosmatau3d.models import combinations
from kosmatau3d.models import constants
from kosmatau3d.models import masspoints

'''
This is a class to handle a combination of fractal masses in an ensemble.
It will have its associated probability, which will scale its intrinsic
intensity and optical depth. It returns a tuple of the combination's
probability, intensity, optical depth, and FUV field.
'''


def initialise(clumpcombination=[], interclumpcombination=[]):
    combinations.clumpCombination = clumpcombination
    return


def setClumpCombination(clumpcombination=[]):
    combinations.clumpCombination = clumpcombination
    return


def getAfuv(verbose=False):
    # The input is the density in the voxel. The probability should be included outside of this module.
    #  (velocity-independent)
    afuv = masspoints.getAfuv()
    for ens in range(len(afuv)):
        afuv[ens] = afuv[ens] * combinations.clumpCombination[ens]
    return [np.exp(-afuv[ens].sum(1)) for ens in range(len(afuv))]


def calculateEmission(test_calc=False, test_opacity=False, test_fv=False, f_v=None, probability=1, debug=False, test=False):
    '''
    This retrieves the emission from the masspoints, with dimensions (masspoints,species,velocity,velocity). It
    sums over the masspoint dimension. The probability will remain dormant until it is further-developed.
    '''
    intensitylist = [[] for _ in range(len(constants.clumpMassNumber))]
    opticaldepthlist = [[] for _ in range(len(constants.clumpMassNumber))]
  
    for ens in range(len(constants.clumpMassNumber)):

        if test_fv:
            f_ds = np.maximum(f_v[ens], 1)
        else:
            f_ds = 1
      
        for c in combinations.clumpCombination[ens]:
            if test_calc:
                intensitylist[ens].append((c * (masspoints.clumpIntensity[ens]*(masspoints.clumpOpticalDepth[ens]
                                                *np.pi*masspoints.clumpRadius[ens].T**2)
                                                /(1-np.exp(-masspoints.clumpOpticalDepth[ens]))
                                                /constants.voxel_size**3/f_ds).T).T)
            else:
                intensitylist[ens].append((c*masspoints.clumpIntensity[ens].T).T)
            if test_opacity:
                opticaldepthlist[ens].append((c * (masspoints.clumpOpticalDepth[ens]*np.pi
                                                   *masspoints.clumpRadius[ens].T**2/constants.voxel_size**3/f_ds).T).T)
            else:
                opticaldepthlist[ens].append((c*masspoints.clumpOpticalDepth[ens].T).T)
        
        combinations.clumpIntensity[ens] = np.array(intensitylist[ens]).sum(1)
        combinations.clumpOpticalDepth[ens] = np.array(opticaldepthlist[ens]).sum(1)
  
    return


def reinitialise():
    # Reinitialise all temporary variables to the correct number of clump sets.
  
    combinations.clumpCombination = [[] for _ in range(len(constants.clumpMassNumber))]
  
    combinations.clumpIntensity = [[] for _ in range(len(constants.clumpMassNumber))]
    combinations.clumpOpticalDepth = [[] for _ in range(len(constants.clumpMassNumber))]
  
    return


# jit_module(nopython=False)
