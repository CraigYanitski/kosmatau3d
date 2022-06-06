import numpy as np
from copy import copy
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


def initialise(clumpcombination=[], totalcombination=[]):
    combinations.clumpCombination = clumpcombination
    combinations.CLmaxCombination = totalcombination
    return


def setClumpCombination(clumpcombination=[]):
    combinations.clumpCombination = clumpcombination
    return


def getAfuv(verbose=False):
    # The input is the density in the voxel. The probability should be included outside of this module.
    #  (velocity-independent)
    afuv = masspoints.getAfuv()
    for ens in range(len(afuv)):
        afuv[ens] = afuv[ens] * combinations.CLmaxCombination[ens]
    return [np.exp(-afuv[ens].sum(1)) for ens in range(len(afuv))]


def calculateEmission(test_calc=False, test_opacity=False, test_fv=False, f_v=None, 
                      suggested_calc=True, probability=1, old_dust=False, debug=False, test=False):
    '''
    This retrieves the emission from the masspoints, with dimensions (masspoints,species,velocity,velocity). It
    sums over the masspoint dimension. The probability will remain dormant until it is further-developed.
    '''
    intensitylist = [[] for _ in range(len(constants.clumpMassNumber))]
    opticaldepthlist = [[] for _ in range(len(constants.clumpMassNumber))]
    dustintensitylist = [[] for _ in range(len(constants.clumpMassNumber))]
    dustopticaldepthlist = [[] for _ in range(len(constants.clumpMassNumber))]
  
    for ens in range(len(constants.clumpMassNumber)):

        if test_fv:
            f_ds = np.maximum(f_v[ens], 1)
        else:
            f_ds = 1

        # The intensity and optical depth are arranged as the dust followed by the chemical transitions
        iDust = constants.wavelengths[constants.nDust].size

        for c in combinations.clumpCombination[ens]:
            if suggested_calc:
                intensity = (c * (masspoints.clumpIntensity[ens][:, iDust:]*(masspoints.clumpOpticalDepth[ens][:, iDust:])
                                  /(1-np.exp(-masspoints.clumpOpticalDepth[ens][:, iDust:]))
                                  ).T).T#/constants.voxel_size/f_ds*(4/3*masspoints.clumpRadius[ens].T)
                i_nan = np.isnan(intensity) | np.isinf(intensity)
                intensity[i_nan] = ((c * (masspoints.clumpIntensity[ens][:, iDust:]
                                          ).T).T)[i_nan]#/constants.voxel_size/f_ds*(4/3*masspoints.clumpRadius[ens].T)
                intensitylist[ens].append(copy(intensity))
                # intensitylist[ens].append((c * (masspoints.clumpIntensity[ens]*(masspoints.clumpOpticalDepth[ens])
                #                                 /(1-np.exp(-masspoints.clumpOpticalDepth[ens]))
                #                                 /(4/3*masspoints.clumpRadius[ens].T)).T).T)
            elif test_calc:
                intensity = (c * (masspoints.clumpIntensity[ens][:, iDust:]*masspoints.clumpOpticalDepth[ens][:, iDust:]
                                  /(1-np.exp(-masspoints.clumpOpticalDepth[ens][:, iDust:]))
                                  /constants.voxel_size/f_ds).T).T
                i_nan = np.isnan(intensity) | np.isinf(intensity)
                intensity[i_nan] = ((c * masspoints.clumpIntensity[ens][:, iDust:].T).T)[i_nan]
                intensitylist[ens].append(copy(intensity))
                # intensitylist[ens].append((c * (masspoints.clumpIntensity[ens]*(masspoints.clumpOpticalDepth[ens])
                #                                 /(1-np.exp(-masspoints.clumpOpticalDepth[ens]))
                #                                 /constants.voxel_size/f_ds).T).T)
            else:
                intensitylist[ens].append((c*masspoints.clumpIntensity[ens][:, iDust:].T).T)
            if suggested_calc:
                opticaldepthlist[ens].append((c * (masspoints.clumpOpticalDepth[ens][:, iDust:]
                                                   ).T).T)#/(4/3*masspoints.clumpRadius[ens].T)
            elif test_opacity:
                opticaldepthlist[ens].append((c * (masspoints.clumpOpticalDepth[ens][:, iDust:]
                                                   /constants.voxel_size/f_ds).T).T)
            else:
                opticaldepthlist[ens].append((c*masspoints.clumpOpticalDepth[ens][:, iDust:].T).T)
        combinations.clumpIntensity[ens] = np.array(intensitylist[ens]).sum(1)
        combinations.clumpOpticalDepth[ens] = np.array(opticaldepthlist[ens]).sum(1)

        if constants.dust != '' and constants.dust != None and constants.dust != 'none':
            if old_dust:
                CLcombinations = copy(combinations.clumpCombination[ens])
            else:
                CLcombinations = copy(combinations.CLmaxCombination[ens])
            for c in CLcombinations:
                if suggested_calc:
                    intensity = (c * (masspoints.clumpIntensity[ens][:, :iDust]*(masspoints.clumpOpticalDepth[ens][:, :iDust])
                                      /(1-np.exp(-masspoints.clumpOpticalDepth[ens][:, :iDust]))
                                      ).T).T#/constants.voxel_size/f_ds*(4/3*masspoints.clumpRadius[ens].T)
                    i_nan = np.isnan(intensity) | np.isinf(intensity)
                    intensity[i_nan] = ((c * (masspoints.clumpIntensity[ens][:, :iDust]
                                              ).T).T)[i_nan]#/constants.voxel_size/f_ds*(4/3*masspoints.clumpRadius[ens].T)
                    dustintensitylist[ens].append(copy(intensity))
                elif test_calc:
                    intensity = (c * (masspoints.clumpIntensity[ens][:, :iDust]*masspoints.clumpOpticalDepth[ens][:, :iDust]
                                      /(1-np.exp(-masspoints.clumpOpticalDepth[ens][:, :iDust]))
                                      /constants.voxel_size/f_ds).T).T
                    i_nan = np.isnan(intensity) | np.isinf(intensity)
                    intensity[i_nan] = ((c * masspoints.clumpIntensity[ens][:, :iDust].T).T)[i_nan]
                    dustintensitylist[ens].append(copy(intensity))
                    # intensitylist[ens].append((c * (masspoints.clumpIntensity[ens]*(masspoints.clumpOpticalDepth[ens])
                    #                                 /(1-np.exp(-masspoints.clumpOpticalDepth[ens]))
                    #                                 /constants.voxel_size/f_ds).T).T)
                else:
                    dustintensitylist[ens].append((c*masspoints.clumpIntensity[ens][:, :iDust].T).T)
                if suggested_calc:
                    dustopticaldepthlist[ens].append((c * (masspoints.clumpOpticalDepth[ens][:, :iDust]
                                                       ).T).T)#/(4/3*masspoints.clumpRadius[ens].T)
                elif test_opacity:
                    dustopticaldepthlist[ens].append((c * (masspoints.clumpOpticalDepth[ens][:, :iDust]
                                                       /constants.voxel_size/f_ds).T).T)
                else:
                    dustopticaldepthlist[ens].append((c*masspoints.clumpOpticalDepth[ens][:, :iDust].T).T)
            combinations.clumpDustIntensity[ens] = np.array(dustintensitylist[ens]).sum(1)
            combinations.clumpDustOpticalDepth[ens] = np.array(dustopticaldepthlist[ens]).sum(1)
        
  
    return


def reinitialise():
    # Reinitialise all temporary variables to the correct number of clump sets.
  
    combinations.clumpCombination = [[] for _ in range(len(constants.clumpMassNumber))]
    combinations.CLmaxCombination = [[] for _ in range(len(constants.clumpMassNumber))]
  
    combinations.clumpIntensity = [[] for _ in range(len(constants.clumpMassNumber))]
    combinations.clumpOpticalDepth = [[] for _ in range(len(constants.clumpMassNumber))]
    combinations.clumpDustIntensity = [[] for _ in range(len(constants.clumpMassNumber))]
    combinations.clumpDustOpticalDepth = [[] for _ in range(len(constants.clumpMassNumber))]
  
    return


# jit_module(nopython=False)
