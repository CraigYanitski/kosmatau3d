"""
This is a module to handle a combination of fractal masses in an ensemble.
It will have its associated probability, which will scale its intrinsic
intensity and optical depth. It returns a tuple of the combination's
probability, intensity, optical depth, and far-UV field.
"""


import numpy as np
from copy import copy
#from numba import jit_module
#import importlib as il

from kosmatau3d.models import combinations
from kosmatau3d.models import constants
from kosmatau3d.models import masspoints


def initialise(clump_combination=[], total_combination=[]):
    '''
    This function assigns the combination data from :code:`models.ensemble` to variables in
    this subpackage. It needs to be executed after :code:`models.ensemble.ensemble.initialise()`.
    '''
    combinations.clump_combination = clump_combination
    combinations.clump_max_combination = total_combination
    return


def set_clump_combination(clump_combination=[]):
    combinations.clump_combination = clump_combination
    return


def get_fuv_extinction(verbose=False):
    '''
    The input is the density in the voxel. 
    The probability should be included outside of this module.

    (velocity-independent)
    '''
    afuv = masspoints.get_taufuv()
    for ens in range(len(afuv)):
        afuv[ens] = afuv[ens] * combinations.clump_max_combination[ens]
    return [np.exp(-afuv[ens].sum(1)) for ens in range(len(afuv))]


def calculate_emission(test_calc=False, test_opacity=False, test_fv=False, f_v=None, 
                       suggested_calc=True, probability=1, old_dust=False, debug=False, test=False):
    '''
    This retrieves the emission from the masspoints, with dimensions (masspoints,species,velocity,velocity). It
    sums over the masspoint dimension. The probability will remain dormant until it is further-developed.
    '''
    species_intensity_list = [[] for _ in range(len(constants.clump_mass_number))]
    species_opticaldepth_list = [[] for _ in range(len(constants.clump_mass_number))]
    dust_intensity_list = [[] for _ in range(len(constants.clump_mass_number))]
    dust_opticaldepth_list = [[] for _ in range(len(constants.clump_mass_number))]
    hi_tb_list = [[] for _ in range(len(constants.clump_mass_number))]
    hi_tau_list = [[] for _ in range(len(constants.clump_mass_number))]
  
    for ens in range(len(constants.clump_mass_number)):

        if test_fv:
            f_ds = np.maximum(f_v[ens], 1)
        else:
            f_ds = 1

        # The intensity and optical depth are arranged as the dust followed by the chemical transitions
        i_dust = constants.wavelengths[constants.n_dust].size

        for c in combinations.clump_combination[ens]:
            if suggested_calc:
                intensity = (c * (masspoints.clump_species_intensity[ens]
                                  *(masspoints.clump_species_optical_depth[ens])
                                  /(1-np.exp(-masspoints.clump_species_optical_depth[ens]))
                                  ).T
                             ).T# / (4/3*masspoints.clump_radius[ens].T)
                            #/constants.voxel_size/f_ds*(4/3*masspoints.clump_radius[ens].T)
                hi_tb_list[ens].append((c * (masspoints.clump_hi_tb[ens]*(masspoints.clump_hi_tau[ens])
                                        / (1-np.exp(-masspoints.clump_hi_tau[ens]))).T).T)
                i_nan = np.isnan(intensity) | np.isinf(intensity)
                intensity[i_nan] = ((c * (masspoints.clump_species_intensity[ens]
                                          ).T).T)[i_nan]
                                   #/ (4/3*masspoints.clump_radius[ens].T)
                                   #/constants.voxel_size/f_ds*(4/3*masspoints.clumpRadius[ens].T)
                species_intensity_list[ens].append(copy(intensity))
                # intensitylist[ens].append((c * (masspoints.clumpIntensity[ens]*(masspoints.clumpOpticalDepth[ens])
                #                                 /(1-np.exp(-masspoints.clumpOpticalDepth[ens]))
                #                                 /(4/3*masspoints.clumpRadius[ens].T)).T).T)
            elif test_calc:
                intensity = (c * (masspoints.clump_species_intensity[ens]*masspoints.clump_species_optical_depth[ens]
                                  /(1-np.exp(-masspoints.clump_species_optical_depth[ens]))
                                  /constants.voxel_size/f_ds).T).T
                i_nan = np.isnan(intensity) | np.isinf(intensity)
                intensity[i_nan] = ((c * masspoints.clump_species_intensity[ens].T).T)[i_nan]
                species_intensity_list[ens].append(copy(intensity))
                # intensitylist[ens].append((c * (masspoints.clumpIntensity[ens]*(masspoints.clumpOpticalDepth[ens])
                #                                 /(1-np.exp(-masspoints.clumpOpticalDepth[ens]))
                #                                 /constants.voxel_size/f_ds).T).T)
            else:
                species_intensity_list[ens].append((c*masspoints.clump_species_intensity[ens].T).T)
            if suggested_calc:
                species_opticaldepth_list[ens].append((c * (masspoints.clump_species_optical_depth[ens]).T).T)
                hi_tau_list[ens].append((c*(masspoints.clump_hi_tau[ens]).T).T)

            elif test_opacity:
                species_opticaldepth_list[ens].append((c * (masspoints.clump_species_optical_depth[ens]
                                                            /constants.voxel_size/f_ds).T).T)
            else:
                species_opticaldepth_list[ens].append((c*masspoints.clump_species_optical_depth[ens].T).T)
        combinations.clump_species_intensity[ens] = np.array(species_intensity_list[ens]).sum(1)
        combinations.clump_species_optical_depth[ens] = np.array(species_opticaldepth_list[ens]).sum(1)
        combinations.clump_hi_tb[ens] = np.array(hi_tb_list[ens]).sum(1)
        combinations.clump_hi_tau[ens] = np.array(hi_tau_list[ens]).sum(1)

        if constants.dust != '' and constants.dust != None and constants.dust != 'none':
            if old_dust:
                CLcombinations = copy(combinations.clump_combination[ens])
            else:
                CLcombinations = copy(combinations.clump_max_combination[ens])
            for c in CLcombinations:
                if suggested_calc:
                    intensity = (c * (masspoints.clump_dust_intensity[ens]
                                      *(masspoints.clump_dust_optical_depth[ens])
                                      /(1-np.exp(-masspoints.clump_dust_optical_depth[ens]))
                                      ).T).T# / (4/3*masspoints.clump_radius[ens].T)
                                #/constants.voxel_size/f_ds*(4/3*masspoints.clumpRadius[ens].T)
                    i_nan = np.isnan(intensity) | np.isinf(intensity)
                    intensity[i_nan] = ((c * (masspoints.clump_dust_intensity[ens]
                                              ).T).T)[i_nan]
                                       #/ (4/3*masspoints.clump_radius[ens].T)
                                       #/constants.voxel_size/f_ds*(4/3*masspoints.clumpRadius[ens].T)
                    dust_intensity_list[ens].append(copy(intensity))
                elif test_calc:
                    intensity = (c * (masspoints.clump_dust_intensity[ens]*masspoints.clump_dust_optical_depth[ens]
                                      /(1-np.exp(-masspoints.clump_dust_optical_depth[ens]))
                                      /constants.voxel_size/f_ds).T).T
                    i_nan = np.isnan(intensity) | np.isinf(intensity)
                    intensity[i_nan] = ((c * masspoints.clump_dust_intensity[ens].T).T)[i_nan]
                    dust_intensity_list[ens].append(copy(intensity))
                    # intensitylist[ens].append((c * (masspoints.clumpIntensity[ens]*(masspoints.clumpOpticalDepth[ens])
                    #                                 /(1-np.exp(-masspoints.clumpOpticalDepth[ens]))
                    #                                 /constants.voxel_size/f_ds).T).T)
                else:
                    dust_intensity_list[ens].append((c*masspoints.clump_dust_intensity[ens].T).T)
                if suggested_calc:
                    dust_opticaldepth_list[ens].append((c * 
                        (masspoints.clump_dust_optical_depth[ens]).T).T)
                        #/(4/3*masspoints.clumpRadius[ens].T)
                elif test_opacity:
                    dust_opticaldepth_list[ens].append((c * (masspoints.clump_dust_optical_depth[ens]
                                                       /constants.voxel_size/f_ds).T).T)
                else:
                    dust_opticaldepth_list[ens].append((c*masspoints.clump_dust_optical_depth[ens].T).T)
            combinations.clump_dust_intensity[ens] = np.array(dust_intensity_list[ens]).sum(1)
            combinations.clump_dust_optical_depth[ens] = np.array(dust_opticaldepth_list[ens]).sum(1)
        
  
    return


def reinitialise():
    '''Reinitialise all temporary variables to the correct number of clump sets.'''
  
    combinations.clump_combination = [[] for _ in range(len(constants.clump_mass_number))]
    combinations.clump_max_combination = [[] for _ in range(len(constants.clump_mass_number))]
  
    combinations.clump_species_intensity = [[] for _ in range(len(constants.clump_mass_number))]
    combinations.clump_species_optical_depth = [[] for _ in range(len(constants.clump_mass_number))]
    combinations.clump_dust_intensity = [[] for _ in range(len(constants.clump_mass_number))]
    combinations.clump_dust_optical_depth = [[] for _ in range(len(constants.clump_mass_number))]

    combinations.clump_hi_tb = [[] for _ in range(len(constants.clump_mass_number))]
    combinations.clump_hi_tau = [[] for _ in range(len(constants.clump_mass_number))]
  
    return


# jit_module(nopython=False)
