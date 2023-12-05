"""
A subpackage to used to interpolate the KOSMA-:math:`\\tau` grids and input 
data and store these methods to use in the evaluation of :code:`kosmatau3d`.
This is a subpackage that can be used for the interpolation of the input data.
It will contain functions to interpolate the intensity or optical depth
for any species, given an index.

The method of interpolation is passed as an argument when initialising
this class. The acceptabled values are 'linear', 'cubic', and 'radial'.
The default method is 'linear'. For the large intensity and optical depth
grids, 'cubic' and 'radial' are the same.
"""


import sys
import warnings

import numpy as np
from numba import jit_module
from numba.core.errors import NumbaDeprecationWarning, NumbaPendingDeprecationWarning

from kosmatau3d.models import constants
from .interpolate import *


warnings.simplefilter('ignore', category=NumbaDeprecationWarning)
warnings.simplefilter('ignore', category=NumbaPendingDeprecationWarning)


initialised = False
clump_species_intensity_interpolation = None
interclump_species_intensity_interpolation = None
clump_species_tau_interpolation = None
interclump_species_tau_interpolation = None
clump_dust_intensity_interpolation = None
interclump_dust_intensity_interpolation = None
clump_dust_tau_interpolation = None
interclump_dust_tau_interpolation = None
clump_column_density_interpolation = None
interclump_column_density_interpolation = None
galaxy_rotation_interpolation = None
velocity_dispersion_interpolation = None
number_density_interpolation = None
h2_scale_height = None
h2_mass_full = None
h2_mass_interpolation = None
clump_h2_column_density_interpolation = None
interclump_h2_column_density_interpolation = None
hi_scale_height = None
hi_mass_full = None
hi_mass_interpolation = None
clump_hi_column_density_interpolation = None
interclump_hi_column_density_interpolation = None
clump_taufuv_interpolation = None
interclump_taufuv_interpolation = None
clump_tg_interpolation = None
interclump_tg_interpolation = None
clump_td_interpolation = None
interclump_td_interpolation = None
fuv_interpolation = None
e_tilde_real = None
e_tilde_imaginary = None


def reset():
    species_intensity_interpolation = None
    species_tau_interpolation = None
    dust_intensity_interpolation = None
    dust_tau_interpolation = None
    galaxy_rotation_interpolation = None
    velocity_dispersion_interpolation = None
    number_density_interpolation = None
    h2_mass_interpolation = None
    hi_mass_interpolation = None
    h2_column_density_interpolation = None
    hi_column_density_interpolation = None
    taufuv_interpolation = None
    tg_interpolation = None
    td_interpolation = None
    fuv_interpolation = None
    e_tilde_real = None
    e_tilde_imaginary = None
    initialised = False
    return


def interpolate_clump_species_intensity(points, verbose=False):
    '''
    This is converted from brightness temperature to Jansky units.
    '''
    # Fully 'encode' the interpolation points to the fortran standard
    if constants.log_encoded:
        points = np.asarray(points)*10
    verbose = verbose or verbose
    if len(species.clump_transitions):
        intensity = np.zeros(len(species.clump_transitions))
        intensity_xi = 0
        for i, index in enumerate(species.clump_transition_indeces):
            if constants.log_encoded:
                if constants.interpolation == 'linear':
                    intensity[i] = (10**(clump_species_intensity_interpolation[i](points)/10))
                elif constants.interpolation == 'radial' or constants.interpolation == 'cubic':
                    intensity[i] = (10**(clump_species_intensity_interpolation[i](points[0], points[1], points[2], points[3])/10))
            
            else:
                if constants.interpolation == 'linear':
                    intensity[i] = (10**clump_species_intensity_interpolation[i](points))
                elif constants.interpolation == 'radial' or constants.interpolation == 'cubic':
                    intensity[i] = (10**clump_species_intensity_interpolation[i](points[0], points[1], points[2], points[3]))
            # if (np.isnan(intensity[i]) or intensity[i]==0):
            #   intensity[i] = 10**-100
            # intensity[i] *= 2*constants.kB/4/np.pi/species.moleculeWavelengths[i]**2/10**-26
        if verbose:
            print('Calculated the intensity for {} species.'.format(len(species.clump_transitions)))
        return intensity
    else: return
    # else:
    #   if verbose:
    #     print('There are no species of this type adding to the intensity.')
    #   intensity = 0


def interpolate_interclump_species_intensity(points, verbose=False):
    '''
    This is converted from brightness temperature to Jansky units.
    '''
    # Fully 'encode' the interpolation points to the fortran standard
    if constants.log_encoded:
        points = np.asarray(points)*10
    verbose = verbose or verbose
    if len(species.interclump_transitions):
        intensity = np.zeros(len(species.interclump_transitions))
        intensity_xi = 0
        for i, index in enumerate(species.interclump_transition_indeces):
            if constants.log_encoded:
                if constants.interpolation == 'linear':
                    intensity[i] = (10**(interclump_species_intensity_interpolation[i](points)/10))
                elif constants.interpolation == 'radial' or constants.interpolation == 'cubic':
                    intensity[i] = (10**(interclump_species_intensity_interpolation[i](points[0], points[1], points[2], points[3])/10))
            
            else:
                if constants.interpolation == 'linear':
                    intensity[i] = (10**interclump_species_intensity_interpolation[i](points))
                elif constants.interpolation == 'radial' or constants.interpolation == 'cubic':
                    intensity[i] = (10**interclump_species_intensity_interpolation[i](points[0], points[1], points[2], points[3]))
            # if (np.isnan(intensity[i]) or intensity[i]==0):
            #   intensity[i] = 10**-100
            # intensity[i] *= 2*constants.kB/4/np.pi/species.moleculeWavelengths[i]**2/10**-26
        if verbose:
            print('Calculated the intensity for {} species.'.format(len(species.interclump_transitions)))
        return intensity
    else: return
    # else:
    #   if verbose:
    #     print('There are no species of this type adding to the intensity.')
    #   intensity = 0


def interpolate_clump_species_tau(points, verbose=False):
    # Fully 'encode' the interpolation points to the fortran standard
    if constants.log_encoded:
        points = np.asarray(points)*10
    verbose = verbose or verbose
    # points = np.log10(points)
    if len(species.clump_transitions):
        tau = np.zeros(len(species.clump_transitions))
        for i, index in enumerate(species.clump_transition_indeces):
            if constants.log_encoded:
                if constants.interpolation == 'linear':
                    tau[i] = (10**(clump_species_tau_interpolation[i](points)/10))
                elif constants.interpolation == 'radial' or constants.interpolation == 'cubic':
                    tau[i] = (10**(clump_species_tau_interpolation[i](points[0], points[1], points[2], points[3])/10))
            
            else:
                if constants.interpolation == 'linear':
                    tau[i] = (10**clump_species_tau_interpolation[i](points))
                elif constants.interpolation == 'radial' or constants.interpolation == 'cubic':
                    tau[i] = (10**clump_species_tau_interpolation[i](points[0], points[1], points[2], points[3]))
            # if np.isnan(tau[i]): tau[i] = 10**-100
            if (tau[i] <= 0):
                # temp = tau[-1]
                # tau[i] = 10**-100
                input('\n<<ERROR>> Negative opacity {} found.\n'.format(temp))
        if verbose:
            print('Calculated the optical depth for {} species.'.format(len(species.clump_transitions)))
        return np.array(tau)
    else:
      return
    # else:
    #   if verbose:
    #     print('There are no species adding to the optical depth.')
    #   tau = 0


def interpolate_interclump_species_tau(points, verbose=False):
    # Fully 'encode' the interpolation points to the fortran standard
    if constants.log_encoded:
        points = np.asarray(points)*10
    verbose = verbose or verbose
    # points = np.log10(points)
    if len(species.interclump_transitions):
        tau = np.zeros(len(species.interclump_transitions))
        for i, index in enumerate(species.interclump_transition_indeces):
            if constants.log_encoded:
                if constants.interpolation == 'linear':
                    tau[i] = (10**(interclump_species_tau_interpolation[i](points)/10))
                elif constants.interpolation == 'radial' or constants.interpolation == 'cubic':
                    tau[i] = (10**(interclump_species_tau_interpolation[i](points[0], points[1], points[2], points[3])/10))
            
            else:
                if constants.interpolation == 'linear':
                    tau[i] = (10**interclump_species_tau_interpolation[i](points))
                elif constants.interpolation == 'radial' or constants.interpolation == 'cubic':
                    tau[i] = (10**interclump_species_tau_interpolation[i](points[0], 
                        points[1], points[2], points[3]))
            # if np.isnan(tau[i]): tau[i] = 10**-100
            if (tau[i] <= 0):
                # temp = tau[-1]
                # tau[i] = 10**-100
                input('\n<<ERROR>> Negative opacity {} found.\n'.format(temp))
        if verbose:
            print('Calculated the optical depth for {} species.'.format(
                len(species.interclump_transitions)))
        return np.array(tau)
    else:
      return
    # else:
    #   if verbose:
    #     print('There are no species adding to the optical depth.')
    #   tau = 0


def interpolate_clump_dust_intensity(points, verbose=False):
    '''
    This will calculate the intensity in Jansky units.
    '''
    # Fully 'encode' the interpolation points to the fortran standard
    if constants.log_encoded:
        points = np.asarray(points)*10
  
        if constants.interpolation == 'linear':
            intensity = []
            for dust in clump_dust_intensity_interpolation:
                intensity.append(10**(dust(points)[0]/10))
        
        elif constants.interpolation == 'radial' or constants.interpolation == 'cubic':
            intensity = (10**(clump_dust_intensity_interpolation(points[0], points[1], points[2], points[3])/10))
    
    else:
        if constants.interpolation == 'linear':
            intensity = []
            for dust in clump_dust_intensity_interpolation:
                intensity.append(10**(dust(points)[0]))
    
        elif constants.interpolation == 'radial' or constants.interpolation == 'cubic':
            intensity = (10**clump_dust_intensity_interpolation(points[0], points[1], points[2], points[3]))
    
    # Convert specific flux Jansky units to brightness temperature Kelvin units
    intensity = np.asarray(intensity) * 10**-26 * constants.wavelengths[constants.n_dust]**2/2/constants.kB
  
    return intensity


def interpolate_interclump_dust_intensity(points, verbose=False):
    '''
    This will calculate the intensity in Jansky units.
    '''
    # Fully 'encode' the interpolation points to the fortran standard
    if constants.log_encoded:
        points = np.asarray(points)*10
  
        if constants.interpolation == 'linear':
            intensity = []
            for dust in interclump_dust_intensity_interpolation:
                intensity.append(10**(dust(points)[0]/10))
        
        elif constants.interpolation == 'radial' or constants.interpolation == 'cubic':
            intensity = (10**(interclump_dust_intensity_interpolation(points[0], points[1], points[2], points[3])/10))
    
    else:
        if constants.interpolation == 'linear':
            intensity = []
            for dust in interclump_dust_intensity_interpolation:
                intensity.append(10**(dust(points)[0]))
    
        elif constants.interpolation == 'radial' or constants.interpolation == 'cubic':
            intensity = (10**interclump_dust_intensity_interpolation(points[0], points[1], points[2], points[3]))
    
    # Convert specific flux Jansky units to brightness temperature Kelvin units
    intensity = np.asarray(intensity) * 10**-26 * constants.wavelengths[constants.n_dust]**2/2/constants.kB
  
    return intensity


def interpolate_clump_dust_tau(points, verbose=False):
    # Fully 'encode' the interpolation points to the fortran standard
    if constants.log_encoded:
        points = np.asarray(points)*10
        
        if constants.interpolation == 'linear':
            tau = []
            for dust in clump_dust_tau_interpolation:
                tau.append(10**(dust(points)[0]/10))
  
        elif constants.interpolation == 'radial' or constants.interpolation == 'cubic':
            tau = (10**(clump_dust_tau_interpolation(points[0], points[1], points[2], points[3])/10))
    
    else:
        if constants.interpolation == 'linear':
            tau = []
            for dust in clump_dust_tau_interpolation:
                tau.append(10**(dust(points)[0]))
    
        elif constants.interpolation == 'radial' or constants.interpolation == 'cubic':
            tau = (10**clump_dust_tau_interpolation(points[0], points[1], points[2], points[3]))
    # if np.isnan(tau[-1]): tau[-1] = 10**-100
    # elif tau[-1]<=0:
    #   temp = tau[-1]
    #   tau[-1] = 10**-100
    #   input('\n<<ERROR>> Negative opacity {} found.\n'.format(temp))
    tau = np.asarray(tau)
  
    if verbose:
      print('Calculated the optical depth for {} species.'.format(len(speciesNumber)))
    return tau


def interpolate_interclump_dust_tau(points, verbose=False):
    # Fully 'encode' the interpolation points to the fortran standard
    if constants.log_encoded:
        points = np.asarray(points)*10
        
        if constants.interpolation == 'linear':
            tau = []
            for dust in interclump_dust_tau_interpolation:
                tau.append(10**(dust(points)[0]/10))
  
        elif constants.interpolation == 'radial' or constants.interpolation == 'cubic':
            tau = (10**(interclump_dust_tau_interpolation(points[0], points[1], points[2], points[3])/10))
    
    else:
        if constants.interpolation == 'linear':
            tau = []
            for dust in interclump_dust_tau_interpolation:
                tau.append(10**(dust(points)[0]))
    
        elif constants.interpolation == 'radial' or constants.interpolation == 'cubic':
            tau = (10**interclump_dust_tau_interpolation(points[0], points[1], points[2], points[3]))
    # if np.isnan(tau[-1]): tau[-1] = 10**-100
    # elif tau[-1]<=0:
    #   temp = tau[-1]
    #   tau[-1] = 10**-100
    #   input('\n<<ERROR>> Negative opacity {} found.\n'.format(temp))
    tau = np.asarray(tau)
  
    if verbose:
      print('Calculated the optical depth for {} species.'.format(len(speciesNumber)))
    return tau


def interpolate_clump_column_density(points, sp):
    idx = constants.abundances.index(sp)
    if constants.interpolation == 'linear':
        return 10**clump_column_density_interpolation[idx](points)[0]
    elif constants.interpolation == 'radial' or constants.interpolation == 'cubic':
        return 10**clump_column_density_interpolation[idx](points[0], points[1], points[2])
    else:
        sys.exit('<<ERROR>>: There is no such method as ' + 
                 '{} to interpolate the '.format(constants.interpolation) +
                 'column density from the KOSMA-tau grid.\n\nExitting...\n\n')


def interpolate_interclump_column_density(points, sp):
    idx = constants.abundances.index(sp)
    if constants.interpolation == 'linear':
        return 10**interclump_column_density_interpolation[idx](points)[0]
    elif constants.interpolation == 'radial' or constants.interpolation == 'cubic':
        return 10**interclump_column_density_interpolation[idx](points[0], points[1], points[2])
    else:
        sys.exit('<<ERROR>>: There is no such method as ' + 
                 '{} to interpolate the '.format(constants.interpolation) +
                 'column density from the KOSMA-tau grid.\n\nExitting...\n\n')


def interpolate_clump_hi_column_density(points):
    if constants.interpolation == 'linear':
        return 10**clump_hi_column_density_interpolation(points)
    elif constants.interpolation == 'radial' or constants.interpolation == 'cubic':
        return 10**clump_hi_column_density_interpolation(points[0], points[1], points[2])
    else:
        sys.exit('<<ERROR>>: There is no such method as ' + 
                 '{} to interpolate the '.format(constants.interpolation) +
                 'column density from the KOSMA-tau grid.\n\nExitting...\n\n')


def interpolate_interclump_hi_column_density(points):
    if constants.interpolation == 'linear':
        return 10**interclump_hi_column_density_interpolation(points)
    elif constants.interpolation == 'radial' or constants.interpolation == 'cubic':
        return 10**interclump_hi_column_density_interpolation(points[0], points[1], points[2])
    else:
        sys.exit('<<ERROR>>: There is no such method as ' + 
                 '{} to interpolate the '.format(constants.interpolation) +
                 'column density from the KOSMA-tau grid.\n\nExitting...\n\n')


def interpolate_clump_h2_column_density(points):
    if constants.interpolation == 'linear':
        return 10**clump_h2_column_density_interpolation(points)
    elif constants.interpolation == 'radial' or constants.interpolation == 'cubic':
        return 10**clump_h2_column_density_interpolation(points[0], points[1], points[2])
    else:
        sys.exit('<<ERROR>>: There is no such method as ' + 
                 '{} to interpolate the '.format(constants.interpolation) +
                 'column density from the KOSMA-tau grid.\n\nExitting...\n\n')


def interpolate_interclump_h2_column_density(points):
    if constants.interpolation == 'linear':
        return 10**interclump_h2_column_density_interpolation(points)
    elif constants.interpolation == 'radial' or constants.interpolation == 'cubic':
        return 10**interclump_h2_column_density_interpolation(points[0], points[1], points[2])
    else:
        sys.exit('<<ERROR>>: There is no such method as ' + 
                 '{} to interpolate the '.format(constants.interpolation) +
                 'column density from the KOSMA-tau grid.\n\nExitting...\n\n')


def interpolate_clump_tg(points):
    if constants.interpolation == 'linear':
        return 10**clump_tg_interpolation(points)
    elif constants.interpolation == 'radial' or constants.interpolation == 'cubic':
        return 10**clump_tg_interpolation(points[0], points[1], points[2])
    else:
        sys.exit('<<ERROR>>: There is no such method as ' + 
                 '{} to interpolate the '.format(constants.interpolation) +
                 'temperature from the KOSMA-tau grid.\n\nExitting...\n\n')


def interpolate_interclump_tg(points):
    if constants.interpolation == 'linear':
        return 10**interclump_tg_interpolation(points)
    elif constants.interpolation == 'radial' or constants.interpolation == 'cubic':
        return 10**interclump_tg_interpolation(points[0], points[1], points[2])
    else:
        sys.exit('<<ERROR>>: There is no such method as ' + 
                 '{} to interpolate the '.format(constants.interpolation) +
                 'temperature from the KOSMA-tau grid.\n\nExitting...\n\n')


def interpolate_clump_td(points):
    if constants.interpolation == 'linear':
        return 10**clump_td_interpolation(points)
    elif constants.interpolation == 'radial' or constants.interpolation == 'cubic':
        return 10**clump_td_interpolation(points[0], points[1], points[2])
    else:
        sys.exit('<<ERROR>>: There is no such method as ' + 
                 '{} to interpolate the '.format(constants.interpolation) +
                 'temperature from the KOSMA-tau grid.\n\nExitting...\n\n')


def interpolate_interclump_td(points):
    if constants.interpolation == 'linear':
        return 10**interclump_td_interpolation(points)
    elif constants.interpolation == 'radial' or constants.interpolation == 'cubic':
        return 10**interclump_td_interpolation(points[0], points[1], points[2])
    else:
        sys.exit('<<ERROR>>: There is no such method as ' + 
                 '{} to interpolate the '.format(constants.interpolation) +
                 'temperature from the KOSMA-tau grid.\n\nExitting...\n\n')


def interpolate_galaxy_rotation(radius):
    return galaxy_rotation_interpolation(radius)


def interpolate_velocity_dispersion(radius):
    return velocity_dispersion_interpolation(radius)


def interpolate_number_density(radius):
    density = number_density_interpolation(radius)
    if (density < 0).any():
        input('<<ERROR>> density {} at radius {} pc!'.format(density, radius))
        sys.exit()
    return density


def interpolate_h2_voxel_filling_factor(radius, height):
    # Calculate voxel-filling factor
    h = h2_scale_height(radius).mean()
    h_min = np.max((-h, height-constants.voxel_size/2))
    h_max = np.min((h, height+constants.voxel_size/2))
    f_vox = (h_max-h_min)/constants.voxel_size
    return f_vox


def interpolate_h2_mass(radius, height):
    # Calculate mass
    mass = h2_mass_interpolation(radius, np.abs(height))
    if (mass < 0).any():
        input('<<ERROR>> clump mass {} at radius {} pc!'.format(mass, radius))
        sys.exit()
    return mass


def interpolate_hi_voxel_filling_factor(radius, height):
    # Calculate voxel-filling factor
    h = hi_scale_height(radius).mean()
    h_min = np.max((-h, height-constants.voxel_size/2))
    h_max = np.min((h, height+constants.voxel_size/2))
    f_vox = (h_max-h_min)/constants.voxel_size
    return f_vox


def interpolate_hi_mass(radius, height):
    # Calculate mass
    mass = hi_mass_interpolation(radius, np.abs(height))
    if (mass < 0).any():
        input('<<ERROR>> interclump mass {} at radius {} pc!'.format(mass, radius))
        sys.exit()
    return mass


def interpolate_clump_taufuv(density, mass):
    return 10**clump_taufuv_interpolation(density, mass)


def interpolate_interclump_taufuv(density, mass):
    return 10**interclump_taufuv_interpolation(density, mass, 0)


def interpolate_fuv(radius, height):
    fuv_temp = fuv_interpolation(radius, abs(height))
    # if np.mean(radius) < 4500:
    #     fuv_temp *= constants.fuv_scale_gc
    return fuv_temp


def __str__():
    return 'Available Interpolations:\n ' \
           '-H2 intensity\n ' \
           '-H2 optical depth\n ' \
           '-H2 mass (galactic)\n ' \
           '-H2 density (galactic)\n ' \
           '-HI intensity\n ' \
           '-HI optical depth\n ' \
           '-HI mass (galactic)\n ' \
           '-HI density (galactic)\n ' \
           '-Voxel rotation (galactic)\n ' \
           '-FUV extinction\n ' \
           '-FUV radiation (galactic)'


# jit_module(nopython=False)
