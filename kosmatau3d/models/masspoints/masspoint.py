"""
This is a module to handle one fractal mass in a combination.
It will have the associated emission and extinction information from the 
KOSMA-:math:`\\tau` simulations, which will be used to contribute to the 
combination intrinsic intensity and optical depth.
At the moment, the emission interpolation is performed for each individual 
species.
It will be an future update to perform the interpolation for all species at the 
same time.
"""


import importlib as il
from time import time

import numpy as np
from copy import copy, deepcopy
from numba import jit_module
import matplotlib.pyplot as plt

from kosmatau3d.models import masspoints
from kosmatau3d.models import constants
from kosmatau3d.models import interpolations
from kosmatau3d.models import species


def set_masspoint_data(density=[], fuv=[0], crir=2e-16):
    '''
    This sets the information for the masspoints used in a given voxel. The density should be in units of
    cm^-3, and the FUV field should be in units of the Draine field (2.7 * 10^-3 erg cm^-2)
    '''

    masspoints.log_crir = np.log10(crir)

    for ens in range(constants.ensembles):
        masspoints.clump_log_density[ens] = np.log10(10.**(constants.clump_log_mass[ens]*(1-3./constants.gamma)) *
                                                     (10.**(constants.clump_log_mass[ens] *
                                                      (1+3./constants.gamma-constants.alpha))).sum() /
                                                     (10.**(constants.clump_log_mass[ens]*(2-constants.alpha))).sum() *
                                                     density[ens]/1.91)
        masspoints.clump_log_density_orig[ens] = deepcopy(masspoints.clump_log_density[ens])
        masspoints.clump_radius[ens] = ((3./(4.*np.pi)*(10.**constants.clump_log_mass[ens]*constants.mass_solar) /
                                        (10.**masspoints.clump_log_density[ens]*constants.mass_h*1.91))**(1./3.) /
                                        constants.pc/100.)
        masspoints.log_fuv[ens] = np.log10(fuv[ens])
        
        i = (masspoints.clump_log_density[ens] < 3) & (np.around(masspoints.clump_log_density[ens]) == 3)
        if i.any():
            masspoints.clump_log_density[ens][i] = 3
        i = (masspoints.clump_log_density[ens] > 7) & (np.around(masspoints.clump_log_density[ens]) == 7)
        if i.any():
            masspoints.clump_log_density[ens][i] = 7
    
    return


def get_taufuv(debug=False):
    '''
    Return the far-UV optical depth for each clump in each ensemble.
    '''
    taufuv = [[] for _ in range(len(constants.clump_mass_number))]
    for ens in range(len(constants.clump_mass_number)):
        if constants.interclump_idx[ens]:
            taufuv[ens] = interpolations.interpolate_interclump_taufuv(masspoints.clump_log_density[ens],
                                                                       constants.clump_log_mass[ens])[0]
        else:
            taufuv[ens] = interpolations.interpolate_clump_taufuv(masspoints.clump_log_density[ens],
                                                                  constants.clump_log_mass[ens])[0]
    return taufuv


def masspoint_emission(interpolation_point, ens, masspoint, velocity=0, verbose=False, debug=False, test=False):
    '''
    This function calculates the emission of a single clump of the given mass/density/radius.
    '''
    
    radius = masspoints.clump_radius[ens][0, masspoint]
  
    if debug:
        print('\n', interpolation_point)
  
    # indeces = species.molecules.getInterpolationIndeces()
  
    intensity_xi = []
    opticaldepth_xi = []
  
    if len(species.clump_transitions):
        # Intensity currently in converted to Jansky, to coinside with the dust continuum
        if constants.interclump_idx[ens]:
            intensity = interpolations.interpolate_interclump_species_intensity(interpolation_point[1:])
            tau = interpolations.interpolate_interclump_species_tau(interpolation_point[1:])
        else:
            intensity = interpolations.interpolate_clump_species_intensity(interpolation_point)
            tau = interpolations.interpolate_clump_species_tau(interpolation_point)
    
        intensity_xi.append(deepcopy(intensity))
        opticaldepth_xi.append(deepcopy(tau))
  
    else:
        intensity_xi.append([])
        opticaldepth_xi.append([])
  
    if debug:
        print(intensity)
    
    if constants.dust != '' and constants.dust != None and constants.dust != 'none':
        # Must convert Janskys in dust intensity file using I/2/constants.kB*(constants.wavelengths)**2*10**-26) -- now calculated in interpolation function
        if constants.interclump_idx[ens]:
            intensity = (interpolations.interpolate_interclump_dust_intensity(interpolation_point[1:]))
            tau = interpolations.interpolate_interclump_dust_tau(interpolation_point[1:])
        else:
            intensity = (interpolations.interpolate_clump_dust_intensity(interpolation_point))
            tau = interpolations.interpolate_clump_dust_tau(interpolation_point)
    
        # Append clump-corrected dust emission
        intensity_xi.append(intensity/np.pi/(np.arcsin(radius/10.)**2))
        opticaldepth_xi.append(deepcopy(tau))
  
    else:
        intensity_xi.append([])
        opticaldepth_xi.append([])
  
    masspoints.clump_species_intensity[ens][masspoint, :] = deepcopy(intensity_xi[0])
    masspoints.clump_species_optical_depth[ens][masspoint, :] = deepcopy(opticaldepth_xi[0])
    masspoints.clump_dust_intensity[ens][masspoint, :] = deepcopy(intensity_xi[1])
    masspoints.clump_dust_optical_depth[ens][masspoint, :] = deepcopy(opticaldepth_xi[1])

    if constants.interclump_idx[ens]:
        t_g = interpolations.interpolate_interclump_tg(interpolation_point[1:])
        t_d = interpolations.interpolate_interclump_td(interpolation_point[1:])
    else:
        t_g = interpolations.interpolate_clump_tg(interpolation_point[1:])
        t_d = interpolations.interpolate_clump_td(interpolation_point[1:])
    masspoints.clump_t_gas[ens][0, masspoint] = copy(t_g)
    masspoints.clump_t_dust[ens][0, masspoint] = copy(t_d)

    if constants.interclump_idx[ens]:
        n_hi = interpolations.interpolate_interclump_hi_column_density(interpolation_point[1:])
        n_h2 = interpolations.interpolate_interclump_h2_column_density(interpolation_point[1:])
        n_species = []
        for sp in constants.abundances:
            n_species.append(interpolations.interpolate_interclump_column_density(interpolation_point[1:], sp))
    else:
        n_hi = interpolations.interpolate_clump_hi_column_density(interpolation_point[1:])
        n_h2 = interpolations.interpolate_clump_h2_column_density(interpolation_point[1:])
        n_species = []
        for sp in constants.abundances:
            n_species.append(interpolations.interpolate_clump_column_density(interpolation_point[1:], sp))
    masspoints.clump_hi_col_dens[ens][0, masspoint] = copy(n_hi)
    masspoints.clump_h2_col_dens[ens][0, masspoint] = copy(n_h2)

    rcl = masspoints.clump_radius[ens][0, masspoint] \
            * 10**((masspoints.clump_log_density_orig[ens][0, masspoint]-masspoints.clump_log_density[ens][0, masspoint])/3.)
    kappa = 7.5419439e-18 / t_g * n_hi/1e20 / rcl
    tb = t_g * (1-np.exp(-4.1146667e18*kappa*rcl))
    masspoints.clump_hi_tb[ens][masspoint, 0] = copy(tb)
    masspoints.clump_hi_tau[ens][masspoint, 0] = kappa * 4/3*rcl*constants.pc*100
    masspoints.clump_N_species[ens][masspoint, :] = np.asarray([n_sp * np.pi*(rcl*constants.pc*100)**2 for n_sp in n_species])
    masspoints.clump_hi_mass[ens][0, masspoint] = n_hi * np.pi*(rcl*constants.pc*100)**2 \
            * constants.mass_h/constants.mass_solar
    masspoints.clump_h2_mass[ens][0, masspoint] = n_h2 * np.pi*(rcl*constants.pc*100)**2 \
            * 2*constants.mass_h/constants.mass_solar

    return  # (intensity,opticalDepth)


def calculate_emission(taufuv=0, timed=False):
    '''
    This function can be called to set-up the clump emissions in the masspoints module. It calls masspointEmission() for
    each clump.
    '''
    # clumpIntensity = [[] for _ in range(len(constants.clumpMassNumber))]
    # clumpOpticalDepth = [[] for _ in range(len(constants.clumpMassNumber))]
  
    for ens in range(constants.ensembles):
        for i in range(constants.clump_log_mass[ens].size):
            t0 = time()
            gridpoint = [masspoints.log_crir,
                         masspoints.clump_log_density[ens][0, i],
                         constants.clump_log_mass[ens][0, i],
                         masspoints.log_fuv[ens]]
            # print(f'interpolation: {gridpoint}')
            masspoint_emission(gridpoint, ens, i)
            if timed:
                print(time()-t0)
    
        # masspoints.clumpIntensity[ens] = np.array(clumpIntensity[ens])
        # masspoints.clumpOpticalDepth[ens] = np.array(clumpOpticalDepth[ens])
  
    return


def plot_intensity(transitions='all', quantity='intensity', n_cl=[], title='', velocity=None, logscale=None,
                   test_calc=False):
    '''
    Plot the emissivity or absorption for each clump in each ensemble.
    It is assumed that the line profile is Gaussian.
    '''

    if isinstance(transitions, str) and transitions in species.clump_transitions:
        transitions = [transitions]
  
    elif isinstance(transitions, list) or isinstance(transitions, np.ndarray):
        pass
  
    else:
        transitions = species.clump_transitions
  
    n_dust = constants.wavelengths[constants.n_dust].size
    if not velocity:
        velocity = np.linspace(-5, 5, num=1000)
    profile = np.exp(-velocity**2/2/constants.clump_dispersion**2)
  
    # ylabel = r'$I_\nu \ (K)$'
  
    for ens in range(constants.ensembles):
  
        if not len(n_cl):
            n_cl = np.ones(masspoints.clump_species_intensity[ens].shape[0])
    
        fig, axes = plt.subplots(masspoints.clump_intensity[ens].shape[0],
                                 figsize=(10, 5*masspoints.clump_intensity[ens].shape[0]))
    
        if isinstance(axes, np.ndarray):
            ax = axes
    
        else:
            ax = [axes]
    
        if quantity == 'emissivity':
            ylabel = r'$\epsilon_{\lambda} \ \left( \frac{K}{pc} \right)$'
    
        elif quantity == 'absorption':
            ylabel = r'$\kappa_{\lambda} \ \left( \frac{1}{pc} \right)$'
    
        elif quantity == 'intensity':
            ylabel = r'$I_{\lambda} \ \left( K \right)$'
    
        for clump in range(masspoints.clump_intensity[ens].shape[0]):
    
            labels = []
      
            for transition in transitions:
      
                if transition not in species.clump_transitions:
                    print('Species {} not in model.'.format(transition))
                    continue
        
                i = n_dust + np.where(transition == np.asarray(species.clump_transitions))[0][0]
                Icl = masspoints.clump_intensity[ens][clump, i]*profile
                tcl = masspoints.clump_optical_depth[ens][clump, i]*profile
                intensity = Icl/tcl*(1-np.exp(-tcl))
                if test_calc:
                    ds = constants.voxel_size
                    rcl = masspoints.clump_radius[ens][clump, 0]
                    eps = Icl * tcl/(1-np.exp(-tcl)) / ds
                    kap = tcl / ds
        
                if quantity == 'emissivity':
                    if test_calc:
                        value = eps
                    else:
                        value = n_cl[clump]*Icl/masspoints.clump_radius[ens][0, clump]/2
                elif quantity == 'absorption':
                    if test_calc:
                        value = kap
                    else:
                        value = n_cl[clump]*tcl/masspoints.clump_radius[ens][0, clump]/2
                elif quantity == 'intensity':
                    if test_calc:
                        value = eps/kap * (1-np.exp(-kap*ds))
                    else:
                        value = Icl/tcl * (1-np.exp(-n_cl[clump]*tcl))
                else:
                    raise ValueError('quantity should be "emissivity", "absorption", or "intensity"')
        
                if logscale:
                    ax[clump].semilogy(velocity, value, ls='-', lw=2)
                else:
                    ax[clump].plot(velocity, value, ls='-', lw=2)
        
                labels.append(transition)
      
            ax[clump].legend(labels=labels, fontsize=14, bbox_to_anchor=(1.05, 1))
            ax[clump].set_title(r'{mass} $M_\odot$ clump {q}'.format(mass=10**constants.clump_log_mass[ens][0, clump],
                                                                     q=quantity), fontsize=20)
            ax[clump].set_xlabel(r'$V_r \ (\frac{km}{s})$', fontsize=20)
            ax[clump].set_ylabel(ylabel, fontsize=20)
            plt.setp(ax[clump].get_xticklabels(), fontsize=14)
            plt.setp(ax[clump].get_yticklabels(), fontsize=14)
    
        # if title:
        #   fig.set_title(title, fontsize=20)
        # else:
        #   fig.set_title('Clump set {} {}'.format(ens+1, quantity), fontsize=20)
        
        fig.tight_layout()
    
        plt.show()
  
    return


def reinitialise():
    '''Reinitialise all temporary variable to the correct number of clump sets.'''
  
    # Properties
  
    masspoints.log_crir = [0 for _ in range(constants.ensembles)]
    masspoints.log_fuv = [0 for _ in range(constants.ensembles)]
  
    masspoints.clump_log_density = [np.zeros((1, constants.clump_mass_number[_])) 
                                    for _ in range(constants.ensembles)]
    masspoints.clump_log_density_orig = [np.zeros((1, constants.clump_mass_number[_])) 
                                         for _ in range(constants.ensembles)]
    masspoints.clump_radius = [np.zeros((1, constants.clump_mass_number[_])) 
                               for _ in range(constants.ensembles)]
    masspoints.clump_t_gas = [np.zeros((1, constants.clump_mass_number[_])) 
                              for _ in range(constants.ensembles)]
    masspoints.clump_t_dust = [np.zeros((1, constants.clump_mass_number[_])) 
                               for _ in range(constants.ensembles)]
    masspoints.clump_hi_col_dens = [np.zeros((1, constants.clump_mass_number[_])) 
                                    for _ in range(constants.ensembles)]
    masspoints.clump_h2_col_dens = [np.zeros((1, constants.clump_mass_number[_])) 
                                    for _ in range(constants.ensembles)]
    masspoints.clump_hi_mass = [np.zeros((1, constants.clump_mass_number[_])) 
                                    for _ in range(constants.ensembles)]
    masspoints.clump_h2_mass = [np.zeros((1, constants.clump_mass_number[_])) 
                                    for _ in range(constants.ensembles)]
    masspoints.clump_N_species = [np.zeros((constants.clump_mass_number[_], 
                                                 len(constants.abundances)))
                                       for _ in range(len(constants.clump_mass_number))]
  
    # Emission
  
    masspoints.clump_intensity = [np.zeros((constants.clump_mass_number[_],
                                            len(species.clump_transitions) 
                                                + constants.wavelengths[constants.n_dust].size))
                                  for _ in range(constants.ensembles)]
    masspoints.clump_optical_depth = [np.zeros((constants.clump_mass_number[_],
                                                len(species.clump_transitions) 
                                                    + constants.wavelengths[constants.n_dust].size))
                                      for _ in range(constants.ensembles)]
    masspoints.clump_species_intensity = [np.zeros((constants.clump_mass_number[_], 
                                                   len(species.clump_transitions)))
                                          for _ in range(len(constants.clump_mass_number))]
    masspoints.clump_species_optical_depth = [np.zeros((constants.clump_mass_number[_], 
                                                       len(species.clump_transitions)))
                                              for _ in range(len(constants.clump_mass_number))]
    masspoints.clump_dust_intensity = [np.zeros((constants.clump_mass_number[_],
                                                 constants.wavelengths[constants.n_dust].size))
                                       for _ in range(len(constants.clump_mass_number))]
    masspoints.clump_dust_optical_depth = [np.zeros((constants.clump_mass_number[_],
                                                     constants.wavelengths[constants.n_dust].size))
                                           for _ in range(len(constants.clump_mass_number))]
    masspoints.clump_hi_tb = [np.zeros((constants.clump_mass_number[_], 1)) 
                              for _ in range(constants.ensembles)]
    masspoints.clump_hi_tau = [np.zeros((constants.clump_mass_number[_], 1)) 
                               for _ in range(constants.ensembles)]
  
    return


# jit_module(nopython=False)
