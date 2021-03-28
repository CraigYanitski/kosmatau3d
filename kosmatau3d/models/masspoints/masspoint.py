import importlib as il
from time import time

import numpy as np
from numba import jit_module
import matplotlib.pyplot as plt

from kosmatau3d.models import masspoints
from kosmatau3d.models import constants
from kosmatau3d.models import interpolations
from kosmatau3d.models import species


'''
This is a module to handle one fractal mass in a combination.
It will have the associated emission and extinction information from the KOSMA-tau simulations,
which will be used to contribute to the Combination class' intrinsic intensity and optical depth.
At the moment, the emission interpolation is performed for each individual species. It will be
an future update to perform the interpolation for all species at the same time.
'''


def setMasspointData(density=[], FUV=[0]):
    '''
    This sets the information for the masspoints used in a given voxel. The density should be in units of
    cm^-3, and the FUV field should be in units of the Draine field (2.7 * 10^-3 erg cm^-2)
    '''
    for ens in range(constants.ensembles):
        masspoints.clumpLogDensity[ens] = np.log10(10.**(constants.clumpLogMass[ens]*(1-3./constants.gamma)) *
                                                   (10.**(constants.clumpLogMass[ens] *
                                                    (1+3./constants.gamma-constants.alpha))).sum() /
                                                   (10.**(constants.clumpLogMass[ens]*(2-constants.alpha))).sum() *
                                                   density[ens]/1.91)
        masspoints.clumpRadius[ens] = ((3./(4.*np.pi)*(10.**constants.clumpLogMass[ens]*constants.massSolar) /
                                       (10.**masspoints.clumpLogDensity[ens]*constants.massH*1.91))**(1./3.) /
                                       constants.pc/100.)
        masspoints.logFUV[ens] = np.log10(FUV[ens])
        
        i = (masspoints.clumpLogDensity[ens] < 3) & (np.around(masspoints.clumpLogDensity[ens]) == 3)
        if i.any():
            masspoints.clumpLogDensity[ens][i] = 3
        i = (masspoints.clumpLogDensity[ens] > 7) & (np.around(masspoints.clumpLogDensity[ens]) == 7)
        if i.any():
            masspoints.clumpLogDensity[ens][i] = 7
    
    return


def getAfuv(debug=False):
    clumpAfuv = [[] for _ in range(len(constants.clumpMassNumber))]
    for ens in range(len(constants.clumpMassNumber)):
        clumpAfuv[ens] = interpolations.interpolateFUVextinction(masspoints.clumpLogDensity[ens],
                                                                 constants.clumpLogMass[ens])[0]
    return clumpAfuv


def masspointEmission(interpolationPoint, ens, masspoint, velocity=0, verbose=False, debug=False, test=False):
    '''
    This function calculates the emission of a single clump of the given mass/density/radius.
    '''
    
    radius = masspoints.clumpRadius[ens][0, masspoint]
  
    if debug:
        print('\n', interpolationPoint)
  
    # indeces = species.molecules.getInterpolationIndeces()
  
    intensity_xi = []
    opticalDepth_xi = []
  
    if len(species.molecules):
        # Intensity currently in converted to Jansky, to coinside with the dust continuum
        intensity = interpolations.interpolateIntensity(interpolationPoint)
        tau = interpolations.interpolateTau(interpolationPoint)
    
        intensity_xi.append(intensity)
        opticalDepth_xi.append(tau)
  
    else:
        intensity_xi.append([])
        opticalDepth_xi.append([])
  
    if debug:
        print(intensity)
    
    if constants.dust:
        # Must convert Janskys in dust intensity file using I/2/constants.kB*(constants.wavelengths)**2*10**-26)
        intensity = (interpolations.interpolateDustIntensity(interpolationPoint))
        tau = interpolations.interpolateDustTau(interpolationPoint)
    
        # Append clump-corrected dust emission
        intensity_xi.append(intensity/np.pi/(np.arcsin(radius/10.)**2))
        opticalDepth_xi.append(tau)
  
    else:
        intensity_xi.append([])
        opticalDepth_xi.append([])
  
    masspoints.clumpIntensity[ens][masspoint, :] = np.append(intensity_xi[1], intensity_xi[0])
    masspoints.clumpOpticalDepth[ens][masspoint, :] = np.append(opticalDepth_xi[1], opticalDepth_xi[0])
  
    return  # (intensity,opticalDepth)


def calculateEmission(tauFUV=0, timed=False):
    '''
    This function can be called to set-up the clump emissions in the masspoints module. It calls masspointEmission() for
    each clump.
    '''
    # clumpIntensity = [[] for _ in range(len(constants.clumpMassNumber))]
    # clumpOpticalDepth = [[] for _ in range(len(constants.clumpMassNumber))]
  
    for ens in range(constants.ensembles):
        for i in range(constants.clumpLogMass[ens].size):
            t0 = time()
            gridpoint = [masspoints.clumpLogDensity[ens][0, i],
                         constants.clumpLogMass[ens][0, i],
                         masspoints.logFUV[ens]]
            masspointEmission(gridpoint, ens, i)
            # clumpIntensity[ens].append(emission[0])
            # clumpOpticalDepth[ens].append(emission[1])
            if timed:
                print(time()-t0)
    
        # masspoints.clumpIntensity[ens] = np.array(clumpIntensity[ens])
        # masspoints.clumpOpticalDepth[ens] = np.array(clumpOpticalDepth[ens])
  
    return


def plotIntensity(molecule='all', quantity='intensity', n_cl=[], title='', velocity=None, logscale=None):

    if isinstance(molecule, str) and molecule in species.molecules:
        molecule = [molecule]
  
    elif isinstance(molecule, list) or isinstance(molecule, np.ndarray):
        pass
  
    else:
        molecule = species.molecules
  
    nDust = constants.wavelengths[constants.nDust].size
    if not velocity:
        velocity = np.linspace(-5, 5, num=1000)
    profile = np.exp(-velocity**2/2/constants.clumpDispersion**2)
  
    # ylabel = r'$I_\nu \ (K)$'
  
    for ens in range(constants.ensembles):
  
        if not len(n_cl):
            n_cl = np.ones(masspoints.clumpIntensity[ens].shape[0])
    
        fig, axes = plt.subplots(masspoints.clumpIntensity[ens].shape[0],
                                 figsize=(10, 5*masspoints.clumpIntensity[ens].shape[0]))
    
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
    
        for clump in range(masspoints.clumpIntensity[ens].shape[0]):
    
            labels = []
      
            for mol in molecule:
      
                if mol not in species.molecules:
                    print('Species {} not in model.'.format(mol))
                    continue
        
                i = nDust + np.where(mol == np.asarray(species.molecules))[0][0]
                Icl = masspoints.clumpIntensity[ens][clump, i]*profile
                tcl = masspoints.clumpOpticalDepth[ens][clump, i]*profile
                intensity = Icl/tcl*(1-np.exp(-tcl))
        
                if quantity == 'emissivity':
                    value = n_cl[clump]*Icl/masspoints.clumpRadius[ens][0, clump]/2
                elif quantity == 'absorption':
                    value = n_cl[clump]*tcl/masspoints.clumpRadius[ens][0, clump]/2
                elif quantity == 'intensity':
                    value = Icl/tcl * (1-np.exp(-n_cl[clump]*tcl))
        
                if logscale:
                    ax[clump].semilogy(velocity, value, ls='-', lw=2)
                else:
                    ax[clump].plot(velocity, value, ls='-', lw=2)
        
                labels.append(mol)
      
            ax[clump].legend(labels=labels, fontsize=14, bbox_to_anchor=(1.05, 1))
            ax[clump].set_title(r'{mass} $M_\odot$ clump {q}'.format(mass=10**constants.clumpLogMass[ens][0, clump],
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
    # Reinitialise all temporary variable to the correct number of clump sets.
  
    # Properties
  
    masspoints.logFUV = [0 for _ in range(constants.ensembles)]
  
    masspoints.clumpLogDensity = [np.zeros(constants.clumpMassNumber[_]) for _ in range(constants.ensembles)]
    masspoints.clumpRadius = [np.zeros(constants.clumpMassNumber[_]) for _ in range(constants.ensembles)]
  
    # Emission
  
    masspoints.clumpIntensity = [np.zeros((constants.clumpMassNumber[_],
                                           len(species.molecules) + constants.wavelengths[constants.nDust].size))
                                 for _ in range(constants.ensembles)]
    masspoints.clumpOpticalDepth = [np.zeros((constants.clumpMassNumber[_],
                                              len(species.molecules) + constants.wavelengths[constants.nDust].size))
                                    for _ in range(constants.ensembles)]
  
    return


# jit_module(nopython=False)
