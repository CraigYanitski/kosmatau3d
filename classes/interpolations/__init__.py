from numba import jit_module
import numpy as np

from numba.errors import NumbaDeprecationWarning, NumbaPendingDeprecationWarning
import warnings
warnings.simplefilter('ignore', category=NumbaDeprecationWarning)
warnings.simplefilter('ignore', category=NumbaPendingDeprecationWarning)

import constants
from .interpolate import *
'''
This is a module that can be used for the interpolation of the input data.
It will contain functions to interpolate the intensity or optical depth
for any species, given an index.

The method of interpolation is passed as an argument when initialising
this class. The acceptabled values are 'linear', 'cubic', and 'radial'.
The default method is 'linear'. For the large intensity and optical depth
grids, 'cubic' and 'radial' are the same.
'''
# PRIVATE
#def __init__(self, species, observations, directory='MilkyWay', interpolate='linear', verbose=False):
#species = species
#indeces = np.append(species[0].getFileIndeces(), species[1].getFileIndeces())
# for element in species:
#   for i in element.getFileIndeces(): indeces.append(i)
#observations = observations
#interpolation = interpolate
#verbose = verbose
intensityInterpolation = None
tauInterpolation = None
dustIntensityInterpolation = None
dustTauInterpolation = None
rotationInterpolation = None
dispersionInterpolation = None
densityInterpolation = None
clumpMassInterpolation = None
interclumpMassInterpolation = None
FUVextinctionInterpolation = None
FUVfieldInterpolation = None
eTildeReal = None
eTildeImaginary = None

# PUBLIC

def reset():
  intensityInterpolation = None
  tauInterpolation = None
  dustIntensityInterpolation = None
  dustTauInterpolation = None
  rotationInterpolation = None
  dispersionInterpolation = None
  densityInterpolation = None
  clumpMassInterpolation = None
  interclumpMassInterpolation = None
  FUVextinctionInterpolation = None
  FUVfieldInterpolation = None
  eTildeReal = None
  eTildeImaginary = None
  return

def interpolateIntensity(points, speciesNumber, verbose=False):
  '''
  This is converted from brightness temperature to Jansky units.
  '''
  verbose = verbose or verbose
  #points = np.log10(points)
  if isinstance(speciesNumber,np.ndarray):
    intensity = np.zeros(len(speciesNumber))
    intensity_xi = 0
    for i,index in enumerate(speciesNumber):
      if constants.interpolation=='linear': intensity[i] = (10**intensityInterpolation[index](points))
      elif constants.interpolation=='radial' or interpolation=='cubic': intensity[i] = (10**intensityInterpolation[index](points[0], points[1], points[2]))
      # if (np.isnan(intensity[i]) or intensity[i]==0):
      #   intensity[i] = 10**-100
      #intensity[i] *= 2*constants.kB/4/np.pi/species.moleculeWavelengths[i]**2/10**-26
    if verbose:
      print('Calculated the intensity for {} species.'.format(len(speciesNumber)))
    return intensity
  else: return
  # else:
  #   if verbose:
  #     print('There are no species of this type adding to the intensity.')
  #   intensity = 0

def interpolateTau(points, speciesNumber, verbose=False):
  verbose = verbose or verbose
  #points = np.log10(points)
  if isinstance(speciesNumber, np.ndarray):
    tau = np.zeros(len(speciesNumber))
    for i,index in enumerate(speciesNumber):
      if constants.interpolation=='linear': tau[i] = (10**tauInterpolation[index](points))
      elif constants.interpolation=='radial' or interpolation=='cubic': tau[i] = (10**tauInterpolation[index](points[0], points[1], points[2]))
      # if np.isnan(tau[i]): tau[i] = 10**-100
      elif (tau[i]<=0):
        #temp = tau[-1]
        # tau[i] = 10**-100
        input('\n<<ERROR>> Negative opacity {} found.\n'.format(temp))
    if verbose:
      print('Calculated the optical depth for {} species.'.format(len(speciesNumber)))
    return np.array(tau)
  else: return
  # else:
  #   if verbose:
  #     print('There are no species adding to the optical depth.')
  #   tau = 0

def interpolateDustIntensity(points, verbose=False):
  '''
  This will calculate the intensity in Jansky units.
  '''
  #verbose = verbose or verbose
  # intensity = np.zeros(333, dtype=np.float64)
  if constants.interpolation=='linear':
    #print(intensity)
    intensity = 10**dustIntensityInterpolation(points)[0]
    #print(intensity)
    # constants.hclambda/constants.kB / \
    #             np.log(1+2*constants.hclambda/constants.wavelengths**2/(10**(dustIntensityInterpolation(points)[0])*10**-26))
    #/2/constants.kB*(constants.wavelengths)**2*10**-26
    intensity *= 10**-26*constants.wavelengths[constants.nDust]**2/2/constants.kB
  elif constants.interpolation=='radial' or interpolation=='cubic': intensity = (10**dustIntensityInterpolation(points[0], points[1], points[2]))
  #if np.isnan(intensity[-1]) or intensity[-1]==0: intensity[-1] = 10**-100
  if verbose:
    print('Calculated the dust intensity.')
  return intensity

def interpolateDustTau(points, verbose=False):
  verbose = verbose or verbose
  #tau = []
  if constants.interpolation=='linear': tau = (10**dustTauInterpolation(points))[0]
  elif constants.interpolation=='radial' or interpolation=='cubic': tau = (10**dustTauInterpolation(points[0], points[1], points[2]))
  # if np.isnan(tau[-1]): tau[-1] = 10**-100
  # elif tau[-1]<=0:
  #   temp = tau[-1]
  #   tau[-1] = 10**-100
  #   input('\n<<ERROR>> Negative opacity {} found.\n'.format(temp))
  if verbose:
    print('Calculated the optical depth for {} species.'.format(len(speciesNumber)))
  return tau

def interpolateRotationalVelocity(radius):
  return rotationInterpolation(radius)

def interpolateVelocityDispersion(radius):
  return dispersionInterpolation(radius)

def interpolateDensity(radius):
  density = densityInterpolation(radius)
  if (density<0).any():
    input('<<ERROR>> density {} at radius {} pc!'.format(density, radius))
    sys.exit()
  return density

def interpolateClumpMass(radius):
  mass = clumpMassInterpolation(radius)
  if (mass<0).any():
    input('<<ERROR>> clump mass {} at radius {} pc!'.format(mass, radius))
    sys.exit()
  return mass

def interpolateInterclumpMass(radius):
  mass = interclumpMassInterpolation(radius)
  if (mass<0).any():
    input('<<ERROR>> interclump mass {} at radius {} pc!'.format(mass, radius))
    sys.exit()
  return mass

def interpolateFUVextinction(density, mass):
  return 10**FUVextinctionInterpolation(density, mass)

def interpolateFUVfield(radius, height):
  return FUVfieldInterpolation(radius, abs(height))

def __str__():
  return 'Available Interpolations:\n -Clump intensity\n -Clump optical depth\n -Clump mass (galactic)\n -Clump density (galactic)\n -Voxel rotation (galactic)\n -UV extinction\n -FUV field (galactic)'

jit_module(nopython=False)
