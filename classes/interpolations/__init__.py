import constants
import numpy as np
from interpolations import interpolate
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
rotationInterpolation = None
dispersionInterpolation = None
densityInterpolation = None
clumpMassInterpolation = None
interclumpMassInterpolation = None
FUVextinctionInterpolation = None
FUVfieldInterpolation = None

# PUBLIC

def interpolateIntensity(points, speciesNumber, verbose=False):
  verbose = verbose or verbose
  #points = np.log10(points)
  if len(speciesNumber):
    intensity = []
    for i in speciesNumber:
      if constants.interpolation=='linear': intensity.append(10**intensityInterpolation[i](points))
      elif constants.interpolation=='radial' or interpolation=='cubic': intensity.append(10**intensityInterpolation[i](points[0], points[1], points[2]))
      if np.isnan(intensity[-1]) or intensity[-1]==0: intensity[-1] = 10**-100
    if verbose:
      print('Calculated the intensity for {} species.'.format(len(speciesNumber)))
  else:
    if verbose:
      print('There are no species of this type adding to the intensity.')
    intensity = 0
  return np.array(intensity)

def interpolateTau(points, speciesNumber, verbose=False):
  verbose = verbose or verbose
  #points = np.log10(points)
  if len(speciesNumber):
    tau = []
    for i in speciesNumber:
      if constants.interpolation=='linear': tau.append(10**tauInterpolation[i](points))
      elif constants.interpolation=='radial' or interpolation=='cubic': tau.append(10**tauInterpolation[i](points[0], points[1], points[2]))
      if np.isnan(tau[-1]): tau[-1] = 10**-100
      elif tau[-1]<=0:
        temp = tau[-1]
        tau[-1] = 10**-100
        input('\n<<ERROR>> Negative opacity {} found.\n'.format(temp))
    if verbose:
      print('Calculated the optical depth for {} species.'.format(len(speciesNumber)))
  else:
    if verbose:
      print('There are no species adding to the optical depth.')
    tau = 0
  return np.array(tau)

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

def interpolateFUVfield(radius):
  return FUVfieldInterpolation(radius)

def __str__():
  return 'Available Interpolations:\n -Clump intensity\n -Clump optical depth\n -Clump mass (galactic)\n -Clump density (galactic)\n -Voxel rotation (galactic)\n -UV extinction\n -FUV field (galactic)'
