import numpy as np
from numba import jit_module
import importlib as il

import masspoints

import constants
import interpolations
import species

# class Masspoint(object):
#   '''
#   This is a class to handle one fractal mass in a combination.
#   It will have the associated emission and extinction information from the KOSMA-tau simulations,
#   which will be used to contribute to the Combination class' intrinsic intensity and optical depth.
#   At the moment, the emission interpolation is performed for each individual species. It will be
#   an future update to perform the interpolation for all species at the same time.
#   '''
#   # PRIVATE

#   def __init__(self, density=0, mass=0, fuv=0, number=1, debugging=False):
#     # self.__species = species     #list of both moleculular and dust species
#     self.__density = density
#     self.__mass = mass 	#mass of this KOSMA-tau simulation
#     #self.__Afuv = self.__interpolations.interpolateFUVextinction(density, mass)
#     self.__FUV = fuv
#     #input('{}: {}'.format(mass, number))
#     #self.__FUV = FUVfield()     #the FUV field for this combination of mass points
#     self.__debugging = debugging
#     if debugging:
#       self.number = number    #number of this KOSMA-tau model contributing to the combination (OBSOLETE)
#       self.intensity_xi = []       #velocity-averaged intensity of this combination of masspoints
#       self.opticalDepth_xi = []    #velocity-averaged optical depth of this combination of masspoints
#     return

#   def __str__(self):
#     return 'Simulated KOSMA-tau clump of mass {}'.format(10**float(self.__mass))

  # PUBLIC
def setMasspointData(density=0, FUV=0):
  '''
  This sets the information for the masspoints used in a given voxel. The density should be in units of
  cm^-3, and the FUV field should be in units of the Draine field (2.7 * 10^-3 erg cm^-2)
  '''
  masspoints.clumpLogDensity = np.log10(10.**(constants.clumpLogMass*(1-3./constants.gamma))*sum(10.**(constants.clumpLogMass*(1+3./constants.gamma-constants.alpha))) / \
                                    sum(10.**(constants.clumpLogMass*(2-constants.alpha)))*density/1.91)
  masspoints.interclumpLogDensity = np.log10(10.**(constants.interclumpLogMass*(1-3./constants.gamma))*sum(10.**(constants.interclumpLogMass*(1+3./constants.gamma-constants.alpha))) / \
                                        sum(10.**(constants.interclumpLogMass*(2-constants.alpha)))*density/1.91)
  masspoints.logFUV = np.log10(FUV)
  masspoints.clumpRadius = ((3./(4.*np.pi)*(10.**constants.clumpLogMass*constants.massSolar)/ \
                              (10.**masspoints.clumpLogDensity*constants.massH*1.91))**(1./3.)/constants.pc/100.)
  masspoints.interclumpRadius = ((3./(4.*np.pi)*(10.**constants.interclumpLogMass*constants.massSolar)/ \
                              (10.**masspoints.interclumpLogDensity*constants.massH*1.91))**(1./3.)/constants.pc/100.)
  return

def getAfuv(debug=False):
  clumpAfuv = interpolations.interpolateFUVextinction(masspoints.clumpLogDensity, constants.clumpLogMass)
  interclumpAfuv = interpolations.interpolateFUVextinction(masspoints.interclumpLogDensity, constants.interclumpLogMass)
  if debug and self.__mass<0:
    print('\n', masspoints.clumpLogDensity, constants.clumpLogMass, clumpAfuv)
    print('\n', masspoints.interclumpLogDensity, constants.interclumpLogMass, interclumpAfuv)
  return (clumpAfuv,interclumpAfuv)

def calculateEmission(Afuv=0):
  '''
  This function can be called to set-up the clump emissions in the masspoints module. It calls masspointEmission() for
  each clump. 
  '''
  clumpIntensity = []
  clumpOpticalDepth = []

  interclumpIntensity = []
  interclumpOpticalDepth = []

  for i in range(constants.clumpLogMass.size):
    gridpoint = [masspoints.clumpLogDensity[0,i], constants.clumpLogMass[0,i], masspoints.logFUV] #include Afuv
    emission = masspointEmission(gridpoint)
    clumpIntensity.append(emission[0])
    clumpOpticalDepth.append(emission[1])

  masspoints.clumpIntensity = np.array(clumpIntensity)
  masspoints.clumpOpticalDepth = np.array(clumpOpticalDepth)

  for i in range(constants.interclumpLogMass.size):
    gridpoint = [masspoints.interclumpLogDensity[0,i], constants.interclumpLogMass[0,i], masspoints.logFUV]
    emission = masspointEmission(gridpoint)
    interclumpIntensity.append(emission[0])
    interclumpOpticalDepth.append(emission[1])

  masspoints.interclumpIntensity = np.array(interclumpIntensity)
  masspoints.interclumpOpticalDepth = np.array(interclumpOpticalDepth)

  return

#@jit(forceobj=False)
def masspointEmission(interpolationPoint, velocity=0, verbose=False, debug=False, test=False):
  '''
  This function calculates the emission of the given clump of the given mass.
  '''
  if debug:
    print('\n', interpolationPoint)

  indeces = species.molecules.getInterpolationIndeces()

  intensity_xi = []
  opticalDepth_xi = []

  if len(indeces):
    # Intensity currently in converted to Jansky, to coinside with the dust continuum
    intensity = interpolations.interpolateIntensity(interpolationPoint, indeces)
    tau = interpolations.interpolateTau(interpolationPoint, indeces)

    intensity_xi.append(intensity)
    opticalDepth_xi.append(tau)

  else:
    intensity_xi.append([])
    opticalDepth_xi.append([])

    # for i in range(len(indeces)):
    #   # Intensity currently in K
    #   if debug: print(np.exp(-(constants.wavelengths-species.moleculeWavelengths[i])**2/ \
    #                  2/(np.log(2))**2/species.moleculeWavelengths[i]**4))
    #   intensity_xi.append(intensity[i] * np.exp(-(constants.wavelengths-species.moleculeWavelengths[i])**2*velocity**2/ \
    #                                             2/constants.clumpDispersion**2/species.moleculeWavelengths[i]**2))
    #                                             #*np.exp(-1/2.*((velocityRange-velocityRange.T)/(constants.clumpDispersion))**2)))
    #   opticalDepth_xi.append(tau[i] * np.exp(-(constants.wavelengths-species.moleculeWavelengths[i])**2*velocity**2/ \
    #                                          2/constants.clumpDispersion**2/species.moleculeWavelengths[i]**2))
    #                                          #*np.exp(-1/2.*((velocityRange-velocityRange.T)/(constants.clumpDispersion))**2)))
    
  if debug: print(intensity)
  #tau = []
  if constants.dust:
    intensity = (interpolations.interpolateDustIntensity(interpolationPoint))#/2/constants.kB*(constants.wavelengths)**2*10**-26)
    tau = interpolations.interpolateDustTau(interpolationPoint)

    intensity_xi.append(intensity)#/2/constants.kB*(constants.wavelengths)**2*10**-26)
    opticalDepth_xi.append(tau)
  #intensity_xi.append(intensity)
  #opticalDepth_xi.append(tau)
  intensity = np.append(intensity_xi[1], intensity_xi[0])
  opticalDepth = np.append(opticalDepth_xi[1], opticalDepth_xi[0])

  return (intensity,opticalDepth)
    
  # def getEmission(self, verbose=False):
  #  if verbose:
  #    print('\nMasspoint intensity, optical depth\n{}\n{}'.format(self.__intensity, self.__opticalDepth))
  #    input()
  #  return (self.__intensity,self.__opticalDepth)
  # def getSpeciesEmission(self, number=-1, debug=False):
  #  if debug:
  #    print('\nMasspoint intensity, optical depth - xi:\n{}\n{}\n'.format(self.__intensity_xi, self.__opticalDepth_xi))
  #    input()
  #  if number==-1: return [self.__intensity_xi,self.__opticalDepth_xi]
  #  else: return [self.__intensity_xi[number],self.__opticalDepth_xi[number]]