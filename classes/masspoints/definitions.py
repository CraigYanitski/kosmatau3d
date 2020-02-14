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
                                    sum(10.**(constants.clumpLogMass*(2-constants.alpha)))*density[0]/1.91)
  masspoints.interclumpLogDensity = np.log10(10.**(constants.interclumpLogMass*(1-3./constants.gamma))*sum(10.**(constants.interclumpLogMass*(1+3./constants.gamma-constants.alpha))) / \
                                        sum(10.**(constants.interclumpLogMass*(2-constants.alpha)))*density[1]/1.91)
  masspoints.logFUV = np.log10(FUV)
  masspoint.clumpRadii = ((3./(4.*np.pi)*(10.**constants.clumplogMass*constants.massSolar)/ \
                              (10.**masspoints.clumpLogDensity*constants.massH*1.91))**(1./3.)/constants.pc)
  return

def getAfuv(debug=False):
  clumpAfuv = interpolations.interpolateFUVextinction(masspoints.clumpLogDensity, constants.clumpLogMass)
  interclumpAfuv = interpolations.interpolateFUVextinction(masspoints.interclumpLogDensity, constants.interclumpLogMass)
  if debug and self.__mass<0:
    print('\n', masspoints.clumpLogDensity, constants.clumpLogMass, clumpAfuv)
    print('\n', masspoints.interclumpLogDensity, constants.interclumpLogMass, interclumpAfuv)
  return (clumpAfuv,interclumpAfuv)

#@jit(forceobj=False)
def calculateEmission(Afuv=0, verbose=False, debug=False, test=False):
  #velocity.resize((len(velocity), 1))
  velocityRange = constants.velocityRange
  #velocityRange = np.linspace(velocity-3*vDispersion, velocity+3*vDispersion, num=7)      #a range of 7 is used to account for the observed velocity +/- 3 sigma
  velocityRange = np.resize(velocityRange, (1, velocityRange.size))
  #print(velocity)
  if debug:
    print('Masspoint velocity argument:\n{}'.format(velocity))
    print('Masspoint velocity range variable:\n{}'.format(constants.velocityRange))
    print('Masspoint velocity difference result:\n{}'.format(velocityRange-velocity))
    #input()
  speciesNumber = len(species.molecules.getInterpolationIndeces()) + len(species.dust.getInterpolationIndeces())
  # if self.__number==0:
  #   intensity_xi = np.full((speciesNumber, velocityRange.size, velocityRange.size), 10**-100)
  #   opticalDepth_xi = np.full((speciesNumber, velocityRange.size, velocityRange.size), 10**-100)
  # else:
  #interpolationPoint = [masspoints.clumpDensity, constants.clumpLogMass, max(np.log10(self.__FUV.getFUV())-Afuv/2.5, 0)]
  if debug==False:
    if test: print('\n', interpolationPoint)
    #input()
  CLinterpolationPoints = [masspoints.clumpLogdensity, constants.clumpLogMass, masspoints.logFUV]
  CLintensity_xi = []
  CLopticalDepth_xi = []

  ICinterpolationPoints = [masspoints.interclumpLogdensity, constants.interclumpLogMass, masspoints.logFUV]
  ICintensity_xi = []
  ICopticalDepth_xi = []
  #for i,element in enumerate(self.__species):      #This is commented-out since numba cannot support the isinstance() function
  #  if verbose: print(element)
  #  if isinstance(element, Molecules):
  for index in species.molecules.getInterpolationIndeces():
    # Intensity currently in K
    CLintensity_xi.append((interpolations.interpolateIntensity(CLinterpolationPoints, [index])*np.exp(-1/2.*((velocityRange-velocityRange.T)/(constants.clumpDispersion))**2)))
    CLopticalDepth_xi.append((interpolations.interpolateTau(CLinterpolationPoints, [index])*np.exp(-1/2.*((velocityRange-velocityRange.T)/(constants.clumpDispersion))**2)))
    
    ICintensity_xi.append((interpolations.interpolateIntensity(ICinterpolationPoints, [index])*np.exp(-1/2.*((velocityRange-velocityRange.T)/(constants.clumpDispersion))**2)))
    ICopticalDepth_xi.append((interpolations.interpolateTau(ICinterpolationPoints, [index])*np.exp(-1/2.*((velocityRange-velocityRange.T)/(constants.clumpDispersion))**2)))
    #self.__intensity_xi[-1] = self.__intensity_xi[-1].sum(1)
    #self.__opticalDepth_xi[-1] = self.__opticalDepth_xi[-1].sum(1)
    if test:
      print('\n{}\n'.format(self.__intensity_xi[-1].max()))
  if debug:
    print('intensity_xi:\n{}\n'.format(self.__intensity_xi[-1]))
  #  elif isinstance(element, Dust):
  for i,index in enumerate(species.dust.getInterpolationIndeces()):
    # Intensity currently in converted to K, to coinside with the molecule emission
    CLintensity_xi.append(np.full((velocityRange.size, velocityRange.size), interpolations.interpolateIntensity(CLinterpolationPoints, [index]))/2/constants.kB*constants.c**2/(species.dust.getFrequencies()[i]*10**9)**2*10**-26)
    CLopticalDepth_xi.append(np.full((velocityRange.size, velocityRange.size), interpolations.interpolateTau(CLinterpolationPoints, [index])))
    
    ICintensity_xi.append(np.full((velocityRange.size, velocityRange.size), interpolations.interpolateIntensity(ICinterpolationPoints, [index]))/2/constants.kB*constants.c**2/(species.dust.getFrequencies()[i]*10**9)**2*10**-26)
    ICopticalDepth_xi.append(np.full((velocityRange.size, velocityRange.size), interpolations.interpolateTau(ICinterpolationPoints, [index])))
  # del interpolationPoint
  #if speciesNumber>1:
    # if self.__debugging:
    #   self.intensity_xi = self.number*np.array(intensity_xi)
    #   self.opticalDepth_xi = self.number*np.array(opticalDepth_xi)
  CLintensity_xi = np.array(CLintensity_xi)
  CLopticalDepth_xi = np.array(CLopticalDepth_xi)

  ICintensity_xi = np.array(ICintensity_xi)
  ICopticalDepth_xi = np.array(ICopticalDepth_xi)
  #else:
    # if self.__debugging:
    #   self.intensity_xi = np.array([intensity_xi])
    #   self.opticalDepth_xi = np.array([opticalDepth_xi])
    #   if np.isnan(self.__intensity).any():
    #     print('\nThere is an invalid intensity:\n', interpolationPoint)
    #     #input()
    # else:
    #   intensity_xi = np.array([intensity_xi])
    #   opticalDepth_xi = np.array([opticalDepth_xi])
  if np.isnan(intensity).any():
    print('\nThere is an invalid intensity:\n', interpolationPoint)
    #input()
  if debug:
    np.set_printoptions(threshold=100000)
    print('\nIntensity xi, optical depth xi:\n{}\n{}\n{}\n{}\n'.format(self.__intensity_xi.shape, self.__intensity_xi[:3,:,:], self.__opticalDepth_xi.shape, self.__opticalDepth_xi[:3,:,:]))
    #input()
  return ((CLintensity_xi,CLopticalDepth_xi), (ICintensity_xi,ICopticalDepth_xi))
    
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