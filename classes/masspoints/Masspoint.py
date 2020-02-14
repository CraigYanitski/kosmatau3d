import numpy as np
from numba import jit
import importlib as il

import constants
import interpolations
import species

class Masspoint(object):
  '''
  This is a class to handle one fractal mass in a combination.
  It will have the associated emission and extinction information from the KOSMA-tau simulations,
  which will be used to contribute to the Combination class' intrinsic intensity and optical depth.
  At the moment, the emission interpolation is performed for each individual species. It will be
  an future update to perform the interpolation for all species at the same time.
  '''
  # PRIVATE

  def __init__(self, density=0, mass=0, fuv=0, number=1, debugging=False):
    # self.__species = species     #list of both moleculular and dust species
    self.__density = density
    self.__mass = mass 	#mass of this KOSMA-tau simulation
    #self.__Afuv = self.__interpolations.interpolateFUVextinction(density, mass)
    self.__FUV = fuv
    #input('{}: {}'.format(mass, number))
    #self.__FUV = FUVfield()     #the FUV field for this combination of mass points
    self.__debugging = debugging
    if debugging:
      self.number = number    #number of this KOSMA-tau model contributing to the combination (OBSOLETE)
      self.intensity_xi = []       #velocity-averaged intensity of this combination of masspoints
      self.opticalDepth_xi = []    #velocity-averaged optical depth of this combination of masspoints
    return

  def __str__(self):
    return 'Simulated KOSMA-tau clump of mass {}'.format(10**float(self.__mass))

  # PUBLIC

  def getAfuv(self, debug=False):
    Afuv = interpolations.interpolateFUVextinction(self.__density, self.__mass)
    if debug and self.__mass<0:
      print('\n', self.__density, self.__mass, Afuv)
    return Afuv

  #@jit(forceobj=False)
  def calculateEmission(self, velocity, vDispersion, Afuv, verbose=False, debug=False, test=False):
    #velocity.resize((len(velocity), 1))
    velocityRange = constants.velocityRange
    #velocityRange = np.linspace(velocity-3*vDispersion, velocity+3*vDispersion, num=7)      #a range of 7 is used to account for the observed velocity +/- 3 sigma
    velocityRange = np.resize(velocityRange, (1, velocityRange.size))
    #print(velocity)
    if debug:
      print('Masspoint velocity argument:\n{}'.format(velocity))
      print('Masspoint velocity range variable:\n{}'.format(velocityRange))
      print('Masspoint velocity difference result:\n{}'.format(velocityRange-velocity))
      #input()
    speciesNumber = len(species.molecules.getInterpolationIndeces()) + len(species.dust.getInterpolationIndeces())
    # if self.__number==0:
    #   intensity_xi = np.full((speciesNumber, velocityRange.size, velocityRange.size), 10**-100)
    #   opticalDepth_xi = np.full((speciesNumber, velocityRange.size, velocityRange.size), 10**-100)
    # else:
    interpolationPoint = [self.__density, self.__mass, max(np.log10(self.__FUV.getFUV())-Afuv/2.5, 0)]
    if debug==False:
      if test: print('\n', interpolationPoint)
      #input()
    intensity_xi = []
    opticalDepth_xi = []
    #for i,element in enumerate(self.__species):      #This is commented-out since numba cannot support the isinstance() function
    #  if verbose: print(element)
    #  if isinstance(element, Molecules):
    for index in species.molecules.getInterpolationIndeces():
      # Intensity currently in K
      intensity_xi.append((interpolations.interpolateIntensity(interpolationPoint, [index])*np.exp(-1/2.*((velocityRange-velocityRange.T)/(constants.clumpDispersion))**2)))
      opticalDepth_xi.append((interpolations.interpolateTau(interpolationPoint, [index])*np.exp(-1/2.*((velocityRange-velocityRange.T)/(constants.clumpDispersion))**2)))
      #self.__intensity_xi[-1] = self.__intensity_xi[-1].sum(1)
      #self.__opticalDepth_xi[-1] = self.__opticalDepth_xi[-1].sum(1)
      if test:
        print('\n{}\n'.format(self.__intensity_xi[-1].max()))
    if debug:
      print('intensity_xi:\n{}\n'.format(self.__intensity_xi[-1]))
    #  elif isinstance(element, Dust):
    for i,index in enumerate(species.dust.getInterpolationIndeces()):
      # Intensity currently in converted to K, to coinside with the molecule emission
      intensity_xi.append(np.full((velocityRange.size, velocityRange.size), interpolations.interpolateIntensity(interpolationPoint, [index]))/2/constants.kB*constants.c**2/(species.dust.getFrequencies()[i]*10**9)**2*10**-26)
      opticalDepth_xi.append(np.full((velocityRange.size, velocityRange.size), interpolations.interpolateTau(interpolationPoint, [index])))
    # del interpolationPoint
    if speciesNumber>1:
      if self.__debugging:
        self.intensity_xi = self.number*np.array(intensity_xi)
        self.opticalDepth_xi = self.number*np.array(opticalDepth_xi)
      else:
        intensity_xi = np.array(intensity_xi)
        opticalDepth_xi = np.array(opticalDepth_xi)
    else:
      if self.__debugging:
        self.intensity_xi = np.array([intensity_xi])
        self.opticalDepth_xi = np.array([opticalDepth_xi])
        if np.isnan(self.__intensity).any():
          print('\nThere is an invalid intensity:\n', interpolationPoint)
          #input()
      else:
        intensity_xi = np.array([intensity_xi])
        opticalDepth_xi = np.array([opticalDepth_xi])
        if np.isnan(intensity).any():
          print('\nThere is an invalid intensity:\n', interpolationPoint)
          #input()
    if debug:
      np.set_printoptions(threshold=100000)
      print('\nIntensity xi, optical depth xi:\n{}\n{}\n{}\n{}\n'.format(self.__intensity_xi.shape, self.__intensity_xi[:3,:,:], self.__opticalDepth_xi.shape, self.__opticalDepth_xi[:3,:,:]))
      #input()
    return (intensity_xi,opticalDepth_xi)
    
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