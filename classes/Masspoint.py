import numpy as np
import importlib as il
from Constants import *
from Molecules import *
from Dust import *
class Masspoint(object):
  '''
  This is a class to handle one fractal mass in a combination.
  It will have the associated emission and extinction information from the KOSMA-tau simulations,
  which will be used to contribute to the Combination class' intrinsic intensity and optical depth.
  At the moment, the emission interpolation is performed for each individual species. It will be
  an future update to perform the interpolation for all species at the same time.
  '''
  # PRIVATE
  def __init__(self, species, interpolations, density=0, mass=0, fuv=0, number=0):
    self.__species = species     #list of both moleculular and dust species
    self.__interpolations = interpolations
    self.__density = density
    self.__mass = mass 	#mass of this KOSMA-tau simulation
    self.__FUV = fuv
    self.__number = number    #number of this KOSMA-tau model contributing to the combination
    #self.__FUV = FUVfield()     #the FUV field for this combination of mass points
    self.__constants = Constants()
    self.__intensity_xi = []       #velocity-averaged intensity of this combination of masspoints
    self.__opticalDepth_xi = []    #velocity-averaged optical depth of this combination of masspoints
    return
  def __str__(self):
    return 'Simulated KOSMA-tau clump of mass {}'.format(10**self.mass)

  # PUBLIC
  def reloadModule(self):
    il.reload(Molecules)
    il.reload(Dust)
    return
  def calculateEmission(self, velocity, vDispersion, verbose=False, debug=False):
    velocity.resize((len(velocity), 1))
    speciesNumber = len(self.__species[0].getInterpolationIndeces()) + len(self.__species[1].getInterpolationIndeces())
    if self.__number==0:
      self.__intensity_xi = np.full((speciesNumber, len(velocity), len(velocity)), 10**-100)
      self.__opticalDepth_xi = np.full((speciesNumber, len(velocity), len(velocity)), 10**-100)
    else:
      for i,element in enumerate(self.__species):
        interpolationPoint = [self.__density, self.__mass, np.log10(self.__FUV.getFUV())]
        if debug:
          print(interpolationPoint)
          input()
        if verbose: print(element)
        if isinstance(element, Molecules):
          for index in element.getInterpolationIndeces():
            self.__intensity_xi.append(self.__interpolations.interpolateIntensity(interpolationPoint, [index])*self.__number*np.exp(-1/2.*((velocity-velocity.T)/(self.__constants.clumpDispersion))**2))
            self.__opticalDepth_xi.append(self.__interpolations.interpolateTau(interpolationPoint, [index])*self.__number*np.exp(-1/2.*((velocity-velocity.T)/(self.__constants.clumpDispersion))**2))
            #self.__intensity_xi[-1] = self.__intensity_xi[-1].sum(1)
            #self.__opticalDepth_xi[-1] = self.__opticalDepth_xi[-1].sum(1)
          if debug: input('intensity_xi:\n{}\n'.format(self.__intensity_xi[-1]))
        elif isinstance(element, Dust):
          for index in element.getInterpolationIndeces():
            self.__intensity_xi.append(np.full(len(velocity), self.__interpolations.interpolateIntensity(interpolationPoint, [index])*self.__number))
            self.__opticalDepth_xi.append(np.full(len(velocity), self.__interpolations.interpolateTau(interpolationPoint, [index])*self.__number))
      if speciesNumber>1:
        self.__intensity_xi = np.array(self.__intensity_xi)
        self.__opticalDepth_xi = np.array(self.__opticalDepth_xi)
      else:
        self.__intensity_xi = np.array([self.__intensity_xi])
        self.__opticalDepth_xi = np.array([self.__opticalDepth_xi])
    self.__intensity = (self.__intensity_xi).sum(0)
    self.__opticalDepth = (self.__opticalDepth_xi).sum(0)
    #self.__opticalDepth = -np.log((np.exp(-np.array(self.__opticalDepth_xi))).sum(0))
    if debug:
      print('\nIntensity, optical depth:\n{}\n{}\n'.format(self.__intensity, self.__opticalDepth))
      input()
    if np.isnan(self.__intensity).any():
      print('\nThere is an invalid intensity:\n', interpolationPoint)
      input()
    return
  def getEmission(self, verbose=False):
    if verbose:
      print('\nMasspoint intensity, optical depth\n{}\n{}'.format(self.__intensity, self.__opticalDepth))
      input()
    return (self.__intensity,self.__opticalDepth)
  def getSpeciesEmission(self, number=-1, debug=False):
    if debug:
      print('\nMasspoint intensity, optical depth - xi:\n{}\n{}\n'.format(self.__intensity_xi, self.__opticalDepth_xi))
      input()
    if number==-1: return [self.__intensity_xi,self.__opticalDepth_xi]
    else: return [self.__intensity_xi[number],self.__opticalDepth_xi[number]]