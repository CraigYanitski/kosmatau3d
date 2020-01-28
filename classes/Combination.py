import numpy as np
import importlib as il
from Masspoint import *
class Combination(object):
  '''
  This is a class to handle a combination of fractal masses in an ensemble.
  It will have its associated probability, which will scale its intrinsic
  intensity and optical depth. It returns a tuple of the combination's
  probability, intensity, optical depth, and FUV field.
  '''
  # PRIVATE
  def __init__(self, species, interpolations, combination=[], masses=[], density=[], fuv=0, probability=0, debugging=False):
    self.__species = species     #list of both moleculular and dust species
    self.__interpolations = interpolations
    self.__combination = combination 	#list of the number of each masspoint
    self.__FUV = fuv                    #the FUV field for this combination of mass points
    self.__density = density
    self.__probability = probability            #the probability of this combination of masspoints
    self.__masspoints = []
    for i,mass in enumerate(masses):
      masspoint = Masspoint(self.__species, self.__interpolations, density=self.__density[i], mass=mass, fuv=fuv, debugging=debugging)
      self.__masspoints.append(masspoint)    #list of masses in the combination
    self.__debugging = debugging
    if debugging:
      self.intensityList = []       #the list of masspoint intensities (OBSOLETE)
      self.opticalDepthList = []    #the list of masspoint optical depths (OBSOLETE)
      self.intensity = []             #velocity-averaged intensity of this combination of masspoints
      self.opticalDepth = []          #velocity-averaged optical depth of this combination of masspoints
    return
  def __setFUV(self, fuvField):
    self.__FUV = self.__probability*fuvField
    return
  def __str__(self):
    return 'Combination {}:\n  ->probability {}\n  ->intensity {}\n  ->optical depth {}\n  ->FUV field {}'\
            .format(self.__combination, self.__probability, sum(self.__intensity), sum(self.__opticalDepth), self.__FUV)

  # PUBLIC
  #def addMolecule(self, element):
  #  self.__listMolecules.append(MolecularEmission(element))
  #  return
  #def addDust(self, element):
  #  self.__listDust.append(DustEmission(element))
  #  return
  def reloadModules(self):
    il.reload(Masspoint)
    for masspoint in self.__masspoints: masspoint.reloadModules()
  def getMasspoints(self):
    return self.__masspoints
  def getProbability(self):
    return self.__probability
  def addFUV(self, fuvField):
    self.__setFUV(fuvField)
    return
  def getAfuv(self):
    Afuv = 0.
    for i in range(len(self.__masspoints)):
      Afuv += self.__combination[i]*self.__masspoints[i].getAfuv()
    #print('{}, {}'.format(self.__probability[self.__probability.nonzero()], Afuv))
    return self.__probability[self.__probability.nonzero()]*np.exp(-Afuv)
  def addMasspoint(self, mass, number):
    self.__masspoints.append(Masspoint(self.__species, self.__interpolations, mass, number))
    self.__combination.append(number)
    return
  def calculateEmission(self, velocity, vDispersion, Afuv, debug=False, test=False):
    #print('Calculating combination emission')
    # if isinstance(self.__intensity,np.ndarray):
    #   print('\nThe emission has already been calculated for this combination.\n')
    #   return
    if debug:
      print(self.__combination)
      input()
    intensityList = []
    opticalDepthList = []
    for i,masspoint in enumerate(self.__masspoints):
      if test:
        print(masspoint)
        (intensity,opticalDepth) = masspoint.calculateEmission(velocity, vDispersion, Afuv, test=False)
        input()
      else:
        (intensity,opticalDepth) = masspoint.calculateEmission(velocity, vDispersion, Afuv)
      intensityList.append(self.__combination[i]*intensity)
      opticalDepthList.append(self.__combination[i]*opticalDepth)
    intensityList = np.array(intensityList)
    opticalDepthList = np.array(opticalDepthList)
    if debug:
      print('\nProbability:', self.__probability)
      print('\n', intensity, '\n\n', opticalDepth)
      input()
    if self.__probability.size==1:
      intensity = (self.__probability*intensityList.sum(0))
      opticalDepth = (self.__probability*np.exp(-opticalDepthList.sum(0)))
    else:
      #print('Probability: {}\n'.format(self.__probability.shape))
      #print('Intensity: {}\n'.format(self.__intensityList.shape))
      #self.__probability.resize(self.__probability.size, 1)
      intensity = []
      opticalDepth = []
      for element in range(len(self.__species[0].getInterpolationIndeces()) + len(self.__species[1].getInterpolationIndeces())):
        intensity.append((self.__probability*intensityList[:,element,:,:].sum(0)))
        opticalDepth.append((self.__probability*np.exp(-opticalDepthList[:,element,:,:].sum(0))))
      if self.__debugging:
        self.__intensity = np.array(intensity)
        self.__opticalDepth = np.array(opticalDepth)
        self.__intensity[self.__intensity==0] = 0
        self.__opticalDepth[self.__opticalDepth==0] = 1
      else:
        intensity = np.array(intensity)
        opticalDepth = np.array(opticalDepth)
        intensity[intensity<=0] = 0
        opticalDepth[opticalDepth<=0] = 1
    if debug:
      print(self.__intensity, self.__opticalDepth)
      input()
    # del intensityList
    # del opticalDepthList
    return (intensity,opticalDepth,self.__FUV.getFUV())
  def getScaledCombinationEmission(self, verbose=False):
    emission = (self.__intensity,self.__opticalDepth,self.__FUV.getFUV())
    if verbose:
      print('\nCombination emission:\n', emission)
      input()
    return emission