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
  def __init__(self, species, interpolations, combination=[], masses=[], density=[], fuv=0, probability=0):
    self.__species = species     #list of both moleculular and dust species
    self.__interpolations = interpolations
    self.__combination = combination 	#list of the number of each masspoint
    self.__FUV = fuv                    #the FUV field for this combination of mass points
    self.__density = density
    self.__masspoints = []
    for i,mass in enumerate(masses):
      masspoint = Masspoint(self.__species, self.__interpolations, density=self.__density[i], mass=mass, fuv=fuv, number=self.__combination[i])
      self.__masspoints.append(masspoint)    #list of masses in the combination
    self.__intensity = []             #velocity-averaged intensity of this combination of masspoints
    self.__opticalDepth = []          #velocity-averaged optical depth of this combination of masspoints
    self.__probability = probability            #the probability of this combination of masspoints
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
  def addMasspoint(self, mass, number):
    self.__masspoints.append(Masspoint(self.__species, self.__interpolations, mass, number))
    self.__combination.append(number)
    return
  def calculateEmission(self, velocity, vDispersion, emission='element', debug=False, test=False):
    #print('Calculating combination emission')
    if isinstance(self.__intensity,np.ndarray):
      print('\nThe emission has already been calculated for this combination.\n')
      return
    if debug:
      print(self.__combination)
      input()
    intensityList = []
    tauList = []
    for i,masspoint in enumerate(self.__masspoints):
      if test:
        print(masspoint)
        masspoint.calculateEmission(velocity, vDispersion, test=True)
        input()
      else: masspoint.calculateEmission(velocity, vDispersion)
      (intensity,opticalDepth) = masspoint.getSpeciesEmission()
      intensityList.append(intensity)
      tauList.append(opticalDepth)
    intensityList = np.array(intensityList)
    tauList = np.array(tauList)
    if debug:
      print('\nProbability:', self.__probability)
      print('\n', intensity, '\n\n', opticalDepth)
      input()
    if emission=='element':    #sum masspoints and prepare to average over combinations
      self.__intensity = (self.__probability*intensityList.sum(0))
      self.__opticalDepth = (self.__probability*np.exp(-tauList.sum(0)))
    elif emission=='all':
      self.__intensity = (self.__probability.T*iList).sum(0)
      self.__opticalDepth = -np.log(self.__probability.T*(np.exp(-tauList))).sum(0) # WRONG
    #self.__intensity = np.array(self.__intensity)
    #self.__opticalDepth = np.array(self.__opticalDepth)
    if debug:
      print(self.__intensity, self.__opticalDepth)
      input()
    return
  def getScaledCombinationEmission(self, verbose=False):
    emission = (self.__intensity,self.__opticalDepth,self.__FUV.getFUV())
    if verbose:
      print('\nCombination emission:\n', emission)
      input()
    return emission