import numpy as np
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
    self.__probability = 0            #the probability of this combination of masspoints
    return
  def __setFUV(self, fuvField):
    self.__FUV = fuvField
    return
  def __str__(self):
    return 'Combination {}:\n  ->probability {}\n  ->intensity {}\n  ->optical depth {}\n  ->FUV field {}'\
            .format(self.__combination, self.__probability, 10**sum(self.__intensity), 10**sum(self.__opticalDepth), self.__FUV)

  # PUBLIC
  #def addMolecule(self, element):
  #  self.__listMolecules.append(MolecularEmission(element))
  #  return
  #def addDust(self, element):
  #  self.__listDust.append(DustEmission(element))
  #  return
  def getMasspoints(self):
    return self.__masspoints
  def addFUV(self, fuvField):
    self.__setFUV(fuvField)
    return
  def addMasspoint(self, mass, number):
    self.__masspoints.append(Masspoint(self.__species, self.__interpolations, mass, number))
    self.__combination.append(number)
    return
  def calculateEmission(self, vrange, vDispersion):
    #print('Calculating combination emission')
    for i,masspoint in enumerate(self.__masspoints):
      (intensity,opticalDepth) = self.__combination[i]*masspoint.calculateEmission(vrange, vDispersion)
      self.__intensity.append(intensity)
      self.__opticalDepth.append(opticalDepth)
    self.__intensity = np.array(self.__intensity)
    self.__opticalDepth = np.array(self.__opticalDepth)
    return
  def getScaledCombinationEmission(self):
    return (self.__probability*np.stack((10**(self.__intensity), 10**(self.__opticalDepth))),self.__probability*self.__FUV.getFUV())