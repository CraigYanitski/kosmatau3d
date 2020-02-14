import numpy as np
import importlib as il

import combinations

import masspoints
# class Combination(object):
#   '''
#   This is a class to handle a combination of fractal masses in an ensemble.
#   It will have its associated probability, which will scale its intrinsic
#   intensity and optical depth. It returns a tuple of the combination's
#   probability, intensity, optical depth, and FUV field.
#   '''
#   # PRIVATE
#   def __init__(self, combination=[], masses=[], density=[], fuv=0, probability=0, maxProbability=1, debugging=False):
#     # self.__species = species     #list of both moleculular and dust species
#     self.__combination = combination 	#list of the number of each masspoint
#     self.__FUV = fuv                    #the FUV field for this combination of mass points
#     self.__density = density
#     self.__maxProbability = maxProbability      #the maximum probability (over all velocities) of this combination of masspoints
#     self.__probability = probability            #the probability of this combination of masspoints
#     self.__masspoints = []
#     for i,mass in enumerate(masses):
#       masspoint = Masspoint(density=self.__density[i], mass=mass, fuv=fuv, debugging=debugging)
#       self.__masspoints.append(masspoint)    #list of masses in the combination
#     self.__debugging = debugging
#     if debugging:
#       self.intensityList = []       #the list of masspoint intensities (OBSOLETE)
#       self.opticalDepthList = []    #the list of masspoint optical depths (OBSOLETE)
#       self.intensity = []             #velocity-averaged intensity of this combination of masspoints
#       self.opticalDepth = []          #velocity-averaged optical depth of this combination of masspoints
#     return

#   def __setFUV(self, fuvField):
#     self.__FUV = self.__probability*fuvField
#     return

#   def __str__(self):
#     return 'Combination {}:\n  ->probability {}\n  ->intensity {}\n  ->optical depth {}\n  ->FUV field {}'\
#             .format(self.__combination, self.__probability, sum(self.__intensity), sum(self.__opticalDepth), self.__FUV)

  # PUBLIC

  #def addMolecule(self, element):
  #  self.__listMolecules.append(MolecularEmission(element))
  #  return
  #def addDust(self, element):
  #  self.__listDust.append(DustEmission(element))
  #  return

# def reloadModules(self):
#   il.reload(Masspoint)
#   for masspoint in self.__masspoints: masspoint.reloadModules()

# def getMasspoints(self):
#   return self.__masspoints

# def getProbability(self):
#   return self.__probability

# def addFUV(self, fuvField):
#   self.__setFUV(fuvField)
#   return

def initialise(clumpCombination=[], interclumpCombination=[]):
  combinations.clumpCombination = clumpCombination
  combinations.interclumpCombination = interclumpCombination
  return

def setClumpCombination(clumpCombination=[]):
  combinations.clumpCombination = clumpCombination
  return

def setInterlumpCombination(interclumpCombination=[]):
  combinations.interclumpCombination = interclumpCombination
  return

def getAfuv(verbose=False):
  # The input is the density in the voxel. The probability should be included outside of this module.
  #  (velocity-independent)
  #Afuv = 0.
  #for i in range(len(constants.masspoints)):
  clAfuv,icAfuv = masspoints.getAfuv()
  clAfuv *= combinations.clumpCombination
  icAfuv *= combinations.interclumpCombination
  if verbose:
    print('clump: {}'.format(clAfuv))
    print('interclump: {}'.format(icAfuv))
  return (np.exp(-clAfuv.sum(1)), np.exp(-icAfuv.sum(1)))#(self.__probability[self.__probability.nonzero()].prod(),self.__probability[self.__probability.nonzero()].prod()*np.exp(-Afuv))

  # def addMasspoint(self, mass, number):
  #   self.__masspoints.append(Masspoint(mass, number))
  #   self.__combination.append(number)
  #   return

def calculateEmission(Afuv=0, probability=1, debug=False, test=False):
  '''
  This retrieves the emission from the masspoints, with dimensions (masspoints,species,velocity,velocity). It
  sums over the masspoint dimension. The probability will remain dormant until it is further-developed.
  '''
  #print('Calculating combination emission')
  # if isinstance(self.__intensity,np.ndarray):
  #   print('\nThe emission has already been calculated for this combination.\n')
  #   return
  if debug:
    print('combinations')
    print(combinations.clumpCombination)
    print(combinations.interclumpCombination)
    input()
  CLintensityList = []
  CLopticalDepthList = []
  ICintensityList = []
  ICopticalDepthList = []
  #for i,masspoint in enumerate(self.__masspoints):
  # if test:
  #   print(masspoint)
  #   CLemission,ICemission = masspoints.calculateEmission(Afuv=Afuv, test=False)
  #   input()
  # else:
  #   CLemission,ICemission = masspoints.calculateEmission(Afuv=Afuv)
  # intensityList.append(self.__combination[i]*intensity)
  # opticalDepthList.append(self.__combination[i]*opticalDepth)
  # intensityList = np.array(intensityList)
  # opticalDepthList = np.array(opticalDepthList)
  for c in combinations.clumpCombination:
    CLintensityList.append((c*masspoints.clumpIntensity.T).T)
    CLopticalDepthList.append((c*masspoints.clumpOpticalDepth.T).T)
  for c in combinations.interclumpCombination:
    ICintensityList.append(c*masspoints.interclumpIntensity)
    ICopticalDepthList.append(c*masspoints.interclumpOpticalDepth)
  if debug:
    print('Combination emissions')
    print('\n', CLintensityList[0], '\n\n', CLopticalDepthList[1])
    print('\n', ICintensityList[0], '\n\n', ICopticalDepthList[1])
    input()

  # if self.__probability.size==1:
  #   intensity = (self.__probability*intensityList.sum(0))
  #   opticalDepth = (self.__probability*np.exp(-opticalDepthList.sum(0)))
  # else:
    #print('Probability: {}\n'.format(self.__probability.shape))
    #print('Intensity: {}\n'.format(self.__intensityList.shape))
    #self.__probability.resize(self.__probability.size, 1)
  # intensity = []
  # opticalDepth = []
  # for element in range(len(species.molecules.getInterpolationIndeces()) + len(species.dust.getInterpolationIndeces())):
  #   intensity.append((probability*intensityList[:,element,:,:].sum(0)))
  #   opticalDepth.append((probability*np.exp(-opticalDepthList[:,element,:,:].sum(0))))
  # if self.__debugging:
  #   self.__intensity = np.array(intensity)
  #   self.__opticalDepth = np.array(opticalDepth)
  #   self.__intensity[self.__intensity==0] = 0
  #   self.__opticalDepth[self.__opticalDepth==0] = 1
  #else:
  combinations.clumpIntensity = np.array(CLintensityList).sum(1)
  combinations.clumpOpticalDepth = np.array(CLopticalDepthList).sum(1)
  combinations.interclumpIntensity = np.array(ICintensityList).sum(1)
  combinations.interclumpOpticalDepth = np.array(ICopticalDepthList).sum(1)

  combinations.clumpIntensity[combinations.clumpIntensity<=0] = 0
  combinations.clumpOpticalDepth[combinations.clumpOpticalDepth<=0] = 1
  combinations.interclumpIntensity[combinations.interclumpIntensity<=0] = 0
  combinations.interclumpOpticalDepth[combinations.interclumpOpticalDepth<=0] = 1

  # if debug:
  #   print(self.__intensity, self.__opticalDepth)
  #   input()
  # del intensityList
  # del opticalDepthList

  return

  # def getScaledCombinationEmission(self, verbose=False):
  # emission = (self.__intensity,self.__opticalDepth,self.__FUV.getFUV())
  # if verbose:
  #   print('\nCombination emission:\n', emission)
  #   input()
  # return emission
