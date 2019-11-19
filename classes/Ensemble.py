import importlib as il
import numpy as np
import sys
from Combination import *
from Binomial import *
from Gauss import *
from Poisson import *
from Constants import *
from FUVfield import *
class Ensemble(object):
  '''
  This class owes itself largely to the work  done by Silke Andree-Labsch and
  Christoph Bruckmann. It contains a number of combinations of fractal mass
  configurations, and therefore contains the ensemble-averaged intensity and
  optical depth.
  It can be either a clump or an interclump ensemble.
  '''
  # PRIVATE
  def __init__(self, clumpType, species, interpolations, combination='binomial', verbose=False):
    self.__species = species     #list of both moleculular and dust species
    self.__interpolations = interpolations
    self.__clumpType = clumpType     #type of mass in ensemble (clump or interclump medium)
    self.__constants = Constants()
    self.__flagCombination = combination
    self.__verbose = verbose
    self.__combinations = []    #list of combinations
    self.__combinationObjects = []    #list of instances of Combination()
    self.__probability = []
    self.__intensity = 0        #ensemble-averaged intensity for each velocity
    self.__opticalDepth = 0     #ensemble-averaged optical depth for each velocity
    self.__FUV = 0
    self.__massObserved = 0     #'observed' mass
    self.__massEnsemble = 0     #ensemble mass
    self.__radiusEnsemble = 0   #ensemble radius
    self.__volumeEnsemble = 0   #ensemble volume
    self.__densityObserved = 0
    self.__densityEnsemble = 0
    if self.__clumpType=='clump':
      self.__massLimits = self.__constants.clumpMassLimits   #maximum mass in ensemble
    elif self.__clumpType=='interclump':
      self.__massLimits = self.__constants.interclumpMassLimits   #maximum mass in ensemble
    self.__masspoints = []
    self.__masspointNumberRange = []
    self.__masspointRadii = []
    self.__masspointDensity = []
    self.__interpolationPoints = []
    self.__Nj = []    #this is 'deltaNji' in the original code, and comepletely confusing regarding the inconsistency with Silke's thesis
    return
  def __setMass(self, mass):
    '''Set the mass.'''
    self.__massObserved = mass
    return
  def __setDensity(self, density):
    '''Set the mass.'''
    self.__densityObserved = density
    return
  def __setFUV(self, fuv):
    self.__FUV = fuv
    return
  def __setExtinction(self, afuv):
    self.__extinction = afuv
    return
  def __setVelocity(self, velocity):
    self.__velocity = velocity
    self.__velocityStep = abs(velocity[-1]-velocity[-2])
    return
  def __setVelocityDispersion(self, velocityDispersion):
    self.__velocityDispersion = velocityDispersion
    return
  def __str__(self):
    return 'The {} ensemble to simulate the fractal structure with {} instances of KOSMA-tau giving\n  {} possible combinations in the line-of-sight of an observer:'\
            .format(self.__clumpType, self.__masspoints.size, len(self.__combinations), self.__combinations)

  # PUBLIC
  def reloadModules(self):
    il.reload(Combination)
    il.reload(Binomial)
    il.reload(Gauss)
    il.reload(Poisson)
    il.reload(FUVfield)
    il.reload(Constants)
    for combination in self.__combinationObjects:
      combination.reloadModules()
    return
  def initialise(self, mass=0, density=0, velocity=0, velocityDispersion=0, FUV=0, extinction=0):
    self.__setMass(mass)
    self.__setDensity(density)
    self.__setVelocity(velocity)
    self.__setVelocityDispersion(velocityDispersion)
    self.__setExtinction(extinction)
    self.__setFUV(FUV)
    self.initialiseEnsemble()
    return
  def setMass(self, mass):
    '''Set the mass.'''
    self.__setMass(mass)
    return
  def getMass(self):
    return self.__massObserved
  #def initialiseEnsemble
  def calculateMasspoints(self):
    '''This is a function to get the clump masses in this ensemble. It will soon be depreciated
       as I will change the reading of the KOSMA-tau files to convert from the fortran 10*log(values).'''
    self.__masspoints = np.arange(self.__massLimits[0], self.__massLimits[1])
    return
  def calculateMasspointDensity(self):
    '''I am not entirely sure this is needed...'''
    return np.log10(self.__masspointDensity)
  def getInterpolationPoints(self):
    '''This might be removed. I'm not certain the interpolation points need to be accessed out-of-class.'''
    return
  def getMasspoints(self):
    '''This is pretty self-explanatory. Return the list of clump masses when called...'''
    return self.__masspoints
  def getCombinations(self):
    return self.__combinations
  def calculateRadii(self, verbose=False):
    '''This function calculates the interpolation points necessary for reading the KOSMA-tau files.'''
    self.__masspointDensity = np.log10(10.**(self.__masspoints*(1-3./self.__constants.gamma))*sum(10.**(self.__masspoints*(1+3./self.__constants.gamma-self.__constants.alpha))) / \
                                       sum(10.**(self.__masspoints*(2-self.__constants.alpha)))*self.__densityObserved/1.91)
    if verbose: print(self.__masspointDensity)
    #if (self.__masspointDensity>self.__constants.densityLimits[1]).any() or (self.__masspointDensity<self.__constants.densityLimits[0]).any(): sys.exit('WARNING: surface density {} outside of KOSMA-tau grid!\n\n'.format(self.__masspointDensity))
    self.__interpolationPoints = np.stack((self.__masspointDensity, self.__masspoints, np.full(self.__masspoints.size, self.__FUV)))
    self.__masspointRadii = (3./(4.*np.pi)*(10.**self.__masspoints*self.__constants.massSolar)/(10.**self.__masspointDensity*self.__constants.massH*1.91))**(1./3.)/self.__constants.pc
    #print(self.__masspointRadii, 'cm')
    return
  def getRadii(self):
    return self.__masspointRadii
  def getDensity(self):
    return self.__masspointDensity
  def calculateCombinations(self, verbose=False):
    '''This function calculates all of the different combinations of clump masses that may be in a line-of-sight.
      It is basically the essence of the probabilistic approach used to create the superposition of clumps.'''
    dimension = len(self.__masspointNumberRange[0])
    ranges = self.__masspointNumberRange
    if verbose: print('\nMasspoint Ranges:\n', ranges)
    ranges[:,:,1] += 1
    combinations = []
    for i in range(len(ranges)):
      if dimension==1:
        grid = np.mgrid[ranges[i,0,0]:ranges[i,0,1]]
        combinations.append(np.array([grid[0].flatten()], dtype=np.int))
      elif dimension==2:
        grid = np.mgrid[ranges[i,0,0]:ranges[i,0,1], ranges[i,1,0]:ranges[i,1,1]]
        combinations.append(np.array([grid[0].flatten(), grid[1].flatten()], dtype=np.int))
      elif dimension==3:
        grid = np.mgrid[ranges[i,0,0]:ranges[i,0,1], ranges[i,1,0]:ranges[i,1,1], ranges[i,2,0]:ranges[i,2,1]]
        combinations.append(np.array([grid[0].flatten(), grid[1].flatten(), grid[2].flatten()], dtype=np.int))
      elif dimension==4:
        grid = np.mgrid[ranges[i,0,0]:ranges[i,0,1], ranges[i,1,0]:ranges[i,1,1], ranges[i,2,0]:ranges[i,2,1], ranges[i,3,0]:ranges[i,3,1]]
        combinations.append(np.array([grid[0].flatten(), grid[1].flatten(), grid[2].flatten(), grid[3].flatten()], dtype=np.int))
      elif dimension==5:
        grid = np.mgrid[ranges[i,0,0]:ranges[i,0,1], ranges[i,1,0]:ranges[i,1,1], ranges[i,2,0]:ranges[i,2,1], ranges[i,3,0]:ranges[i,3,1], \
                        ranges[i,4,0]:ranges[i,4,1]]
        combinations.append(np.array([grid[0].flatten(), grid[1].flatten(), grid[2].flatten(), grid[3].flatten(), grid[4].flatten()], dtype=np.int))
      elif dimension==6:
        grid = np.mgrid[ranges[i,0,0]:ranges[i,0,1], ranges[i,1,0]:ranges[i,1,1], ranges[i,2,0]:ranges[i,2,1], ranges[i,3,0]:ranges[i,3,1], \
                        ranges[i,4,0]:ranges[i,4,1], ranges[i,5,0]:ranges[i,5,1]]
        combinations.append(np.array([grid[0].flatten(), grid[1].flatten(), grid[2].flatten(), grid[3].flatten(), grid[4].flatten(), grid[5].flatten()], dtype=np.int))
      elif dimension==7:
        grid = np.mgrid[ranges[i,0,0]:ranges[i,0,1], ranges[i,1,0]:ranges[i,1,1], ranges[i,2,0]:ranges[i,2,1], ranges[i,3,0]:ranges[i,3,1], \
                        ranges[i,4,0]:ranges[i,4,1], ranges[i,5,0]:ranges[i,5,1], ranges[i,6,0]:ranges[i,6,1]]
        combinations.append(np.array([grid[0].flatten(), grid[1].flatten(), grid[2].flatten(), grid[3].flatten(), grid[4].flatten(), grid[5].flatten(), grid[6].flatten()], dtype=np.int))
      else: sys.exit('\nThere are too many masses for the current grid ({}).\nExitting. . .\n\n'.format(dimension))
    #combinations = np.array(combinations)
    if verbose: print('\nCalculated combinations:\n', combinations)
    return combinations
  def createCombinationObjects(self, verbose=False):
    '''This function removes all of the unnecessary degenerate looping during this calculation.
       Of course it is possible because of the wonders of numpy.ndarray(). . .'''
    verbose = self.__verbose or verbose
    self.__Nj = (self.__massObserved*10.**(self.__masspoints*(1-self.__constants.alpha))) / sum(10.**(self.__masspoints*(2-self.__constants.alpha)))
    #print('\nNj:\n', self.__Nj)
    self.__massEnsemble = sum(self.__Nj*10.**self.__masspoints)
    self.__radiusEnsemble = sum(self.__Nj*self.__masspointRadii)
    self.__volumeEnsemble = sum(self.__Nj*np.pi*4./3.*self.__masspointRadii**3)
    self.__densityEnsemble = self.__massEnsemble/self.__volumeEnsemble
    if verbose: print('velocity, mean, dispersion', self.__velocity, self.__velocity.mean(), self.__velocityDispersion)
    self.__deltaNji = np.array([self.__Nj]).T/np.sqrt(2*np.pi)/self.__velocityDispersion*(np.exp(-0.5*((self.__velocity-self.__velocity.mean())/self.__velocityDispersion)**2)).T*self.__velocityStep
    surfaceProbability = np.array([np.pi*self.__masspointRadii**2/self.__constants.pixelWidth**2])    #this is 'pTab' in the original code
    probableNumber = (self.__deltaNji*surfaceProbability.T)   #this is 'expectedValTab' in the original code
    try: standardDeviation = np.sqrt(self.__deltaNji*surfaceProbability.T*(1-surfaceProbability.T))    #this is 'standardDeriTab' in the original code
    except ValueError:
      input('\nObserved mass, sufaceProbability, standardDeviation**2:\n', self.__massObserved, surfaceProbability, '\n', self.__deltaNji*surfaceProbability.T*(1-surfaceProbability.T))
    if verbose: print('\nsuface probability, expected number, standard deviation:\n', surfaceProbability, '\n', probableNumber, '\n', standardDeviation)
    #print(surfaceProbability, probableNumber, standardDeviation)
    lower = np.maximum(np.zeros([self.__masspoints.size, self.__velocity.size]), np.floor(probableNumber-self.__constants.nSigma*standardDeviation))
    upper = np.minimum(self.__deltaNji, np.ceil(probableNumber+self.__constants.nSigma*standardDeviation))
    if verbose: print('\nupper,lower:\n',upper,'\n',lower)
    self.__masspointNumberRange = np.array([lower, upper]).T
    if verbose: print('\nMasspoint number range:\n', self.__masspointNumberRange)
    self.__combinations = self.calculateCombinations()
    #self.__probability = np.zeros((len(self.__combinations),2))
    if verbose: print('\nCombinations:\n', self.__combinations)
    for i,combinations in enumerate(self.__combinations):
      if verbose: print(combinations)
      self.__combinations[i] = np.array(combinations).T
      probability = []
      if verbose: print('\nEnsemble combinations:\n', combinations)
      for combination in self.__combinations[i]:
        combination = np.array([combination])
        if verbose:
          print('\nCombination:\n', combination)
          input()
        if self.__flagCombination=='binomial':
          if np.any(probableNumber>self.__constants.pnGauss) and np.any(self.__deltaNji>self.__constants.nGauss):
            # use gauss!
            g = Gauss(probableNumber, standardDeviation)
            probability.append(g.gaussfunc(combination))
            #print('gauss!!...')
          else:
            # use binomial 
            # <<This will likely print an error when there are more masspoints>>
            b = Binomial(self.__deltaNji, surfaceProbability, debug=verbose) # n and p for binominal 
            probability.append(b.binomfunc(combination))
        elif self.__flagCombination=='poisson':
          if np.any(probableNumber>self.__constants.pnGauss) and np.any(self.__deltaNji>self.__constants.nGauss):
            # use gauss
            g = Gauss(probableNumber, standardDeviation)
            probability.append(g.gaussfunc(combination))
            #print('gauss!!...')
          else:
            # use poisson
            po = Poisson(probableNumber)
            probability.append(po.poissonfunc(combination))
        else: probability.append([0, 0])
        if verbose:
          print(probability[-1])
          input()
        if (probability[-1]==np.nan).any():
          print('\nThere is an invalid probability:\n')
          input()
        self.__combinationObjects.append(Combination(self.__species, self.__interpolations, combination=combination.flatten(), masses=self.__masspoints, density=self.__masspointDensity, fuv=self.__FUV, probability=probability[-1]))
      self.__probability.append(probability)
    #for i,combination in enumerate(self.__combinations): self.__probability[i] = self.__probability[i](combination)
    if verbose: print('Probability:', self.__probability)
    return
  def initialiseEnsemble(self):

    self.calculateMasspoints()
    self.calculateRadii()
    self.createCombinationObjects()
    if self.__verbose: print(self.__clumpType)
    return
  def calculate(self, debug=False):
    '''Maybe <<PARALLELISE>> this??

       This is a function to cycle through the Combination instances to create a large numpy.ndarray,
       which is used to calculate the final sums needed for the voxel.'''
    emissionResult = []
    FUVresult = []
    if debug:
      print('Calculating {} ensemble emission'.format(self.__clumpType))
      input()
    for combination in self.__combinationObjects:
      combination.calculateEmission(self.__velocity, self.__velocityDispersion)
      result = combination.getScaledCombinationEmission() #<<this needs to be altered>>
      emissionResult = result[0]
      self.__FUV = result[1]
    emissionResult = np.array(emissionResult)
    if debug:
      print(emissionResult)
      input()
    self.__intensity = emissionResult[0]#.sum(3)[0]
    self.__opticalDepth = emissionResult[1]#.sum(3)[1]
    if debug:
      print('\nIntensity\n', self.__intensity, '\nOptical depth\n', self.__opticalDepth)
      input()
    return
  def getEnsembleEmission(self):
    '''This returns the ensemble emission...nothing more.'''
    if (self.__intensity==np.nan).any():
      print('\nInvalid intensity:\n', self.__intensity)
      input()
    return (self.__intensity,self.__opticalDepth,self.__FUV)
  def getCombinations(self):
    return self.__combinations
  def getCombinationObjects(self):
    return self.__combinationObjects