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
  def __init__(self, clumpType, species, interpolations, combination='binomial'):
    self.__species = species     #list of both moleculular and dust species
    self.__interpolations = interpolations
    self.__clumpType = clumpType     #type of mass in ensemble (clump or interclump medium)
    self.__constants = Constants()
    self.__flagCombination = combination
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
    return
  def __setVelocityDispersion(self, velocityDispersion):
    self.__velocityDispersion = velocityDispersion
    return
  def __str__(self):
    return 'The {} ensemble to simulate the fractal structure with {} instances of KOSMA-tau giving\n  {} possible combinations in the line-of-sight of an observer'\
            .format(self.__clumpType, self.__masspoints.size, self.__combinations.shape[0])

  # PUBLIC
  def initialise(self, mass=0, density=0, velocity=0, velocityDispersion=0, FUV=0, extinction=0):
    self.__setMass(mass)
    self.__setDensity(density)
    self.__setVelocity(velocity)
    self.__setVelocityDispersion(velocityDispersion)
    self.__setExtinction(extinction)
    self.__setFUV(FUV)
    self.calculate()
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
  def calculateRadii(self):
    '''This function calculates the interpolation points necessary for reading the KOSMA-tau files.'''
    self.__masspointDensity = np.log10(self.__masspoints**(1-3/self.__constants.gamma)*sum(self.__masspoints**(1+3/self.__constants.gamma-self.__constants.alpha)) / \
                                       sum(self.__masspoints**(2-self.__constants.alpha))*self.__densityObserved/1.91)
    if (self.__masspointDensity>self.__constants.densityLimits[1]).any() or (self.__masspointDensity<self.__constants.densityLimits[0]).any(): sys.exit('WARNING: surface density outside of KOSMA-tau grid!\n\n')
    self.__interpolationPoints = np.stack((self.__masspointDensity, self.__masspoints, np.full(self.__masspoints.size, self.__FUV)))
    self.__masspointRadii = ((3/(4.0*np.pi)*(self.__masspoints*self.__constants.massSolar)/(self.__interpolationPoints[0]*self.__constants.massH*1.91))**(1./3.))/self.__constants.pc
    return
  def calculateCombinations(self):
    '''This function calculates all of the different combinations of clump masses that may be in a line-of-sight.
      It is basically the essence of the probabilistic approach used to create the superposition of clumps.'''
    dimension = len(self.__masspointNumberRange)
    ranges = self.__masspointNumberRange
    ranges[:,1] += 1
    if dimension==1:
      grid = np.mgrid[ranges[0][0]:ranges[0][1]]
      combinations = np.array([grid[0].flatten()], dtype=np.int)
    elif dimension==2:
      grid = np.mgrid[ranges[0][0]:ranges[0][1], ranges[1][0]:ranges[1][1]]
      combinations = np.array([grid[0].flatten(), grid[1].flatten()], dtype=np.int)
    elif dimension==3:
      grid = np.mgrid[ranges[0][0]:ranges[0][1], ranges[1][0]:ranges[1][1], ranges[2][0]:ranges[2][1]]
      combinations = np.array([grid[0].flatten(), grid[1].flatten(), grid[2].flatten()], dtype=np.int)
    elif dimension==4:
      grid = np.mgrid[ranges[0][0]:ranges[0][1], ranges[1][0]:ranges[1][1], ranges[2][0]:ranges[2][1], ranges[3][0]:ranges[3][1]]
      combinations = np.array([grid[0].flatten(), grid[1].flatten(), grid[2].flatten(), grid[3].flatten()], dtype=np.int)
    elif dimension==5:
      grid = np.mgrid[ranges[0][0]:ranges[0][1], ranges[1][0]:ranges[1][1], ranges[2][0]:ranges[2][1], ranges[3][0]:ranges[3][1], \
                      ranges[4][0]:ranges[4][1]]
      combinations = np.array([grid[0].flatten(), grid[1].flatten(), grid[2].flatten(), grid[3].flatten(), grid[4].flatten()], dtype=np.int)
    elif dimension==6:
      grid = np.mgrid[ranges[0][0]:ranges[0][1], ranges[1][0]:ranges[1][1], ranges[2][0]:ranges[2][1], ranges[3][0]:ranges[3][1], \
                      ranges[4][0]:ranges[4][1], ranges[5][0]:ranges[5][1]]
      combinations = np.array([grid[0].flatten(), grid[1].flatten(), grid[2].flatten(), grid[3].flatten(), grid[4].flatten(), grid[5].flatten()], dtype=np.int)
    elif dimension==7:
      grid = np.mgrid[ranges[0][0]:ranges[0][1], ranges[1][0]:ranges[1][1], ranges[2][0]:ranges[2][1], ranges[3][0]:ranges[3][1], \
                      ranges[4][0]:ranges[4][1], ranges[5][0]:ranges[5][1], ranges[6][0]:ranges[6][1]]
      combinations = np.array([grid[0].flatten(), grid[1].flatten(), grid[2].flatten(), grid[3].flatten(), grid[4].flatten(), grid[5].flatten(), grid[6].flatten()], dtype=np.int)
    else: sys.exit('\nThere are too many masses for the current grid ({}).\nExitting. . .\n\n'.format(dimension))
    return combinations.T
  def createCombinationObjects(self):
    '''This function removes all of the unnecessary degenerate looping during this calculation.
       Of course it is possible because of the wonders of numpy.ndarray(). . .'''
    self.__Nj = (self.__massObserved*10**(self.__masspoints.astype(float)*(1-self.__constants.alpha))) / sum(10**(self.__masspoints.astype(float)*(2-self.__constants.alpha)))
    #print('\nNj:\n', self.__Nj)
    self.__massEnsemble = sum(self.__Nj*self.__masspoints)
    self.__radiusEnsemble = sum(self.__Nj*self.__masspointRadii)
    self.__volumeEnsemble = sum(self.__Nj*np.pi*4./3.*self.__masspointRadii**3)
    self.__densityEnsemble = self.__massEnsemble/self.__volumeEnsemble
    apparentSurfaceDensity = (np.pi*self.__masspointRadii**2/self.__constants.pixelWidth**2)    #this is 'pTab' in the original code
    expectedMass = (self.__Nj*apparentSurfaceDensity)   #this is 'expectedValTab' in the original code
    standardDeviation = ((self.__Nj*apparentSurfaceDensity*(1-apparentSurfaceDensity))**0.5)    #this is 'standardDeriTab' in the original code
    #print('\napparent suface density, expected mass, standard deviation:\n', apparentSurfaceDensity, expectedMass, standardDeviation)
    lower = np.maximum(np.zeros(self.__masspoints.size), np.floor(expectedMass-self.__constants.nSigma*standardDeviation))
    upper = np.minimum(self.__Nj, np.ceil(expectedMass+self.__constants.nSigma*standardDeviation))
    self.__masspointNumberRange = np.array([lower, upper]).T
    #print('\nMasspoint number range:\n', self.__masspointNumberRange)
    self.__combinations = self.calculateCombinations()
    self.__probability = np.zeros((len(self.__combinations),2))
    #print('\nCombinations:\n', self.__combinations)
    for i,combination in enumerate(self.__combinations):
      print(combination)
      if self.__flagCombination=='Gauss':
        if np.any(expectedMass>self.__constants.pnGauss) and np.any(self.__Nj>self.__constants.nGauss):
          # use gauss!
          g = Gauss(expectedMass, standardDeviation)
          self.__probability[i] = g.gaussfunc(combination)
          #print('gauss!!...')
        else:
          # use binomial 
          # <<This will likely print an error when there are more masspoints>>
          b = Binomial(self.__Nj, apparentSurfaceDensity) # n and p for binominal 
          self.__probability[i] = b.binomfunc(combination)
      else:
        if np.any(expectedMass>self.__constants.pnGauss) and np.any(self.__Nj>self.__constants.nGauss):
          # use gauss
          g = Gauss(expectedMass, standardDeviation)
          self.__probability[i] = g.gaussfunc(combination)
          #print('gauss!!...')
        else:
          # use poisson
          po = Poisson(expectedMass)
          self.__probability[i] = po.poissonfunc(combination)
    #for i,combination in enumerate(self.__combinations): self.__probability[i] = self.__probability[i](combination)
    # print(self.__probability)
    return
  def calculate(self):
    '''Maybe <<PARALLELISE>> this??

       This is a function to cycle through the Combination instances to create a large numpy.ndarray,
       which is used to calculate the final sums needed for the voxel.'''
    self.calculateMasspoints()
    self.calculateRadii()
    self.createCombinationObjects()
    print('Calculating ensemble emission')
    result = []
    #print(self.__combinations)
    for combination in self.__combinations:
      self.__combinationObjects.append(Combination(self.__species, self.__interpolations, combination=combination, masses=self.__masspoints, density=self.__masspointDensity, fuv=self.__FUV, probability=self.__probability))
      self.__combinationObjects[-1].calculateEmission(self.__velocity, self.__velocityDispersion)
      result.append(self.__combinationObjects[-1].getScaledCombinationEmission()) #<<this needs to be altered>>
    result = np.array(result)
    self.__intensity = result#.sum(3)[0]
    self.__opticalDepth = result#.sum(3)[1]
    if not isinstance(self.__FUV, FUVfield): self.__FUV = result#.average(3)[2]
    return
  def getEnsembleEmission(self):
    '''This returns the ensemble emission...nothing more.'''
    return (self.__intensity,self.__opticalDepth,self.__FUV)
  def getCombinations(self):
    return self.__combinations
  def getCombinationObjects(self):
    return self.__combinationObjects