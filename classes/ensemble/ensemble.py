import importlib as il
import pprint
import numpy as np
from numba import jit_module
import sys
import gc

import constants
import ensemble
import combinations
import masspoints

import ensembleStatistics as stat

'''
This class owes itself largely to the work  done by Silke Andree-Labsch and
Christoph Bruckmann. It contains a number of combinations of fractal mass
configurations, and therefore contains the ensemble-averaged intensity and
optical depth.
It can be either a clump or an interclump ensemble.
'''

def initialise(velocity=0, ensembleDispersion=0, clumpMass=0, interclumpMass=0, verbose=False):
  if verbose:
    print('Ensemble instance initialised\n')
  ensemble.clumpMass = clumpMass
  ensemble.interclumpMass = interclumpMass
  createCombinationObjects(velocity, ensembleDispersion)
  return

def calculateCombinations(clumpN, test=True, verbose=False):
  '''
  This function calculates all of the different combinations of clump masses that may be in a line-of-sight.
  It is basically the essence of the probabilistic approach used to create the superposition of clumps.
  '''
  dimension = len(clumpN[0])
  ranges = clumpN
  if verbose:
    print('\nMasspoint Ranges:\n{}\n'.format(ranges))
    input()
  ranges[:,:,1] += 1
  combinations = []
  for i in range(len(ranges)):

    if test:
      # grid = np.meshgrid(*[np.arange(ranges[i,j,0],ranges[i,j,1]) for j in range(dimension)])
      grid = np.meshgrid(*[np.arange(0,ranges[i,j,1]) for j in range(dimension)])
      combinations.append(np.array([grid[j].flatten() for j in range(len(grid))]))
    
    else:
      if dimension==1:
        grid = np.arange(ranges[i,0,0], ranges[i,0,1])
        combinations.append(np.array([grid.flatten()], dtype=np.int))
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
      else:
        sys.exit('\nThere are too many masses for the current grid ({}).\nExitting. . .\n\n'.format(dimension))
  
  #combinations = np.array(combinations)
  if verbose:
    print('\nCalculated combinations:\n', combinations)
  return combinations

def createCombinationObjects(velocity, ensembleDispersion, verbose=False, debug=False):
  '''
  This function removes all of the unnecessary degenerate looping during this calculation.
  Of course it is possible because of the wonders of numpy.ndarray(). . .
  '''
  #verbose = self.__verbose or verbose
  #if verbose: print(self.__clumpType)
  ensemble.clumpNj = (ensemble.clumpMass*10.**(constants.clumpLogMass*(1-constants.alpha))) / (10.**(constants.clumpLogMass*(2-constants.alpha))).sum()
  ensemble.interclumpNj = (ensemble.interclumpMass*10.**(constants.interclumpLogMass*(1-constants.alpha))) / (10.**(constants.interclumpLogMass*(2-constants.alpha))).sum()
  
  if verbose:
    print('\nClump Nj:\n', ensemble.clumpNj)
    print('\nInterclump Nj:\n', ensemble.interclumpNj)

  clumpMassEnsemble = (ensemble.clumpNj*10.**constants.clumpLogMass).sum()
  clumpRadiusEnsemble = (ensemble.clumpNj*masspoints.clumpRadius).sum()
  clumpVolumeEnsemble = (ensemble.clumpNj*np.pi*4./3.*masspoints.clumpRadius.T**3).sum()

  interclumpMassEnsemble = (ensemble.interclumpNj*10.**constants.interclumpLogMass).sum()
  interclumpRadiusEnsemble = (ensemble.interclumpNj*masspoints.interclumpRadius).sum()
  interclumpVolumeEnsemble = (ensemble.interclumpNj*np.pi*4./3.*masspoints.interclumpRadius.T**3).sum()
  
  if verbose:
    print(clumpMassEnsemble, clumpVolumeEnsemble)
    print(interclumpMassEnsemble, interclumpVolumeEnsemble)
  
  clumpDensityEnsemble = clumpMassEnsemble/clumpVolumeEnsemble
  interclumpDensityEnsemble = interclumpMassEnsemble/interclumpVolumeEnsemble
  
  if verbose:
    print('velocity, mean, dispersion', self.__velocity, self.__velocity.mean(), self.__velocityDispersion)
  
  ensembleDispersion = np.sqrt(constants.ensembleDispersion**2+ensembleDispersion**2)
  ensemble.clumpDeltaNji = (np.array(ensemble.clumpNj).T/np.sqrt(2*np.pi)/ensembleDispersion*\
                           (np.exp(-0.5*((constants.velocityRange-velocity)/ensembleDispersion)**2)).T*\
                           constants.velocityStep)
  ensemble.interclumpDeltaNji = (np.array(ensemble.interclumpNj).T/np.sqrt(2*np.pi)/ensembleDispersion*\
                                (np.exp(-0.5*((constants.velocityRange-velocity)/ensembleDispersion)**2)).T*\
                                constants.velocityStep)
  
  #self.__deltaNji = np.array([self.__deltaNji]).T

  # print('\n\nC L U M P S\n\n')
  resize = constants.clumpNmax/ensemble.clumpDeltaNji[-1,:]   #scaling factor to set the maximum number of the largest clump
  i_finite = ~(np.isnan(resize)|np.isinf(resize))
  resize =  resize[i_finite]
  # print('resize:\n', resize)
  # print('number:\n', ensemble.clumpDeltaNji)
  # print('size:\n', constants.resolution)
  ensemble.clumpNormalisedDeltaNji = ensemble.clumpDeltaNji[:,i_finite]*resize
  # print('normalised number:\n', ensemble.clumpNormalisedDeltaNji)
  # print('normalised size:\n', constants.resolution*np.sqrt(resize))
  clumpSurfaceProbability = np.array(np.pi*masspoints.clumpRadius.T**2/constants.resolution**2/resize)    #this is 'pTab' in the original code
  if (clumpSurfaceProbability>=1).any():    #increase voxel size if the clumps are too large
    supersize = np.ceil(np.pi*masspoints.clumpRadius.max()**2/constants.resolution**2/resize)
    # print('supersize:\n', supersize)
    ensemble.clumpNormalisedDeltaNji *= supersize
    # print('renormalised number:\n', ensemble.clumpNormalisedDeltaNji)
    # print('renormalised size:\n', constants.resolution*np.sqrt(resize*supersize))
    clumpSurfaceProbability = np.array(np.pi*masspoints.clumpRadius.T**2/constants.resolution**2/supersize/resize)
  clumpProbableNumber = (ensemble.clumpNormalisedDeltaNji*clumpSurfaceProbability)   #this is 'expectedValTab' in the original code
  clumpStandardDeviation = np.sqrt(ensemble.clumpNormalisedDeltaNji*clumpSurfaceProbability*(1-clumpSurfaceProbability))    #this is 'standardDeriTab' in the original code
  
  # print('\nClump surface P, Pn, Sn\n', clumpSurfaceProbability, '\n\n', clumpProbableNumber, '\n\n', clumpStandardDeviation, '\n')

  resizeMax = constants.clumpNmax/ensemble.clumpNj[0,-1]   #scaling factor to set the maximum number of the largest clump
  ensemble.clumpNormalisedNj = ensemble.clumpNj*resizeMax
  CLmaxSurfaceProbability = np.array(np.pi*masspoints.clumpRadius.T**2/constants.resolution**2/resizeMax)
  if (CLmaxSurfaceProbability>=1).any():
    supersize = np.ceil(np.pi*masspoints.clumpRadius.max()**2/constants.resolution**2/resizeMax)
    ensemble.clumpNormalisedNj *= supersize
    CLmaxSurfaceProbability = np.array(np.pi*masspoints.clumpRadius.T**2/constants.resolution**2/resizeMax/supersize)
  CLmaxProbableNumber = (ensemble.clumpNormalisedNj*CLmaxSurfaceProbability.flatten())
  CLmaxStandardDeviation = np.sqrt(ensemble.clumpNormalisedNj*(CLmaxSurfaceProbability*(1-CLmaxSurfaceProbability)).flatten())

  # print('\nMax clump surface P, Pn, Sn\n', CLmaxSurfaceProbability, '\n\n', CLmaxProbableNumber, '\n\n', CLmaxStandardDeviation, '\n')

  # print('\n\nI N T E R C L U M P S\n\n')
  resize = constants.interclumpNmax/ensemble.interclumpDeltaNji[-1,:]   #scaling factor to set the maximum number of the largestclump
  i_finite = ~(np.isnan(resize)|np.isinf(resize))
  resize =  resize[i_finite]
  # print('resize:\n', resize)
  # print('number:\n', ensemble.interclumpDeltaNji)
  # print('size:\n', constants.resolution)
  ensemble.interclumpNormalisedDeltaNji = ensemble.interclumpDeltaNji[:,i_finite]*resize
  # print('normalised number:\n', ensemble.interclumpDeltaNji)
  # print('normalised size:\n', constants.resolution*np.sqrt(resize))
  interclumpSurfaceProbability = np.array(np.pi*masspoints.interclumpRadius.T**2/constants.resolution**2/resize)    #this is 'pTab' in the original code
  if (interclumpSurfaceProbability>=1).any():    #increase voxel size if the interclumps are too large
    supersize = np.ceil(np.pi*masspoints.interclumpRadius.max()**2/constants.resolution**2/resize)
    # print('supersize:\n', supersize)
    ensemble.interclumpNormalisedDeltaNji *= supersize
    # print('normalised number:\n', ensemble.interclumpNormalisedDeltaNji)
    # print('normalised size:\n', constants.resolution*np.sqrt(resize*supersize))
    interclumpSurfaceProbability = np.array(np.pi*masspoints.interclumpRadius.T**2/constants.resolution**2/resize/supersize)
  interclumpProbableNumber = (ensemble.interclumpNormalisedDeltaNji*interclumpSurfaceProbability)   #this is 'expectedValTab' in the original code
  interclumpStandardDeviation = np.sqrt(ensemble.interclumpNormalisedDeltaNji*interclumpSurfaceProbability*(1-interclumpSurfaceProbability))    #this is 'standardDeriTab' in the original code
  
  # print('\nInterclump surface P, Pn, Sn\n', interclumpSurfaceProbability, '\n\n', interclumpProbableNumber, '\n\n', interclumpStandardDeviation, '\n')

  resizeMax = constants.clumpNmax/ensemble.interclumpNj[0,-1]   #scaling factor to set the maximum number of the largest clump
  ensemble.interclumpNormalisedNj = ensemble.interclumpNj*resizeMax
  ICmaxSurfaceProbability = np.array(np.pi*masspoints.interclumpRadius.T**2/(constants.resolution*np.sqrt(resizeMax))**2)
  if (ICmaxSurfaceProbability>=1).any():
    supersize = np.ceil(np.pi*masspoints.interclumpRadius.max()**2/constants.resolution**2/resizeMax)
    ensemble.interclumpNormalisedNj *= supersize
    ICmaxSurfaceProbability = np.array(np.pi*masspoints.interclumpRadius.T**2/constants.resolution**2/resizeMax/supersize)
  ICmaxProbableNumber = (ensemble.interclumpNormalisedNj*ICmaxSurfaceProbability.flatten())
  ICmaxStandardDeviation = np.sqrt(ensemble.interclumpNormalisedNj*(ICmaxSurfaceProbability*(1-ICmaxSurfaceProbability)).flatten())

  # print('\nMax interclump surface P, Pn, Sn\n', ICmaxSurfaceProbability, '\n\n', ICmaxProbableNumber, '\n\n', ICmaxStandardDeviation, '\n')
  
  if verbose:
    print('\nNj\n', self.__Nj)
    print('\nsurface probability {}, expected number {}, standard deviation {}\n'.format(surfaceProbability, maxProbableNumber, maxStandardDeviation))
    print('\nDelta Nji\n', self.__deltaNji)
    print('\nsurface probability, expected number, standard deviation:\n', surfaceProbability, '\n', probableNumber, '\n', standardDeviation)
    input()

  lower = np.zeros([constants.clumpLogMass.size, 1])
  clumpLower = np.maximum(lower, np.floor(clumpProbableNumber-constants.nSigma*clumpStandardDeviation))
  clumpUpper = np.minimum(ensemble.clumpNormalisedDeltaNji, np.ceil(clumpProbableNumber+constants.nSigma*clumpStandardDeviation))
  CLmaxLower = np.maximum(np.zeros([constants.clumpLogMass.size, 1]), np.floor(CLmaxProbableNumber-constants.nSigma*CLmaxStandardDeviation))
  CLmaxUpper = np.minimum(ensemble.clumpNormalisedNj, np.ceil(CLmaxProbableNumber+constants.nSigma*CLmaxStandardDeviation))
  if verbose:
    print('\nupper,lower:\n',clumpUpper,'\n',clumpLower)
  
  lower = np.zeros([constants.interclumpLogMass.size, 1])
  interclumpLower = np.maximum(lower, np.floor(interclumpProbableNumber-constants.nSigma*interclumpStandardDeviation))
  interclumpUpper = np.minimum(ensemble.interclumpNormalisedDeltaNji, np.ceil(interclumpProbableNumber+constants.nSigma*interclumpStandardDeviation))
  ICmaxLower = np.maximum(np.zeros([constants.interclumpLogMass.size, 1]), np.floor(ICmaxProbableNumber-constants.nSigma*ICmaxStandardDeviation))
  ICmaxUpper = np.minimum(ensemble.interclumpNormalisedNj, np.ceil(ICmaxProbableNumber+constants.nSigma*ICmaxStandardDeviation))
  if verbose:
    print('\nupper,lower:\n',interclumpUpper,'\n',interclumpLower)
  
  ensemble.clumpNumberRange = np.array([clumpLower, clumpUpper], dtype=np.int).T
  ensemble.interclumpNumberRange = np.array([interclumpLower, interclumpUpper], dtype=np.int).T
  
  if verbose:
    print('\nMasspoint number range:\n', ensemble.clumpNumberRange.round())
    print('\nMasspoint number range:\n', ensemble.interclumpNumberRange.round())
  
  ensemble.clumpCombinations = calculateCombinations(ensemble.clumpNumberRange)
  ensemble.interclumpCombinations = calculateCombinations(ensemble.interclumpNumberRange)
  
  #########################################################################################

  if verbose:
    input(clumpCombinations)
    input(interclumpCombinations)
  #self.__deltaNji.resize(1,2,51)
  #probableNumber.resize(1,2,51)
  #standardDeviation.resize(1,2,51)
  #input(self.__velocity)
  #input(self.__velocityBins)
  #self.__probability = np.zeros((len(self.__combinations),2))
  largestCombination = ((ensemble.clumpNumberRange[:,:,1]-0).prod(1))
  # largestCombination = ((ensemble.clumpNumberRange[:,:,1]-ensemble.clumpNumberRange[:,:,0]).prod(1))
  ensemble.clumpLargestIndex = largestCombination.argmax()#np.where((ensemble.clumpNumberRange[:,:,1]-ensemble.clumpNumberRange[:,:,0]).prod(1)==largestCombination)[0][0]
  ensemble.clumpLargestCombination = int(largestCombination.max())
  largestCombination = ((ensemble.interclumpNumberRange[:,:,1]-0).prod(1))
  # largestCombination = ((ensemble.interclumpNumberRange[:,:,1]-ensemble.interclumpNumberRange[:,:,0]).prod(1))
  ensemble.interclumpLargestIndex = largestCombination.argmax()#np.where((ensemble.interclumpNumberRange[:,:,1]-ensemble.interclumpNumberRange[:,:,0]).prod(1)==largestCombination)[0][0]
  ensemble.interclumpLargestCombination = int(largestCombination.max())
  
  if verbose:
    print('\nCombinations:\n{}'.format(ensemble.clumpCombinations))
    input('\nLargest number of combinations: {}\n'.format(ensemble.clumpCombinations[ensemble.clumpLargestCombination]))
    print('\nCombinations:\n{}'.format(ensemble.interclumpCombinations))
    input('\nLargest number of combinations: {}\n'.format(ensemble.interclumpCombinations[ensemble.interclumpLargestCombination]))
  
  # Probability (clean-up soon)
  # Clump
  probabilityList = []
  combinationIndeces = []
  clumpLargestCombination = ensemble.clumpCombinations[ensemble.clumpLargestIndex].T
  # print(clumpLargestCombination)
  for i,combinations in enumerate(ensemble.clumpCombinations):   #loop over combinations of masspoints in each velocity bin
    ensemble.clumpCombinations[i] = np.array(combinations).T
    probability = np.zeros((ensemble.clumpLargestCombination, constants.clumpLogMass.size))
    maxProbability = np.zeros((ensemble.clumpLargestCombination, constants.clumpLogMass.size))
    combinationIndeces.append(i)
    if verbose:
      print('\nEnsemble combinations:\n', combinations)
    if combinations.any() and ~(np.isinf(combinations)|np.isnan(combinations)).any():    #calculate the probability if there are any masspoints in this velocity bin
      for combination in ensemble.clumpCombinations[i]:
        combination = np.array([combination])
        index = np.where((clumpLargestCombination==combination).all(1))[0][0]
        # print(index)
        if verbose:
          print('\nCombination:\n', combination)
          input()
        if constants.probability=='binomial':
          if verbose:
            print(clumpProbableNumber.shape, combination.shape)
          iGauss = (clumpProbableNumber[:,i]>constants.pnGauss) & (ensemble.clumpNormalisedDeltaNji[:,i]>constants.nGauss)
          # print(iGauss.shape, iGauss)
          if iGauss.any():
            # use gauss!
            if verbose:
              print('Gauss')
            g = stat.Gauss(clumpProbableNumber[:,i], clumpStandardDeviation[:,i], debug=debug)
            probability[index,iGauss] = (g.gaussfunc(combination))[iGauss]
          if (~iGauss).any():
            # use binomial
            if verbose:
              print('Binomial')
            # <<This will likely print an error when there are more masspoints>>
            b = stat.Binomial(ensemble.clumpNormalisedDeltaNji[:,i], clumpSurfaceProbability[:,i], debug=debug) # n and p for binominal 
            probability[index,~iGauss] = (b.binomfunc(combination))[~iGauss]
          iGauss = ((CLmaxProbableNumber>constants.pnGauss) & (ensemble.clumpNormalisedNj>constants.nGauss))[0]
          # print(iGauss.shape, iGauss)
          if iGauss.any:
            # print(combination)
            g = stat.Gauss(CLmaxProbableNumber, CLmaxStandardDeviation, debug=debug)
            maxProbability[index,iGauss] = (g.gaussfunc(combination))[iGauss]
          if (~iGauss).any():
            b = stat.Binomial(ensemble.clumpNormalisedNj, CLmaxSurfaceProbability, debug=debug) # n and p for binominal 
            maxProbability[index,~iGauss] = (b.binomfunc(combination))[~iGauss]

        elif constants.probability=='poisson':
          if np.any(clumpProbableNumber[:,i]>constants.pnGauss) and np.any(ensemble.clumpNormalisedDeltaNji[:,i]>constants.nGauss):
            # use gauss
            if verbose:
              print('Gauss')
            g = stat.Gauss(clumpProbableNumber[:,i], clumpStandardDeviation[:,i])
            probability.append(g.gaussfunc(combination))
          else:
            # use poisson
            if verbose:
              print('Poisson')
            po = stat.Poisson(clumpProbableNumber[:,i])
            probability.append(po.poissonfunc(combination))
          if np.any(CLmaxProbableNumber>constants.pnGauss) and np.any(ensemble.clumpNormalisedNj>constants.nGauss):
            g = stat.Gauss(CLmaxProbableNumber, CLmaxStandardDeviation, debug=debug)
            maxProbability[index,:] = (g.gaussfunc(combination))
          else:
            po = stat.Poisson(CLmaxProbableNumber) # n and p for binominal 
            maxProbability[index,:] = (po.binomfunc(combination))

    if np.shape(probability)!=np.shape(maxProbability) and debug:
      for i in range(len(probability)):
        print(np.array(probability[i]).shape)
        print(np.array(maxProbability[i]).shape)
      input()
    #print((probability))
    while(len(probability) < ensemble.clumpLargestCombination):
      probability.append(np.zeros((constants.clumpLogMass.size)))
      maxProbability.append(np.zeros((constants.clumpLogMass.size)))
      if debug:
        print('Probability length:', len(probability), ', last element shape:', probability[-1].shape)
        input()
    if (np.array(probability[-1])==np.nan).any():
      print('\nThere is an invalid probability:', probability[-1], '\n')
      input()
    #print(probabilityList, probability)
    probabilityList.append(np.array(probability))
    if i==ensemble.clumpLargestIndex:
      maxProbabilityList = np.array(maxProbability)
  clumpCombinationIndeces = np.array(combinationIndeces)
  ensemble.clumpProbability = np.array(probabilityList)
  ensemble.CLmaxProbability = np.array(maxProbabilityList)

  # Interclump
  probabilityList = []
  combinationIndeces = []
  interclumpLargestCombination = ensemble.interclumpCombinations[ensemble.interclumpLargestIndex].T
  # print(interclumpLargestCombination)
  for i,combinations in enumerate(ensemble.interclumpCombinations):   #loop over combinations of masspoints in each velocity bin
    ensemble.interclumpCombinations[i] = np.array(combinations).T
    probability = np.zeros((ensemble.interclumpLargestCombination, constants.interclumpLogMass.size))
    maxProbability = np.zeros((ensemble.interclumpLargestCombination, constants.interclumpLogMass.size))
    combinationIndeces.append(i)
    if verbose:
      print(combinations)
      print('\nEnsemble combinations:\n', combinations)
    if combinations.any():    #calculate the probability if there are any masspoints in this velocity bin
      for combination in ensemble.interclumpCombinations[i]:
        combination = np.array([combination])
        index = np.where((interclumpLargestCombination==combination).all(1))[0][0]
        # print(combination)
        if verbose:
          print('\nCombination:\n', combination)
          input()
        if constants.probability=='binomial':
          if verbose:
            print(interclumpProbableNumber.shape, combination.shape)
          iGauss = (interclumpProbableNumber[:,i]>constants.pnGauss) & (ensemble.interclumpNormalisedDeltaNji[:,i]>constants.nGauss)
          if iGauss.any():
            # use gauss!
            if verbose:
              print('Gauss')
            # print(interclumpSurfaceProbability, i)
            g = stat.Gauss(interclumpProbableNumber[:,i], interclumpStandardDeviation[:,i], debug=debug, verbose=False)
            probability[index,iGauss] = (g.gaussfunc(combination, verbose=False))[iGauss]
            #print('gauss!!...')
          if (~iGauss).any():
            # use binomial
            if verbose:
              print('Binomial')
            # <<This will likely print an error when there are more masspoints>>
            # print(interclumpProbableNumber[:,i])
            # print(ensemble.interclumpDeltaNji[:,i], interclumpSurfaceProbability)
            # print(ensemble.interclumpNj, interclumpSurfaceProbability)
            b = stat.Binomial(ensemble.interclumpNormalisedDeltaNji[:,i], interclumpSurfaceProbability[:,i], debug=debug) # n and p for binominal 
            probability[index,~iGauss] = (b.binomfunc(combination))[~iGauss]
            # print(combination, b.binomfunc(combination))
            # input()
          iGauss = ((ICmaxProbableNumber>constants.pnGauss) & (ensemble.interclumpNormalisedNj>constants.nGauss))[0]
          if iGauss.any():
            g = stat.Gauss(ICmaxProbableNumber, ICmaxStandardDeviation, debug=debug, verbose=False)
            maxProbability[index,iGauss] = (g.gaussfunc(combination, verbose=False))[iGauss]
          if (~iGauss).any():
            b = stat.Binomial(ensemble.interclumpNormalisedNj, ICmaxSurfaceProbability, debug=debug) # n and p for binominal 
            maxProbability[index,~iGauss] = (b.binomfunc(combination))[~iGauss]

        elif constants.probability=='poisson':
          if np.any(interclumpProbableNumber[:,i]>constants.pnGauss) and np.any(ensemble.interclumpNormalisedDeltaNji[:,i]>constants.nGauss):
            # use gauss
            if verbose:
              print('Gauss')
            g = stat.Gauss(interclumpProbableNumber[:,i], interclumpStandardDeviation[:,i])
            probability.append(g.gaussfunc(combination))
          else:
            # use poisson
            if verbose:
              print('Poisson')
            po = stat.Poisson(interclumpProbableNumber[:,i])
            probability.append(po.poissonfunc(combination))
          if np.any(ICmaxProbableNumber>constants.pnGauss) and np.any(ensemble.interclumpNormalisedNj>constants.nGauss):
            g = stat.Gauss(ICmaxProbableNumber, ICmaxStandardDeviation, debug=debug)
            maxProbability[index,:] = (g.gaussfunc(combination))
          else:
            po = stat.Poisson(ICmaxProbableNumber) # n and p for binominal 
            maxProbability[index,:] = (po.binomfunc(combination))
    if np.shape(probability)!=np.shape(maxProbability) and debug:
      for i in range(len(probability)):
        print(np.array(probability[i]).shape)
        print(np.array(maxProbability[i]).shape)
      input()
    while(len(probability) < ensemble.interclumpLargestCombination):
      probability.append(np.zeros((constants.interclumpLogMass.size)))
      maxProbability.append(np.zeros((constants.interclumpLogMass.size)))
      if debug:
        print('Probability length:', len(probability), ', last element shape:', probability[-1].shape)
        input()
    if (np.array(probability[-1])==np.nan).any():
      print('\nThere is an invalid probability:', probability[-1], '\n')
      input()
    probabilityList.append(np.array(probability))
    if i==ensemble.interclumpLargestIndex:
      maxProbabilityList = np.array(maxProbability)
  interclumpCombinationIndeces = np.array(combinationIndeces)
  ensemble.interclumpProbability = np.array(probabilityList)
  ensemble.ICmaxProbability = np.array(maxProbabilityList)

  #print(self.__probability.prod(2)[:,2].shape, '\n', self.__maxProbability.max(0))
  #input()
  if debug:
    print('Clump')
    for i in range(len(ensemble.clumpProbability)): input(np.array(ensemble.clumpProbability[i]).shape)
  if verbose:
    print('\nProbability ({}):\n{}\n'.format(np.shape(ensemble.clumpProbability), ensemble.clumpProbability))
    if debug:
      for i in range(len(ensemble.clumpProbability)): print('Probability shapes:\n{}\n'.format(np.array(ensemble.clumpProbability[i].shape)))
      for i in clumpCombinationIndeces: print('Combination sizes:\n{}\n'.format(np.array(ensemble.clumpCombinations[i].size)))
    input()
  if debug:
    print('Interclump')
    for i in range(len(ensemble.interclumpProbability)): input(np.array(ensemble.interclumpProbability[i]).shape)
  if verbose:
    print('\nProbability ({}):\n{}\n'.format(np.shape(ensemble.clumpProbability), ensemble.clumpProbability))
    if debug:
      for i in range(len(ensemble.interclumpProbability)): print('Probability shapes:\n{}\n'.format(np.array(ensemble.interclumpProbability[i].shape)))
      for i in interclumpCombinationIndeces: print('Combination sizes:\n{}\n'.format(np.array(ensemble.interclumpCombinations[i].size)))
    input()

  ensemble.clumpIndeces = ensemble.clumpProbability.prod(2).sum(1).nonzero()[0]
  if ensemble.clumpIndeces.size>constants.clumpMaxIndeces: constants.clumpMaxIndeces = ensemble.clumpIndeces.size
  ensemble.interclumpIndeces = ensemble.interclumpProbability.prod(2).sum(1).nonzero()[0]
  if ensemble.interclumpIndeces.size>constants.interclumpMaxIndeces: constants.interclumpMaxIndeces = ensemble.interclumpIndeces.size

  ensemble.clumpCombinations = np.array(ensemble.clumpCombinations[ensemble.clumpIndeces[0]:(ensemble.clumpIndeces[-1]+1)])
  ensemble.interclumpCombinations = np.array(ensemble.interclumpCombinations[ensemble.interclumpIndeces[0]:(ensemble.interclumpIndeces[-1]+1)])

  ensemble.clumpProbability = ensemble.clumpProbability[ensemble.clumpIndeces,:,:]
  ensemble.interclumpProbability = ensemble.interclumpProbability[ensemble.interclumpIndeces,:,:]

  ensemble.clumpLargestIndex = np.where(ensemble.clumpIndeces==ensemble.clumpLargestIndex)[0][0]
  ensemble.interclumpLargestIndex = np.where(ensemble.interclumpIndeces==ensemble.interclumpLargestIndex)[0][0]

  return

#@jit(nopython=False)

def calculate(Afuv, debug=False, test=False):
  '''

  !!!!  <<DEGENERATE>>

  Maybe <<PARALLELISE>> this??

  This is a function to cycle through the Combination instances to create a large numpy.ndarray,
  which is used to calculate the final sums needed for the voxel.
  '''
  intensityResult = []
  opticalDepthResult = []
  
  for combination in ensemble.clumpCombinationObjects:
    if (combination in ensemble.clumpCombinationObjects[-10:]) and test:
      result = combination.calculateEmission(Afuv=Afuv, test=True)
    else:
      result = combination.calculateEmission(Afuv=Afuv)

    intensityResult.append(result[0])
    opticalDepthResult.append(result[1])

  intensityResult = np.array(intensityResult)
  opticalDepthResult = np.array(opticalDepthResult)

  if debug:
    print('\nIntensity\n{}\nOptical depth\n{}\n'.format(intensityResult, opticalDepthResult))
  if debug:
    print(intensityResult)
    print(opticalDepthResult)
    input()

  intensity = intensityResult.sum(0)
  opticalDepth = -np.log((opticalDepthResult.sum(0)).astype(np.float))
  gc.collect()
  if debug:
    print('\nIntensity\n{}\nOptical depth\n{}\n'.format(self.__intensity, self.__opticalDepth))
    input()

  return (intensity,opticalDepth)

def printEnsembleParameters():
  
  np.set_printoptions(precision=4, suppress=True)

  print('\nC L U M P\n')
  print('Nj:\n{}'.format(ensemble.clumpNj.astype(np.float)))
  print('delta Nji:\n{}\n'.format(ensemble.clumpDeltaNji.astype(np.float).T))
  print('Normalised Nj:\n{}'.format(ensemble.clumpNormalisedNj.astype(np.float)))
  print('Normalised delta Nj:\n{}\n'.format(ensemble.clumpNormalisedDeltaNji.astype(np.float).T))
  
  print('\nI N T E R C L U M P\n')
  print('Nj:\n{}'.format(ensemble.interclumpNj.astype(np.float)))
  print('delta Nji:\n{}\n'.format(ensemble.interclumpDeltaNji.astype(np.float).T))
  print('Normalised Nj:\n{}'.format(ensemble.interclumpNormalisedNj.astype(np.float)))
  print('Normalised Delta Nji:\n{}\n'.format(ensemble.interclumpNormalisedDeltaNji.astype(np.float).T))
  
  return

# jit_module(nopython=False)
