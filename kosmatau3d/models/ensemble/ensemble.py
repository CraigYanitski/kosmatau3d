import sys
import gc
import importlib as il
import pprint
import warnings

import numpy as np
from scipy import stats
from numba import jit_module

from kosmatau3d.models import constants
from kosmatau3d.models import ensemble
from kosmatau3d.models import combinations
from kosmatau3d.models import masspoints
# from kosmatau3d.models import ensembleStatistics as stat

'''
This class owes itself largely to the work  done by Silke Andree-Labsch and
Christoph Bruckmann. It contains a number of combinations of fractal mass
configurations, and therefore contains the ensemble-averaged intensity and
optical depth.
It can be either a clump or an interclump ensemble.
'''


def initialise(ensembledispersion=0, ensemblemass=0, suggested_calc=True, verbose=False):
    if verbose:
        print('Ensemble instance initialised\n')
    ensemble.clumpMass = ensemblemass
    createCombinationObjects(ensembledispersion, dtype=constants.dtype, suggested_calc=suggested_calc)
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
    # ranges[:,:,1] += 1
    combinations = []
    for i in range(len(ranges)):
  
        if test:
            # grid = np.meshgrid(*[np.arange(ranges[i,j,0],ranges[i,j,1]) for j in range(dimension)])
            grid = np.meshgrid(*[np.arange(0, ranges[i, j, 1]) for j in range(dimension)])
            combinations.append(np.array([grid[j].flatten() for j in range(len(grid))]))
        
        else:
            if dimension == 1:
                grid = np.arange(ranges[i, 0, 0], ranges[i, 0, 1])
                combinations.append(np.array([grid.flatten()], dtype=np.int))
            elif dimension == 2:
                grid = np.mgrid[ranges[i, 0, 0]:ranges[i, 0, 1], ranges[i, 1, 0]:ranges[i, 1, 1]]
                combinations.append(np.array([grid[0].flatten(), grid[1].flatten()], dtype=np.int))
            elif dimension == 3:
                grid = np.mgrid[ranges[i, 0, 0]:ranges[i, 0, 1], ranges[i, 1, 0]:ranges[i, 1, 1],
                                ranges[i, 2, 0]:ranges[i, 2, 1]]
                combinations.append(np.array([grid[0].flatten(), grid[1].flatten(), grid[2].flatten()], dtype=np.int))
            elif dimension == 4:
                grid = np.mgrid[ranges[i, 0, 0]:ranges[i, 0, 1], ranges[i, 1, 0]:ranges[i, 1, 1],
                                ranges[i, 2, 0]:ranges[i, 2, 1], ranges[i, 3, 0]:ranges[i, 3, 1]]
                combinations.append(np.array([grid[0].flatten(), grid[1].flatten(), grid[2].flatten(),
                                              grid[3].flatten()],
                                             dtype=np.int))
            elif dimension == 5:
                grid = np.mgrid[ranges[i, 0, 0]:ranges[i, 0, 1], ranges[i, 1, 0]:ranges[i, 1, 1],
                                ranges[i, 2, 0]:ranges[i, 2, 1], ranges[i, 3, 0]:ranges[i, 3, 1],
                                ranges[i, 4, 0]:ranges[i, 4, 1]]
                combinations.append(np.array([grid[0].flatten(), grid[1].flatten(), grid[2].flatten(),
                                              grid[3].flatten(), grid[4].flatten()],
                                             dtype=np.int))
            elif dimension == 6:
                grid = np.mgrid[ranges[i, 0, 0]:ranges[i, 0, 1], ranges[i, 1, 0]:ranges[i, 1, 1],
                                ranges[i, 2, 0]:ranges[i, 2, 1], ranges[i, 3, 0]:ranges[i, 3, 1],
                                ranges[i, 4, 0]:ranges[i, 4, 1], ranges[i, 5, 0]:ranges[i, 5, 1]]
                combinations.append(np.array([grid[0].flatten(), grid[1].flatten(), grid[2].flatten(),
                                              grid[3].flatten(), grid[4].flatten(), grid[5].flatten()],
                                             dtype=np.int))
            elif dimension == 7:
                grid = np.mgrid[ranges[i, 0, 0]:ranges[i, 0, 1], ranges[i, 1, 0]:ranges[i, 1, 1],
                                ranges[i, 2, 0]:ranges[i, 2, 1], ranges[i, 3, 0]:ranges[i, 3, 1],
                                ranges[i, 4, 0]:ranges[i, 4, 1], ranges[i, 5, 0]:ranges[i, 5, 1],
                                ranges[i, 6, 0]:ranges[i, 6, 1]]
                combinations.append(np.array([grid[0].flatten(), grid[1].flatten(), grid[2].flatten(),
                                              grid[3].flatten(), grid[4].flatten(), grid[5].flatten(),
                                              grid[6].flatten()],
                                             dtype=np.int))
            else:
                sys.exit('\nThere are too many masses for the current grid ({}).\nExitting. . .\n\n'.format(dimension))
    
    # combinations = np.array(combinations)
    if verbose:
        print('\nCalculated combinations:\n', combinations)
    return combinations


def createCombinationObjects(ensembleDispersion, dtype=np.float64, suggested_calc=True, verbose=False, debug=False):
    '''
    This function removes all of the unnecessary degenerate looping during this calculation.
    Of course it is possible because of the wonders of numpy.ndarray(). . .
    '''
    # verbose = self.__verbose or verbose
    # if verbose: print(self.__clumpType)
    
    for ens in range(len(constants.clumpMassNumber)):
  
        # This can be used to add more dispersion if needed. It is left over from Silke's Orion Bar application.
        # It has moved to the VoxelGrid setup.
        # ensembleDispersion[ens] = np.sqrt(constants.ensembleDispersion**2+ensembleDispersion[ens]**2)
        
        ensemble.clumpNj[ens] = (ensemble.clumpMass[ens]*10.**(constants.clumpLogMass[ens]*(1-constants.alpha))) / \
                                (10.**(constants.clumpLogMass[ens]*(2-constants.alpha))).sum()
    
        if verbose:
            print('\nClump Nj:\n', ensemble.clumpNj[ens])
            print('\nInterclump Nj:\n', ensemble.interclumpNj)
    
        clumpMassEnsemble = (ensemble.clumpNj[ens]*10.**constants.clumpLogMass[ens]).sum()
        clumpRadiusEnsemble = (ensemble.clumpNj[ens]*masspoints.clumpRadius[ens]).sum()
        clumpVolumeEnsemble = (ensemble.clumpNj[ens]*np.pi*4./3.*masspoints.clumpRadius[ens].T**3).sum()
        
        if verbose:
            print(clumpMassEnsemble, clumpVolumeEnsemble)
            print(interclumpMassEnsemble, interclumpVolumeEnsemble)
        
        clumpDensityEnsemble = clumpMassEnsemble/clumpVolumeEnsemble
        
        if verbose:
            print('velocity, mean, dispersion', self.__velocity, self.__velocity.mean(), self.__velocityDispersion)

        system_dispersion = np.sqrt(np.abs(ensembleDispersion[ens]**2 - constants.clumpDispersion**2))
        if suggested_calc:
            dispersion = np.maximum(ensembleDispersion[ens], constants.clumpDispersion)
            velocityStep = dispersion/constants.velocity_resolution
            ensemble.clumpVelocities[ens] = np.linspace(-constants.n_sigma*dispersion,
                                                        constants.n_sigma*dispersion,
                                                        num=np.round(2*constants.n_sigma*dispersion /
                                                                     velocityStep).astype(np.int)+1)
            velocityStep = ensemble.clumpVelocities[ens][1]-ensemble.clumpVelocities[ens][0]
            ensemble.clumpDeltaNji[ens] = (np.array(ensemble.clumpNj[ens]).T
                                           * constants.clumpDispersion/dispersion *
                                           (np.exp(-0.5*(ensemble.clumpVelocities[ens]/dispersion)**2)).T)
        else:
            velocityStep = np.minimum(system_dispersion/constants.velocity_resolution,
                                      constants.clumpDispersion/constants.velocity_resolution)
            ensemble.clumpVelocities[ens] = np.linspace(-constants.n_sigma*system_dispersion,
                                                        constants.n_sigma*system_dispersion,
                                                        num=np.round(2*constants.n_sigma*system_dispersion /
                                                                     velocityStep).astype(np.int)+1)
            velocityStep = ensemble.clumpVelocities[ens][1]-ensemble.clumpVelocities[ens][0]
            ensemble.clumpDeltaNji[ens] = (np.array(ensemble.clumpNj[ens]).T/np.sqrt(2*np.pi)/system_dispersion *
                                           (np.exp(-0.5*(ensemble.clumpVelocities[ens]/system_dispersion)**2)).T *
                                           velocityStep)
        
        warnings.filterwarnings('ignore', category=RuntimeWarning)
    
        # print('\n\nC L U M P S\n\n')
        pj = np.pi*masspoints.clumpRadius[ens].T**2 / constants.voxel_size**2

        if constants.clumpNmax[ens]:
            # scaling factor to set the maximum number of the largest clump
            normalise = constants.clumpNmax[ens]/ensemble.clumpDeltaNji[ens][-1, :]
        else:
            normalise = np.clip(np.around(ensemble.clumpDeltaNji[ens][-1, :]), 1, None)\
                        / ensemble.clumpDeltaNji[ens][-1, :]
        ensemble.clumpNormalisedDeltaNji[ens] = np.around(ensemble.clumpDeltaNji[ens]*normalise)
        clumpSurfaceProbability = np.array(pj/normalise)
        i_change = (clumpSurfaceProbability >= 1).any(0)
        if i_change.any():  # increase voxel size if the clumps are too large
            n_resized = np.ceil(ensemble.clumpDeltaNji[ens][-1, i_change] * pj.max())
            normalise[i_change] = n_resized / ensemble.clumpDeltaNji[ens][-1, i_change]
            ensemble.clumpNormalisedDeltaNji[ens][:, i_change] = np.around(ensemble.clumpDeltaNji[ens][:, i_change]
                                                                           *normalise[i_change])
            clumpSurfaceProbability = np.asarray(pj / normalise)
        clumpProbableNumber = (ensemble.clumpNormalisedDeltaNji[ens]*clumpSurfaceProbability)
        clumpStandardDeviation = np.sqrt(ensemble.clumpNormalisedDeltaNji[ens]*clumpSurfaceProbability *
                                         (1-clumpSurfaceProbability))
        
        ensemble.clumpSurfaceProbability[ens] = clumpSurfaceProbability
        ensemble.clumpProbableNumber[ens] = clumpProbableNumber
        ensemble.clumpStandardDeviation[ens] = clumpStandardDeviation
    
        if constants.clumpNmax[ens]:
            # scaling factor to set the maximum number of the largest clump
            normaliseMax = constants.clumpNmax[ens]/ensemble.clumpNj[ens][0,-1]
        else:
            normaliseMax = np.clip(ensemble.clumpNj[ens][0, -1], 1, None) / ensemble.clumpNj[ens][0, -1]
        ensemble.clumpNormalisedNj[ens] = np.around(ensemble.clumpNj[ens]*normaliseMax)
        CLmaxSurfaceProbability = np.array(pj / normaliseMax)
        i_change = (CLmaxSurfaceProbability >= 1).any(0)
        if i_change.any():
            n_resized = np.ceil(ensemble.clumpNj[ens][0, -1] * pj.max())
            normaliseMax = n_resized / ensemble.clumpNj[ens][0, -1]
            ensemble.clumpNormalisedNj[ens] = np.around(ensemble.clumpNj[ens]*normaliseMax)
            CLmaxSurfaceProbability = np.asarray(pj / normaliseMax)
        CLmaxProbableNumber = (ensemble.clumpNormalisedNj[ens].T*CLmaxSurfaceProbability)
        CLmaxStandardDeviation = np.sqrt(ensemble.clumpNormalisedNj[ens].T*(CLmaxSurfaceProbability *
                                                                          (1-CLmaxSurfaceProbability)))

        ensemble.CLmaxSurfaceProbability[ens] = CLmaxSurfaceProbability
        ensemble.CLmaxProbableNumber[ens] = CLmaxProbableNumber
        ensemble.CLmaxStandardDeviation[ens] = CLmaxStandardDeviation
    
        warnings.filterwarnings('default', category=RuntimeWarning)
        
        if verbose:
            print('\nNj\n', self.__Nj)
            print('\nsurface probability {}, expected number {}, standard deviation {}\n'.format(surfaceProbability,
                                                                                                 maxProbableNumber,
                                                                                                 maxStandardDeviation))
            print('\nDelta Nji\n', self.__deltaNji)
            print('\nsurface probability, expected number, standard deviation:\n', surfaceProbability,
                  '\n', probableNumber, '\n', standardDeviation)
            input()
    
        lower = np.zeros([constants.clumpLogMass[ens].size, 1])
        clumpLower = np.maximum(lower,
                                np.floor(clumpProbableNumber-constants.n_sigma*clumpStandardDeviation))
        clumpUpper = np.minimum(ensemble.clumpNormalisedDeltaNji[ens],
                                np.ceil(clumpProbableNumber+constants.n_sigma*clumpStandardDeviation))
        CLmaxLower = np.maximum(lower,
                                np.floor(CLmaxProbableNumber-constants.n_sigma*CLmaxStandardDeviation))
        CLmaxUpper = np.minimum(ensemble.clumpNormalisedNj[ens].T,
                                np.ceil(CLmaxProbableNumber+constants.n_sigma*CLmaxStandardDeviation))
        if verbose:
            print('\nupper,lower:\n', clumpUpper, '\n', clumpLower)
        
        ensemble.clumpNumberRange[ens] = np.array([clumpLower, clumpUpper+1], dtype=np.int).T
        ensemble.CLmaxNumberRange[ens] = np.array([CLmaxLower, CLmaxUpper+1], dtype=np.int).T

        if verbose:
            print('\nMasspoint number range:\n', ensemble.clumpNumberRange[ens].round())
            print('\nMasspoint number range:\n', ensemble.interclumpNumberRange.round())
        
        ensemble.clumpCombinations[ens] = calculateCombinations(ensemble.clumpNumberRange[ens])
        ensemble.CLmaxCombinations[ens] = calculateCombinations(ensemble.CLmaxNumberRange[ens])

        #########################################################################################
    
        if verbose:
            input(clumpCombinations)
            input(interclumpCombinations)
        
        largestCombination = ((ensemble.clumpNumberRange[ens][:,:,1]-0).prod(1))
        ensemble.clumpLargestIndex[ens] = largestCombination.argmax()
        ensemble.clumpLargestCombination[ens] = int(largestCombination.max())
        
        if verbose:
            print('\nCombinations:\n{}'.format(ensemble.clumpCombinations[ens]))
            input('\nLargest number of combinations: ' +
                  '{}\n'.format(ensemble.clumpCombinations[ensemble.clumpLargestCombination[ens]]))
            print('\nCombinations:\n{}'.format(ensemble.interclumpCombinations))
            input('\nLargest number of combinations: ' +
                  '{}\n'.format(ensemble.interclumpCombinations[ensemble.interclumpLargestCombination]))
        
        # Probability (clean-up soon)
        # Clump
        probabilityList = []
        combinationIndeces = []
        clumpLargestCombination = ensemble.clumpCombinations[ens][ensemble.clumpLargestIndex[ens]].T
        # print(clumpLargestCombination)
        # loop over combinations of masspoints in each velocity bin
        for i,combinations in enumerate(ensemble.clumpCombinations[ens]):
            ensemble.clumpCombinations[ens][i] = combinations.T
            probability = np.zeros((ensemble.clumpLargestCombination[ens], constants.clumpLogMass[ens].size),
                                   dtype=dtype)
            # maxProbability = np.zeros((ensemble.clumpLargestCombination[ens], constants.clumpLogMass[ens].size),
            #                           dtype=dtype)
            combinationIndeces.append(i)
            if verbose:
                print('\nEnsemble combinations:\n', combinations)
            # calculate the probability if there are any masspoints in this velocity bin
            if combinations.any() and ~(np.isinf(combinations) | np.isnan(combinations)).any():
                # print(i)
                for combination in ensemble.clumpCombinations[ens][i]:
                    # combination = np.array([combination])
                    index = np.where((clumpLargestCombination == combination).all(1))[0][0]
                    # print(index)
                    if verbose:
                        print('\nCombination:\n', combination)
                        input()
                    if constants.probability == 'binomial':
                        if verbose:
                            print(clumpProbableNumber.shape, combination.shape)
                        iGauss = (clumpProbableNumber[:, i] > constants.pnGauss) & \
                                 (ensemble.clumpNormalisedDeltaNji[ens][:, i] > constants.nGauss)
                        # print(iGauss.shape, iGauss)
                        if iGauss.any():
                            # use gauss!
                            if verbose:
                                print('Gauss')
              
                            if not constants.gauss:
                                constants.gauss = True
              
                            # if constants.scipyProbability:
                            probability[index, iGauss] = stats.norm.pdf(combination, loc=clumpProbableNumber[:, i],
                                                                        scale=clumpStandardDeviation[:, i])[iGauss]
                            
                            # else:
                            #   g = stat.Gauss(clumpProbableNumber[:,i], clumpStandardDeviation[:,i], debug=debug)
                            #   probability[index,iGauss] = (g.gaussfunc(combination))[iGauss]
            
                        if (~iGauss).any():
                            # use binomial
                            if verbose:
                                print('Binomial')
                            # <<This will likely print an error when there are more masspoints>>
                            
                            # if constants.scipyProbability:
                            probability[index, ~iGauss] = stats.binom.pmf(combination,
                                                                          ensemble.clumpNormalisedDeltaNji[ens][:, i],
                                                                          clumpSurfaceProbability[:, i])[~iGauss]
                            
                            # else:
                            #   b = stat.Binomial(ensemble.clumpNormalisedDeltaNji[ens][:,i],
                            #                     clumpSurfaceProbability[:,i], debug=debug)  # n and p for binominal
                            #   probability[index,~iGauss] = (b.binomfunc(combination))[~iGauss]
                        
                        # iGauss = ((CLmaxProbableNumber > constants.pnGauss) &
                        #           (ensemble.clumpNormalisedNj[ens] > constants.nGauss))[0]
                        # # print(iGauss.shape, iGauss)
                        # if iGauss.any():
                        #     # print(combination)
                        #
                        #     # if constants.scipyProbability:
                        #     maxProbability[index, iGauss] = stats.norm.pdf(combination,
                        #                                                    loc=CLmaxProbableNumber.flatten(),
                        #                                                    scale=CLmaxStandardDeviation.flatten())[iGauss]
                        #
                        #     # else:
                        #     #   g = stat.Gauss(CLmaxProbableNumber, CLmaxStandardDeviation, debug=debug)
                        #     #   maxProbability[index,iGauss] = (g.gaussfunc(combination))[iGauss]
                        #
                        # if (~iGauss).any():
                        #
                        #     # if constants.scipyProbability:
                        #     maxProbability[index, ~iGauss] = stats.binom.pmf(combination,
                        #                                                      ensemble.clumpNormalisedNj[ens].flatten(),
                        #                                                      CLmaxSurfaceProbability.flatten())[~iGauss]
                        #
                        #     # else:
                        #     #   b = stat.Binomial(ensemble.clumpNormalisedNj[ens],
                        #     #                     CLmaxSurfaceProbability, debug=debug)  # n and p for binominal
                        #     #   maxProbability[index,~iGauss] = (b.binomfunc(combination))[~iGauss]
          
                    elif constants.probability == 'poisson':
                        iGauss = (clumpProbableNumber[:, i] > constants.pnGauss) & \
                                 (ensemble.clumpNormalisedDeltaNji[ens][:, i] > constants.nGauss)
                        if iGauss.any():
                            # use gauss
                            if verbose:
                                print('Gauss')
                            probability[index, iGauss] = stats.norm.pdf(combination, loc=clumpProbableNumber[:, i],
                                                                        scale=clumpStandardDeviation[:, i])[iGauss]
                        if (~iGauss).any():
                            # use poisson
                            if verbose:
                                print('Poisson')
                            probability[index, ~iGauss] = stats.pmf(combination, clumpProbableNumber[:, i])[~iGauss]
                        # iGauss = ((CLmaxProbableNumber > constants.pnGauss) &
                        #           (ensemble.clumpNormalisedNj[ens] > constants.nGauss))[0]
                        # if iGauss.any():
                        #     maxProbability[index, iGauss] = stats.norm.pdf(combination,
                        #                                                    loc=CLmaxProbableNumber.flatten(),
                        #                                                    scale=CLmaxStandardDeviation.flatten())[iGauss]
                        # if (~iGauss).any():
                        #     maxProbability[index, ~iGauss] = stats.pmf(combination, clumpProbableNumber[:, i])[~iGauss]
      
            # if np.shape(probability) != np.shape(maxProbability) and debug:
            #     for i in range(len(probability)):
            #         print(np.array(probability[i]).shape)
            #         print(np.array(maxProbability[i]).shape)
            #     input()
            #print((probability))
            while(len(probability) < ensemble.clumpLargestCombination[ens]):
                probability.append(np.zeros((constants.clumpLogMass[ens].size)))
                maxProbability.append(np.zeros((constants.clumpLogMass[ens].size)))
                if debug:
                    print('Probability length:', len(probability), ', last element shape:', probability[-1].shape)
                    input()
            if (np.array(probability[-1]) == np.nan).any():
                print('\n\nThere is an invalid probability:', probability[-1], '\n\n')
                input()
            # print(probabilityList, probability)
            probabilityList.append(np.array(probability))
            # if i == ensemble.clumpLargestIndex[ens]:
            #     maxProbabilityList = np.array(maxProbability)
        clumpCombinationIndeces = np.array(combinationIndeces)
        # maxProbability = np.array(maxProbabilityList, dtype=dtype)
        ensemble.clumpProbability[ens] = np.array(probabilityList)
        # ensemble.CLmaxProbability[ens] = maxProbability
        for i, combinations in enumerate(ensemble.CLmaxCombinations[ens]):
            ensemble.CLmaxCombinations[ens][i] = combinations.T
            maxProbability = np.zeros_like(combinations.T, dtype=dtype)
            for index, combination in enumerate(ensemble.CLmaxCombinations[ens][i]):
                    print(combination)
                    if constants.probability == 'binomial':
                        # if verbose:
                        print(CLmaxProbableNumber.shape, combination.shape)

                        iGauss = ((CLmaxProbableNumber > constants.pnGauss) &
                                  (ensemble.clumpNormalisedNj[ens] > constants.nGauss))[0]
                        print(iGauss.shape, iGauss)
                        if iGauss.any():
                            # print(combination)

                            # if constants.scipyProbability:
                            maxProbability[index, iGauss] = stats.norm.pdf(combination,
                                                                           loc=CLmaxProbableNumber.flatten(),
                                                                           scale=CLmaxStandardDeviation.flatten())[iGauss]

                            # else:
                            #   g = stat.Gauss(CLmaxProbableNumber, CLmaxStandardDeviation, debug=debug)
                            #   maxProbability[index,iGauss] = (g.gaussfunc(combination))[iGauss]

                        if (~iGauss).any():

                            # if constants.scipyProbability:
                            maxProbability[index, ~iGauss] = stats.binom.pmf(combination,
                                                                             ensemble.clumpNormalisedNj[ens].flatten(),
                                                                             CLmaxSurfaceProbability.flatten())[~iGauss]
                    elif constants.probability == 'poisson':
                        iGauss = ((CLmaxProbableNumber > constants.pnGauss) &
                                  (ensemble.clumpNormalisedNj[ens] > constants.nGauss))[0]
                        if iGauss.any():
                            maxProbability[index, iGauss] = stats.norm.pdf(combination,
                                                                           loc=CLmaxProbableNumber.flatten(),
                                                                           scale=CLmaxStandardDeviation.flatten())[iGauss]
                        if (~iGauss).any():
                            maxProbability[index, ~iGauss] = stats.pmf(combination, clumpProbableNumber[:, i])[~iGauss]
        maxProbability = np.array(maxProbability, dtype=dtype)
        ensemble.CLmaxProbability[ens] = maxProbability

        # ensemble.CLmaxProbability[ens] = np.array([maxProbability[i]/maxProbability.prod(1).sum(0)
        #                                            for i in range(maxProbability.shape[0])], dtype=np.float128)
  
    for ens in range(len(constants.clumpMassNumber)):
        if debug:
            print('Clump')
            for i in range(len(ensemble.clumpProbability[ens])):
                input(np.array(ensemble.clumpProbability[ens][i]).shape)
        if verbose:
            print('\nProbability ({}):\n{}\n'.format(np.shape(ensemble.clumpProbability[ens]), ensemble.clumpProbability[ens]))
            if debug:
                for i in range(len(ensemble.clumpProbability[ens])):
                    print('Probability shapes:\n{}\n'.format(np.array(ensemble.clumpProbability[ens][i].shape)))
                for i in clumpCombinationIndeces:
                    print('Combination sizes:\n{}\n'.format(np.array(ensemble.clumpCombinations[ens][i].size)))
            input()
        ensemble.clumpIndeces[ens] = ensemble.clumpProbability[ens].prod(2).sum(1).nonzero()[0]
        if ensemble.clumpIndeces[ens].size > constants.clumpMaxIndeces[ens]:
            constants.clumpMaxIndeces[ens] = ensemble.clumpIndeces[ens].size
        # print(type(ensemble.clumpCombinations[ens]))
        ensemble.clumpCombinations[ens] = ensemble.clumpCombinations[ens][ensemble.clumpIndeces[ens][0]:
                                                                          (ensemble.clumpIndeces[ens][-1]+1)]
        # ensemble.clumpCombinations[ens] = np.array(ensemble.clumpCombinations[ens][ensemble.clumpIndeces[ens][0]:
        #                                                                            (ensemble.clumpIndeces[ens][-1]+1)])
        # print(type(ensemble.clumpCombinations[ens]))
        clumpProbability = ensemble.clumpProbability[ens][ensemble.clumpIndeces[ens], :, :]
        ensemble.clumpProbability[ens] = np.array([clumpProbability[i]/clumpProbability.prod(2).sum(1)[i]
                                                   for i in range(clumpProbability.shape[0])])
        ensemble.clumpLargestIndex[ens] = np.where(ensemble.clumpIndeces[ens] == ensemble.clumpLargestIndex[ens])[0][0]
    # print(constants.clumpMaxIndeces)
    
    return


def calculate(afuv, debug=False, test=False):
    '''
    !!!!  <<DEGENERATE>>
  
    Maybe <<PARALLELISE>> this??
  
    This is a function to cycle through the Combination instances to create a large numpy.ndarray,
    which is used to calculate the final sums needed for the voxel.
    '''
  
    intensityresult = []
    opticaldepthresult = []
    
    for combination in ensemble.clumpCombinationObjects:
        if (combination in ensemble.clumpCombinationObjects[-10:]) and test:
            result = combination.calculateEmission(Afuv=afuv, test=True)
        else:
            result = combination.calculateEmission(Afuv=afuv)
    
        intensityresult.append(result[0])
        opticaldepthresult.append(result[1])
  
    intensityresult = np.array(intensityresult)
    opticaldepthresult = np.array(opticaldepthresult)
  
    intensity = intensityresult.sum(0)
    opticaldepth = -np.log((opticaldepthresult.sum(0)).astype(np.float))
    gc.collect()
  
    return (intensity,opticaldepth)


def reinitialise():
    # Reinitialise all temporary variables to the correct number of clump sets.
  
    ensemble.clumpMass = 0
  
    ensemble.clumpNj = [[] for _ in range(len(constants.clumpMassNumber))]
    ensemble.clumpDeltaNji = [[] for _ in range(len(constants.clumpMassNumber))]
    ensemble.clumpNormalisedNj = [[] for _ in range(len(constants.clumpMassNumber))]
    ensemble.clumpNormalisedDeltaNji = [[] for _ in range(len(constants.clumpMassNumber))]
    ensemble.clumpSurfaceProbability = [[] for _ in range(len(constants.clumpMassNumber))]
    ensemble.clumpProbableNumber = [[] for _ in range(len(constants.clumpMassNumber))]
    ensemble.clumpStandardDeviation = [[] for _ in range(len(constants.clumpMassNumber))]
    ensemble.CLmaxSurfaceProbability = [[] for _ in range(len(constants.clumpMassNumber))]
    ensemble.CLmaxProbableNumber = [[] for _ in range(len(constants.clumpMassNumber))]
    ensemble.CLmaxStandardDeviation = [[] for _ in range(len(constants.clumpMassNumber))]
  
    ensemble.clumpNumberRange = [[] for _ in range(len(constants.clumpMassNumber))]
    ensemble.CLmaxNumberRange = [[] for _ in range(len(constants.clumpMassNumber))]

    ensemble.clumpCombinations = [[] for _ in range(len(constants.clumpMassNumber))]
    ensemble.CLmaxCombinations = [[] for _ in range(len(constants.clumpMassNumber))]

    ensemble.clumpLargestCombination = [0 for _ in range(len(constants.clumpMassNumber))]
    ensemble.clumpLargestIndex = [0 for _ in range(len(constants.clumpMassNumber))]
  
    ensemble.clumpProbability = [[] for _ in range(len(constants.clumpMassNumber))]
    ensemble.CLmaxProbability = [[] for _ in range(len(constants.clumpMassNumber))]
  
    ensemble.clumpIndeces = [[] for _ in range(len(constants.clumpMassNumber))]
  
    return


def print_ensembleparameters():
  
    np.set_printoptions(precision=4, suppress=True)
  
    for i in range(len(ensemble.clumpNj)):
        print('\nC L U M P   S E T   {}\n'.format(i+1))
        print('Nj:\n{}'.format(ensemble.clumpNj[i].astype(np.float)))
        print('delta Nji:\n{}\n'.format(ensemble.clumpDeltaNji[i].astype(np.float).T))
        print('Normalised Nj:\n{}'.format(ensemble.clumpNormalisedNj[i].astype(np.float)))
        print('Normalised delta Nj:\n{}\n'.format(ensemble.clumpNormalisedDeltaNji[i].astype(np.float).T))
    
    return


# jit_module(nopython=False)
