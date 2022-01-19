import importlib as il
from copy import copy, deepcopy
from time import time

import numpy as np
from numba import jit
from scipy.interpolate import interp1d
from logging import getLogger

from kosmatau3d.models import (
  constants,
  interpolations,
  observations,
  species,
  ensemble,
  combinations,
  masspoints)


class Voxel(object):
    '''
    This is a class to handle each voxel in KOSMA-tau^3. It contains ensembles
    of spherical PDR simulations (to mimic the fractal structure of PDRs), a FUV
    field element, and element absorption, separated for dust and molecules. It
    should have one clump ensemble and one interclump ensemble. This is to account
    for the diffuse surroundings of a clump.
    '''
    
    # PRIVATE
    
    def __init__(self, index=0, debugging=False):
      
        self.__index = index  # index of voxel in VoxelGrid, sort of like its ID in the overall model
        
        self.__debugging = debugging
        
        self.__mass = 0
        self.__ensembleMass = 0
        self.__velocity = 0  # velocity of mass at voxel point
        self.__ensembleDispersion = 0  # dispersion of velocity at voxel point
        
        self.__clumpVelocityIndeces = []
        self.__modelMass = []
        self.__volumeFactor = []
    
        self.__intensitySpecies = []   # intensity of species transitions in voxel
        self.__opticalDepthSpecies = []  # optical depth of species transitions in voxel
        self.__intensityDust = []  # intensity of dust continuum in voxel
        self.__opticalDepthDust = []  # optical depth of dust continuum in voxel
        
        self.__FUV = 0.
        self.__Afuv = 0.
        
        self.__x = 0
        self.__y = 0
        self.__z = 0
        
        self.__r = 0
        self.__phi = 0

        # testing flags
        self.test_calc = False
        self.test_opacity = False
        self.test_pexp = False
        self.test_fv = False

        return
  
    def __setMass(self):
        self.__mass = self.__clumpMass.sum()
  
    def __setClumpMass(self, r):
        mass = [interpolations.interpolateClumpMass(r), interpolations.interpolateInterclumpMass(r)]
        self.__clumpMass = constants.clumpMassFactor*np.asarray(mass).mean(1)
        return
  
    def __setVelocity(self, r):
        velocity = interpolations.interpolateRotationalVelocity(r)
        
        if constants.fromEarth:
    
            # Calculate the correction to the voxel velocity vectors
            relativeRpol = np.sqrt((self.__x-constants.rGalEarth)**2+self.__y**2)
            relativePhi = np.arctan2(self.__y, self.__x-constants.rGalEarth)
            relativeSigma = np.arccos((self.__r**2+relativeRpol**2-constants.rGalEarth**2)/(2*self.__r*relativeRpol))
            sigma = np.arctan2(self.__z, abs(self.__x-constants.rGalEarth))
      
            # Correct the relative velocity of the voxel
            velocityEarth = interpolations.interpolateRotationalVelocity(constants.rGalEarth)
            velocityCirc = velocity.mean() - velocityEarth*self.__r/constants.rGalEarth
      
            self.__velocity = np.sign(relativePhi) * velocityCirc * np.sin(relativeSigma) * np.cos(sigma)
      
            if self.__r == 0:
                self.__velocity = 0
            # self.__velocity = (velocity.mean()) * np.sin(self.__phi)
        
        else:
            self.__velocity = np.array(velocity)
    
        # Use this to check the evaluation of the velocity field. It is still not working correctly...
        # print(self.__velocity)
    
        ensembleDispersion = interpolations.interpolateVelocityDispersion(r)
        
        self.__ensembleDispersion = ensembleDispersion.mean()
        
        return
  
    def __setEnsembleDensity(self, r):
        density = interpolations.interpolateDensity(r)
        self.__ensembleDensity = constants.densityFactor*density.mean()
        return
  
    def __setExtinction(self):
        self.__UVextinction = interpolations.interpolateFUVextinction(self.__ensembleDensity,
                                                                      self.__clumpMass+self.__interclumpMass)
        return
  
    def __setFUV(self, r, z):
        # This is in units of the Draine field
        fuv = interpolations.interpolateFUVfield(r, z)/constants.normUV*constants.globalUV
        self.__FUV = np.clip(fuv, 1, None)
        #self.__FUV = FUVfield(fuv)
        return
  
    def __str__(self):
        return 'Voxel {}'.format(self.__index)

    # PUBLIC
  
    # def reloadModules(self):
    #   il.reload(Ensemble)
    #   il.reload(FUVfield)
    #   self.__clump.reloadModules()
    #   self.__interclump.reloadModules()
    #   return
  
    def setIndex(self, index):
        self.__index = index
        return
  
    def getIndex(self):
        return self.__index
  
    def setPosition(self, x, y, z, r, phi):
        self.__x = x
        self.__y = y
        self.__z = z
        self.__r = r
        self.__phi = phi
        return
  
    def getFUV(self):
        return self.__FUV
  
    def setProperties(self, voxel_size=1, molecules='all', dust='PAH', alpha=1.84, gamma=2.31, clumpMassNumber=[3,1],
                      clumpMassRange=[[0,2],[-2]], clumpNmax=[1, 100], velocityRange=[-10,10], velocityNumber=51,
                      ensembleMass=100, velocity=0., ensembleDispersion=1, velocity_resolution=3, volumeFactor=None,
                      ensembleDensity=[15000, 1911], FUV=[20000, 1], crir=2e-16,
                      from_grid=False, new_grid=False, change_interpolation=False,
                      timed=False, verbose=False, debug=False):
        '''
          This method calculates the radii assuming an origin of (0,0). It then averages
         over a subgrid of 3x3. It might be improved later by having functionality to
         change the subgrid dimensions.
    
         There are a few kwargs that can be used to initialise a Voxel instance:
    
         MODEL PARAMETERS
    
         alpha: The power law distribution required for initial mass function dN(M)/dM = M^-alpha.
                The default is 1.84 as determined in Heithausen et al. (1998).
         gamma: The power law mass-size relation satisfying M = C R^gamma. The default is 2.31 as
                determined in Heithausen et al. (1998).
         velocityRange: The range (list) of the observed radial velocities. By default it
                        is [-10, 10].
         velocityNumber: The number of observed velocities in the specified range (including
                         endpoints). The default is 51.
         clumpMassNumber: The number of clump masses for each set. This should be a list of
                          integers with the same length as clumpMassRange. The default is
                          [3, 1].
         clumpMassRange: The range (list) of clump masses for each set. This is a list of
                         ranges with the same length as clumpMassNumber. The ranges can have
                         length 1, in which case it only evaluates that clump mass. The
                         masses are in units of dex(Msol). By default it is [[0, 2], [-2]].
         clumpDensity: This changes the default clump density in the constants module. This
                       overrides the voxel density kwarg, and is useful for creating interclumps
                       medium that is not scaled according to the clumps. the default is
                       [None, 1911], so the first clump set takes the voxel density. This is in
                       units of cm^-3.
         clumpFUV: This changes the default clump FUV in the constants module. This overrides
                   the voxel FUV kwarg. It is a list of length of the number of clump sets.
                   Use a value of None to use the voxel FUV field. The default is [None, 1], in
                   units of the Draine field.
         clumpNmax: The maximum number of the largest clump in the ensemble. This helps to limit
                    the calculation of a large combination of clumps by rescaling the voxel (thus
                    eliminating any change to the brightness temperature). Larger combinations are used
                    if the rescaled voxel is smaller than the largest clump. The default is [1, 100].
         voxel_size: The voxel size in parsecs. The default is 1 pc.
         molecules: The molecules included in the model. This is a list of strings, where each
                    string has the element name (as in the grid) followed by the transition
                    number (Which is taken from KOSMA-tau). It is set to the first transition
                    of C+, C, and CO by default.
         dust: This is a string to limit how much of the dust continuum is included. The
               available continuum ranges from 1 nm to 3.1 mm. This can be limited to the range
               of molecular lines by setting dust='molecular' (22 dust lines), or to the range of
               PAH dust features by setting dust='PAH' (201 dust lines). All dust lines (333)
               will be calculated if dust=''
    
         VOXEL PROPERTIES
    
         ensembleMass: This can be either a float or list of floats, depending on how many clump
                      sets you have. The default is 100 Msol.
         velocity: The radial velocity of the observed Voxel instance. This is used with the
                   model velocities to determine the number of clumps seen at a given velocity.
                   The default is 0 km/s.
         ensembleDispersion: The velocity dispersion of the ensemble. This is used with the
                             model velocities to determine the number of clumps at a given
                             velocity. The default is 1 km/s.
         ensembleDensity: The observed hydrogen density in the voxel. This is different from and overridden by
                  any default density specified in the constants module. the default is
                  15000 cm^-3.
         FUV: The FUV field of the voxel. All clumps are isotropically radiated with this
              field. It is different from and overridden by any default FUV field specified
              in the constants module. The default is 20000 Draine units.
         crir: The primary cosmic ray ionisation rate with respect to molecular hydrogen (zeta_H2). By default,
               the local rate is used (2e-16).
    
         GENERAL FLAGS
    
         fromFile: This flag can be set to retrieve the voxel properties from a file. It is
                   False by default.
         verbose: This is mainly used for debugging purposes. It prints various statements
                  as the code is evaluated
    
        '''
        # print('Voxel instance initialised')

        self.__logger = getLogger(__name__)
        if debug:
            self.__logger.setLevel('DEBUG')
        elif verbose:
            self.__logger.setLevel('INFO')
        else:
            self.__logger.setLevel('WARNING')
        
        if timed:
            t0 = time()
    
        x, y = np.meshgrid(np.linspace(self.__x-.5*constants.voxel_size, self.__x+.5*constants.voxel_size, 3),
                           np.linspace(self.__y-.5*constants.voxel_size, self.__y+.5*constants.voxel_size, 3))
        r = np.array([x.flatten(), y.flatten()]).T
        r = np.linalg.norm(r, axis=1)

        if constants.voxel_size != voxel_size:
            constants.voxel_size = voxel_size

        if constants.alpha != alpha or constants.gamma != gamma:
            constants.changeMassFunctionParameters(alpha=alpha, gamma=gamma)

        if constants.velocityBin != velocityRange or constants.velocityNumber != velocityNumber:
            constants.changeVelocityRange(velocityRange)
            constants.changeVelocityNumber(velocityNumber)

        if constants.velocity_resolution != velocity_resolution:
            constants.velocity_resolution = velocity_resolution

        if constants.clumpLogMassRange != clumpMassRange or constants.clumpMassNumber != clumpMassNumber\
                or constants.clumpNmax != clumpNmax:
            constants.addClumps(massRange=clumpMassRange, num=clumpMassNumber, Nmax=clumpNmax, reset=True)

        if new_grid or change_interpolation or \
                not interpolations.initialised or not observations.initialised:
            constants.changeDustWavelengths(dust)
            observations.methods.initialise()
            species.addMolecules(molecules)
            interpolations.initialise()

        if timed:
            t1 = time()
            self.__logger.info('Model setup: {}'.format(t1-t0))

        masspoints.reinitialise()
        combinations.reinitialise()
        ensemble.reinitialise()

        self.__velocity = velocity

        if isinstance(ensembleMass, list) or isinstance(ensembleMass, np.ndarray):
            self.__ensembleMass = ensembleMass
        else:
            self.__ensembleMass = [ensembleMass] * len(constants.clumpMassNumber)

        if isinstance(ensembleDispersion, list) or isinstance(ensembleDispersion, np.ndarray):
            self.__ensembleDispersion = ensembleDispersion
        else:
            self.__ensembleDispersion = [ensembleDispersion] * len(constants.clumpMassNumber)

        if isinstance(volumeFactor, float) or isinstance(volumeFactor, int):
            volumeFactor = [volumeFactor] * len(constants.clumpMassNumber)
        if volumeFactor:
            ensembleDensity = [ensembleMass[ens]*constants.massSolar/constants.massH/volumeFactor[ens] /
                               constants.voxel_size**3/constants.pc**3/100**3
                               for ens in range(len(constants.clumpMassNumber))]
        if isinstance(ensembleDensity, list) or isinstance(ensembleDensity, np.ndarray):
            self.__ensembleDensity = ensembleDensity
        else:
            self.__ensembleDensity = [ensembleDensity] * len(constants.clumpMassNumber)

        if isinstance(FUV, list) or isinstance(FUV, np.ndarray):
            self.__FUV = FUV
        else:
            self.__FUV = [FUV] * len(clumpMassNumber)

        self.__crir = crir

        velocity = self.__velocity
    
        # This will allow the code to reuse the standard clump density constants for voxel sets (ie. do not alter the model constants)
        # density = copy(constants.clumpDensity)
        # fuv = copy(constants.clumpFUV)
        # for i in range(len(density)):
            # if density[i]=='auto': density[i] = self.__ensembleDensity
            # if fuv[i]==None: fuv[i] = self.__FUV
    
        masspoints.setMasspointData(density=self.__ensembleDensity, FUV=self.__FUV, crir=self.__crir)
        ensemble.initialise(ensembledispersion=self.__ensembleDispersion, ensemblemass=self.__ensembleMass)
        combinations.initialise(clumpcombination=[ensemble.clumpCombinations[ens][ensemble.clumpLargestIndex[ens]]
                                                  for ens in range(len(constants.clumpMassNumber))])
        
        if timed:
            t2 = time()
            self.__logger.info('Modules initialised:'.format(t2-t1))
    
        for ens in range(len(constants.clumpMassNumber)):
            self.__modelMass.append((ensemble.clumpDeltaNji[ens].sum(1)*10**constants.clumpLogMass[ens]).sum())
            self.__volumeFactor.append((ensemble.clumpNj[ens]*(masspoints.clumpRadius[ens]**3*4*np.pi/3)).sum() /
                                       constants.voxel_size**3)
            if abs(self.__modelMass[ens]-self.__ensembleMass[ens]) > 0.1*self.__ensembleMass[ens]:
                self.__logger.error('ERROR: Voxel {} mass difference for clump set {} greater than 10%'
                                    .format(self.__index, ens+1))
    
        # clumpAfuv,interclumpAfuv = combinations.getAfuv()
        # clumpAfuv *= ensemble.CLmaxProbability
        # interclumpAfuv *= ensemble.ICmaxProbability
    
        self.__clumpVelocities = copy(ensemble.clumpVelocities)
        self.__clumpVelocityIndeces = copy(ensemble.clumpIndeces)   # all list entries should be the same
        
        self.__intensitySpecies = [np.zeros((constants.velocityRange.size,
                                             len(species.molecules)), dtype=np.float64)
                                   for _ in range(constants.ensembles)]
        self.__opticalDepthSpecies = [np.zeros((constants.velocityRange.size,
                                                len(species.molecules)), dtype=np.float64)
                                      for _ in range(constants.ensembles)]
        self.__intensityDust = [np.zeros((constants.velocityRange.size,
                                          constants.wavelengths[constants.nDust].size), dtype=np.float64)
                                for _ in range(constants.ensembles)]
        self.__opticalDepthDust = [np.zeros((constants.velocityRange.size,
                                             constants.wavelengths[constants.nDust].size), dtype=np.float64)
                                   for _ in range(constants.ensembles)]

        self.__emissivity_species = deepcopy(self.__intensitySpecies)
        self.__absorption_species = deepcopy(self.__opticalDepthSpecies)
        self.__emissivity_dust = deepcopy(self.__intensityDust)
        self.__absorption_dust = deepcopy(self.__opticalDepthDust)
    
        # This gives an error if there are too many clumps in a line-of-sight; tau_FUV is too large for this equation...
        Afuv = combinations.getAfuv()
        self.__tauFUV = [-np.log((ensemble.CLmaxProbability[ens].prod(1)*Afuv[ens]).sum())
                         for ens in range(len(constants.clumpMassNumber))]
        
        if timed:
            self.__logger.info('setProperties() time of execution:'.format(time()-t0))

        self.__logger.info('voxel initialised')

        return
  
    def getPosition(self):
        return np.array([self.__x, self.__y, self.__z])
  
    def getDensity(self):
        return self.__ensembleDensity
  
    def getEnsembleMass(self):
        return self.__ensembleMass
  
    def getModelMass(self):
        return self.__modelMass
  
    def getVolumeFillingFactor(self):
        return self.__volumeFactor
  
    def getVelocity(self):
        return self.__velocity, self.__ensembleDispersion
  
    def getClumpVelocity(self):
        return self.__clumpVelocities, self.__clumpVelocityIndeces
  
    def getFUVabsorption(self):
        return self.__tauFUV
  
    # @jit
    def calculateEmission(self, test_calc=False, test_opacity=False, test_pexp=False, test_fv=False,
                          verbose=False, timed=False):

        self.test_calc = test_calc
        self.test_opacity = test_opacity
        self.test_pexp = test_pexp
        self.test_fv = test_fv

        if timed:
            t0 = time()
        
        masspoints.calculateEmission(tauFUV=self.__tauFUV, timed=timed)
        
        if timed:
            t1 = time()
            self.__logger.info('Masspoint emission calculated:',format(t1-t0))
        
        combinations.calculateEmission(test_calc=test_calc, test_opacity=test_opacity,
                                       test_fv=test_fv, f_v=self.__volumeFactor)
        
        if timed:
            t2 = time()
            self.__logger.info('Combination emission calculated:'.format(t2-t1))
    
        clumpIntensity = [[] for _ in range(len(constants.clumpMassNumber))]
        clumpOpticalDepth = [[] for _ in range(len(constants.clumpMassNumber))]
        clumpIntensityDust = [[] for _ in range(len(constants.clumpMassNumber))]
        clumpOpticalDepthDust = [[] for _ in range(len(constants.clumpMassNumber))]
    
        # interclumpIntensity = []
        # interclumpOpticalDepth = []
    
        iDust = constants.wavelengths[constants.nDust].size
        f_ds = [np.maximum(1, self.__volumeFactor[ens]) for ens in range(len(constants.clumpMassNumber))]
    
        # Clump
        for ens in range(constants.ensembles):
    
            vel = constants.velocityRange   # v_obs
            clumpVel = ensemble.clumpVelocities[ens]   # v_sys
            # print(vel.shape, clumpVel.shape)
            nv = np.abs(vel-self.__velocity) <= 4*np.maximum(self.__ensembleDispersion[ens],
                                                             constants.clumpDispersion)
            factor = 1/np.sqrt(2*np.pi*constants.clumpDispersion**2) \
                     * np.exp(-(vel[nv].reshape(1, -1)-clumpVel.reshape(-1, 1)-self.__velocity)**2
                              /2/constants.clumpDispersion**2)
            
            # clumpIntensity[ens] = np.zeros((ensemble.clumpVelocities[ens].size,
            #                                 constants.velocityRange.size,
            #                                 combinations.clumpIntensity[ens].shape[0],
            #                                 combinations.clumpIntensity[ens].shape[1]))
            # clumpIntensityDust[ens] = np.zeros((constants.velocityRange.size,
            #                                     combinations.clumpIntensity[ens].shape[0],
            #                                     combinations.clumpIntensity[ens].shape[1]))
            # clumpOpticalDepth[ens] = np.zeros((ensemble.clumpVelocities[ens].size,
            #                                    constants.velocityRange.size,
            #                                    combinations.clumpIntensity[ens].shape[0],
            #                                    combinations.clumpIntensity[ens].shape[1]))
            # clumpOpticalDepthDust[ens] = np.zeros((constants.velocityRange.size,
            #                                        combinations.clumpIntensity[ens].shape[0],
            #                                        combinations.clumpIntensity[ens].shape[1]))
      
            if nv.any():
            
                for i, probability in enumerate(ensemble.clumpProbability[ens]):
          
                    if timed:
                        t3 = time()
                        self.__logger.info('Start I_xi calculation:'.format(t3-t2))
            
                    intensity = copy(combinations.clumpIntensity[ens][:, iDust:])
                    # shape (v_obs, combination, wavelength)
                    intensity = np.array([intensity*factor[ensemble.clumpIndeces[ens][i], j]
                                          for j in range(factor.shape[1])])
                    # shape (v_sys, v_obs, combination, wavelength)
                    clumpIntensity[ens].append(np.array([(probability.prod(1)/probability.prod(1).sum(0) *
                                                          intensity[j].T).T
                                                         for j in range(factor.shape[1])]))
                    # shape (combination, wavelength)
                    intensityDust = copy(combinations.clumpIntensity[ens][:, :iDust])
                    # shape (v_obs, combination, wavelength)
                    clumpIntensityDust[ens].append((probability.prod(1)/probability.prod(1).sum(0)
                                                    * intensityDust.T).T)
            
                    opticalDepth = copy(combinations.clumpOpticalDepth[ens][:, iDust:])
                    opticalDepth = np.array([opticalDepth*factor[ensemble.clumpIndeces[ens][i], j]
                                             for j in range(factor.shape[1])])
                    opticalDepthDust = copy(combinations.clumpOpticalDepth[ens][:, :iDust])
                    if self.test_pexp:
                        clumpOpticalDepthDust[ens].append((probability.prod(1)/probability.prod(1).sum(0)
                                                           * opticalDepthDust.T).T)
                        clumpOpticalDepth[ens].append(np.array([(probability.prod(1)/probability.prod(1).sum(0)
                                                                 * opticalDepth[j].T).T
                                                                for j in range(factor.shape[1])]))
                    elif self.test_opacity:
                        clumpOpticalDepthDust[ens].append((probability.prod(1)/probability.prod(1).sum(0)
                                                           * np.exp(-opticalDepthDust.T
                                                                    *constants.voxel_size*f_ds[ens])).T)
                        clumpOpticalDepth[ens].append(np.array([(probability.prod(1)/probability.prod(1).sum(0)
                                                                 * np.exp(-opticalDepth[j].T
                                                                          *constants.voxel_size*f_ds[ens])).T
                                                                for j in range(factor.shape[1])]))
                    else:
                        clumpOpticalDepthDust[ens].append((probability.prod(1)/probability.prod(1).sum(0)
                                                           * np.exp(-opticalDepthDust.T)).T)
                        clumpOpticalDepth[ens].append(np.array([(probability.prod(1)/probability.prod(1).sum(0)
                                                                 * np.exp(-opticalDepth[j].T)).T
                                                                for j in range(factor.shape[1])]))
            
                    if timed:
                        t4 = time()
                        self.__logger.info('End I_xi calculation:'.format(t4-t3))
          
                # All of these have shape (velocity, wavelength)
                self.__intensitySpecies[ens][nv, :] = self.__intensitySpecies[ens][nv, :] + \
                                                      (np.array(clumpIntensity[ens]).sum(2)
                                                       ).sum(0).astype(constants.dtype)
                # self.__intensityDust[ens][:,:] = self.__intensityDust[ens][:,:] +
                #                                  np.array([np.array(clumpIntensityDust[ens]).sum(1).sum(0)
                #                                  for _ in range(factor.shape[1])]).astype(constants.dtype)
                # self.__opticalDepthDust[ens][:,:] = self.__opticalDepthDust[ens][:,:] +
                #                                     np.array([-np.log(np.array(clumpOpticalDepthDust[ens]).sum(1)).sum(0)
                #                                     for _ in range(factor.shape[1])]).astype(constants.dtype)
                #                                     (factor.shape[1])]).astype(constants.dtype)
                self.__intensityDust[ens][:, :] = (self.__intensityDust[ens][:, :] +
                                                   np.array([np.array(clumpIntensityDust[ens]).sum(1).sum(0)
                                                             for _ in range(self.__intensityDust[ens].shape[0])
                                                             ]).astype(constants.dtype))
                if self.test_pexp:
                    self.__opticalDepthSpecies[ens][nv, :] = self.__opticalDepthSpecies[ens][nv, :] + \
                                                             (np.array(clumpOpticalDepth[ens]).sum(2)
                                                              ).sum(0).astype(constants.dtype)
                    self.__opticalDepthDust[ens][:, :] = (self.__opticalDepthDust[ens][:, :] +
                                                          np.array([np.array(clumpOpticalDepthDust[ens]).sum(1).sum(0)
                                                                    for _ in range(self.__opticalDepthDust[ens].shape[0])
                                                                    ]).astype(constants.dtype))
                else:
                    self.__opticalDepthSpecies[ens][nv, :] = self.__opticalDepthSpecies[ens][nv, :] + \
                                                             (-np.log(np.array(clumpOpticalDepth[ens]).sum(2))
                                                              ).sum(0).astype(constants.dtype)
                    self.__opticalDepthDust[ens][:, :] = (self.__opticalDepthDust[ens][:, :] +
                                                          np.array([-np.log(np.array(clumpOpticalDepthDust[ens]).sum(1)
                                                                            ).sum(0)
                                                                    for _ in range(self.__opticalDepthDust[ens].shape[0])
                                                                    ]).astype(constants.dtype))
          
                if iDust > 1:
                    intensityDustInterp = interp1d(constants.wavelengths[constants.nDust],
                                                   self.__intensityDust[ens].max(0), fill_value='extrapolate')
                    opticalDepthDustInterp = interp1d(constants.wavelengths[constants.nDust],
                                                      self.__opticalDepthDust[ens].max(0), fill_value='extrapolate')
          
                for i, transition in enumerate(species.moleculeWavelengths):
                    if iDust > 1:
                        self.__intensitySpecies[ens][:, i] += intensityDustInterp(transition)
                        self.__opticalDepthSpecies[ens][:, i] += opticalDepthDustInterp(transition)
                    else:
                        self.__intensitySpecies[ens][:, i] += self.__intensityDust[ens].max()
                        self.__opticalDepthSpecies[ens][:, i] += self.__opticalDepthDust[ens].max()
            
            else:
                self.__logger.warning('Voxel with velocity {} not within given observing velocity range.'
                                      .format(self.__velocity))

        if self.test_calc:
            self.__emissivity_species = deepcopy(self.__intensitySpecies)
            self.__emissivity_dust = deepcopy(self.__intensityDust)
        elif self.test_fv:
            f_ds = [np.maximum(1, self.__volumeFactor[ens]) for ens in range(len(constants.clumpMassNumber))]
            self.__emissivity_species = [self.__intensitySpecies[ens]/constants.voxel_size/f_ds[ens]
                                         for ens in range(len(constants.clumpMassNumber))]
            self.__emissivity_dust = [self.__intensityDust[ens]/constants.voxel_size/f_ds[ens]
                                      for ens in range(len(constants.clumpMassNumber))]
        else:
            self.__emissivity_species = [self.__intensitySpecies[ens]/constants.voxel_size
                                         for ens in range(len(constants.clumpMassNumber))]
            self.__emissivity_dust = [self.__intensityDust[ens]/constants.voxel_size
                                      for ens in range(len(constants.clumpMassNumber))]
        if self.test_opacity:
            self.__absorption_species = deepcopy(self.__opticalDepthSpecies)
            self.__absorption_dust = deepcopy(self.__opticalDepthDust)
        elif self.test_fv:
            f_ds = [np.maximum(1, self.__volumeFactor[ens]) for ens in range(len(constants.clumpMassNumber))]
            self.__absorption_species = [self.__opticalDepthSpecies[ens]/constants.voxel_size/f_ds[ens]
                                         for ens in range(len(constants.clumpMassNumber))]
            self.__absorption_dust = [self.__opticalDepthDust[ens]/constants.voxel_size/f_ds[ens]
                                      for ens in range(len(constants.clumpMassNumber))]
        else:
            self.__absorption_species = [self.__opticalDepthSpecies[ens]/constants.voxel_size
                                         for ens in range(len(constants.clumpMassNumber))]
            self.__absorption_dust = [self.__opticalDepthDust[ens]/constants.voxel_size
                                      for ens in range(len(constants.clumpMassNumber))]
    
        self.__logger.info('Voxel emission calculated.')
          
        if timed:
            self.__logger.info('calculateEmission() time of execution:'.format(time()-t0))
        
        return
  
    def getSpeciesEmissivity(self, include_dust=False):
        epsilon_species = deepcopy(self.__emissivity_species)
        if include_dust:
            epsilon_dust = self.__emissivity_dust
            for ens in range(len(constants.clumpMassNumber)):
                if epsilon_dust[ens].shape[1] > 1:
                    epsilon_dust_interp = interp1d(constants.wavelengths[constants.nDust],
                                                   epsilon_dust[ens].max(0), fill_value='extrapolate')
                for i, transition in enumerate(species.moleculeWavelengths):
                    if epsilon_dust[ens].shape[1] > 1:
                        epsilon_species[ens][:, i] -= epsilon_dust_interp(transition)
                else:
                    epsilon_species[ens] -= epsilon_dust[ens]
        return np.asarray(epsilon_species)
  
    def getSpeciesAbsorption(self, include_dust=False):
        kappa_species = deepcopy(self.__absorption_species)
        if include_dust:
            kappa_dust = self.__absorption_dust
            for ens in range(len(constants.clumpMassNumber)):
                if kappa_dust[ens].shape[1] > 1:
                    kappa_dust_interp = interp1d(constants.wavelengths[constants.nDust],
                                                   kappa_dust[ens].max(0), fill_value='extrapolate')
                    for i, transition in enumerate(species.moleculeWavelengths):
                        if kappa_dust[ens].shape[1] > 1:
                            kappa_species[ens][:, i] -= kappa_dust_interp(transition)
                else:
                    kappa_species[ens] -= kappa_dust[ens]
        return np.asarray(kappa_species)
  
    def getSpeciesIntensity(self, integrated=False, include_dust=False):
        epsilon = self.getSpeciesEmissivity(include_dust=include_dust)
        kappa = self.getSpeciesAbsorption(include_dust=include_dust)
        intensity = np.zeros_like(epsilon)
        for ens in range(len(constants.clumpMassNumber)):
            if self.__volumeFactor[ens] > 1 and self.test_fv:
                ds = self.__volumeFactor[ens]*constants.voxel_size
            else:
                ds = constants.voxel_size
            intensity[ens] = epsilon[ens]/kappa[ens]*(1-np.exp(-kappa[ens]*ds))
            i_nan = np.isnan(intensity[ens])
            intensity[ens][i_nan] = epsilon[ens][i_nan]
        if integrated:
            intensity_final = np.zeros((intensity.shape[0], intensity.shape[2]))
            for ens in range(len(constants.clumpMassNumber)):
                intensity_final[ens] = np.trapz(intensity[ens], constants.velocityRange, axis=0)
        else:
            intensity_final = intensity
        return intensity_final
  
    def getDustEmissivity(self):
        return np.asarray(self.__emissivity_dust)
  
    def getDustAbsorption(self):
        return np.asarray(self.__absorption_dust)
  
    def getDustIntensity(self):
        epsilon = self.getDustEmissivity()
        kappa = self.getDustAbsorption()
        intensity = np.zeros_like(epsilon)
        for ens in range(len(constants.clumpMassNumber)):
            if self.__volumeFactor[ens] > 1 and self.test_fv:
                ds = self.__volumeFactor[ens]*constants.voxel_size
            else:
                ds = constants.voxel_size
            intensity[ens] = epsilon[ens]/kappa[ens]*(1-np.exp(-kappa[ens]*ds))
            i_nan = np.isnan(intensity[ens])
            intensity[ens][i_nan] = epsilon[ens][i_nan]
        return intensity
  
    def plotMolecule(self, molecule='', quantity='intensity', moleculeName='', title='', logscale=False):
        # Plot the molecule emission with respect to the voxel velocity structure
    
        import matplotlib.pyplot as plt
    
        if molecule in species.molecules and isinstance(molecule, str):
            molecule = [molecule]
    
        elif isinstance(molecule, list):
            pass
    
        else:
            molecule = species.molecules
    
        if quantity == 'emissivity':
            value = self.getSpeciesEmissivity()
            ylabel = r'$\epsilon_{\lambda} \ \left( \frac{K}{pc} \right)$'
    
        elif quantity == 'absorption':
            value = self.getSpeciesAbsorption()
            ylabel = r'$\kappa_{\lambda} \ \left( \frac{1}{pc} \right)$'
    
        elif quantity == 'intensity':
            value = self.getSpeciesIntensity()
            ylabel = r'$I_{\lambda} \ \left( K \right)$'
        else:
            self.__logger.warning('Quantity {} not available.'.format(quantity))
            return
    
        fig, axes = plt.subplots(len(value), figsize=(10, 5*len(value)))
    
        for ens in range(constants.ensembles):
    
            if isinstance(axes, np.ndarray):
                ax = axes[ens]
      
            else:
                ax = axes
      
            vel = constants.velocityRange  # [self.__clumpVelocityIndeces[ens]]
            labels = []
      
            for n, mol in enumerate(molecule):
      
                if mol not in species.molecules:
                    self.__logger.warning('Species {} not in model.'.format(mol))
                    continue
        
                i = np.where(mol == np.asarray(species.molecules))[0][0]
                
                if logscale:
                    ax.semilogy(vel, value[ens][:, i], marker='o', ms=2, ls='-', lw=1)
                else:
                    ax.plot(vel, value[ens][:, i], marker='o', ms=2, ls='-', lw=1)
                
                if isinstance(moleculeName, list):
                    labels.append(moleculeName[n])
                elif moleculeName:
                    labels.append(moleculeName)
                else:
                    labels.append(mol)
      
            ax.legend(labels=labels, fontsize=14, bbox_to_anchor=(1.05, 1))
            ax.set_xlabel(r'$V_r \ (\frac{km}{s})$', fontsize=20)
            ax.set_ylabel(ylabel, fontsize=20)
            plt.setp(ax.get_xticklabels(), fontsize=14)
            plt.setp(ax.get_yticklabels(), fontsize=14)
      
            if title:
                ax.set_title(title, fontsize=20)
            else:
                ax.set_title('Clump set {} {}'.format(ens+1, quantity), fontsize=20)
        
        fig.tight_layout()
    
        plt.show()
    
        return
  
    def plotSpectrum(self, quantity='intensity', title=''):
        # Plot the either the intesity or optical depth spectrum at the voxel velocity (the velocity with the largest
        #  ensemble).

        import matplotlib as mpl
        import matplotlib.pyplot as plt
        from matplotlib.lines import Line2D

        mpl.rcParams['text.usetex'] = True
    
        moleculeColour = {
                           'C+': 'xkcd:orange',
                           'C': 'xkcd:red',
                           'O': 'xkcd:violet',
                           'CO': 'xkcd:maroon',
                           '13C+': 'xkcd:light blue',
                           '13C': 'xkcd:blue',
                           '13CO': 'xkcd:sapphire',
                           'HCO+': 'xkcd:neon green',
                           'H13CO+': 'xkcd:forrest green',
                           'H3O+': 'xkcd:vomit'
                         }
    
        constants.resortWavelengths()
        nDust = constants.wavelengths[constants.nDust].size
    
        if quantity == 'emissivity':
            valueDust = self.getDustEmissivity()
            valueSpecies = self.getSpeciesEmissivity()
            ylabel = r'$\epsilon_{\lambda} \ \left( \frac{K}{pc} \right)$'
    
        elif quantity == 'absorption':
            valueDust = self.getDustAbsorption()
            valueSpecies = self.getSpeciesAbsorption()
            ylabel = r'$\kappa_{\lambda} \ \left( \frac{1}{pc} \right)$'
    
        elif quantity == 'intensity':
            valueDust = self.getDustIntensity()
            valueSpecies = self.getSpeciesIntensity()
            ylabel = r'$I_{\lambda} \ \left( K \right)$'
        
        else:
            self.__logger.warning('Quantity {} not available.'.format(quantity))
            return
    
        fig, axes = plt.subplots(constants.ensembles, figsize=(10, 5*constants.ensembles))
    
        for ens in range(constants.ensembles):
    
            if isinstance(axes, np.ndarray):
                ax = axes[ens]
      
            else:
                ax = axes
      
            vel = int(valueSpecies[ens].shape[0]/2)
            
            if valueDust[ens].shape[1] > 1:
                dustInterp = interp1d(constants.wavelengths[constants.nDust], valueDust[ens][vel, :],
                                      fill_value='extrapolate')
                ax.loglog(1e6*constants.wavelengths[constants.nDust], valueDust[ens][vel, :],
                          ls='--', lw=0.5, c='xkcd:black')
            else:
                dust_spectrum = np.full(2, valueDust[ens][vel, 0])
                ax.loglog(1e6*constants.wavelengths[[0, 20]], dust_spectrum, ls='--', lw=0.5, c='xkcd:black')
            
            molecules = []
            colour = []
      
            for i, molecule in enumerate(species.molecules):
                if valueDust[ens].shape[1] > 1:
                    valTemp = dustInterp(species.moleculeWavelengths[i])
                else:
                    valTemp = valueDust[ens][vel, 0]
                if not molecule.split()[0] in molecules:
                    colour.append(moleculeColour[molecule.split()[0]])
                    molecules.append(molecule.split()[0])
                ax.loglog([1e6*species.moleculeWavelengths[i]]*2, [valTemp, valueSpecies[ens][vel, i]],
                          ls='-', lw=1, c=colour[-1])
      
            lines = [Line2D([0], [0], color=moleculeColour[molecule], lw=1) for molecule in molecules]
            # mols  = ['dust', ]
            leg = ax.legend(labels=molecules, handles=lines, fontsize=16, bbox_to_anchor=(1.05, 1))
            
            ax.set_xlabel(r'$\lambda \ (\mu m)$', fontsize=20)
            ax.set_ylabel(ylabel, fontsize=20)
            ax.set_ylim([valueDust[ens][vel, :].min()*10e-3, valueDust[ens][vel, :].max()*10e3])
            plt.setp(ax.get_xticklabels(), fontsize=14)
            plt.setp(ax.get_yticklabels(), fontsize=14)
      
            if title:
                ax.set_title(title, fontsize=20)
            else:
                ax.set_title('Clump Set {} {} spectrum'.format(ens+1, quantity), fontsize=20)
        
        # fig.tight_layout()
    
        plt.show()
    
        return
      
    def printVoxel(self):
        # Coming soon...
        return
