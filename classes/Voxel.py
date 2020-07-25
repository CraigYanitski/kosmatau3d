import importlib as il
import numpy as np

import constants
import interpolations

import ensemble
import combinations
import masspoints

#from ensemble import Ensemble

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
    
    #self.__species = species     #list of both moleculular and dust species
    
    self.__index = index         #index of voxel in VoxelGrid, sort of like its ID in the overall model
    
    self.__debugging = debugging
    
    self.__velocity = 0     #velocity of mass at voxel point
    self.__ensembleDispersion = 0   #dispersion of velocity at voxel point
    
    self.__clumpVelocityIndeces = []
    self.__interclumpVelocityIndeces = []

    self.__intensity = 0      #intensity of emissions at voxel point
    self.__opticalDepth = 0   #optical depth at voxel point
    
    self.__FUV = 0.
    self.__Afuv = 0.
    
    # self.__clump = Ensemble('clump', debugging=debugging)    #clumpy ensemble at voxel point
    # self.__interclump = Ensemble('interclump', debugging=debugging)    #diffuse interclump ensemble at voxel point
    
    self.__x = 0
    self.__y = 0
    self.__z = 0
    
    self.__r = 0
    self.__phi = 0
    
    return

  def __setMass(self):
    self.__mass = self.__clumpMass+self.__interclumpMass

  def __setClumpMass(self, r):
    mass = interpolations.interpolateClumpMass(r)
    self.__clumpMass = constants.clumpMassFactor*mass.mean()
    #self.__clump.setMass(self.__clumpMass)
    return

  def __setInterclumpMass(self, r):
    mass = interpolations.interpolateInterclumpMass(r)
    self.__interclumpMass = constants.interclumpMassFactor*mass.mean()
    #self.__interclump.setMass(self.__interclumpMass)
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

      if self.__r==0: self.__velocity = 0
      #self.__velocity = (velocity.mean()) * np.sin(self.__phi)
    
    else:
      self.__velocity = np.array(velocity)

    # Use this to check the evaluation of the velocity field. It is still not working correctly...
    #print(self.__velocity)

    ensembleDispersion = interpolations.interpolateVelocityDispersion(r)
    
    self.__ensembleDispersion = ensembleDispersion.mean()
    
    return

  def __setDensity(self, r):
    density = interpolations.interpolateDensity(r)
    self.__density = constants.densityFactor*density.mean()
    return

  def __setExtinction(self):
    self.__UVextinction = interpolations.interpolateFUVextinction(self.__density, self.__clumpMass+self.__interclumpMass)
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

  def setProperties(self, clumpMass=1387835, interclumpMass=1387835, mass=2775651, velocity=0., ensembleDispersion=1, density=14885, FUV=21591, debug=False):
    ''' This method calculates the radii assuming an origin of (0,0). It then averages
       over a subgrid of 3x3. It might be improved later by having functionality to
       change the subgrid dimensions.'''
    #print('Voxel instance initialised')

    x,y = np.meshgrid(np.linspace(self.__x-.5*constants.resolution, self.__x+.5*constants.resolution,3), \
                      np.linspace(self.__y-.5*constants.resolution, self.__y+.5*constants.resolution,3))
    r = np.array([x.flatten(), y.flatten()]).T
    r = np.linalg.norm(r, axis=1)

    if not debug:

      self.__setClumpMass(r)
      self.__setInterclumpMass(r)
      self.__setMass()
      self.__setVelocity(r)
      self.__setDensity(r)
      #self.__setExtinction()
      self.__setFUV(self.__r, self.__z)

      # This can convert the velocity from Cartesian to radial. It is assumed a scalar value is a radial velocity.
      if isinstance(self.__velocity, float):
        velocity = self.__velocity
      else:
        velocity = np.linalg.norm(self.__velocity)

    else:
      
      self.__clumpMass = clumpMass
      self.__interclumpMass = interclumpMass
      self.__mass = mass
      self.__velocity = velocity
      self.__ensembleDispersion = ensembleDispersion
      self.__density = density
      self.__FUV = FUV

      velocity = self.__velocity

    masspoints.setMasspointData(density=self.__density, FUV=self.__FUV)
    ensemble.initialise(velocity=velocity, ensembleDispersion=self.__ensembleDispersion, clumpMass=self.__clumpMass, interclumpMass=self.__interclumpMass)
    combinations.initialise(clumpCombination=ensemble.clumpCombinations[ensemble.clumpLargestIndex], \
                            interclumpCombination=ensemble.interclumpCombinations[ensemble.interclumpLargestIndex])

    # clumpAfuv,interclumpAfuv = combinations.getAfuv()
    # clumpAfuv *= ensemble.CLmaxProbability
    # interclumpAfuv *= ensemble.ICmaxProbability

    self.__clumpVelocityIndeces = ensemble.clumpIndeces
    self.__interclumpVelocityIndeces = ensemble.interclumpIndeces

    Afuv = combinations.getAfuv()
    self.__Afuv = -(np.log10((ensemble.CLmaxProbability.prod(1)*Afuv[0]).sum()) + np.log10((ensemble.ICmaxProbability.prod(1)*Afuv[1]).sum()))
    
    return

  def getPosition(self):
    return np.array([self.__x, self.__y, self.__z])

  # def getClumps(self):
  #   return (self.__clump, self.__interclump)

  def getDensity(self):
    return self.__density

  def getMass(self):
    return self.__mass

  def getClumpMass(self):
    return self.__clumpMass

  def getInterclumpMass(self):
    return self.__interclumpMass

  def getVelocity(self):
    return self.__velocity,self.ensembleDispersion

  def getClumpVelocity(self):
    return constants.velocityRange[self.__clumpVelocityIndeces]

  def getInterclumpVelocity(self):
    return constants.velocityRange[self.__interclumpVelocityIndeces]

  def getAfuv(self):
    return np.float64(self.__Afuv)

  def calculateEmission(self, verbose=False):
    if verbose:
      print('\nCalculating voxel V{} emission\nFUV extinction: {}'.format(self.__index, self.__Afuv))
    # iClump,tauClump = self.__clump.calculate(self.__Afuv, test=False)
    # iInterclump,tauInterclump = self.__interclump.calculate(self.__Afuv, test=False)
    
    masspoints.calculateEmission(Afuv=self.__Afuv)
    combinations.calculateEmission()

    clumpIntensity = []
    clumpOpticalDepth = []
    interclumpIntensity = []
    interclumpOpticalDepth = []
    iDust = constants.wavelengths[constants.nDust].size

    # Clump
    for probability in ensemble.clumpProbability:
      clumpIntensity.append((probability.prod(1)*combinations.clumpIntensity.T).T)
      clumpOpticalDepth.append((probability.prod(1)*np.exp(-combinations.clumpOpticalDepth.T)).T)

    vel = constants.velocityRange[ensemble.clumpIndeces]
    factor = np.exp(-(vel.reshape(1,-1)-vel.reshape(-1,1))**2/2/constants.clumpDispersion**2).sum(1)
    # print('clump factor:', factor)
    shape = np.array(clumpIntensity).shape
    clumpIntensityAvg = np.zeros((shape[0],shape[2]))
    clumpIntensityAvg[:,:iDust] = np.array(clumpIntensity)[:,:,:iDust].sum(1)
    clumpIntensityAvg[:,iDust:] = factor.reshape(-1,1)*np.array(clumpIntensity)[:,:,iDust:].sum(1)
    clumpOpticalDepthAvg = np.zeros((shape[0],shape[2]))
    clumpOpticalDepthAvg[:,:iDust] = -np.log(np.array(clumpOpticalDepth)[:,:,:iDust].sum(1))
    clumpOpticalDepthAvg[:,iDust:] = -np.log(factor.reshape(-1,1)*np.array(clumpOpticalDepth)[:,:,iDust:].sum(1))

    # Interclump
    for probability in ensemble.interclumpProbability:
      interclumpIntensity.append((probability.prod(1)*combinations.interclumpIntensity.T).T)
      interclumpOpticalDepth.append((probability.prod(1)*np.exp(-combinations.interclumpOpticalDepth.T)).T)

    vel = constants.velocityRange[ensemble.interclumpIndeces]
    factor = np.exp(-(vel.reshape(1,-1)-vel.reshape(-1,1))**2/2/constants.clumpDispersion**2).sum(1)
    # print('interclump factor:', factor)
    shape = np.array(interclumpIntensity).shape
    interclumpIntensityAvg = np.zeros((shape[0],shape[2]))
    interclumpIntensityAvg[:,:iDust] = np.array(interclumpIntensity)[:,:,:iDust].sum(1)
    interclumpIntensityAvg[:,iDust:] = factor.reshape(-1,1)*np.array(interclumpIntensity)[:,:,iDust:].sum(1)
    interclumpOpticalDepthAvg = np.zeros((shape[0],shape[2]))
    interclumpOpticalDepthAvg[:,:iDust] = -np.log(np.array(interclumpOpticalDepth)[:,:,:iDust].sum(1))
    interclumpOpticalDepthAvg[:,iDust:] = -np.log(factor.reshape(-1,1)*np.array(interclumpOpticalDepth)[:,:,iDust:].sum(1))
    
    if verbose:
      print('\nClump and interclump intensity:\n', clumpIntensity, '\n', interclumpIntensity)
      print('\nClump and interclump optical depth:\n', clumpOpticalDepth, '\n', interclumpOpticalDepth)
      input()

    # Sum over ensembles
    self.__intensity = (clumpIntensityAvg,interclumpIntensityAvg)
    self.__opticalDepth = (clumpOpticalDepthAvg,interclumpOpticalDepthAvg)

    if verbose: print('\nShape: ', self.__intensity.shape, '\n\n')
    #if isinstance(FUVclump, FUVfield): self.__FUV = FUVfield(np.average(FUVclump.getFUV()+FUVinterclump.getFUV()))
    return

  def getEmission(self, verbose=False):
    if verbose:
      print(self.__x, self.__y, self.__z)
      print(self.__Afuv, '\n')
    emission = (self.__intensity, self.__opticalDepth)
    if verbose: print(emission)
    return emission
    
  def printVoxel(self):
    # OUTDATED
    if self.__debugging:
      iClump,tauClump,FUVclump = self.__clump.getEnsembleEmission()
      iInterclump,tauInterclump,FUVinterclump = self.__interclump.getEnsembleEmission()
    names = self.__species[0].getMolecules() + self.__species[1].getDust()
    transitions = self.__species[0].getTransitions() + self.__species[1].getTransitions()
    totalIntensity = np.array([names, transitions, self.__intensity.sum(2).max(1)]).T
    clumpIntensity = np.array([names, transitions, iClump.sum(2).max(1)]).T
    interclumpIntensity = np.array([names, transitions, iInterclump.sum(2).max(1)]).T
    totalTau = np.array([names, transitions, self.__opticalDepth.sum(2).max(1)]).T
    clumpTau = np.array([names, transitions, tauClump.sum(2).max(1)]).T
    interclumpTau = np.array([names, transitions, tauInterclump.sum(2).max(1)]).T
    print('\nVoxel {}\n  ->Cartesian position: ({}, {}, {})'.format(self.__index, self.__x, self.__y, self.__z))
    print('  ->mass {} M_sol'.format(self.__mass))
    print('  ->intensity')
    for i in range(len(names)):
      print('    {} {}: {}'.format(totalIntensity[i][0],totalIntensity[i][1],totalIntensity[i][2]))
    if self.__debugging:
      print('    -clump')
      for i in range(len(names)):
        print('      {} {}: {}'.format(clumpIntensity[i][0],clumpIntensity[i][1],clumpIntensity[i][2]))
      print('    -interclump')
      for i in range(len(names)):
        print('      {} {}: {}'.format(interclumpIntensity[i][0],interclumpIntensity[i][1],interclumpIntensity[i][2]))
    print('  ->optical depth')
    for i in range(len(names)):
      print('    {} {}: {}'.format(totalTau[i][0],totalTau[i][1],totalTau[i][2]))
    if self.__debugging:
      print('    -clump')
      for i in range(len(names)):
        print('      {} {}: {}'.format(clumpTau[i][0],clumpTau[i][1],clumpTau[i][2]))
      print('    -interclump')
      for i in range(len(names)):
        print('      {} {}: {}'.format(interclumpTau[i][0],interclumpTau[i][1],interclumpTau[i][2]))
    print('  ->FUV field {}'.format(self.__FUV.getFUV()))
    return