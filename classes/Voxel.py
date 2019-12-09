import importlib as il
from Ensemble import *
from FUVfield import *
class Voxel(object):
  '''
  This is a class to handle each voxel in KOSMA-tau^3. It contains ensembles
  of spherical PDR simulations (to mimic the fractal structure of PDRs), a FUV
  field element, and element absorption, separated for dust and molecules. It
  should have one clump ensemble and one interclump ensemble. This is to account
  for the diffuse surroundings of a clump.
  '''
  # PRIVATE
  def __init__(self, species, interpolations, index, vNumber=51):
    self.__species = species     #list of both moleculular and dust species
    self.__interpolations = interpolations
    self.__index = index    #index of voxel in VoxelGrid, sort of like its ID in the overall model
    self.__constants = Constants()
    self.__velocity = 0     #velocity of mass at voxel point
    self.__velocityDispersion = 0   #dispersion of velocity at voxel point
    self.__velocityNumber = vNumber
    self.__velocityRange = 0
    self.__intensity = 0      #intensity of emissions at voxel point
    self.__opticalDepth = 0   #optical depth at voxel point
    self.__FUV = 0
    self.__clump = Ensemble('clump', self.__species, self.__interpolations)    #clumpy ensemble at voxel point
    self.__interclump = Ensemble('interclump', self.__species, self.__interpolations)    #diffuse interclump ensemble at voxel point
    self.__x = 0
    self.__y = 0
    self.__z = 0
    self.__r = 0
    self.__phi = 0
    self.__scale = 0
    return
  def __setMass(self):
    self.__mass = self.__clumpMass+self.__interclumpMass
  def __setClumpMass(self):
    self.__clumpMass = self.__interpolations.interpolateClumpMass(self.__r)
    self.__clump.setMass(self.__clumpMass)
    return
  def __setInterclumpMass(self):
    self.__interclumpMass = self.__interpolations.interpolateInterclumpMass(self.__r)
    self.__interclump.setMass(self.__interclumpMass)
    return
  def __setVelocity(self):
    self.__velocity = self.__interpolations.interpolateRotationalVelocity(self.__r)
    self.__velocityDispersion = self.__interpolations.interpolateVelocityDispersion(self.__r)
    self.__velocityRange = np.linspace(self.__velocity-self.__velocityDispersion/2., self.__velocity+self.__velocityDispersion/2, num=self.__velocityNumber)
    return
  def __setDensity(self, densityFactor=2):
    self.__density = densityFactor*self.__interpolations.interpolateDensity(self.__r)
    return
  def __setExtinction(self):
    self.__UVextinction = self.__interpolations.interpolateFUVextinction(self.__density, self.__clumpMass+self.__interclumpMass)
    return
  def __setFUV(self):
    fuv = self.__interpolations.interpolateFUVfield(self.__r)/self.__constants.normUV*self.__constants.globalUV
    fuv = np.clip(fuv, 1, None)
    self.__FUV = FUVfield(fuv)
    return
  def __str__(self):
    return 'Voxel {}'.format(self.__index)

  # PUBLIC
  def reloadModules(self):
    il.reload(Ensemble)
    il.reload(FUVfield)
    self.__clump.reloadModules()
    self.__interclump.reloadModules()
    return
  def setIndex(self, index):
    self.__index = index
    return
  def setPosition(self, x, y, z, r, phi, scale):
    self.__x = x
    self.__y = y
    self.__z = z
    self.__r = r
    self.__phi = phi
    self.__scale = scale
    return
  def setProperties(self):
    #print('Voxel instance initialised')
    self.__setClumpMass()
    self.__setInterclumpMass()
    self.__setMass()
    self.__setVelocity()
    self.__setDensity()
    self.__setExtinction()
    self.__setFUV()
    self.__clump.initialise(mass=self.__clumpMass, density=self.__density, velocity=self.__velocity, velocityDispersion=self.__velocityDispersion, FUV=self.__FUV, extinction=self.__UVextinction)
    self.__interclump.initialise(mass=self.__interclumpMass, density=1911, velocity=self.__velocity, velocityDispersion=self.__velocityDispersion, FUV=self.__FUV, extinction=self.__UVextinction)
    return
  def getPosition(self):
    return (self.__x, self.__y, self.__z)
  def getClumps(self):
    return (self.__clump, self.__interclump)
  def getVelocity(self):
    return (self.__velocity, self.__velocityDispersion, self.__velocityRange)
  def calculateEmission(self, verbose=False):
    if verbose:
      print('\nCalculating voxel V{} emission'.format(self.__index))
    self.__clump.calculate(test=False)
    self.__interclump.calculate(test=False)
    iClump,tauClump,FUVclump = self.__clump.getEnsembleEmission()
    iInterclump,tauInterclump,FUVinterclump = self.__interclump.getEnsembleEmission()
    if verbose:
      print('\nClump and interclump intensity:', iClump, iInterclump)
      print('\nClump and interclump optical depth:', tauClump, tauInterclump)
      input()
    # Sum over ensembles
    self.__intensity = (iClump+iInterclump)
    self.__opticalDepth = (tauClump+tauInterclump)
    if verbose: print('\nShape: ', self.__intensity.shape, '\n\n')
    if isinstance(FUVclump, FUVfield): self.__FUV = FUVfield(np.average(FUVclump.getFUV()+FUVinterclump.getFUV()))
    return
  def getEmission(self, verbose=False):
    emission = ((self.__intensity.sum(2)), (self.__opticalDepth.sum(2)), self.__FUV)
    if verbose: print(emission)
    return emission
  def printVoxel(self):
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
    print('    -clump')
    for i in range(len(names)):
      print('      {} {}: {}'.format(clumpIntensity[i][0],clumpIntensity[i][1],clumpIntensity[i][2]))
    print('    -interclump')
    for i in range(len(names)):
      print('      {} {}: {}'.format(interclumpIntensity[i][0],interclumpIntensity[i][1],interclumpIntensity[i][2]))
    print('  ->optical depth')
    for i in range(len(names)):
      print('    {} {}: {}'.format(totalTau[i][0],totalTau[i][1],totalTau[i][2]))
    print('    -clump')
    for i in range(len(names)):
      print('      {} {}: {}'.format(clumpTau[i][0],clumpTau[i][1],clumpTau[i][2]))
    print('    -interclump')
    for i in range(len(names)):
      print('      {} {}: {}'.format(interclumpTau[i][0],interclumpTau[i][1],interclumpTau[i][2]))
    print('  ->FUV field {}'.format(self.__FUV.getFUV()))
    return