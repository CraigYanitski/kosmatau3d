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
  def __init__(self, species, interpolations, index, vNumber=10):
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
    self.__velocityRange = np.linspace(self.__velocity-self.__velocityDispersion, self.__velocity+self.__velocityDispersion, self.__velocityNumber)
    return
  def __setDensity(self):
    self.__density = self.__interpolations.interpolateDensity(self.__r)
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
    #try self.__FUV.getFUV():
      return 'Voxel {}\n  ->Cartesian position: ({}, {}, {})\n  ->mass {}\n  ->intensity {}\n  ->optical depth {}\n  ->FUV field {}'.format(self.__index, self.__x, self.__y, self.__z, self.__mass, 10**self.__intensity, 10**self.__opticalDepth, self.__FUV.getFUV())
    #except:
    #  return 'Voxel {}\n  ->Cartesian position: ({}, {}, {})\n  ->mass {}\n  ->intensity {}\n  ->optical depth {}\n  ->FUV field {}'.format(self.__index, self.__x, self.__y, self.__z, self.__mass, 10**self.__intensity, 10**self.__opticalDepth, self.__FUV)

  # PUBLIC
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
    self.__setClumpMass()
    self.__setInterclumpMass()
    self.__setMass()
    self.__setVelocity()
    self.__setDensity()
    self.__setExtinction()
    self.__setFUV()
    self.__clump.initialise(mass=self.__clumpMass, density=self.__density, velocity=self.__velocityRange, velocityDispersion=self.__velocityDispersion, FUV=self.__FUV, extinction=self.__UVextinction)
    self.__interclump.initialise(mass=self.__interclumpMass, density=self.__density, velocity=self.__velocityRange, velocityDispersion=self.__velocityDispersion, FUV=self.__FUV, extinction=self.__UVextinction)
    return
  def getPosition(self):
    return (self.__x, self.__y, self.__z)
  def getClumps(self):
    return (self.__clump, self.__interclump)
  def calculateEmission(self, verbose=False):
    if verbose: print('\nCalculating voxel V{} emission'.format(self.__index))
    self.__clump.calculate()
    self.__interclump.calculate()
    iClump,tauClump,FUVclump = self.__clump.getEnsembleEmission()
    iInterclump,tauInterclump,FUVinterclump = self.__interclump.getEnsembleEmission()
    self.__intensity = iClump+iInterclump
    self.__opticalDepth = tauClump+tauInterclump
    if isinstance(FUVclump, FUVfield): self.__FUV = FUVfield(FUVclump.getFUV()+FUCinterclump.getFUV())
    return
  def getEmission(self, verbose=False):
    emission = (self.__intensity, self.__opticalDepth, self.__FUV)
    if verbose: print(emission)
    return emission