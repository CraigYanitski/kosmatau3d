from Ensemble import *
class Voxel(object):
  '''
  This is a class to handle each voxel in KOSMA-tau^3. It contains ensembles
  of spherical PDR simulations (to mimic the fractal structure of PDRs), a FUV
  field element, and element absorption, separated for dust and molecules. It
  should have one clump ensemble and one interclump ensemble. This is to account
  for the diffuse surroundings of a clump.
  '''
  # PRIVATE
  def __init__(self, species, interpolations, index):
    self.__species = species     #list of both moleculular and dust species
    self.__interpolations = interpolations
    self.__index = index    #index of voxel in VoxelGrid, sort of like its ID in the overall model
    self.__velocity = 0     #velocity of mass at voxel point
    self.__velocityDispersion = 0   #dispersion of velocity at voxel point
    self.__intensity = 0      #intensity of emissions at voxel point
    self.__opticalDepth = 0   #optical depth at voxel point
    self.__clump = Ensemble('clump', self.__species, self.__interpolations)    #clumpy ensemble at voxel point
    self.__interclump = Ensemble('interclump', self.__species, self.__interpolations)    #diffuse interclump ensemble at voxel point
    self.__x = 0
    self.__y = 0
    self.__z = 0
    self.__r = 0
    self.__phi = 0
    self.__scale = 0
    return
  def __setClumpMass(self):
    self.__clumpMass = self.__interpolations.interpolateClumpMass(self.__r)
    self.__clump.setMass(self.__clumpMass)
    return
  def __setInterclumpMass(self):
    self.__interclumpMass = self.__interpolations.interpolateInterclumpMass(self.__r)
    self.__interclump.setMass(self.__interclumpMass)
    return
  def __setVelocity(self):
    self.__velocity = self.__interpolations.interpolateRotationVelocity(self.__r)
    self.__velocityDispersion = 0.5*((self.__interpolations.interpolateRotationVelocity(self.__r+0.5*self.__scale)-self.__velocity)**2 + \
                                     (self.__interpolations.interpolateRotationVelocity(self.__r-0.5*self.__scale)-self.__velocity)**2)
    return
  def __setDensity(self):
    self.__density = self.__interpolations.interpolateDensity(self.__r)
    return
  def __setVelocity(self):
    self.__UVextinction = self.__interpolations.interpolateFUVextinction(self.__r)
    return
  def __setVelocity(self):
    self.__FUV = self.__interpolations.interpolateFUVfield(self.__r)
    return

  # PUBLIC
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
    self.__setVelocity()
    self.__setDensity()
    self.__setExtinction()
    self.__setFUV()
    properties = [self.__velocity, self.__velocityDispersion, self.__FUV, self.__UVextinction]
    self.__interclump.initialise(properties)
    self.__clump.initialise(properties)
    return
  def calculateEmission(self):
    iClump,tauClump,FUVclump = self.__clump.getEmission()
    iInterclump,tauInterclump,FUVinterclump = self.__interclump.getEmission()
    return np.array(iClump+iInterclump, tauClump+tauInterclump, FUVclump)