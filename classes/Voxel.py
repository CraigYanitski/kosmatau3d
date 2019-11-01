class Voxel():
  '''
  This is a class to handle each voxel in KOSMA-tau^3. It contains ensembles
  of spherical PDR simulations (to mimic the fractal structure of PDRs), a FUV
  field element, and element absorption, separated for dust and molecules. It
  should have one clump ensemble and one interclump ensemble. This is to account
  for the diffuse surroundings of a clump.
  '''
  # PRIVATE
  def __init__(self, index):
    self.__index = index    #index of voxel in VoxelGrid, sort of like its ID in the overall model
    self.__velocity = 0     #velocity of mass at voxel point
    self.__velocityDispersion = 0   #dispersion of velocity at voxel point
    self.__intensity = 0      #intensity of emissions at voxel point
    self.__opticalDepth = 0   #optical depth at voxel point
    self.__clump = Ensemble('clump')    #clumpy ensemble at voxel point
    self.__interclump = Ensemble('interclump')    #diffuse interclump ensemble at voxel point
    self.__x = 0
    self.__y = 0
    self.__z = 0
    self.__r = 0
    self.__phi = 0
    return
  def __setClumpMass(self):
    self.__clumpMass = interpolations.interpolateClumpMass(self.__r)
    self.__clump.setMass(self.__clumpMass)
    return
  def __setInterclumpMass(self):
    self.__interclumpMass = interpolations.interpolateInterclumpMass(self.__r)
    self.__interclump.setMass(self.__interclumpMass)
    return
  def __setVelocity(self):
    self.__velocity = interpolations.interpolateRotationVelocity(self.__r)
    self.__velocityDispersion = 
    return
  def __setDensity(self):
    self.__density = interpolations.interpolateDensity(self.__r)
    return
  def __setVelocity(self):
    self.__UVextinction = interpolations.interpolateFUVextinction(self.__r)
    return
  def __setVelocity(self):
    self.__FUV = interpolations.interpolateFUVfield(self.__r)
    return

  # PUBLIC
  def setPosition(self, x, y, z, r, phi):
    self.__x = x
    self.__y = y
    self.__z = z
    self.__r = r
    self.__phi = phi
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
