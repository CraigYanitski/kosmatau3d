class VoxelGrid(Model):
  '''
  This is a class to handle all of the voxels in KOSMA-tau^3. It contains a
  specified arrangement of voxels, and must coordinate with the Dimensions class
  to make the Shape class functional.
  '''
  # PRIVATE
  def __init__(self, species, indeces=1):
    self.__voxelNumber = indeces
    self.__voxels = []
    for i in indeces: self.__voxels.append(Voxel(i))
    self.__map = {}       #dictionary object to map the voxel indeces to the correct location
    self.__species = species
    self.__voxelIntensity = []
    self.__voxelOpticalDepth = []
    self.__interpolations = Interpolate(self.__species)
    return

  # PUBLIC
  #def createGrid(self, indeces):
  #  for i in indeces: self.__voxels.append(Voxel(i))
  #  return
  def initialiseVoxels(self):
    x,y,z,scale = dimensions.voxelCartesianPosition()
    r,phi = dimensions.voxelPolarPoasition()
    for i,voxel in enumerate(self.__voxels):
      voxel.setPosition(x, y, z, r, phi, scale)
      voxel.setProperties(self.__interpolations)
    return
  def calculateEmission(self):
    for i,voxel in enumerate(self.__voxelNumber):
      emission = voxel.calculateEmission()
      self.__voxelIntensity.append(emission[0])
      self.__voxelOpticalDepth.append(emission[1])
    return
  def allVoxels(self):
    # Just in case all of the Voxel() instances need to be retrieved
    return self.__voxels
  def totalEmission(self):
    # Return the emission from all of the voxels, separated by observed velocity
    return np.array(self.__voxelIntensity,self.__voxelOpticalDepth)
