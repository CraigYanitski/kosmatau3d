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
    self.__index = index    #index of voxel in VoxelGrid
    self.__velocity = 0     #velocity of mass at voxel point
    self.__velocityDispersion = 0   #dispersion of velocity at voxel point
    self.__intensity = 0      #intensity of emissions at voxel point
    self.__opticalDepth = 0   #optical depth at voxel point
    self.__clump = Ensemble('clump')    #clumpy ensemble at voxel point
    self.__interclump = Ensemble('interclump')    #diffuse interclump ensemble at voxel point
    self.__x = 0
    self.__y = 0
    self.__z = 0
    return

  # PUBLIC
  def initialiseVelocity(self, size=10):
    super()._Model
    return
  def createClump(self, ):
    self.__clump.initialise()
    return
