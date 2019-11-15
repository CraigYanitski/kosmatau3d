from Dimensions import *
class Shape():
  '''
  This class defines the intrinsic shape of the PDR. It is used to modify the
  PDR structure without having to run a separate simulation (coming soon).
  '''
  # PRIVATE
  def __init__(self, x, y, z, resolution=1000, modelType=''):
    self.__type = modelType
    self.__scale = resolution
    self.__dimensions = Dimensions(x, y, z, resolution=self.__scale)
    #self.__grid = VoxelGrid(self.__dimensions.indeces())
    return

  # PUBLIC
  def getType(self):
    return self.__type
  def getDimensions(self):
    return self.__dimensions
  #def setup(self, x=0, y=0, z=0, r=0, i=0, name='box'):
