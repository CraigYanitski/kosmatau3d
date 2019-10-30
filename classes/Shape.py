class Shape():
  '''
  This class defines the intrinsic shape of the PDR. It is used to modify the
  PDR structure without having to run a separate simulation (coming soon).
  '''
  # PRIVATE
  def __init__(self, x, y, z, type=''):
    self.__type = type
    self.__dimensions = Dimensions(x, y, z)
    #self.__grid = VoxelGrid(self.__dimensions.indeces())
    return

  # PUBLIC
  #def setup(self, x=0, y=0, z=0, r=0, i=0, name='box'):
