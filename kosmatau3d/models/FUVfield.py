class FUVfield(object):
  '''
  This is a class to handle the FUV field for each voxel. It will be added
  in due time.
  '''
  # PRIVATE
  def __init__(self, fuv=0):
    self.__FUV = fuv
    return
  def __setFUV(self, fuv):
    self.__FUV = fuv
    return
  def __str__(self):
    return 'FUV field: {}'.format(self.__FUV)

  # PUBLIC
  def setFUV(self, fuv):
    self.__setFUV(fuv)
    return
  def getFUV(self):
    return self.__FUV