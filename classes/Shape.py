import importlib as il
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
    self.__x = []
    self.__y = []
    self.__z = []
    self.__r = []
    self.__phi = []
    if self.__type=='disk': createDisk()
    if self.__type=='spheroid': createSpheroid()
    return

  # PUBLIC
  def reloadModules(self):
    il.reload(Dimensions)
    return
  def getType(self):
    return self.__type
  def getDimensions(self):
    return self.__dimensions
  def getResolution(self):
    return self.__scale
  def voxelNumber(self):
    return self.__r.size
  #def setup(self, x=0, y=0, z=0, r=0, i=0, name='box'):

  def createDisk(self):
    x,y,z = self.__dimensions.voxelCartesianPosition()
    r,phi = self.__dimensions.voxelPolarPosition()
    X = []
    Y = []
    Z = []
    R = []
    PHI = []
    for i in range(len(r)):
      if np.sqrt(x[i]**2+y[i]**2)<x.max:
        X.append(x[i])
        Y.append(y[i])
        Z.append(z[i])
        R.append(r[i])
        PHI.append(phi[i])
    self.__x = np.array(X)
    self.__y = np.array(Y)
    self.__z = np.array(Z)
    self.__r = np.array(R)
    self.__phi = np.array(PHI)
    return

  def createSpheroid(self):
    x,y,z = self.__dimensions.voxelCartesianPosition()
    r,phi = self.__dimensions.voxelPolarPosition()
    X = []
    Y = []
    Z = []
    R = []
    PHI = []
    for i in range(len(r)):
      if (x[i]/x.max)**2+(y[i]/y.max)**2+(z[i]/z.max)**2<1:
        X.append(x[i])
        Y.append(y[i])
        Z.append(z[i])
        R.append(r[i])
        PHI.append(phi[i])
    self.__x = np.array(X)
    self.__y = np.array(Y)
    self.__z = np.array(Z)
    self.__r = np.array(R)
    self.__phi = np.array(PHI)
    return