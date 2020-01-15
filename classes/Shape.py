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
    if self.__type=='disk': self.__createDisk()
    if self.__type=='spheroid': self.__createSpheroid()
    return
  def __createDisk(self):
    x,y,z = self.__dimensions.voxelCartesianPosition()
    r,phi = self.__dimensions.voxelPolarPosition()
    X = []
    Y = []
    Z = []
    R = []
    PHI = []
    print(x.max)
    for i in range(len(r)):
      if np.sqrt(x[i]**2+y[i]**2)<max(x):
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
  def __createSpheroid(self):
    x,y,z = self.__dimensions.voxelCartesianPosition()
    r,phi = self.__dimensions.voxelPolarPosition()
    X = []
    Y = []
    Z = []
    R = []
    PHI = []
    for i in range(len(r)):
      if (x[i]/(x))**2+(y[i]/(y))**2+(z[i]/(z))**2<1:
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
  def voxelCartesianPositions(self):
    # Return the Cartesian coordinates of each voxel
    return (self.__x, self.__y, self.__z, self.__scale)
  def voxelPolarPositions(self):
    # Return the polar coordinates of each voxel
    return (self.__r, self.__phi)
