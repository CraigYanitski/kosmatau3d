import numpy as np
from tqdm import tqdm
from Voxel import *
from Interpolate import *
class VoxelGrid(object):
  '''
  This is a class to handle all of the voxels in KOSMA-tau^3. It contains a
  specified arrangement of voxels, and must coordinate with the Dimensions class
  to make the Shape class functional.
  '''
  # PRIVATE
  def __init__(self, dimensions):
    self.__dimensions = dimensions
    self.__voxelNumber = self.__dimensions.voxelNumber()
    self.__voxels = []
    self.__map = {}       #dictionary object to map the voxel indeces to the correct location
    self.__species = None
    self.__voxelIntensity = []
    self.__voxelOpticalDepth = []
    self.__voxelFUV = []
    self.__interpolations = None
    return
  def __initialiseGrid(self, species, observations):
    self.__species = species
    self.__interpolations = Interpolate(self.__species, observations)
    for i in range(self.__voxelNumber): self.__voxels.append(Voxel(self.__species, self.__interpolations, i))
  def __str__(self):
    return 'VoxelGrid\n  ->{} voxels\n  ->intensity {}\n  ->optical depth {}'.format(self.__voxelNumber, sum(self.__voxelIntensity), sum(self.__voxelOpticalDepth))

  # PUBLIC
  #def createGrid(self, indeces):
  #  for i in indeces: self.__voxels.append(Voxel(i))
  #  return
  def getDimensions(self):
    return self.__dimensions
  def getInterpolations(self):
    return self.__interpolations
  def initialiseVoxels(self, species, observations, verbose=False):
    self.__initialiseGrid(species, observations)
    print('Initialising Grid...')
    x,y,z,scale = self.__dimensions.voxelCartesianPosition()
    r,phi = self.__dimensions.voxelPolarPosition()
    self.__unusedVoxels = []
    with tqdm(total=len(self.__voxels), desc='Voxels initialised') as progress:
      for i,voxel in enumerate(self.__voxels):
        if r[i]<=2*max(x):
          if verbose: print(r[i])
          voxel.setIndex(i-len(self.__unusedVoxels))
          voxel.setPosition(x[i], y[i], z[i], r[i], phi, scale)
          voxel.setProperties()
        else: self.__unusedVoxels.append(i)
        progress.update(1)
    for i in self.__unusedVoxels[::-1]:
      self.__voxels.remove(self.__voxels[i])
    self.__voxelNumber = len(self.__voxels)
    return
  def calculateEmission(self, verbose=False):
    print('Calculating grid emission')
    with tqdm(total=len(self.__voxels), desc='Calculating voxel emissions') as progress:
      for i,voxel in enumerate(self.__voxels):
        voxel.calculateEmission()
        emission = voxel.getEmission()
        if verbose: print(emission)
        self.__voxelIntensity.append(emission[0])
        self.__voxelOpticalDepth.append(emission[1])
        self.__voxelFUV.append(emission[2].getFUV())
        progress.update(1)
    return
  def getVoxelNumber(self):
    return self.__voxelNumber
  def allVoxels(self):
    # Just in case all of the Voxel() instances need to be retrieved
    return self.__voxels
  def totalEmission(self):
    # Return the emission from all of the voxels, separated by observed velocity
    return np.array([self.__voxelIntensity,self.__voxelOpticalDepth])
