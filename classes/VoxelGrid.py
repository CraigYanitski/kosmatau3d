import numpy as np
from tqdm import tqdm
#import progressbar as pb
import importlib as il
import gc
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
    self.__x = []
    self.__y = []
    self.__z = []
    return
  def __initialiseGrid(self, species, observations):
    self.__species = species
    self.__interpolations = Interpolate(self.__species, observations)
    for i in range(self.__voxelNumber): self.__voxels.append(Voxel(self.__species, self.__interpolations, i))
  def __str__(self):
    return 'VoxelGrid\n  ->{} voxels\n  ->intensity {}\n  ->optical depth {}'.format(self.__voxelNumber, sum(self.__voxelIntensity), -np.log(np.exp(-self.__voxelOpticalDepth)))

  # PUBLIC
  #def createGrid(self, indeces):
  #  for i in indeces: self.__voxels.append(Voxel(i))
  #  return
  def reloadModules(self):
    il.reload(Voxel)
    il.reload(Interpolate)
    for voxel in self.__grid:
      voxel.reloadModules()
    self.__interpolations.reloadModules()
    return
  def getDimensions(self):
    return self.__dimensions
  def getInterpolations(self):
    return self.__interpolations
  def initialiseVoxels(self, species, observations, verbose=False):
    self.__initialiseGrid(species, observations)
    print('\nInitialising Grid...')
    x,y,z,scale = self.__dimensions.voxelCartesianPosition()
    r,phi = self.__dimensions.voxelPolarPosition()
    self.__unusedVoxels = []
    with tqdm(total=len(self.__voxels), desc='Voxels initialised', miniters=1, dynamic_ncols=True) as progress:
      for i,voxel in enumerate(self.__voxels):
        if r[i]<=max(x):
          if verbose:
              print('\nMax X, Radius:', max(x), r[i], '\n')
          self.__x.append(x[i])
          self.__y.append(y[i])
          self.__z.append(z[i])
          voxel.setIndex(i-len(self.__unusedVoxels))
          voxel.setPosition(x[i], y[i], z[i], r[i], phi, scale)
          voxel.setProperties()
        else: self.__unusedVoxels.append(i)
        progress.update()
      progress.close()
    for i in self.__unusedVoxels[::-1]:
      self.__voxels.remove(self.__voxels[i])
    self.__voxelNumber = len(self.__voxels)
    return
  def calculateEmission(self, verbose=False):
    print('\nCalculating grid emission...')
    with tqdm(total=len(self.__voxels), desc='Voxel emissions', miniters=1, dynamic_ncols=True) as progress:
      for i,voxel in enumerate(self.__voxels):
        voxel.calculateEmission()
        gc.collect()
        emission = voxel.getEmission()
        if verbose: print(emission)
        progress.update()
        self.__voxelIntensity.append(emission[0])
        self.__voxelOpticalDepth.append(emission[1])
        self.__voxelFUV.append(emission[2])
      progress.close()
      self.__voxelIntensity = np.array(self.__voxelIntensity)
      self.__voxelOpticalDepth = np.array(self.__voxelOpticalDepth)
    print('\nCalculation complete.\n')
    del emission
    del progress
    del voxel
    return
  def getVoxelNumber(self):
    return self.__voxelNumber
  def getVoxelPositions(self):
    return np.array([self.__x, self.__y, self.__z])
  def allVoxels(self):
    # Just in case all of the Voxel() instances need to be retrieved
    return self.__voxels
  def totalEmission(self):
    # Return the emission from all of the voxels, separated by observed velocity
    return np.array([self.__voxelIntensity,self.__voxelOpticalDepth])
  def printVoxels(self):
    for voxel in self.__voxels: voxel.printVoxel()
    return