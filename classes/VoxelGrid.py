import numpy as np
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
    for i in range(self.__voxelNumber): self.__voxels.append(Voxel(i))
  def __str__(self):
    return 'VoxelGrid\n  -{} voxels'.format(self.__voxelNumber)

  # PUBLIC
  #def createGrid(self, indeces):
  #  for i in indeces: self.__voxels.append(Voxel(i))
  #  return
  def getDimensions(self):
    return self.__dimensions
  def initialiseVoxels(self, species, observations):
    self.__initialiseGrid(species, observations)
    x,y,z,scale = self.__dimensions.voxelCartesianPosition()
    r,phi = self.__dimensions.voxelPolarPosition()
    for i,voxel in enumerate(self.__voxels):
      voxel.setPosition(x, y, z, r, phi, scale)
      voxel.setProperties(self.__interpolations)
    return
  def calculateEmission(self):
    for i,voxel in enumerate(self.__voxelNumber):
      emission = voxel.calculateEmission()
      self.__voxelIntensity.append(emission[0])
      self.__voxelOpticalDepth.append(emission[1])
      self.__voxelFUV.append(emission[2])
    return
  def allVoxels(self):
    # Just in case all of the Voxel() instances need to be retrieved
    return self.__voxels
  def totalEmission(self):
    # Return the emission from all of the voxels, separated by observed velocity
    return np.array(self.__voxelIntensity,self.__voxelOpticalDepth)
