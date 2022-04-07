import numpy as np

from kosmatau3d.models import constants


class Dimensions(object):
    '''
    This is a class to contain the dimensions of VoxelGrid(). It will
    coordinate with class VoxelGrid to arrange the voxels. the input dimensions must be
    given in kpc.
  
    At the moment this will just work in the plane of a disk.
    '''
    # PRIVATE
    def __init__(self, x, y, z, i=0):
        self.__i = i
        self.__x = np.floor_divide(x*1000, constants.voxel_size)*constants.voxel_size
        self.__y = np.floor_divide(y*1000, constants.voxel_size)*constants.voxel_size
        self.__z = np.floor_divide(z*1000, constants.voxel_size)*constants.voxel_size
        self.dimensions = (self.__x, self.__y, self.__z)
        self.__xOffset = self.__x/2.
        self.__yOffset = self.__y/2.
        self.__zOffset = self.__z/2.
        self.__xRange = np.arange(0, self.__x+constants.voxel_size, constants.voxel_size)
        self.__yRange = np.arange(0, self.__y+constants.voxel_size, constants.voxel_size)
        self.__zRange = np.arange(0, self.__z+constants.voxel_size, constants.voxel_size)
        grid = np.meshgrid(self.__xRange-self.__xOffset, self.__yRange-self.__yOffset, self.__zRange-self.__zOffset)
        self.__xPositions = grid[0].flatten()
        self.__yPositions = grid[1].flatten()
        self.__zPositions = grid[2].flatten()
        self.__r, self.__phi = self.__convertToPolar()
        # self.__h = self.hCalc(self.__x,self.__y,self.__z,self.__i)
        return
    
    def __convertToPolar(self):
        # This is a function that calculates the polar coordinates of each voxel
        r = ((self.__xPositions)**2 + (self.__yPositions)**2)**0.5
        phi = np.arctan2(self.__yPositions, self.__xPositions)
        return r, phi
    
    def __str__(self):
        return 'dimensions {}pc x {}pc x {}pc'.format(self.__x, self.__y, self.__z)
  
    # PUBLIC
    def getDimension(self):
        return (self.__x, self.__y, self.__z)

    def voxelNumber(self):
        # Return the number of voxels required
        return len(self.__xPositions)
    
    def voxelCartesianPosition(self):
        # Return the Cartesian coordinates of each voxel
        return (self.__xPositions, self.__yPositions, self.__zPositions)
    
    def voxelPolarPosition(self):
        # Return the polar coordinates of each voxel
        return (self.__r, self.__phi)
    
    def hCalc(self):
        #calculates height h from koordinates relative to disk-plain
        # This will soon be rewritten.
        h = (self.__yPositions - np.tan(self.__i)*self.__zPositions) * np.cos(self.__i)
        return h
