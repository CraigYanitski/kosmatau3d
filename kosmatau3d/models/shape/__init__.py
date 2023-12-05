"""
A subpackage to specify the layout of voxels in a larger three-dimensional
model.
This is left over from an early version of the code that used too many classes,
and as such it is planned to be revised and moved into a module.
"""


import importlib as il

from numba import jit_module

from .dimensions import *
from kosmatau3d.models import constants


class Shape():
    '''
    This class defines the intrinsic shape of the PDR. It is used to modify the
    PDR structure without having to run a separate simulation (coming soon).
    '''
  
    # PRIVATE
  
    def __init__(self, x, y, z, modelType=''):
        self.__type = modelType
        self.dimensions = (x, y, z)
        # self.__dimensions = Dimensions(x, y, z, resolution=self.__scale)
        # self.__x = []
        # self.__y = []
        # self.__z = []
        # self.__r = []
        # self.__phi = []
        if self.__type == 'disk':
            self.__createDisk()
        elif self.__type == 'spheroid':
            self.__createSpheroid()
        elif self.__type == 'shell':
            self.__createShell()
        else:
            self.__createBlock(x, y, z, shape=True)
        return
    
    def __createBlock(self, x, y, z, shape=False):
        '''I would like to have this function absorb the responsibilities of the Dimensions() class soon.'''
        self.__dimensions = Dimensions(x, y, z)
        x, y, z = self.__dimensions.voxelCartesianPosition()
        r, phi = self.__dimensions.voxelPolarPosition()
        if shape:
            self.__x = x
            self.__y = y
            self.__z = z
            self.__r = r
            self.__phi = phi
        return (x, y, z)
    
    def __createDisk(self):
        x, y, z = self.__createBlock(self.dimensions[0], self.dimensions[1], self.dimensions[2])
        # x, y, z = self.__dimensions.voxelCartesianPosition()
        r, phi = self.__dimensions.voxelPolarPosition()
        X = []
        Y = []
        Z = []
        R = []
        PHI = []
        for i in range(len(r)):
            if np.sqrt(x[i]**2+y[i]**2) < x.max():
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
        x, y, z = self.__createBlock(self.dimensions[0], self.dimensions[1], self.dimensions[2])
        # x, y, z = self.__dimensions.voxelCartesianPosition()
        r, phi = self.__dimensions.voxelPolarPosition()
        X = []
        Y = []
        Z = []
        R = []
        PHI = []
        for i in range(len(r)):
            if (x[i]/x.max())**2+(y[i]/y.max())**2+(z[i]/z.max())**2<1:
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
    
    def __createShell(self, radius, theta, phi):
        '''
        This function is far from finished. I would like to create an elegant method of initialising a grid of an
        appropriate size and sculpting the shell. This requires me to first solve the necessary equations.
        '''
        x, y, z = self.__dimensions.voxelCartesianPosition()
        r, phi = self.__dimensions.voxelPolarPosition()
        X = []
        Y = []
        Z = []
        R = []
        PHI = []
        for i in range(len(r)):
            if (x[i]/x.max())**2+(y[i]/y.max())**2+(z[i]/z.max())**2<1:
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

    def getShape(self):
        return self.__dimensions.getDimension()

    def getResolution(self):
        return self.__scale
    
    def voxelNumber(self):
        return self.__r.size
    
    #def setup(self, x=0, y=0, z=0, r=0, i=0, name='box'):
    def voxelCartesianPositions(self):
        # Return the Cartesian coordinates of each voxel
        return (self.__x, self.__y, self.__z)
    
    def voxelPolarPositions(self):
        # Return the polar coordinates of each voxel
        return (self.__r, self.__phi)
  
  # jit_module(nopython=False)
