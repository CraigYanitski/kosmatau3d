import importlib as il
import sys

import numpy as np
import matplotlib.pyplot as plt
# from mpl_toolkits.mplot3d import Axes3D

from kosmatau3d.models import shape  # import Shape
# import statistics
from kosmatau3d.models import observations
from kosmatau3d.models import species
from .VoxelGrid import *


class Model(object):
    '''
    This is the highest class in the hierarchy of the KOSMA-tau^3 simulation.
    It contains all of the information needed to properly model a PDR (I think).
    '''
    # PRIVATE
  
    def __init__(self, history_path='', directory='', x=0, y=0, z=0, modelType='', resolution=1000, molecules=[],
                 dust='', clumpMassRange=[], clumpMassNumber=[], clumpNmax=[], velocityRange=[], velocityNumber=0,
                 clumpMassFactor=[], FUVfactor=1, densityFactor=1, globalUV=10, verbose=False):
      
        if not len(clumpMassRange):
            sys.exit('<<ERROR>> Define mass sets in argument.')
            # sys.exit()
        if not len(velocityRange):
            sys.exit('<<ERROR>> Define observing velocities in argument.')
            # sys.exit()
        
        # this just adds a label to the type of model being created. ie 'disk', 'bar', 'sphere', etc.
        constants.type = modelType
        constants.voxel_size = float(resolution)
        constants.HISTORYPATH = history_path
        constants.changeDirectory(directory)
        
        # Factors
        constants.clumpMassFactor = clumpMassFactor
        constants.FUVFactor = FUVfactor
        constants.densityFactor = densityFactor
        constants.globalUV = globalUV
    
        constants.changeVelocityRange(velocityRange)
        constants.changeVelocityNumber(velocityNumber)
        constants.addClumps(massRange=clumpMassRange, num=clumpMassNumber, Nmax=clumpNmax, reset=True)
        
        observations.methods.initialise()
        self.__addSpecies(molecules)
        constants.changeDustWavelengths(dust)

        # Shape() object to create the parameters for the grid of voxels
        self.__shape = shape.Shape(x, y, z, modelType=modelType)
        # VoxelGrid() object to build the model and calculate the emission
        self.__grid = VoxelGrid(self.__shape)
        # Orientation() object to change the viewing angle and expected spectra
        # self.__orientation = Orientation(self.__shape.getDimensions())
        self.__speciesNames = []  # this is a list of the species names for easy printout
        self.__verbose = verbose
        self.__intensityMap = []
        self.__mapPositions = []
        return
  
    def __addSpecies(self, speciesTransition):
        species.addMolecules(speciesTransition)
        return
  
    def __str__(self):
        printout = 'A {} model of {} voxels'.format(constants.type, self.getGrid().getVoxelNumber())
        if self.__verbose:
            printout += '\n  arranged in {}'.format(self.__shape.getDimensions())
            printout += '\n\nConsidering {} species:\n{}\n{}'.format(len(self.speciesNames),
                                                                     self.__molecules,
                                                                     self.__dust)
        emission = self.__grid.totalEmission()
        printout += '\n\nTotal intensity: {}\nTotal optical depth: {}'.format(emission[0].sum(),
                                                                              np.log(np.exp(emission[1]).sum()))
        return printout
  
    # PUBLIC
  
    def getType(self):
        return constants.type
  
    def getShape(self):
        return self.__shape
  
    def getGrid(self):
        return self.__grid
  
    def getOrientation(self):
        return self.__orientation
  
    # def getObservations(self):
    #   return self.__observations
  
    def getSpecies(self):
        return species.speciesNames
  
    def getSpeciesNames(self):
        return species.speciesNames
  
    # def reloadModules(self):
    #   il.reload(Shape)
    #   il.reload(VoxelGrid)
    #   il.reload(Orientation)
    #   il.reload(Observations)
    #   il.reload(Molecules)
    #   il.reload(Dust)
    #   self.__shape.reloadModules()
    #   self.__grid.reloadModules()
    #   #self.__observations.reloadModules()
    #   #self.__orientation.reloadModules()
    #   self.__molecules.reloadModules()
    #   self.__dust.reloadModules()
    #   return
  
    def calculateModel(self, **kwargs):
        self.__grid.calculateEmission(**kwargs)
        return
  
    def writeEmission(self, debug=False):
        self.__grid.writeEmission(debug=debug)
        return
  
    def getIntensityMap(self):
        return (self.__mapPositions, self.__intensityMap)
  
    def printIntensityMap(self, index=None):
  
        print(self.__species[0].getMolecules(), self.__species[1].getDust())
    
        if not index==None:
    
            position = self.__mapPositions[index]
            intensity = self.__intensityMap[index]
      
            print('At position x={} y={}, the intensity is'.format(position[0], position[1]))
            
            for element in range(intensity.shape[0]):
      
                i = intensity[element].argmax()
                print('{}: {} centered at {} km/s'.format(self.__speciesNames[element],
                                                          intensity[element][intensity[element].nonzero()],
                                                          self.__constants.velocityBins[i]))
            
            print()
    
        else:
          
            for index in range(len(self.__mapPositions)):
            
                position = self.__mapPositions[index]
                intensity = self.__intensityMap[index]
              
                print('At position x={} y={}, the intensity is'.format(position[0], position[1]))
              
                for element in range(intensity.shape[0]):
              
                    i = intensity[element].argmax()
                    print('{}: {} centered at {} km/s'.format(self.__speciesNames[element],
                                                              intensity[element][intensity[element].nonzero()],
                                                              self.__constants.velocityBins[i]))
              
                print()
        
        return
      
    def plotModel(self, plot='total intensity', debug=False):
        positions = self.__grid.getVoxelPositions()
        limits = [positions.min(), positions.max()]
        if plot == 'total intensity':
            if debug:
                print(self.__grid.totalEmission().shape)
            weights = (self.__grid.totalEmission()[0]).max(2).sum(1)
            plot = r'$I \ (\chi)$'
        elif plot == 'total optical depth':
            weights = (self.__grid.totalEmission()[1]).max(2).sum(1)
            plot = r'$\tau$'
        elif plot == 'clump intensity':
            weights = (self.__grid.clumpEmission()[0]).max(2).sum(1)
            plot = r'$I \ (\chi)$'
        elif plot == 'clump optical depth':
            weights = (self.__grid.clumpEmission()[1]).max(2).sum(1)
            plot = r'$\tau$'
        elif plot == 'interclump intensity':
            weights = (self.__grid.interclumpEmission()[0]).max(2).sum(1)
            plot = r'$I \ (\chi)$'
        elif plot == 'interclump optical depth':
            weights = (self.__grid.interclumpEmission()[1]).max(2).sum(1)
            plot = r'$\tau$'
        elif plot == 'FUV':
            weights = (self.__grid.getFUV())
            plot = r'$FUV \ (\chi)$'
        elif plot == 'Afuv':
            weights = (self.__grid.getAfuv())
            plot = r'$\tau_{FUV}$'
        elif plot == 'velocity':
            weights = (self.__grid.getVelocity())
            plot = r'$v_{rot} \ (\frac{km}{s})$'
        fig = plt.figure()
        ax = fig.add_subplot(111, projection='3d')
        model = ax.scatter(positions[0], positions[1], positions[2],
                           c=weights, cmap=plt.cm.hot, marker='s', s=27, alpha=1, linewidths=0)
        ax.set_xlim(limits)
        ax.set_ylim(limits)
        ax.set_zlim(limits)
        cbar = plt.colorbar(model)
        ax.set_title('PDR Emission within the Milky Way')
        ax.set_xlabel('X (pc)')
        ax.set_ylabel('Y (pc)')
        ax.set_zlabel('Z (pc)')
        cbar.set_label(plot, rotation=0)
        plt.show()
        return
  
    def hdu_header(self, name='', units='', filename=None, dim=None, data=None):
        if filename is None:
            return
    
        header = fits.Header()
        
        header['SIMPLE'] = (True, 'conforms to FITS standard')
        header['BITPIX'] = (-64, 'element size')
        header['NAXIS'] = (len(dim), 'number of axes')
        header['EXTEND'] = True
        
        for i in range(len(dim)):
            header['NAXIS{}'.format(i+1)] = dim[i]
        
        # header['NAXIS1'] = self.__constants.velocityBins.size#emission[0].data[0,0,0,:].size
        # header['NAXIS2'] = len(self.__species[0].getMolecules())+len(self.__species[1].getDust())
        # header['NAXIS3'] = 2#emission[0].data[:,0,0,0].size
        # header['NAXIS4'] = self.__voxelNumber#emission[0].data[0,:,0,0].size
      
        header['NAME'] = name
        header['UNITS'] = units
        
        if '.fits' not in filename:
            filename = constants.HISTORYPATH + constants.directory + 'r{}_n{}/'.format(constants.resolution,
                                                                                       self.__shape.voxelNumber()) + \
                                               filename + '.fits'
        else:
            filename = constants.HISTORYPATH + constants.directory + filename
        
        if os.path.exists(filename):
            os.remove(filename)
        
        hdu = fits.PrimaryHDU(data=data, header=header)
        hdu.writeto(filename, overwrite=True)
        
        return
