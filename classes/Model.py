import importlib as il
import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D

import shape #import Shape
import statistics
import observations
import radiativeTransfer
import species

from VoxelGrid import *

class Model(object):
  '''
  This is the highest class in the hierarchy of the KOSMA-tau^3 simulation.
  It contains all of the information needed to properly model a PDR (I think).
  '''
  # PRIVATE

  def __init__(self, x, y, z, modelName='', modelType='', resolution=1000, verbose=False):
    constants.type = modelType   #this just adds a label to the type of model being created. ie 'disk', 'bar', 'sphere', etc.
    constants.resolution = float(resolution)
    constants.change.changeDirectory(modelName)
    observations.methods.initialise()
    self.__shape = shape.Shape(x, y, z, modelType=modelType)      #Shape() object to create the parameters for the grid of voxels
    self.__grid = VoxelGrid(self.__shape)   #VoxelGrid() object to build the model and calculate the emission
    #self.__orientation = Orientation(self.__shape.getDimensions())      #Orientation() object to change the viewing angle and expected spectra
    # self.__molecules = Molecules()   #Molecules() object to centralise the molecules in model
    # self.__dust = Dust()             #Dust() object to centralise the dust in the model
    #self.__species = [self.__molecules, self.__dust]    #this is a list of the species objects being considered
    self.__speciesNames = []          #this is a list of the species names for easy printout
    self.__verbose = verbose
    self.__intensityMap = []
    self.__mapPositions = []
    return

  def __str__(self):
    printout = 'A {} model of {} voxels'.format(constants.type, self.getGrid().getVoxelNumber())
    if self.__verbose:
      printout += '\n  arranged in {}'.format(self.__shape.getDimensions())
      printout += '\n\nConsidering {} species:\n{}\n{}'.format(len(self.speciesNames), self.__molecules, self.__dust)
    emission = self.__grid.totalEmission()
    printout += '\n\nTotal intensity: {}\nTotal optical depth: {}'.format(emission[0].sum(), np.log(np.exp(emission[1]).sum()))
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

  def calculateModel(self):
    self.__grid.calculateVoxels()
    return

  def addSpecies(self, speciesTransition):
    for element in speciesTransition:
      molecule = element.split(' ')[0]
      transition = int(element.split(' ')[1])
      species.addMolecule(molecule, transition)

    species.moleculeWavelengths = species.molecules.getWavelengths()
    constants.resortWavelengths()
    return

  def writeEmission(self):
    self.__grid.writeEmission()
    return

  def setLOS(self, x=0, y=0, z=0, dim='xy'):
    self.__orientation.setLOS(self.__grid, x=x, y=y, z=z, dim=dim)
    return

  def calculateObservation(self, velocity=[0], dim='xy'):
    xPositions,yPositions,zPositions = self.__grid.getVoxelPositions()
    #print('\nx\n', np.unique(xArray), '\ny\n', np.unique(yArray), '\nz\n', np.unique(zArray))
    position = []
    intensityMap = []
    if dim=='xy':
      Array = np.unique([xPositions,yPositions], axis=1).T
      if self.__verbose:
        print('\nx\n{}\n\ny\n{}\n'.format(Array[:,0],Array[:,1]))
      for x,y in Array:
        #for y in np.unique(yArray):
          radiativeTransfer.orientation.setLOS(self.__grid, x=x, y=y, dim=dim)
          position.append([x,y])
          intensity = radiativeTransfer.orientation.calculateRadiativeTransfer(velocity)
          intensityMap.append(intensity)
    if dim=='xz':
      xArray,zArray = np.unique([xPositions,zPositions], axis=1)
      for x in np.unique(xArray):
        for z in np.unique(zArray):
          radiativeTransfer.orientation.setLOS(self.__grid, x=x, z=z, dim=dim)
          position.append([x,z])
          intensity = radiativeTransfer.orientation.calculateRadiativeTransfer(velocity)
          intensityMap.append(intensity)
    if dim=='yz':
      yArray,zArray = np.unique([yPositions,zPositions], axis=1)
      for y in np.unique(yArray):
        for z in np.unique(zArray):
          radiativeTransfer.orientation.setLOS(self.__grid, y=y, z=z, dim=dim)
          position.append([y,z])
          intensity = radiativeTransfer.orientation.calculateRadiativeTransfer(velocity)
          intensityMap.append(intensity)
    self.__mapPositions = np.array(position)
    self.__intensityMap = np.array(intensityMap)
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
        print('{}: {} centered at {} km/s'.format(self.__speciesNames[element], intensity[element][intensity[element].nonzero()], self.__constants.velocityBins[i]))
      print()
    else:
      for index in range(len(self.__mapPositions)):
        position = self.__mapPositions[index]
        intensity = self.__intensityMap[index]
        print('At position x={} y={}, the intensity is'.format(position[0], position[1]))
        for element in range(intensity.shape[0]):
          i = intensity[element].argmax()
          print('{}: {} centered at {} km/s'.format(self.__speciesNames[element], intensity[element][intensity[element].nonzero()], self.__constants.velocityBins[i]))
        print()
    return
    
  def plotModel(self, plot='total intensity', debug=False):
    positions = self.__grid.getVoxelPositions()
    limits = [positions.min(), positions.max()]
    if plot=='total intensity':
      if debug: print(self.__grid.totalEmission().shape)
      weights = (self.__grid.totalEmission()[0]).max(2).sum(1)
      plot = r'$I \ (\chi)$'
    elif plot=='total optical depth':
      weights = (self.__grid.totalEmission()[1]).max(2).sum(1)
      plot = r'$\tau$'
    elif plot=='clump intensity':
      weights = (self.__grid.clumpEmission()[0]).max(2).sum(1)
      plot = r'$I \ (\chi)$'
    elif plot=='clump optical depth':
      weights = (self.__grid.clumpEmission()[1]).max(2).sum(1)
      plot = r'$\tau$'
    elif plot=='interclump intensity':
      weights = (self.__grid.interclumpEmission()[0]).max(2).sum(1)
      plot = r'$I \ (\chi)$'
    elif plot=='interclump optical depth':
      weights = (self.__grid.interclumpEmission()[1]).max(2).sum(1)
      plot = r'$\tau$'
    elif plot=='FUV':
      weights = (self.__grid.getFUV())
      plot = r'$FUV \ (\chi)$'
    elif plot=='Afuv':
      weights = (self.__grid.getAfuv())
      plot = r'$\tau_{FUV}$'
    elif plot=='velocity':
      weights = (self.__grid.getVelocity())
      plot = r'$v_{rot} \ (\frac{km}{s})$'
    fig = plt.figure()
    ax = fig.add_subplot(111, projection='3d')
    model = ax.scatter(positions[0], positions[1], positions[2], c=weights, cmap=plt.cm.hot, marker='s', s=27, alpha=1, linewidths=0)
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