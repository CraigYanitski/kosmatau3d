import importlib as il
import numpy as np
from Shape import *
from VoxelGrid import *
from Orientation import *
from Observations import *
from Molecules import *
from Dust import *
class Model(object):
  '''
  This is the highest class in the hierarchy of the KOSMA-tau^3 simulation.
  It contains all of the information needed to properly model a PDR (I think).
  '''
  # PRIVATE
  def __init__(self, x, y, z, modelType='', resolution=1000, verbose=False):
    self.__type = modelType   #this just adds a label to the type of model being created. ie 'disk', 'bar', 'sphere', etc.
    self.__scale = float(resolution)
    self.__constants = Constants()
    self.__shape = Shape(x, y, z, modelType=modelType, resolution=self.__scale)      #Shape() object to create the parameters for the grid of voxels
    self.__grid = VoxelGrid(self.__shape)   #VoxelGrid() object to build the model and calculate the emission
    self.__observations = Observations(self.__shape.getResolution())    #Observations() object to centralise the required data for the program
    self.__orientation = Orientation(self.__shape.getDimensions(), self.__observations)      #Orientation() object to change the viewing angle and expected spectra
    self.__molecules = Molecules()   #Molecules() object to centralise the molecules in model
    self.__dust = Dust()             #Dust() object to centralise the dust in the model
    self.__species = [self.__molecules, self.__dust]    #this is a list of the species objects being considered
    self.__speciesNames = []          #this is a list of the species names for easy printout
    self.__verbose = verbose
    self.__intensityMap = []
    self.__mapPositions = []
    return
  def __str__(self):
    printout = 'A {} model of {} voxels'.format(self.__type, self.getGrid().getVoxelNumber())
    if self.__verbose:
      printout += '\n  arranged in {}'.format(self.__shape.getDimensions())
      printout += '\n\nConsidering {} species:\n{}\n{}'.format(len(self.speciesNames), self.__molecules, self.__dust)
    emission = self.__grid.totalEmission()
    printout += '\n\nTotal intensity: {}\nTotal optical depth: {}'.format(emission[0].sum(), np.log(np.exp(emission[1]).sum()))
    return printout

  # PUBLIC
  def getType(self):
    return self.__type
  def getShape(self):
    return self.__shape
  def getGrid(self):
    return self.__grid
  def getOrientation(self):
    return self.__orientation
  def getObservations(self):
    return self.__observations
  def getSpecies(self):
    return self.__species
  def getSpeciesNames(self):
    return self.__molecules.getMolecules() + self.__dust.getDust()
  def reloadModules(self):
    il.reload(Shape)
    il.reload(VoxelGrid)
    il.reload(Orientation)
    il.reload(Observations)
    il.reload(Molecules)
    il.reload(Dust)
    self.__shape.reloadModules()
    self.__grid.reloadModules()
    #self.__observations.reloadModules()
    #self.__orientation.reloadModules()
    self.__molecules.reloadModules()
    self.__dust.reloadModules()
    return
  def initialiseModel(self):
    self.__grid.initialiseVoxels(self.__species, self.__observations)
    return
  def addDust(self, dust, transition):
    (numbers,species,transitions,frequencies) = self.__observations.speciesData
    i = (species==dust)&(transitions==transition)
    if transition in self.__dust.getTransitions():
      self.__dust.addTransition(dust, transition, frequencies[i], numbers[i])
      # self.__dustNumber.append(numbers[species=='dust' and transitions==transition])
      # self.__dustTransitions['dust'].append(transition)
      # self.__dustFrequencies['dust'].append(frequencies[species=='dust' and transitions==transition])
    else:
      self.__dust.addDust(dust, transition, frequencies[i], numbers[i])
      # self.__dustNames.append('dust')
      # self.__dustNumber.append(numbers[species=='dust' and transitions==transition])
      # self.__dustTransitions['dust'].append(transition)
      # self.__dustFrequencies['dust'].append(frequencies[species=='dust' and transitions==transition])
    #self.__species = [self.__molecules, self.__dust]
    self.__speciesNames = np.append(self.__molecules.getMolecules(), self.__dust.getDust())
    return
  def addMolecule(self, molecule, transition):
    (numbers,species,transitions,frequencies) = self.__observations.speciesData
    i = (species==molecule)&(transitions==transition)
    if molecule in self.__molecules.getMolecules():
      self.__molecules.addTransition(molecule, transition, frequencies[i], numbers[i])
      #self.__moleculeNumber.append(numbers[species==molecule and transitions==transition])
      #self.__moleculeTransitions[molecule].append(transition)
      #self.__moleculeFrequencies[molecule].append(frequencies[species=='dust' and transitions==transition])
    else:
      self.__molecules.addMolecule(molecule, transition, frequencies[i], numbers[i])
      #self.__moleculeNames.append(molecule)
      #self.__moleculeNumber.append(numbers[species==molecule and transitions==transition])
      #self.__moleculeTransitions[molecule].append(transition)
      #self.__moleculeFrequencies[molecule].append(frequencies[species==molecule and transitions==transition])
    #self.__species = [self.__molecules, self.__dust]
    self.__speciesNames = np.append(self.__molecules.getMolecules(), self.__dust.getDust())
    return
  def addSpecies(self, speciesTransition):
    # find transition number as defined in
    # SetUpOrionBarModelEnvironment.v.1.1.nb
    #from PDR import _globals  
    #(number,species,transition,frequency) = self.__observations.__speciesData
    #gbl._globals['compound']['number'] = []
    #gbl._globals['compound']['frequency'] = []
    for species in speciesTransition:
      element = species.split(' ')[0]
      transition = int(species.split(' ')[1])
      if element=='Dust':
        self.addDust(element, transition)
      else:
        self.addMolecule(element, transition)
    return
  def calculateEmission(self):
    self.__grid.calculateEmission()
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
          self.__orientation.setLOS(self.__grid, x=x, y=y, dim=dim)
          position.append([x,y])
          intensity = self.__orientation.calculateRadiativeTransfer(velocity)
          intensityMap.append(intensity)
    if dim=='xz':
      xArray,zArray = np.unique([xPositions,zPositions], axis=1)
      for x in np.unique(xArray):
        for z in np.unique(zArray):
          self.__orientation.setLOS(self.__grid, x=x, z=z, dim=dim)
          position.append([x,z])
          intensity = self.__orientation.calculateRadiativeTransfer(velocity)
          intensityMap.append(intensity)
    if dim=='yz':
      yArray,zArray = np.unique([yPositions,zPositions], axis=1)
      for y in np.unique(yArray):
        for z in np.unique(zArray):
          self.__orientation.setLOS(self.__grid, y=y, z=z, dim=dim)
          position.append([y,z])
          intensity = self.__orientation.calculateRadiativeTransfer(velocity)
          intensityMap.append(intensity)
    self.__mapPositions = np.array(position)
    self.__intensityMap = np.array(intensityMap)
    return
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