import importlib as il
import sys

import numpy as np
import matplotlib.pyplot as plt
# from mpl_toolkits.mplot3d import Axes3D
from logging import getLogger, basicConfig, FileHandler, Formatter

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
  
    def __init__(self, history_path='', directory='', folder='', 
                 tau_grid_file='clump_tau_LineCenter.dat', 
                 tb_grid_file='clump_Tmb_LineCenter', 
                 tau_fuv_grid_file='RhoMassAFUV.dat',
                 h2_mass_file='h2_mass_profile.dat', 
                 hi_mass_file='hi_mass_profile.dat', 
                 density_file='densities_clouds.dat', 
                 fuv_file='galactic_FUV_complete.dat', 
                 l_range=(912, 2066),
                 velocity_file='rot_milki2018_14.dat',
                 x=0, y=0, z=0, modelType='', resolution=1000,
                 molecules='all', dust='molecular', velocityRange=(), velocityNumber=0,
                 clumpMassRange=((0, 2), (-2)), clumpMassNumber=(3, 1), clumpNmax=(1, 100), 
                 ensemble_mass_factor=(1, 1), interclump_hi_ratio=1,
                 interclump_fillingfactor=None, interclumpLogFUV=None, clumpLogFUV=None,
                 hi_mass_factor=1, h2_mass_factor=1, fuv_factor=1, density_factor=1, globalUV=10, 
                 r_cmz=0, zeta_cmz=1e-14, zeta_sol=2e-16, suggested_calc=True,
                 dilled=False, timed=False, verbose=False, debug=False):
      
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
        constants.history = folder
        constants.changeDirectory(directory)
        
        # Factors
        constants.ensemble_mass_factor = ensemble_mass_factor
        constants.hi_mass_factor = hi_mass_factor
        constants.h2_mass_factor = h2_mass_factor
        constants.interclump_hi_ratio= interclump_hi_ratio
        constants.interclump_fillingfactor = interclump_fillingfactor
        constants.density_factor = density_factor
        constants.fuv_factor = fuv_factor
        # constants.globalUV = globalUV
        constants.clumpLogFUV = clumpLogFUV
        constants.interclumpLogFUV = interclumpLogFUV
        constants.r_cmz = r_cmz
        constants.zeta_cmz = zeta_cmz
        constants.zeta_sol = zeta_sol

        # Clump properties
        constants.changeVelocityRange(velocityRange)
        constants.changeVelocityNumber(velocityNumber)
        constants.addClumps(massRange=clumpMassRange, num=clumpMassNumber, Nmax=clumpNmax, reset=True)
        
        # Read grid & input data, specify transitions, and interpolate
        observations.methods.initialise_grid(tau_grid_file=tau_grid_file, 
                                             tb_grid_file=tb_grid_file, 
                                             tau_fuv_file=tau_fuv_file)
        observations.methods.initialise_input(h2_mass_file=h2_mass_file, 
                                              hi_mass_file=hi_mass_file, 
                                              density_file=density_file, 
                                              fuv_file=fuv_file,
                                              velocity_file=velocity_file)
        constants.dust = dust
        self.__addSpecies(molecules)
        constants.changeDustWavelengths(constants.dust)
        if not interpolations.initialised:
            interpolations.initialise_grid(dilled=dilled)
            interpolations.initialise_model(l_range=l_range)

        # Initialise logger
        self.__logger = getLogger()

        # Shape() object to create the parameters for the grid of voxels
        self.__shape = shape.Shape(x, y, z, modelType=modelType)
        # VoxelGrid() object to build the model and calculate the emission
        self.__grid = VoxelGrid(self.__shape, suggested_calc=suggested_calc)
        # Orientation() object to change the viewing angle and expected spectra
        # self.__orientation = Orientation(self.__shape.getDimensions())
        self.__speciesNames = []  # this is a list of the species names for easy printout
        self.__timed = timed
        self.__verbose = verbose
        self.__debug = debug
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
        # Point logger to file
        format_str = '\n\n%(levelname)s [%(name)s]: %(message)s\n\n'
        filename = constants.HISTORYPATH + constants.directory + constants.history + 'log.txt'
        filehandler = FileHandler(filename, mode='w')
        if self.__debug:
            basicConfig(format=format_str, level='DEBUG')
            #self.__logger.setLevel('DEBUG')
            #filehandler.setLevel('DEBUG')
        elif self.__verbose:
            basicConfig(format=format_str, level='INFO')
            #self.__logger.setLevel('INFO')
            #filehandler.setLevel('INFO')
        else:
            basicConfig(format=format_str, level='WARNING')
            #self.__logger.setLevel('WARNING')
            #filehandler.setLevel('WARNING')
        #filehandler.setFormatter(Formatter(format_str))
        #self.__logger.addHandler(filehandler)

        # Calculate emission
        self.__grid.calculateEmission(**kwargs)

        return
  
    def writeEmission(self):
        self.__grid.writeEmission(verbose=self.__verbose, debug=self.__debug)
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
      
    def plotModel(self, plot='total intensity'):
        positions = self.__grid.getVoxelPositions()
        limits = [positions.min(), positions.max()]
        if plot == 'total intensity':
            if self.__debug:
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
