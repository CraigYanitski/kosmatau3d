import astropy.units as u
import importlib as il
import logging
import matplotlib.pyplot as plt
import numpy as np
import os
import sys

from astropy.io import fits
from logging import getLogger, basicConfig, FileHandler, Formatter

from kosmatau3d.models import shape  # import Shape
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
                 tb_grid_file='clump_Tmb_LineCenter.dat', 
                 tau_fuv_grid_file='RhoMassAFUV.dat',
                 h2_mass_file='h2_mass_profile.dat', 
                 hi_mass_file='hi_mass_profile.dat', 
                 density_file='densities_clouds.dat', 
                 fuv_file='galactic_FUV_complete.dat', 
                 l_range=(912, 2066), average_fuv=False, like_clumps=False,
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
        constants.tb_grid_file = tb_grid_file
        constants.tau_grid_file = tau_grid_file
        constants.tau_fuv_grid_file = tau_fuv_grid_file
        observations.methods.initialise_grid(tau_grid_file=tau_grid_file, 
                                             tb_grid_file=tb_grid_file, 
                                             tau_fuv_grid_file=tau_fuv_grid_file)
        constants.h2_mass_file = h2_mass_file
        constants.hi_mass_file = hi_mass_file
        constants.density_file = density_file
        constants.fuv_file = fuv_file
        constants.velocity_file = velocity_file
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
            constants.average_fuv = average_fuv
            constants.l_range = l_range
            constants.like_clumps = like_clumps
            interpolations.initialise_model(l_range=l_range, average_fuv=average_fuv, 
                                            like_clumps=like_clumps)

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


class SyntheticModel(object):
    '''
    This is an object to load individual `kosmatau3d` models. This in merely for
    the convenience of examining the model information in a consistent manner.
    There is an optional argument when initialising to set a base directory, which
    makes it easier to load multiple models in succession. Due to the complexity
    of the `kosmatau3d` models, it is not recommended to load multiple models at
    the same time.
    '''
    
    def __init__(self, base_dir=''):
        '''
        This initialises the object along with the base directory.
        The owned objects of `base_dir` and `files` are created.
        `files` can be modified again when loading a model, but for now it
        has the default filenames created with `kosmatau3d`.

        :param base_dir: the base directory to use when loading models. Default: `''`.


        '''

        self.base_dir = base_dir
        self.files = {'intensity': 'synthetic_intensity', 
                      'optical_depth': 'synthetic_optical_depth', 
                      'dust_absorption': 'dust_absorption', 
                      'dust_emissivity': 'dust_emissivity', 
                      'species_absorption': 'species_absorption', 
                      'species_emissivity': 'species_emissivity', 
                      'clump_number': 'clump_number',
                      'clump_radius': 'clump_radius',
                      'density': 'voxel_density', 
                      'ensemble_dispersion': 'voxel_ensemble_dispersion', 
                      'ensemble_mass': 'voxel_ensemble_mass', 
                      'fuv_absorption': 'voxel_FUVabsorption', 
                      'fuv': 'voxel_fuv', 
                      'position': 'voxel_position', 
                      'velocity': 'voxel_velocity', 
                      'los_count': 'sightlines', 
                      'log': 'log', }
        
        return

    def test_kwargs(model, **kwargs):
        
        def change_files(self, **kwargs):
            
            # # Ensure model path exists
            # if os.path.exists(self.base_dir + kwargs['directory']):
            #     pass
            # else:
            #     print(f'Directory {self.base_dir + kwargs["directory"]} is not valid!! Check to see '
            #           + 'if `base_dir` was set when initialising or if the specified directory was correct.')
            #     return
            
            # Update model file information
            kwarg_keys = kwargs.keys()
            # self.files['directory'] = kwargs['directory']
            if 'intensity' in kwarg_keys:
                self.files['intensity'] = kwargs['intensity']
            if 'optical_depth' in kwarg_keys:
                self.files['optical_depth'] = kwargs['optical_depth']
            if 'dust_absorption' in kwarg_keys:
                self.files['dust_absorption'] = kwargs['dust_absorption']
            if 'dust_emissivity' in kwarg_keys:
                self.files['dust_emissivity'] = kwargs['dust_emissivity']
            if 'species_absorption' in kwarg_keys:
                self.files['species_absorption'] = kwargs['species_absorption']
            if 'species_emissivity' in kwarg_keys:
                self.files['species_emissivity'] = kwargs['species_emissivity']
            if 'clump_number' in kwarg_keys:
                self.files['clump_number'] = kwargs['clump_number']
            if 'clump_radius' in kwarg_keys:
                self.files['clump_radius'] = kwargs['clump_radius']
            if 'density' in kwarg_keys:
                self.files['density'] = kwargs['density']
            if 'ensemble_dispersion' in kwarg_keys:
                self.files['ensemble_dispersion'] = kwargs['ensemble_dispersion']
            if 'ensemble_mass' in kwarg_keys:
                self.files['ensemble_mass'] = kwargs['ensemble_mass']
            if 'fuv_absorption' in kwarg_keys:
                self.files['fuv_absorption'] = kwargs['fuv_absorption']
            if 'fuv' in kwarg_keys:
                self.files['fuv'] = kwargs['fuv']
            if 'position' in kwarg_keys:
                self.files['position'] = kwargs['position']
            if 'velocity' in kwarg_keys:
                self.files['velocity'] = kwargs['velocity']
            if 'los_count' in kwarg_keys:
                self.files['los_count'] = kwargs['los_count']
            if 'log' in kwarg_keys:
                self.files['log'] = kwargs['log']

            return model(self, **kwargs)

        return change_files

    @test_kwargs
    def load_model(self, directory='', map_units='deg', **kwargs):
        '''
        Load all of the data for one model. Any additional information such as 
        observing velocities, latitude, and longitude are computed as well.

        **Note** that this can be quite computationally expensive. It is not recommended 
        to load multiple models at the same time.

        :param directory: The directory of all of the model. Note that this is 
                          appended to `self.base_dir`.
        :param kwargs: optional kwargs to modify the model files used (specified in 
                       `self.files`).


        '''

        # Ensure model path exists
        if os.path.exists(self.base_dir + directory + self.files['intensity'] + '.fits'):
            self.files['directory'] = directory
        else:
            print(f'Directory {self.base_dir + directory} is not valid!! Check to see '
                  +'if `base_dir` was set when initialising or if the specified directory was correct.')
            return

        # # Update model file information
        # kwarg_keys = kwargs.keys()
        # if 'intensity' in kwarg_keys:
        #     self.files['intensity'] = kwargs['intensity']
        # if 'optical_depth' in kwarg_keys:
        #     self.files['optical_depth'] = kwargs['optical_depth']
        # if 'dust_absorption' in kwarg_keys:
        #     self.files['dust_absorption'] = kwargs['dust_absorption']
        # if 'dust_emissivity' in kwarg_keys:
        #     self.files['dust_emissivity'] = kwargs['dust_emissivity']
        # if 'species_absorption' in kwarg_keys:
        #     self.files['species_absorption'] = kwargs['species_absorption']
        # if 'species_emissivity' in kwarg_keys:
        #     self.files['species_emissivity'] = kwargs['species_emissivity']
        # if 'density' in kwarg_keys:
        #     self.files['density'] = kwargs['density']
        # if 'ensemble_dispersion' in kwarg_keys:
        #     self.files['ensemble_dispersion'] = kwargs['ensemble_dispersion']
        # if 'ensemble_mass' in kwarg_keys:
        #     self.files['ensemble_mass'] = kwargs['ensemble_mass']
        # if 'fuv_absorption' in kwarg_keys:
        #     self.files['fuv_absorption'] = kwargs['fuv_absorption']
        # if 'fuv' in kwarg_keys:
        #     self.files['fuv'] = kwargs['fuv']
        # if 'position' in kwarg_keys:
        #     self.files['position'] = kwargs['position']
        # if 'velocity' in kwarg_keys:
        #     self.files['velocity'] = kwargs['velocity']
        # if 'los_count' in kwarg_keys:
        #     self.files['los_count'] = kwargs['los_count']
        # if 'log' in kwarg_keys:
        #     self.files['log'] = kwargs['log']
        
        # Load all model data (can be expensive for memory)
        self.intensity_file = fits.open(self.base_dir + directory 
                                        + self.files['intensity'] + '.fits')
        self.map_positions = self.intensity_file[0].data
        self.intensity_species = self.intensity_file[1].data
        self.intensity_dust = self.intensity_file[2].data
        self.optical_depth_file = fits.open(self.base_dir + directory 
                                            + self.files['optical_depth'] + '.fits')
        self.optical_depth_species = self.optical_depth_file[1].data
        self.optical_depth_dust = self.optical_depth_file[2].data
        self.dust_absorption_file = fits.open(self.base_dir + directory 
                                              + self.files['dust_absorption'] + '.fits')
        self.dust_absorption = self.dust_absorption_file[0].data
        self.dust_emissivity_file = fits.open(self.base_dir + directory 
                                              + self.files['dust_emissivity'] + '.fits')
        self.dust_emissivity = self.dust_emissivity_file[0].data
        self.species_absorption_file = fits.open(self.base_dir + directory 
                                                 + self.files['species_absorption'] + '.fits')
        self.species_absorption = self.species_absorption_file[0].data
        self.species_emissivity_file = fits.open(self.base_dir + directory 
                                                 + self.files['species_emissivity'] + '.fits')
        self.species_emissivity = self.species_emissivity_file[0].data
        if os.path.exists(self.base_dir + directory + self.files['clump_number'] + '.fits'):
            self.clump_number_file = fits.open(self.base_dir + directory 
                                               + self.files['clump_number'] + '.fits')
            self.clump_number = self.clump_number_file[0].data
        if os.path.exists(self.base_dir + directory + self.files['clump_radius'] + '.fits'):
            self.clump_radius_file = fits.open(self.base_dir + directory 
                                               + self.files['clump_radius'] + '.fits')
            self.clump_radius = self.clump_radius_file[0].data
        self.density_file = fits.open(self.base_dir + directory 
                                      + self.files['density'] + '.fits')
        self.density = self.density_file[0].data
        self.ensemble_dispersion_file = fits.open(self.base_dir + directory 
                                                  + self.files['ensemble_dispersion'] + '.fits')
        self.ensemble_dispersion = self.ensemble_dispersion_file[0].data
        self.ensemble_mass_file = fits.open(self.base_dir + directory 
                                            + self.files['ensemble_mass'] + '.fits')
        self.ensemble_mass = self.ensemble_mass_file[0].data
        self.fuv_absorption_file = fits.open(self.base_dir + directory 
                                             + self.files['fuv_absorption'] + '.fits')
        self.fuv_absorption = self.fuv_absorption_file[0].data
        self.fuv_file = fits.open(self.base_dir + directory 
                                  + self.files['fuv'] + '.fits')
        self.fuv = self.fuv_file[0].data
        self.position_file = fits.open(self.base_dir + directory 
                                       + self.files['position'] + '.fits')
        self.position = self.position_file[0].data
        self.velocity_file = fits.open(self.base_dir + directory 
                                       + self.files['velocity'] + '.fits')
        self.velocity = self.velocity_file[0].data
        self.los_count = np.loadtxt(self.base_dir + directory + self.files['los_count'] + '.csv')
        with open(self.base_dir + directory + self.files['log'] + '.txt') as f:
            self.log = f.readlines()
        
        # Extract headers and create additional axes
        self.info = self.species_absorption_file[0].header['COMMENT']
        self.species_header = self.intensity_file[1].header
        self.species_header['BUNIT'] = self.intensity_file[1].header['BUNIT'] + '/' \
                                       + self.optical_depth_file[1].header['BUNIT']
        self.dust_header = self.intensity_file[2].header
        self.dust_header['BUNIT'] = self.intensity_file[2].header['BUNIT'] + '/' \
                                    + self.optical_depth_file[2].header['BUNIT']
        self.species = self.species_header['SPECIES'].split(', ')
        self.dust = self.dust_header['DUST'].split(', ')
        self.dust_header['BUNIT'] = self.intensity_file[2].header['BUNIT'] + '/' \
                                    + self.optical_depth_file[2].header['BUNIT']
        self.map_lon = np.linspace(self.species_header['CRVAL2'] 
                                   - self.species_header['CDELT2']*(self.species_header['CRPIX2']-0.5), 
                                   self.species_header['CRVAL2'] 
                                   + self.species_header['CDELT2']*(self.species_header['NAXIS2']
                                                                    -self.species_header['CRPIX2']-0.5), 
                                   num=self.species_header['NAXIS2'])
        self.map_lat = np.linspace(self.species_header['CRVAL3'] 
                                   - self.species_header['CDELT3']*(self.species_header['CRPIX3']-0.5), 
                                   self.species_header['CRVAL3'] 
                                   + self.species_header['CDELT3']*(self.species_header['NAXIS3']
                                                                    -self.species_header['CRPIX3']-0.5), 
                                   num=self.species_header['NAXIS3'])
        self.map_vel = np.linspace(self.species_header['CRVAL4'] 
                                   - self.species_header['CDELT4']*(self.species_header['CRPIX4']), 
                                   self.species_header['CRVAL4'] 
                                   + self.species_header['CDELT4']*(self.species_header['NAXIS4']
                                                                    -self.species_header['CRPIX4']-1), 
                                   num=self.species_header['NAXIS4'])
        
        # convert from radians to degrees if specified
        if map_units == 'deg' or map_units == 'degrees':
            self.map_lon *= 180/np.pi
            self.map_lat *= 180/np.pi

        return

    def get_dust_wavelengths(self):
        wav = list(map(u.Unit, self.dust))
        return list([w.decompose() for w in wav])

    def get_species_intensity(self, transition=None, idx=None, include_dust=False, integrated=False):

        if transition is None and idx is None:
            transition = self.species
            idx = range(len(self.species))
        elif isinstance(transition, (tuple, list, np.ndarray)):
            idx = (self.species.index(t) for t in transition)
        elif isinstance(idx, (tuple, list, np.ndarray)):
            transition = (self.species[i] for i in idx)
        elif isinstance(transition, str) and transition in self.species:
            transition = (transition, )
            idx = (self.species.index(transition[0]), )
        elif isinstance(idx, int):
            transition = (self.species[idx], )
            idx = (idx, )
        else:
            print('Please enter a valid value for transition or idx')
            return

        intensity = []

        for i in idx:
            if include_dust:
                intensity_temp = self.intensity_species[:, :, :, i]
            else:
                if len(self.dust) > 10:
                    wav_dust = self.get_dust_wavelengths()
                    f = interp1d(self.intensity_dust, wav_dust, axis=2, kind='slinear', 
                                 fill_value='extrapolate')
                    intensity_dust = f(wav_species[i])
                else:
                    intensity_dust = self.intensity_dust.mean(2)
                intensity_temp = self.intensity_species[:, :, :, i] - intensity_dust
            if integrated:
                intensity.append(np.trapz(intensity_temp, self.map_vel, axis=0))
            else:
                intensity.append(intensity_temp)

        if len(intensity) > 1:
            return deepcopy(intensity)
        else:
            return deepcopy(intensity[0])

    def get_species_optical_depth(self, transition=None, idx=None, include_dust=False, integrated=False):

        if transition is None and idx is None:
            transition = self.species
            idx = range(len(self.species))
        elif isinstance(transition, (tuple, list, np.ndarray)):
            idx = (self.species.index(t) for t in transition)
        elif isinstance(idx, (tuple, list, np.ndarray)):
            transition = (self.species[i] for i in idx)
        elif isinstance(transition, str) and transition in self.species:
            transition = (transition, )
            idx = (self.species.index(transition[0]), )
        elif isinstance(idx, int):
            transition = (self.species[idx], )
            idx = (idx, )
        else:
            print('Please enter a valid value for transition or idx')
            return

        optical_depth = []

        for i in idx:
            if include_dust:
                optical_depth_temp = self.optical_depth_species[:, :, :, i]
            else:
                if len(self.dust) > 10:
                    wav_dust = self.get_dust_wavelengths()
                    f = interp1d(self.optical_depth_dust, wav_dust, axis=2, kind='slinear', 
                                 fill_value='extrapolate')
                    optical_depth_dust = f(wav_species[i])
                else:
                    optical_depth_dust = self.optical_depth_dust.mean(2)
                optical_depth_temp = self.optical_depth_species[:, :, :, i] - optical_depth_dust
            if integrated:
                optical_depth.append(np.trapz(optical_depth_temp, self.map_vel, axis=0))
            else:
                optical_depth.append(optical_depth_temp)

        if len(optical_depth) > 1:
            return deepcopy(optical_depth)
        else:
            return deepcopy(optical_depth[0])

    def get_dust_intensity(self, wavelength=None, idx=None):

        if wavelength is None and idx is None:
            wavelength = self.dust
            idx = range(len(self.dust))
        elif isinstance(wavelength, (tuple, list, np.ndarray)):
            idx = (self.dust.index(t) for t in wavelength)
        elif isinstance(idx, (tuple, list, np.ndarray)):
            wavelength = (self.dust[i] for i in idx)
        elif isinstance(wavelength, str) and wavelength in self.dust:
            wavelength = (wavelength, )
            idx = (self.dust.index(wavelength[0]), )
        elif isinstance(idx, int):
            wavelength = (self.dust[idx], )
            idx = (idx, )
        else:
            print('Please enter a valid value for wavelength or idx')
            return

        intensity = []

        for i in idx:
            intensity.append(self.intensity_dust[:, :, i])
        
        if len(intensity) > 1:
            return deepcopy(intensity)
        else:
            return deepcopy(intensity[0])

    def get_dust_optical_depth(self, wavelength=None, idx=None):

        if wavelength is None and idx is None:
            wavelength = self.dust
            idx = range(len(self.dust))
        elif isinstance(wavelength, (tuple, list, np.ndarray)):
            idx = (self.dust.index(t) for t in wavelength)
        elif isinstance(idx, (tuple, list, np.ndarray)):
            wavelength = (self.dust[i] for i in idx)
        elif isinstance(wavelength, str) and wavelength in self.dust:
            wavelength = (wavelength, )
            idx = (self.dust.index(wavelength[0]), )
        elif isinstance(idx, int):
            wavelength = (self.dust[idx], )
            idx = (idx, )
        else:
            print('Please enter a valid value for wavelength or idx')
            return

        optical_depth = []

        for i in idx:
            optical_depth.append(self.optical_depth_dust[:, :, i])
        
        if len(optical_depth) > 1:
            return deepcopy(optical_depth)
        else:
            return deepcopy(optical_depth[0])


