import astropy.units as u
import importlib as il
import logging
import matplotlib as mpl
import matplotlib.pyplot as plt
import numpy as np
import os
import sys

from astropy.io import fits
from copy import copy, deepcopy
from logging import getLogger, basicConfig, FileHandler, Formatter
from scipy.interpolate import interp1d
from scipy.stats import binned_statistic

from kosmatau3d.models import constants
from kosmatau3d.models import interpolations
from kosmatau3d.models import observations
from kosmatau3d.models import shape  # import Shape
from kosmatau3d.models import species
from .VoxelGrid import VoxelGrid


class Model(object):
    '''
    This is the highest class in the hierarchy of the KOSMA-tau^3 simulation.
    It contains all of the information needed to properly model a PDR (I think).
    '''
    # PRIVATE
  
    def __init__(self, history_path='', directory='', folder='', 
                 clump_tau_grid_file='clump_tau_LineCenter.dat', 
                 clump_tb_grid_file='clump_Tmb_LineCenter.dat', 
                 clump_taufuv_grid_file='RhoMassAFUV.dat',
                 clump_column_density_file='clumpMeanCols.dat',
                 clump_temperature_file='clumpTemperatures_filled.dat',
                 interclump_tau_grid_file='interclumpTauLineCenter.dat', 
                 interclump_dust_tau_grid_file='interclumpDustTau.dat', 
                 interclump_tb_grid_file='interclumpTmbLineCenter.dat', 
                 interclump_dust_tb_grid_file='interclumpDustSED.dat', 
                 interclump_taufuv_grid_file='interclumpTauFUV.dat',
                 interclump_column_density_file='interclumpMeanCols.dat',
                 interclump_temperature_file='interclumpTemperatures_filled.dat',
                 h2_surface_density_file='h2_surface-density.dat', 
                 hi_surface_density_file='hi_surface-density.dat', 
                 h2_scale_height_file='h2_scale-height.dat', 
                 hi_scale_height_file='hi_scale-height.dat', 
                 h_number_density_file='h_number-density.dat', 
                 fuv_file='galactic_FUV_complete.dat', 
                 l_range=(912, 2066), average_fuv=False, scale_gc=1.0, mhi_gc=1.0, mh2_gc=1.0, r_gc=4400, 
                 like_clumps=False, all_full=False, 
                 velocity_file='rot_milki2018_14.dat', disp_core=None, r_core=4400, disp_gmc=None, 
                 x=0, y=0, z=0, model_type='', resolution=1000,
                 transitions='all', dust='molecular', velocity_range=(), velocity_number=0,
                 clump_mass_range=((0, 2), (-2)), clump_mass_number=(3, 1), clump_n_max=(1, 100), 
                 ensemble_mass_factor=(1, 1), interclump_idx=(False, True), interclump_hi_ratio=1,
                 interclump_fillingfactor=None, interclump_density=1911, interclump_log_fuv=None, clump_log_fuv=None,
                 hi_mass_factor=1, h2_mass_factor=1, fuv_factor=1, density_factor=1, global_uv=10, 
                 r_cmz=0, zeta_cmz=1e-14, zeta_sol=2e-16, new_grid=True, suggested_calc=True,
                 dilled=False, timed=False, verbose=False, debug=False):
      
        if not len(clump_mass_range):
            sys.exit('<<ERROR>> Define mass sets in argument.')
            # sys.exit()
        if not len(velocity_range):
            sys.exit('<<ERROR>> Define observing velocities in argument.')
            # sys.exit()
        
        # this just adds a label to the type of model being created. ie 'disk', 'bar', 'sphere', etc.
        constants.type = model_type
        constants.voxel_size = float(resolution)
        constants.HISTORYPATH = history_path
        constants.history = folder
        constants.change_directory(directory)
        
        # Factors
        constants.ensemble_mass_factor = ensemble_mass_factor
        constants.hi_mass_factor = hi_mass_factor
        constants.h2_mass_factor = h2_mass_factor
        constants.interclump_hi_ratio = interclump_hi_ratio
        constants.interclump_fillingfactor = interclump_fillingfactor
        constants.density_factor = density_factor
        constants.interclump_density = interclump_density
        constants.fuv_factor = fuv_factor
        # constants.globalUV = globalUV
        constants.clump_log_fuv = clump_log_fuv
        constants.interclump_log_fuv = interclump_log_fuv
        constants.r_cmz = r_cmz
        constants.zeta_cmz = zeta_cmz
        constants.zeta_sol = zeta_sol

        # Clump properties
        constants.change_velocity_range(velocity_range)
        constants.change_velocity_number(velocity_number)
        constants.add_clumps(mass_range=clump_mass_range, num=clump_mass_number, n_max=clump_n_max, reset=True)
        constants.set_interclump_ensemble(interclump_idx)
        
        # Read grid & input data, specify transitions, and interpolate
        constants.clump_species_tb_grid_file = clump_tb_grid_file
        constants.clump_species_tau_grid_file = clump_tau_grid_file
        constants.clump_dust_tb_grid_file = clump_tb_grid_file
        constants.clump_dust_tau_grid_file = clump_tau_grid_file
        constants.clump_taufuv_grid_file = clump_taufuv_grid_file
        constants.clump_column_density_file = clump_column_density_file
        constants.clump_temperature_file = clump_temperature_file
        constants.interclump_species_tb_grid_file = interclump_tb_grid_file
        constants.interclump_species_tau_grid_file = interclump_tau_grid_file
        constants.interclump_dust_tb_grid_file = interclump_tb_grid_file
        constants.interclump_dust_tau_grid_file = interclump_tau_grid_file
        constants.interclump_taufuv_grid_file = interclump_taufuv_grid_file
        constants.interclump_column_density_file = interclump_column_density_file
        constants.interclump_temperature_file = interclump_temperature_file
        observations.methods.initialise_grid(clump_tau_grid_file=clump_tau_grid_file, 
                                             interclump_tau_grid_file=interclump_tau_grid_file, 
                                             interclump_dust_tau_grid_file=interclump_dust_tau_grid_file, 
                                             clump_tb_grid_file=clump_tb_grid_file, 
                                             interclump_tb_grid_file=interclump_tb_grid_file, 
                                             interclump_dust_tb_grid_file=interclump_dust_tb_grid_file, 
                                             clump_taufuv_grid_file=clump_taufuv_grid_file, 
                                             interclump_taufuv_grid_file=interclump_taufuv_grid_file,
                                             clump_column_density_file=clump_column_density_file,
                                             interclump_column_density_file=interclump_column_density_file,
                                             clump_temperature_file=clump_temperature_file,
                                             interclump_temperature_file=interclump_temperature_file)
        constants.h2_surface_density_file = h2_surface_density_file
        constants.hi_surface_density_file = hi_surface_density_file
        constants.h2_scale_height_file = h2_scale_height_file
        constants.hi_scale_height_file = hi_scale_height_file
        constants.h_number_density_file = h_number_density_file
        constants.fuv_file = fuv_file
        constants.velocity_file = velocity_file
        observations.methods.initialise_input(h2_surface_density_file=h2_surface_density_file, 
                                              hi_surface_density_file=hi_surface_density_file, 
                                              h2_scale_height_file=h2_scale_height_file, 
                                              hi_scale_height_file=hi_scale_height_file, 
                                              h_number_density_file=h_number_density_file, 
                                              fuv_file=fuv_file,
                                              velocity_file=velocity_file)
        constants.dust = dust
        self.__add_transitions(transitions)
        constants.change_dust_wavelengths(constants.dust)
        constants.fuv_scale_gc = scale_gc
        constants.mhi_scale_gc = mhi_gc
        constants.mh2_scale_gc = mh2_gc
        constants.r_gc = r_gc
        constants.disp_gc = disp_core
        constants.disp_r_gc = r_core
        if disp_gmc:
            constants.disp_gmc = disp_gmc
        else:
            constants.disp_gmc = 1.1*constants.voxel_size**0.38
        if not interpolations.initialised or new_grid:
            interpolations.initialise_grid(dilled=dilled)
            constants.average_fuv = average_fuv
            constants.l_range = l_range
            constants.like_clumps = like_clumps
            constants.all_full = all_full
            interpolations.initialise_model(l_range=l_range, average_fuv=average_fuv, 
                                            all_full=all_full, like_clumps=like_clumps)

        # Initialise logger
        self.__logger = getLogger()

        # Shape() object to create the parameters for the grid of voxels
        self.__shape = shape.Shape(x, y, z, modelType=model_type)
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
  
    def __add_transitions(self, species_transition):
        species.add_transitions(species_transition)
        species.add_transitions(species_transition, interclump=True)
        return
  
    def __str__(self):
        printout = 'A {} model of {} voxels'.format(constants.type, self.getGrid().getVoxelNumber())
        if self.__verbose:
            printout += '\n  arranged in {}'.format(self.__shape.getDimensions())
            printout += '\n\nConsidering {} species:\n{}\n{}'.format(len(self.species_names),
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
        return species.species_names
  
    def getSpeciesNames(self):
        return species.species_names
  
    # def reloadModules(self):
    #   il.reload(shape)
    #   il.reload(voxelgrid)
    #   il.reload(orientation)
    #   il.reload(observations)
    #   il.reload(molecules)
    #   il.reload(dust)
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
        #filehandler.setformatter(Formatter(format_str))
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
            filename = constants.HISTORYPATH + constants.directory + \
                       'r{}_n{}/'.format(constants.resolution, self.__shape.voxelNumber()) + \
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
        self.hi_model = False
        self.hi_map = False
        self.files = {'intensity': 'synthetic_intensity', 
                      'optical_depth': 'synthetic_optical_depth', 
                      'hi_intensity': 'synthetic_hi_intensity', 
                      'hi_optical_depth': 'synthetic_hi_optical_depth', 
                      'f_vox': 'voxel-filling_factor', 
                      'dust_absorption': 'dust_absorption', 
                      'dust_emissivity': 'dust_emissivity', 
                      'species_absorption': 'species_absorption', 
                      'hi_emissivity': 'hi_emissivity', 
                      'hi_absorption': 'hi_absorption', 
                      'species_emissivity': 'species_emissivity', 
                      'clump_number': 'clump_number',
                      'clump_radius': 'clump_radius',
                      'density': 'voxel_density', 
                      'ensemble_dispersion': 'voxel_ensemble_dispersion', 
                      'ensemble_mass': 'voxel_ensemble_mass',
                      'hi_mass': 'voxel_hi_mass',
                      'h2_mass': 'voxel_h2_mass',
                      'fuv_absorption': 'voxel_FUVabsorption', 
                      'fuv': 'voxel_fuv', 
                      'position': 'voxel_position', 
                      'velocity': 'voxel_velocity', 
                      'los_count': 'sightlines', 
                      'log': 'log', }
        
        # initialise files
        self.intensity_file = None
        self.hi_intensity_file = None
        self.optical_depth_file = None
        self.hi_optical_depth_file = None

        self.dust_absorption_file = None
        self.dust_optical_depth_file = None
        self.species_absorption_file = None
        self.species_emissivity_file = None
        self.hi_absorption_file = None
        self.hi_optical_depth_file = None

        self.f_vox_file = None
        self.clump_number_file = None
        self.clump_radius_file = None
        self.density_file = None
        self.ensemble_dispersion_file = None
        self.ensemble_mass_file = None
        self.hi_mass_file = None
        self.h2_mass_file = None
        self.fuv_absorption_dispersion_file = None
        self.fuv_file = None
        self.position_file = None
        self.velocity_file = None
        
        return

    def close_files(model, **kwargs):
        '''
        Close any open FITS files.
        '''

        def wrapper(self, **kwargs):

            if isinstance(self.intensity_file, fits.hdu.hdulist.HDUList):
                self.intensity_file.close()
            if isinstance(self.hi_intensity_file, fits.hdu.hdulist.HDUList):
                self.hi_intensity_file.close()
            if isinstance(self.optical_depth_file, fits.hdu.hdulist.HDUList):
                self.optical_depth_file.close()
            if isinstance(self.hi_optical_depth_file, fits.hdu.hdulist.HDUList):
                self.hi_optical_depth_file.close()

            if isinstance(self.dust_absorption_file, fits.hdu.hdulist.HDUList):
                self.dust_absorption_file.close()
            if isinstance(self.dust_optical_depth_file, fits.hdu.hdulist.HDUList):
                self.dust_optical_depth_file.close()
            if isinstance(self.species_absorption_file, fits.hdu.hdulist.HDUList):
                self.species_absorption_file.close()
            if isinstance(self.species_emissivity_file, fits.hdu.hdulist.HDUList):
                self.species_emissivity_file.close()
            if isinstance(self.hi_absorption_file, fits.hdu.hdulist.HDUList):
                self.hi_absorption_file.close()
            if isinstance(self.hi_optical_depth_file, fits.hdu.hdulist.HDUList):
                self.hi_optical_depth_file.close()

            if isinstance(self.f_vox_file, fits.hdu.hdulist.HDUList):
                self.f_vox_file.close()
            if isinstance(self.clump_number_file, fits.hdu.hdulist.HDUList):
                self.clump_number_file.close()
            if isinstance(self.clump_radius_file, fits.hdu.hdulist.HDUList):
                self.clump_radius_file.close()
            if isinstance(self.density_file, fits.hdu.hdulist.HDUList):
                self.density_file.close()
            if isinstance(self.ensemble_dispersion_file, fits.hdu.hdulist.HDUList):
                self.ensemble_dispersion_file.close()
            if isinstance(self.ensemble_mass_file, fits.hdu.hdulist.HDUList):
                self.ensemble_mass_file.close()
            if isinstance(self.hi_mass_file, fits.hdu.hdulist.HDUList):
                self.hi_mass_file.close()
            if isinstance(self.h2_mass_file, fits.hdu.hdulist.HDUList):
                self.h2_mass_file.close()
            if isinstance(self.fuv_absorption_dispersion_file, fits.hdu.hdulist.HDUList):
                self.fuv_absorption_dispersion_file.close()
            if isinstance(self.fuv_file, fits.hdu.hdulist.HDUList):
                self.fuv_file.close()
            if isinstance(self.position_file, fits.hdu.hdulist.HDUList):
                self.position_file.close()
            if isinstance(self.velocity_file, fits.hdu.hdulist.HDUList):
                self.velocity_file.close()

            return model(self, **kwargs)

        return wrapper

    def change_files(model, **kwargs):
        '''
        Change the specified filenames.
        '''

        def wrapper(self, **kwargs):
            
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
            if 'hi_intensity' in kwarg_keys:
                self.files['hi_intensity'] = kwargs['hi_intensity']
            if 'hi_optical_depth' in kwarg_keys:
                self.files['hi_optical_depth'] = kwargs['hi_optical_depth']
            if 'f_vox' in kwarg_keys:
                self.files['f_vox'] = kwargs['f_vox']
            if 'dust_absorption' in kwarg_keys:
                self.files['dust_absorption'] = kwargs['dust_absorption']
            if 'dust_emissivity' in kwarg_keys:
                self.files['dust_emissivity'] = kwargs['dust_emissivity']
            if 'species_absorption' in kwarg_keys:
                self.files['species_absorption'] = kwargs['species_absorption']
            if 'species_emissivity' in kwarg_keys:
                self.files['species_emissivity'] = kwargs['species_emissivity']
            if 'hi_absorption' in kwarg_keys:
                self.files['hi_absorption'] = kwargs['hi_absorption']
            if 'hi_emissivity' in kwarg_keys:
                self.files['hi_emissivity'] = kwargs['hi_emissivity']
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
            if 'hi_mass' in kwarg_keys:
                self.files['hi_mass'] = kwargs['hi_mass']
            if 'h2_mass' in kwarg_keys:
                self.files['h2_mass'] = kwargs['h2_mass']
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

        return wrapper

    @change_files
    @close_files
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
            model_files = os.listdir(self.base_dir + directory)
        else:
            raise FileNotFoundError(f'Directory {self.base_dir + directory} is not valid!! Check '
                                    +'to see if `base_dir` was set when initialising or if the '
                                    +'specified directory was correct.')

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
        if (self.files['f_vox'] + '.fits') in model_files:
            self.f_vox_file = fits.open(self.base_dir + directory + self.files['f_vox'] + '.fits')
            self.f_vox = self.f_vox_file[0].data
        else:
            self.f_vox_file = np.nan
            self.f_vox = np.nan
        if os.path.exists(self.base_dir + directory + self.files['hi_intensity'] + '.fits'):
            self.hi_map = True
            self.hi_intensity_file = fits.open(self.base_dir + directory 
                                               + self.files['hi_intensity'] + '.fits')
            self.hi_map_positions = self.hi_intensity_file[0].data
            self.hi_intensity_species = self.hi_intensity_file[1].data
            self.hi_intensity_dust = self.hi_intensity_file[2].data
            self.hi_optical_depth_file = fits.open(self.base_dir + directory 
                                                   + self.files['hi_optical_depth'] + '.fits')
            self.hi_optical_depth_species = self.hi_optical_depth_file[1].data
            self.hi_optical_depth_dust = self.hi_optical_depth_file[2].data
        else:
            print('HI intensity map not found')
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
        if os.path.exists(self.base_dir + directory + self.files['hi_emissivity'] + '.fits'):
            # print('Loading HI emissivity')
            self.hi_model = True
            self.hi_emissivity_file = fits.open(self.base_dir + directory 
                                                + self.files['hi_emissivity'] + '.fits')
            self.hi_emissivity = self.hi_emissivity_file[0].data
            self.hi_absorption_file = fits.open(self.base_dir + directory 
                                                + self.files['hi_absorption'] + '.fits')
            self.hi_absorption = self.hi_absorption_file[0].data
        else:
            self.hi_model = False
            self.hi_absorption_file = np.nan
            self.hi_absorption = np.nan
            self.hi_emissivity_file = np.nan
            self.hi_emissivity = np.nan
            print('HI emissivity not found')
        if os.path.exists(self.base_dir + directory + self.files['clump_number'] + '.fits'):
            self.clump_number_file = fits.open(self.base_dir + directory 
                                               + self.files['clump_number'] + '.fits')
            self.clump_number = self.clump_number_file[0].data
        else:
            self.clump_number_file = np.nan
            self.clump_number = np.nan
        if os.path.exists(self.base_dir + directory + self.files['clump_radius'] + '.fits'):
            self.clump_radius_file = fits.open(self.base_dir + directory 
                                               + self.files['clump_radius'] + '.fits')
            self.clump_radius = self.clump_radius_file[0].data
        else:
            self.clump_radius_file = np.nan
            self.clump_radius = np.nan
        self.density_file = fits.open(self.base_dir + directory 
                                      + self.files['density'] + '.fits')
        self.density = self.density_file[0].data
        self.ensemble_dispersion_file = fits.open(self.base_dir + directory 
                                                  + self.files['ensemble_dispersion'] + '.fits')
        self.ensemble_dispersion = self.ensemble_dispersion_file[0].data
        self.ensemble_mass_file = fits.open(self.base_dir + directory 
                                            + self.files['ensemble_mass'] + '.fits')
        self.ensemble_mass = self.ensemble_mass_file[0].data
        if (self.files['hi_mass'] + '.fits') in model_files:
            self.hi_mass_file = fits.open(self.base_dir + directory + self.files['hi_mass'] + '.fits')
            self.hi_mass = self.hi_mass_file[0].data
        else:
            self.hi_mass_file = np.nan
            self.hi_mass = np.nan
        if (self.files['h2_mass'] + '.fits') in model_files:
            self.h2_mass_file = fits.open(self.base_dir + directory + self.files['h2_mass'] + '.fits')
            self.h2_mass = self.h2_mass_file[0].data
        else:
            self.h2_mass_file = np.nan
            self.h2_mass = np.nan
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
        if os.path.exists(self.base_dir + directory + self.files['log'] + '.txt'):
            with open(self.base_dir + directory + self.files['log'] + '.txt') as f:
                self.log = f.readlines()
        else:
            self.log = np.nan
        
        # Extract headers and create additional axes
        self.info = self.species_absorption_file[0].header['COMMENT']
        self.species_header = self.intensity_file[1].header
        self.species_header['BUNIT'] = self.intensity_file[1].header['BUNIT'] + '/' \
                                       + self.optical_depth_file[1].header['BUNIT']
        self.dust_header = self.intensity_file[2].header
        self.dust_header['BUNIT'] = self.intensity_file[2].header['BUNIT'] + '/' \
                                    + self.optical_depth_file[2].header['BUNIT']
        self.species = np.array(self.species_header['SPECIES'].split(', '))
        self.dust = np.array(self.dust_header['DUST'].split(', '))
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
            self.map_lon = (self.map_lon * 180/np.pi).round(decimals=8)
            self.map_lat = (self.map_lat * 180/np.pi).round(decimals=8)

        # Extract voxel size from model info
        self.ds = float(self.info[1].split()[1])

        return

    def get_dust_wavelengths(self):
        wav = list(map(u.Unit, self.dust))
        return list([w.decompose() for w in wav])

    def get_volume_filling_factor(self):
        return (4/3*np.pi*self.clump_number*self.clump_radius**3).sum(1)/self.ds**3

    def get_model_species_emissivity(self, transition=None, idx=None, include_dust=False):

        if transition is None and idx is None:
            transition = self.species
            idx = range(len(self.species))
        elif isinstance(transition, (tuple, list, np.ndarray)):
            idx = (list(self.species).index(t) for t in transition)
        elif isinstance(idx, (tuple, list, np.ndarray)):
            transition = (self.species[i] for i in idx)
        elif isinstance(transition, str) and transition in self.species:
            transition = (transition, )
            idx = (list(self.species).index(transition[0]), )
        elif isinstance(idx, int):
            transition = (self.species[idx], )
            idx = (idx, )
        else:
            print('Please enter a valid value for transition or idx')
            return

        emissivity = []

        for i in idx:
            if include_dust:
                emissivity.append(self.species_emissivity[:, :, i])
            else:
                # if len(self.dust) > 10:
                #     wav_dust = self.get_emissivity_wavelengths()
                #     f = interp1d(self.emissivity_dust, wav_dust, axis=2, kind='slinear', 
                #                  fill_value='extrapolate')
                #     emissivity_dust = f(wav_species[i])
                # else:
                #     emissivity_dust = self.emissivity_dust.mean(2)
                emissivity_dust = self.species_emissivity[:, :, i].min(1).reshape(-1, 1)
                emissivity.append(self.species_emissivity[:, :, i] - emissivity_dust)

        if len(emissivity) > 1:
            return np.array(emissivity)
        else:
            return np.array(emissivity[0])

    def get_model_species_absorption(self, transition=None, idx=None, include_dust=False):

        if transition is None and idx is None:
            transition = self.species
            idx = range(len(self.species))
        elif isinstance(transition, (tuple, list, np.ndarray)):
            idx = (list(self.species).index(t) for t in transition)
        elif isinstance(idx, (tuple, list, np.ndarray)):
            transition = (self.species[i] for i in idx)
        elif isinstance(transition, str) and transition in self.species:
            transition = (transition, )
            idx = (list(self.species).index(transition[0]), )
        elif isinstance(idx, int):
            transition = (self.species[idx], )
            idx = (idx, )
        else:
            print('Please enter a valid value for transition or idx')
            return

        absorption = []

        for i in idx:
            if include_dust:
                absorption.append(self.species_absorption[:, :, i])
            else:
                # if len(self.dust) > 10:
                #     wav_dust = self.get_absorption_wavelengths()
                #     f = interp1d(self.absorption_dust, wav_dust, axis=2, kind='slinear', 
                #                  fill_value='extrapolate')
                #     absorption_dust = f(wav_species[i])
                # else:
                #     absorption_dust = self.absorption_dust.mean(2)
                absorption_dust = self.species_absorption[:, :, i].min(1).reshape(-1, 1)
                absorption.append(self.species_absorption[:, :, i] - absorption_dust)

        if len(absorption) > 1:
            return np.array(absorption)
        else:
            return np.array(absorption[0])

    def get_model_species_intensity(self, transition=None, idx=None, include_dust=False, integrated=False):

        eps = self.get_model_species_emissivity(transition=transition, idx=idx, include_dust=include_dust)
        kap = self.get_model_species_absorption(transition=transition, idx=idx, include_dust=include_dust)
        intensity = np.zeros_like(eps)
        i_nan = kap == 0
        intensity[~i_nan] = eps[~i_nan]/kap[~i_nan] * (1 - np.exp(-kap[~i_nan]*self.ds))
        intensity[i_nan] = eps[i_nan]

        return np.array(intensity)

    def get_model_hi_emissivity(self, include_dust=False):

        if include_dust:
            emissivity = self.hi_emissivity[:, :, 0]
        else:
            # if len(self.dust) > 10:
            #     wav_dust = self.get_dust_wavelengths()
            #     f = interp1d(self.dust_emissivity, wav_dust, axis=2, kind='slinear', 
            #                  fill_value='extrapolate')
            #     emissivity_dust = f(wav_species[i])
            # else:
            #     emissivity_dust = self.emissivity_dust.mean(2)
            emissivity_dust = self.hi_emissivity[:, :, 0].min(1).reshape(-1, 1)
            emissivity = self.hi_emissivity[:, :, 0] - emissivity_dust

        return emissivity

    def get_model_hi_absorption(self, include_dust=False):

        if include_dust:
            absorption = self.hi_absorption[:, :, 0]
        else:
            # if len(self.dust) > 10:
            #     wav_dust = self.get_dust_wavelengths()
            #     f = interp1d(self.dust_absorption, wav_dust, axis=2, kind='slinear', 
            #                  fill_value='extrapolate')
            #     absorption_dust = f(wav_species[i])
            # else:
            #     absorption_dust = self.absorption_dust.mean(2)
            absorption_dust = self.hi_absorption[:, :, 0].min(1).reshape(-1, 1)
            absorption = self.hi_absorption[:, :, 0] - absorption_dust

        return absorption

    def get_model_hi_intensity(self, include_dust=False, integrated=False):

        eps = self.get_model_hi_emissivity(include_dust=include_dust)
        kap = self.get_model_hi_absorption(include_dust=include_dust)
        intensity = np.zeros_like(eps)
        i_nan = kap == 0
        intensity[~i_nan] = eps[~i_nan]/kap[~i_nan] * (1 - np.exp(-kap[~i_nan]*self.ds))
        intensity[i_nan] = eps[i_nan]

        return intensity

    def get_model_dust_emissivity(self, wavelength=None, idx=None):

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

        emissivity = []

        for i in idx:
            emissivity.append(self.dust_emissivity[:, :, i])
        
        if len(emissivity) > 1:
            return np.array(emissivity)
        else:
            return np.array(emissivity[0])

    def get_model_dust_absorption(self, wavelength=None, idx=None):

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

        absorption = []

        for i in idx:
            absorption.append(self.dust_absorption[:, :, i])
        
        if len(absorption) > 1:
            return np.array(absorption)
        else:
            return np.array(absorption[0])

    def get_model_dust_intensity(self, wavelength=None, idx=None):

        eps = self.get_model_dust_emissivity(wavelength=wavelength, idx=idx)
        kap = self.get_model_dust_absorption(wavelength=wavelength, idx=idx)
        intensity = np.zeros_like(eps)
        i_nan = kap == 0
        intensity[~i_nan] = eps[~i_nan]/kap[~i_nan] * (1 - np.exp(-kap[~i_nan]*self.ds))
        intensity[i_nan] = eps[i_nan]

        return np.array(intensity)

    def get_species_intensity(self, transition=None, idx=None, include_dust=False, integrated=False):

        if transition is None and idx is None:
            transition = self.species
            idx = range(len(self.species))
        elif isinstance(transition, (tuple, list, np.ndarray)):
            idx = (list(self.species).index(t) for t in transition)
        elif isinstance(idx, (tuple, list, np.ndarray)):
            transition = (self.species[i] for i in idx)
        elif isinstance(transition, str) and transition in self.species:
            transition = (transition, )
            idx = (list(self.species).index(transition[0]), )
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
                # if len(self.dust) > 10:
                #     wav_dust = self.get_dust_wavelengths()
                #     f = interp1d(self.intensity_dust, wav_dust, axis=2, kind='slinear', 
                #                  fill_value='extrapolate')
                #     intensity_dust = f(wav_species[i])
                # else:
                #     intensity_dust = self.intensity_dust.mean(2)
                intensity_dust = self.intensity_species[:, :, :, i].min(0)
                intensity_temp = self.intensity_species[:, :, :, i] - intensity_dust
            if integrated:
                intensity.append(np.trapz(intensity_temp, self.map_vel, axis=0))
            else:
                intensity.append(copy(intensity_temp))

        if len(intensity) > 1:
            return np.array(intensity)
        else:
            return np.array(intensity[0])

    def get_species_optical_depth(self, transition=None, idx=None, include_dust=False):

        if transition is None and idx is None:
            transition = self.species
            idx = range(len(self.species))
        elif isinstance(transition, (tuple, list, np.ndarray)):
            idx = (list(self.species).index(t) for t in transition)
        elif isinstance(idx, (tuple, list, np.ndarray)):
            transition = (self.species[i] for i in idx)
        elif isinstance(transition, str) and transition in self.species:
            transition = (transition, )
            idx = (list(self.species).index(transition[0]), )
        elif isinstance(idx, int):
            transition = (self.species[idx], )
            idx = (idx, )
        else:
            print('Please enter a valid value for transition or idx')
            return

        optical_depth = []

        for i in idx:
            if include_dust:
                optical_depth.append(self.optical_depth_species[:, :, :, i])
            else:
                # if len(self.dust) > 10:
                #     wav_dust = self.get_dust_wavelengths()
                #     f = interp1d(self.optical_depth_dust, wav_dust, axis=2, kind='slinear', 
                #                  fill_value='extrapolate')
                #     optical_depth_dust = f(wav_species[i])
                # else:
                #     optical_depth_dust = self.optical_depth_dust.mean(2)
                optical_depth_dust = self.optical_depth_species[:, :, :, i].min(0)
                optical_depth.append(self.optical_depth_species[:, :, :, i] - optical_depth_dust)

        if len(optical_depth) > 1:
            return np.array(optical_depth)
        else:
            return np.array(optical_depth[0])

    def get_hi_intensity(self, include_dust=False, integrated=False):

        if include_dust:
            intensity_temp = self.hi_intensity_species[:, :, :, 0]
        else:
            # if len(self.dust) > 10:
            #     wav_dust = self.get_dust_wavelengths()
            #     f = interp1d(self.intensity_dust, wav_dust, axis=2, kind='slinear', 
            #                  fill_value='extrapolate')
            #     intensity_dust = f(wav_species[i])
            # else:
            #     intensity_dust = self.intensity_dust.mean(2)
            intensity_dust = self.hi_intensity_species[:, :, :, 0].min(0)
            intensity_temp = self.hi_intensity_species[:, :, :, 0] - intensity_dust
        if integrated:
            intensity = np.trapz(intensity_temp, self.map_vel, axis=0)
        else:
            intensity = deepcopy(intensity_temp)

        return intensity

    def get_hi_optical_depth(self, include_dust=False):

        if include_dust:
            optical_depth = self.hi_optical_depth_species[:, :, :, 0]
        else:
            # if len(self.dust) > 10:
            #     wav_dust = self.get_dust_wavelengths()
            #     f = interp1d(self.intensity_dust, wav_dust, axis=2, kind='slinear', 
            #                  fill_value='extrapolate')
            #     intensity_dust = f(wav_species[i])
            # else:
            #     intensity_dust = self.intensity_dust.mean(2)
            optical_depth_dust = self.hi_optical_depth_species[:, :, :, 0].min(0)
            optical_depth = self.hi_optical_depth_species[:, :, :, 0] - optical_depth_dust

        return optical_depth

    def get_dust_intensity(self, wavelength=None, idx=None):

        if wavelength is None and idx is None:
            wavelength = self.dust
            idx = range(len(self.dust))
        elif isinstance(wavelength, (tuple, list, np.ndarray)):
            idx = (list(self.dust).index(t) for t in wavelength)
        elif isinstance(idx, (tuple, list, np.ndarray)):
            wavelength = (self.dust[i] for i in idx)
        elif isinstance(wavelength, str) and wavelength in self.dust:
            wavelength = (wavelength, )
            idx = (list(self.dust).index(wavelength[0]), )
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
            return np.array(intensity)
        else:
            return np.array(intensity[0])

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
            return np.array(optical_depth)
        else:
            return np.array(optical_depth[0])

    def plot_model_quantity(self, quantity=None, transition=None, transition2=None, ens=None, log=False, 
                            stat='max', include_dust=False, integrated=False,
                            cmap_kwargs={'cmap':'magma', 'marker':'s', 'alpha':0.8, 's':27},
                            cbar_kwargs={}, verbose=False, **kwargs):
        '''
        Return a 3DAxes with the specitied quantity shown in the colour scale.
        '''

        mpl.rcParams['text.usetex'] = False
        mpl.rcParams['font.family'] = 'Nimbus Roman'
        plt.rcParams['mathtext.fontset'] = 'stix'

        if verbose:
            print('cmap')
            print(cmap_kwargs)
            print()
            print('cbar')
            print(cbar_kwargs)
            print()

        if not quantity:
            print('Please specify a property to plot')

        if quantity == 'intensity' and transition in self.species:
            # transition2 = transition[1]
            transition = transition
            value = self.get_model_species_intensity(transition=transition, 
                                                     include_dust=include_dust, 
                                                     integrated=integrated)
            if integrated:
                clabel = r'$\varpi_\nu$ (K km s$^{-1}$)'
            else:
                clabel = r'$I_\nu$ (K)'
        elif quantity == 'emissivity' and transition in self.species:
            value = self.species_emissivity[:, :, self.species==transition]
            clabel = r'$\epsilon_\nu$ (K pc$^{-1}$)'
        elif quantity == 'absorption' and transition in self.species:
            value = self.species_absorption[:, :, self.species==transition]
            clabel = r'$\kappa_\nu$ (pc$^{-1}$)'
        if quantity == 'intensity' and transition in self.dust:
            # transition2 = transition[1]
            transition = transition
            value = self.get_model_dust_intensity(transition=transition)
            clabel = r'$I_\nu$ (K)'
        elif quantity == 'emissivity' and transition in self.dust:
            value = self.species_emissivity[:, :, self.dust==transition]
            clabel = r'$\epsilon_\nu$ (K pc$^{-1}$)'
        elif quantity == 'absorption' and transition in self.dust:
            value = self.species_absorption[:, :, self.dust==transition]
            clabel = r'$\kappa_\nu$ (pc$^{-1}$)'
        if quantity == 'intensity' and transition == 'HI' and self.hi_model:
            # transition2 = transition[1]
            transition = transition
            value = self.get_model_hi_intensity(include_dust=include_dust, 
                                                integrated=integrated)
            if integrated:
                clabel = r'$\varpi_\nu$ (K km s$^{-1}$)'
            else:
                clabel = r'$I_\nu$ (K)'
        elif quantity == 'emissivity' and transition == 'HI' and self.hi_model:
            value = self.hi_emissivity[:, :, 0]
            clabel = r'$\epsilon_\nu$ (K pc$^{-1}$)'
        elif quantity == 'absorption' and transition == 'HI' and self.hi_model:
            value = self.hi_absorption[:, :, 0]
            clabel = r'$\kappa_\nu$ (pc$^{-1}$)'
        elif quantity == 'FUV' or quantity == 'fuv':
            value = self.fuv[:, 0]
            clabel = r'$\chi$ ($\chi_\mathrm{D}$)'
        elif quantity == 'dispersion' or quantity == 'sigma':
            value = self.ensemble_dispersion[:, 0]
            clabel = r'$\sigma_\mathrm{ens}$ (km s$^{-1}$)'
        elif quantity == 'f_vox' or quantity == 'voxel-filling factor':
            value = self.f_vox[:, 0]
            clabel = r'$f_\mathrm{vox}$'
        elif quantity == 'f_vol' or quantity == 'volume-filling factor':
            value = self.get_volume_filling_factor()
            clabel = r'$f_\mathrm{voxel}$'
        elif quantity == 'm_h' or quantity == 'atomic mass':
            if isinstance(ens, int):
                value = self.hi_mass[:, ens]
            else:
                value = self.hi_mass.sum(1)
            clabel = r'$M_\mathrm{H^0}$ ($M_\odot$)'
        elif quantity == 'm_h2' or quantity == 'molecular mass':
            if isinstance(ens, int):
                value = self.h2_mass[:, ens]
            else:
                value = self.h2_mass.sum(1)
            clabel = r'$M_\mathrm{H^2}$ ($M_\odot$)'
        else:
            print('Quantity not available.')
            return

        if (quantity in ['emissivity', 'absorption']) or (quantity == 'intensity' and integrated is False):
            if stat == 'std' or stat == 'sigma':
                value = value.std(1)
            elif stat == 'mean':
                value = value.mean(1)
            elif stat == 'max':
                value = value.max(1)
            elif stat == 'min':
                value = value.min(1)
            else:
                print(f'{stat} not available')
                return

        if log:
            value = np.log10(value)
            clabel = r'log$_{10}$ ' + clabel

        X, Y, Z = self.position.T
        lims = (X.min(), X.max())

        fig = plt.figure(figsize=(10, 8))
        ax = fig.add_subplot(111, projection='3d')

        cm = ax.scatter(X, Y, Z, c=value, **cmap_kwargs)
        cb = fig.colorbar(cm, **cbar_kwargs)

        ax.set_xlabel('X (kpc)', fontsize=32)
        ax.set_ylabel('Y (kpc)', fontsize=32)
        ax.set_zlabel('Z (kpc)', fontsize=32)
        cb.ax.set_ylabel(clabel, fontsize=32)
        ax.tick_params(labelsize=16)
        ax.set_xlim(lims)
        ax.set_ylim(lims)
        ax.set_zlim(lims)
        ax.view_init(elev=90, azim=270)

        return ax

    def radial_plot(self, quantity='intensity', transition=['HI'], transition2=None, idx=0, lat=0, 
                    include_dust=False, integrated=False, log=False, scale=False, 
                    ls='-', lw=2, color='xkcd:maroon', label='', 
                    bins=36, bin_lim=(0, 18000), stat='mean'):

        mpl.rcParams['text.usetex'] = False
        mpl.rcParams['font.family'] = 'Nimbus Roman'
        plt.rcParams['mathtext.fontset'] = 'stix'

        species = self.species
        positions = self.position
        rgal = np.sqrt(positions[:, 0]**2+positions[:, 1]**2+positions[:, 2]**2)

        fig, ax = plt.subplots(figsize=(15, 10))

        i_lat = np.where(self.map_lat == lat)[0][0]
        bins = np.linspace(*bin_lim, num=bins)
        bins_mid = bins[:-1] + self.ds/2

        if label == '':
            label = quantity
        
        if scale:
            f_vox = self.f_vox
        else:
            f_vox = 1

        if quantity in ['intensity', 'emissivity', 'absorption']:
            for i, t in enumerate(transition):
                if t == 'HI':
                    eps = self.get_model_hi_emissivity(include_dust=include_dust) / f_vox
                    kap = self.get_model_hi_absorption(include_dust=include_dust) / f_vox
                    intensity = self.get_model_hi_intensity(include_dust=include_dust) / f_vox
                else:
                    eps = self.get_model_species_emissivity(transition=t, include_dust=include_dust) / f_vox
                    kap = self.get_model_species_absorption(transition=t, include_dust=include_dust) / f_vox
                    intensity = self.get_model_species_intensity(transition=t, include_dust=include_dust) / f_vox
                
                # intensity = eps/kap * (1-np.exp(-kap*self.ds))
                
                if quantity == 'emissivity':
                    value = eps.max()
                    ylabel = r'$\epsilon$ (K pc$^{-1}$ kpc$^{-1}$)'
                elif quantity == 'absorption':
                    value = kap.max()
                    ylabel = r'$\kappa$ (pc$^{-1}$ kpc$^{-1}$)'
                elif integrated:
                    value = np.trapz(intensity, self.map_vel, axis=1)
                    ylabel = r'$W$ (K km s$^{-1}$ kpc$^{-1}$)'
                else:
                    value = intensity.max()
                    ylabel = r'$I$ (K kpc$^{-1}$)'

                value_stat,_,_ = binned_statistic(rgal, value, statistic=stat, bins=bins)
                ax.plot(bins_mid, value_stat/self.ds*1e3, label=t)

            ax.legend(fontsize=16)
            ax.tick_params(labelsize=16)
            ax.set_xlabel(r'$R_\mathrm{gal}$ (kpc)', fontsize=24)
            ax.set_ylabel(ylabel,fontsize=24)
            return ax

        elif quantity == 'mass':
            ylabel = r'$M_\mathrm{ens}$ (M$_\odot$ kpc$^{-1}$)'
            #for idx in range(self.ensemble_mass.shape[1]):
            value = f_vox * self.ensemble_mass[:, 0]
            value_stat,_,_ = binned_statistic(rgal, value, statistic=stat, bins=bins)
            label = f'ens 0 mass'
            ax.semilogy(bins_mid, value_stat, ls='--', lw=3, label=label)
            value = f_vox * self.ensemble_mass[:, 1]
            value_stat,_,_ = binned_statistic(rgal, value, statistic=stat, bins=bins)
            label = f'ens 1 mass'
            ax.semilogy(bins_mid, value_stat, ls='--', lw=3, label=label)
            value = f_vox * self.hi_mass.sum(1)
            value_stat,_,_ = binned_statistic(rgal, value, statistic=stat, bins=bins)
            label = f'H mass'
            ax.semilogy(bins_mid, value_stat, lw=1, label=label)
            value = f_vox * self.h2_mass.sum(1)
            value_stat,_,_ = binned_statistic(rgal, value, statistic=stat, bins=bins)
            label = f'H2 mass'
            ax.semilogy(bins_mid, value_stat, lw=1, label=label)
            ax.legend(fontsize=16)
            ax.tick_params(labelsize=16)
            ax.set_xlabel(r'$R_\mathrm{gal}$ (kpc)', fontsize=24)
            ax.set_ylabel(ylabel, fontsize=24)
            return ax

        elif quantity == 'ensemble mass':
            value = []
            label_suffix = [' total', ' clump', ' interclump']
            value.append(f_vox * self.ensemble_mass.sum(1))
            value.append(f_vox * self.ensemble_mass[:, 0])
            value.append(f_vox * self.ensemble_mass[:, 1])
            ylabel = r'$M_\mathrm{ens}$ (M$_\odot$ kpc$^{-1}$)'

        elif quantity == 'hi mass':
            value = []
            label_suffix = [' total', ' clump', ' interclump']
            value.append(f_vox * self.hi_mass.sum(1))
            value.append(f_vox * self.hi_mass[:, 0])
            value.append(f_vox * self.hi_mass[:, 1])
            ylabel = r'$M_\mathrm{HI}$ (M$_\odot$ kpc$^{-1}$)'

        elif quantity == 'h2 mass':
            value = []
            label_suffix = [' total', ' clump', ' interclump']
            value.append(f_vox * self.h2_mass.sum(1))
            value.append(f_vox * self.h2_mass[:, 0])
            value.append(f_vox * self.h2_mass[:, 1])
            ylabel = r'$M_\mathrm{H_2}$ (M$_\odot$ kpc$^{-1}$)'

        elif quantity == 'ensemble density':
            value = self.density[:, idx]
            ylabel = r'n_\mathrm{ens, '+f'{idx}'+'} (cm$^{-3}$ kpc$^{-1}$)'

        elif quantity == 'ensemble FUV':
            value = self.fuv[:, idx]
            ylabel = r'u_\mathrm{FUV, '+f'{idx}'+'} ($\chi_\mathrm{D}$ kpc$^{-1}$)'

        else:
            print(f'quantity {quantity} not available...')
            exit()

        if log:
            value = np.log10(value)
            ylabel = r'$\mathrm{log}_{10}$ ' + ylabel

        if 'mass' in quantity:
            for i, v in enumerate(value):
                value_stat,_,_ = binned_statistic(rgal, v, statistic=stat, bins=bins)
                ax.plot(bins_mid/1000, value_stat/self.ds*1e3, lw=lw[i], ls=ls[i], color=color[i], label=label+label_suffix[i])
        else:
            value_stat,_,_ = binned_statistic(rgal, value, statistic=stat, bins=bins)
            ax.plot(bins_mid/1000, value_stat/self.ds*1e3, lw=lw, ls=ls, color=color, label=label)
        ax.legend(fontsize=36)
        ax.tick_params(labelsize=36)
        ax.set_xlabel(r'$R_\mathrm{gal}$ (kpc)', fontsize=42)
        ax.set_ylabel(ylabel,fontsize=42)

        return ax

