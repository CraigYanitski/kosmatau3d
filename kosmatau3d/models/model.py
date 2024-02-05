"""
A submodule containing the classes :code:`Model` and :code:`SyntheticModel`,
which are used to create and load large models containing many voxels, 
respectively.

:code:`Model()`
-------------

For creating models, the process is to first create an instance with the 
specified parameters, then calculate the emission in all of the voxels.
All of the model data will be saved into its own directory, so be sure to give 
a unique name to the :code:`folder` kwarg.

.. code-block:: python

   from kosmatau3d import models

   model = models.Model() #create model with default kwargs
   model.calculateModel(kind='linear') #calculate emission in voxels

This method of executing the code will save much of the model data, including 
both the intrinsic voxel data and the computed emission data.
The synthetic observations can then be computed using the 
:code:`radiativeTransfer` submodule:

.. code-block:: python

   models.radiativeTransfer.calculateObservation(**kwargs)

Note that the evaluation of KOSMA-:math:`\tau` line emission and HI line
emission is currently separated, and one can which between these modes of 
operation using the :code:`hi` kwarg.

:code:`SyntheticModel()`
----------------------

For loading the models into memory, one must first create a 
:code:`SyntheticModel` instance, then load the model directory.
From there, one can access all of the information of the model (from voxel
masses, densities, abundances, emission, synthetic observations, etc. in 
addition to methods that can be used to plot the data).

.. code-block:: python

   from kosmatau3d import models

   model = models.SyntheticModel(**kwargs)
   model.load_model(directory="/path/to/model/")
"""

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
from .voxel_grid import VoxelGrid


class Model(object):
    """
    This is the highest class in the hierarchy of the :code:`kosmatau3d`
    simulation.
    It contains all of the information needed to properly model a PDR.
    """

    # PRIVATE

    def __init__(
        self,
        history_path="",
        directory="",
        folder="",
        clump_tau_grid_file="clump_tau_LineCenter.dat",
        clump_tb_grid_file="clump_Tmb_LineCenter.dat",
        clump_taufuv_grid_file="RhoMassAFUV.dat",
        clump_column_density_file="clumpMeanCols.dat",
        clump_temperature_file="clumpTemperatures_filled.dat",
        interclump_tau_grid_file="interclumpTauLineCenter.dat",
        interclump_dust_tau_grid_file="interclumpDustTau.dat",
        interclump_tb_grid_file="interclumpTmbLineCenter.dat",
        interclump_dust_tb_grid_file="interclumpDustSED.dat",
        interclump_taufuv_grid_file="interclumpTauFUV.dat",
        interclump_column_density_file="interclumpMeanCols.dat",
        interclump_temperature_file="interclumpTemperatures_filled.dat",
        h2_surface_density_file="h2_surface-density.dat",
        hi_surface_density_file="hi_surface-density.dat",
        h2_scale_height_file="h2_scale-height.dat",
        hi_scale_height_file="hi_scale-height.dat",
        h_number_density_file="h_number-density.dat",
        fuv_file="galactic_FUV_complete.dat",
        l_range=(912, 2066),
        average_fuv=False,
        scale_gc=1.0,
        mhi_gc=1.0,
        mh2_gc=1.0,
        r_gc=4400,
        like_clumps=False,
        all_full=False,
        velocity_file="rot_milki2018_14.dat",
        disp_core=None,
        r_core=4400,
        disp_gmc=None,
        x=0,
        y=0,
        z=0,
        model_type="",
        resolution=1000,
        abundances=["ELECTR", "H", "H2", "H+"],
        transitions="all",
        dust="molecular",
        velocity_range=(),
        velocity_number=0,
        clump_mass_range=((0, 2), (-2)),
        clump_mass_number=(3, 1),
        clump_n_max=(1, 100),
        ensemble_mass_factor=(1, 1),
        interclump_idx=(False, True),
        interclump_wnm_idx=(False, False),
        interclump_hi_ratio=1,
        interclump_wnm_ratio=0.2,
        interclump_wnm_log_fuv=None,
        interclump_f_fuv_wnm=1e4,
        interclump_fillingfactor=None,
        interclump_density=1911,
        interclump_log_fuv=None,
        clump_log_fuv=None,
        hi_mass_factor=1,
        h2_mass_factor=1,
        fuv_factor=1,
        density_factor=1,
        global_uv=10,
        r_cmz=0,
        zeta_cmz=1e-14,
        zeta_sol=2e-16,
        new_grid=True,
        suggested_calc=True,
        dilled=False,
        timed=False,
        verbose=False,
        debug=False,
    ):
        """
        The `Model()` class, which is used to create large three-dimensional
        PDR models and as well as a synthetic observation. For large models, it
        is best to stream the voxel data to the hard disk and avoid keeping the
        data in memory.
        """

        if not len(clump_mass_range):
            sys.exit("<<ERROR>> Define mass sets in argument.")
            # sys.exit()
        if not len(velocity_range):
            sys.exit("<<ERROR>> Define observing velocities in argument.")
            # sys.exit()

        # this just adds a label to the type of model being created.
        # ie 'disk', 'bar', 'sphere', etc.
        constants.type = model_type
        constants.voxel_size = float(resolution)
        constants.HISTORYPATH = history_path
        constants.history = folder
        constants.change_directory(directory)

        # Clump properties
        constants.change_velocity_range(velocity_range)
        constants.change_velocity_number(velocity_number)
        constants.add_clumps(
            mass_range=clump_mass_range,
            num=clump_mass_number,
            n_max=clump_n_max,
            reset=True,
        )
        constants.set_interclump_ensemble(interclump_idx)
        constants.set_interclump_wnm(interclump_wnm_idx)

        # Factors
        constants.ensemble_mass_factor = ensemble_mass_factor
        constants.hi_mass_factor = hi_mass_factor
        constants.h2_mass_factor = h2_mass_factor
        if np.any(constants.interclump_idx) == False:
            constants.interclump_hi_ratio = 0
        else:
            constants.interclump_hi_ratio = interclump_hi_ratio
        if np.any(constants.interclump_wnm_idx) == False:
            constants.interclump_wnm_ratio = 0
        else:
            constants.interclump_wnm_ratio = interclump_wnm_ratio
        constants.interclump_fillingfactor = interclump_fillingfactor
        constants.density_factor = density_factor
        constants.interclump_density = interclump_density
        constants.fuv_factor = fuv_factor
        # constants.globalUV = globalUV
        constants.clump_log_fuv = clump_log_fuv
        constants.interclump_log_fuv = interclump_log_fuv
        constants.interclump_wnm_log_fuv = interclump_wnm_log_fuv
        constants.interclump_f_fuv_wnm = interclump_f_fuv_wnm
        constants.r_cmz = r_cmz
        constants.zeta_cmz = zeta_cmz
        constants.zeta_sol = zeta_sol

        # Read grid & input data, specify transitions, and interpolate
        if "ELECTR" in abundances:
            abundances.remove("ELECTR")
        if "H" in abundances:
            abundances.remove("H")
        if "H2" in abundances:
            abundances.remove("H2")
        if "H+" in abundances:
            abundances.remove("H+")
        abun = ["ELECTR", "H", "H2", "H+"]
        for sp in abundances:
            abun.append(sp)
        constants.abundances = copy(abun)
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
        observations.methods.initialise_grid(
            clump_tau_grid_file=clump_tau_grid_file,
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
            interclump_temperature_file=interclump_temperature_file,
        )
        constants.h2_surface_density_file = h2_surface_density_file
        constants.hi_surface_density_file = hi_surface_density_file
        constants.h2_scale_height_file = h2_scale_height_file
        constants.hi_scale_height_file = hi_scale_height_file
        constants.h_number_density_file = h_number_density_file
        constants.fuv_file = fuv_file
        constants.velocity_file = velocity_file
        observations.methods.initialise_input(
            h2_surface_density_file=h2_surface_density_file,
            hi_surface_density_file=hi_surface_density_file,
            h2_scale_height_file=h2_scale_height_file,
            hi_scale_height_file=hi_scale_height_file,
            h_number_density_file=h_number_density_file,
            fuv_file=fuv_file,
            velocity_file=velocity_file,
        )
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
            constants.disp_gmc = 1.1 * constants.voxel_size**0.38
        if not interpolations.initialised or new_grid:
            interpolations.initialise_grid(dilled=dilled)
            constants.average_fuv = average_fuv
            constants.l_range = l_range
            constants.like_clumps = like_clumps
            constants.all_full = all_full
            interpolations.initialise_model(
                l_range=l_range,
                average_fuv=average_fuv,
                all_full=all_full,
                like_clumps=like_clumps,
            )

        # Initialise logger
        self.__logger = getLogger()

        # Shape() object to create the parameters for the grid of voxels
        self.__shape = shape.Shape(x, y, z, modelType=model_type)
        # VoxelGrid() object to build the model and calculate the emission
        self.__grid = VoxelGrid(self.__shape, suggested_calc=suggested_calc)
        # Orientation() object to change the viewing angle and expected spectra
        # self.__orientation = Orientation(self.__shape.getDimensions())
        self.__speciesNames = (
            []
        )  # this is a list of the species names for easy printout
        self.__timed = timed
        self.__verbose = verbose
        self.__debug = debug
        self.__intensityMap = []
        self.__mapPositions = []
        return

    def __add_transitions(self, species_transition):
        """
        Add transition(s) to the list of transitions in model output.
        Note that this is for the transitions used in radiative transfer.
        """
        species.add_transitions(species_transition)
        species.add_transitions(species_transition, interclump=True)
        return

    def __str__(self):
        """
        Print information of the model. This will become useful for smaller
        models that can be kept in memory.
        """
        printout = "A {} model of {} voxels".format(
            constants.type, self.getGrid().getVoxelNumber()
        )
        if self.__verbose:
            printout += "\n  arranged in {}".format(self.__shape.getDimensions())
            printout += "\n\nConsidering {} species:\n{}\n{}".format(
                len(self.species_names), self.__molecules, self.__dust
            )
        emission = self.__grid.totalEmission()
        printout += "\n\nTotal intensity: {}\nTotal optical depth: {}".format(
            emission[0].sum(), np.log(np.exp(emission[1]).sum())
        )
        return printout

    # PUBLIC

    def getType(self):
        """
        Return type of model. Currently only the 'disk' models are working.
        """
        return constants.type

    def getShape(self):
        """
        Return `Shape` instance of model.

        This will soon be moved from a class to a subpackage...
        """
        return self.__shape

    def getGrid(self):
        """
        Return `voxelGrid` instance of the model. It contains all of the
        information regarding voxels used in the model.
        """
        return self.__grid

    def getOrientation(self):
        """
        Return orientation of model. This is important for integrating the
        radiative transfer equation.
        """
        return self.__orientation

    # def getObservations(self):
    #   return self.__observations

    def getSpecies(self):
        """
        Return list of species names.
        """
        return species.species_names

    def getSpeciesNames(self):
        """
        Same as `getSpecies()`.
        """
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
        """
        Compute voxel the emissivity and absorption spectra for each voxel.
        """
        # Point logger to file
        format_str = "\n\n%(levelname)s [%(name)s]: %(message)s\n\n"
        filename = (
            constants.HISTORYPATH + constants.directory + constants.history + "log.txt"
        )
        filehandler = FileHandler(filename, mode="w")
        if self.__debug:
            basicConfig(format=format_str, level="DEBUG")
            # self.__logger.setLevel('DEBUG')
            # filehandler.setLevel('DEBUG')
        elif self.__verbose:
            basicConfig(format=format_str, level="INFO")
            # self.__logger.setLevel('INFO')
            # filehandler.setLevel('INFO')
        else:
            basicConfig(format=format_str, level="WARNING")
            # self.__logger.setLevel('WARNING')
            # filehandler.setLevel('WARNING')
        # filehandler.setformatter(Formatter(format_str))
        # self.__logger.addHandler(filehandler)

        # Calculate emission
        self.__grid.calculateEmission(**kwargs)

        return

    def writeEmission(self):
        """
        Write voxel spectra to hard disk. This is only useful when keeping
        the model in memory.
        """
        self.__grid.writeEmission(verbose=self.__verbose, debug=self.__debug)
        return

    def getIntensityMap(self):
        """
        Return synthetic intensity. This is only useful when keeping the
        model in memory.
        """
        return (self.__mapPositions, self.__intensityMap)

    def printIntensityMap(self, index=None):
        """
        Print synthetic intensity data to screen.
        """

        print(self.__species[0].getMolecules(), self.__species[1].getDust())

        if not index == None:
            position = self.__mapPositions[index]
            intensity = self.__intensityMap[index]

            print(
                "At position x={} y={}, the intensity is".format(
                    position[0], position[1]
                )
            )

            for element in range(intensity.shape[0]):
                i = intensity[element].argmax()
                print(
                    "{}: {} centered at {} km/s".format(
                        self.__speciesNames[element],
                        intensity[element][intensity[element].nonzero()],
                        self.__constants.velocityBins[i],
                    )
                )

            print()

        else:
            for index in range(len(self.__mapPositions)):
                position = self.__mapPositions[index]
                intensity = self.__intensityMap[index]

                print(
                    "At position x={} y={}, the intensity is".format(
                        position[0], position[1]
                    )
                )

                for element in range(intensity.shape[0]):
                    i = intensity[element].argmax()
                    print(
                        "{}: {} centered at {} km/s".format(
                            self.__speciesNames[element],
                            intensity[element][intensity[element].nonzero()],
                            self.__constants.velocityBins[i],
                        )
                    )

                print()

        return

    def plotModel(self, plot="total intensity"):
        """
        Outdated method to plot model information. Might be useful when
        keeping model data in memory.
        """
        positions = self.__grid.getVoxelPositions()
        limits = [positions.min(), positions.max()]
        if plot == "total intensity":
            if self.__debug:
                print(self.__grid.totalEmission().shape)
            weights = (self.__grid.totalEmission()[0]).max(2).sum(1)
            plot = r"$I \ (\chi)$"
        elif plot == "total optical depth":
            weights = (self.__grid.totalEmission()[1]).max(2).sum(1)
            plot = r"$\tau$"
        elif plot == "clump intensity":
            weights = (self.__grid.clumpEmission()[0]).max(2).sum(1)
            plot = r"$I \ (\chi)$"
        elif plot == "clump optical depth":
            weights = (self.__grid.clumpEmission()[1]).max(2).sum(1)
            plot = r"$\tau$"
        elif plot == "interclump intensity":
            weights = (self.__grid.interclumpEmission()[0]).max(2).sum(1)
            plot = r"$I \ (\chi)$"
        elif plot == "interclump optical depth":
            weights = (self.__grid.interclumpEmission()[1]).max(2).sum(1)
            plot = r"$\tau$"
        elif plot == "FUV":
            weights = self.__grid.getFUV()
            plot = r"$FUV \ (\chi)$"
        elif plot == "Afuv":
            weights = self.__grid.getAfuv()
            plot = r"$\tau_{FUV}$"
        elif plot == "velocity":
            weights = self.__grid.getVelocity()
            plot = r"$v_{rot} \ (\frac{km}{s})$"
        fig = plt.figure()
        ax = fig.add_subplot(111, projection="3d")
        model = ax.scatter(
            positions[0],
            positions[1],
            positions[2],
            c=weights,
            cmap=plt.cm.hot,
            marker="s",
            s=27,
            alpha=1,
            linewidths=0,
        )
        ax.set_xlim(limits)
        ax.set_ylim(limits)
        ax.set_zlim(limits)
        cbar = plt.colorbar(model)
        ax.set_title("PDR Emission within the Milky Way")
        ax.set_xlabel("X (pc)")
        ax.set_ylabel("Y (pc)")
        ax.set_zlabel("Z (pc)")
        cbar.set_label(plot, rotation=0)
        plt.show()
        return

    def hdu_header(self, name="", units="", filename=None, dim=None, data=None):
        if filename is None:
            return

        header = fits.Header()

        header["SIMPLE"] = (True, "conforms to FITS standard")
        header["BITPIX"] = (-64, "element size")
        header["NAXIS"] = (len(dim), "number of axes")
        header["EXTEND"] = True

        for i in range(len(dim)):
            header["NAXIS{}".format(i + 1)] = dim[i]

        # header['NAXIS1'] = self.__constants.velocityBins.size#emission[0].data[0,0,0,:].size
        # header['NAXIS2'] = len(self.__species[0].getMolecules())+len(self.__species[1].getDust())
        # header['NAXIS3'] = 2#emission[0].data[:,0,0,0].size
        # header['NAXIS4'] = self.__voxelNumber#emission[0].data[0,:,0,0].size

        header["NAME"] = name
        header["UNITS"] = units

        if ".fits" not in filename:
            filename = (
                constants.HISTORYPATH
                + constants.directory
                + "r{}_n{}/".format(constants.resolution, self.__shape.voxelNumber())
                + filename
                + ".fits"
            )
        else:
            filename = constants.HISTORYPATH + constants.directory + filename

        if os.path.exists(filename):
            os.remove(filename)

        hdu = fits.PrimaryHDU(data=data, header=header)
        hdu.writeto(filename, overwrite=True)

        return


class SyntheticModel(object):
    """
    This is a class to load individual :code:`kosmatau3d` models.
    This is merely for the convenience of examining the model information in
    a consistent manner.
    There is an optional argument when initialising to set a base directory,
    which makes it easier to load multiple models in succession.
    Due to the complexity of the `kosmatau3d` models, it is not recommended
    to load multiple models at the same time.
    """

    def __init__(self, base_dir=""):
        """
        This initialises the object along with the base directory.
        The owned objects of :code:`base_dir` and :code:`files` are created.
        :code:`files` can be modified again when loading a model, but for now
        it has the default filenames created with :code:`kosmatau3d`.

        :param base_dir: the base directory to use when loading models.
            Default: `""`.


        """

        self.base_dir = base_dir
        self.dir_path = ""
        self.hi_model = False
        self.hi_map = False
        self.files = {
            "intensity": "synthetic_intensity",
            "optical_depth": "synthetic_optical_depth",
            "hi_intensity": "synthetic_hi_intensity",
            "hi_optical_depth": "synthetic_hi_optical_depth",
            "f_vox": "voxel-filling_factor",
            "species_number": "species_number",
            "t_gas": "voxel_t_gas",
            "t_dust": "voxel_t_dust",
            "dust_absorption": "dust_absorption",
            "dust_emissivity": "dust_emissivity",
            "species_absorption": "species_absorption",
            "hi_emissivity": "hi_emissivity",
            "hi_absorption": "hi_absorption",
            "species_emissivity": "species_emissivity",
            "clump_number": "clump_number",
            "clump_radius": "clump_radius",
            "density": "voxel_density",
            "ensemble_dispersion": "voxel_ensemble_dispersion",
            "ensemble_mass": "voxel_ensemble_mass",
            "hi_mass": "voxel_hi_mass",
            "h2_mass": "voxel_h2_mass",
            "fuv_absorption": "voxel_FUVabsorption",
            "fuv": "voxel_fuv",
            "position": "voxel_position",
            "velocity": "voxel_velocity",
            "los_count": "sightlines",
            "log": "log",
        }

        # initialise files
        self.intensity_file = None
        self.hi_intensity_file = None
        self.optical_depth_file = None
        self.hi_optical_depth_file = None

        self.species_number_file = None
        self.dust_absorption_file = None
        self.dust_optical_depth_file = None
        self.species_absorption_file = None
        self.species_emissivity_file = None
        self.hi_absorption_file = None
        self.hi_optical_depth_file = None

        self.t_gas_file = None
        self.t_dust_file = None
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

    def __open_file(self, filename):
        """
        Open and return model data if file exists, otherwise return NaN.
        """
        if (filename + ".fits") in self.model_files:
            file_data = fits.open(os.path.join(self.dir_path, filename + ".fits"))
        elif (filename + ".csv") in self.model_files:
            file_data = np.loadtxt(os.path.join(self.dir_path, filename + ".csv"))
        elif (filename + ".txt") in self.model_files:
            with open(os.path.join(self.dir_path, filename + ".txt")) as f:
                file_data = f.readlines()
        else:
            file_data = np.nan
        return file_data

    def __hdul_data(self, hdul, idx=0):
        """
        Verify and return data in an HDUList, otherwise return NaN.
        """
        if isinstance(hdul, fits.HDUList):
            if isinstance(idx, list):
                return tuple(hdul[i].data for i in idx)
            else:
                return hdul[idx].data
        elif isinstance(idx, list):
            return tuple(np.nan for i in idx)
        else:
            return np.nan

    def close_files(model, **kwargs):
        """
        Close any open FITS files.
        """

        def wrapper(self, **kwargs):
            if isinstance(self.intensity_file, fits.hdu.hdulist.HDUList):
                self.intensity_file.close()
            if isinstance(self.hi_intensity_file, fits.hdu.hdulist.HDUList):
                self.hi_intensity_file.close()
            if isinstance(self.optical_depth_file, fits.hdu.hdulist.HDUList):
                self.optical_depth_file.close()
            if isinstance(self.hi_optical_depth_file, fits.hdu.hdulist.HDUList):
                self.hi_optical_depth_file.close()

            if isinstance(self.species_number_file, fits.hdu.hdulist.HDUList):
                self.species_number_file.close()
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

            if isinstance(self.t_gas_file, fits.hdu.hdulist.HDUList):
                self.t_gas_file.close()
            if isinstance(self.t_dust_file, fits.hdu.hdulist.HDUList):
                self.t_dust_file.close()
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
            if isinstance(
                self.fuv_absorption_dispersion_file, fits.hdu.hdulist.HDUList
            ):
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
        """
        Change the specified filenames.
        """

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
            if "intensity" in kwarg_keys:
                self.files["intensity"] = kwargs["intensity"]
            if "optical_depth" in kwarg_keys:
                self.files["optical_depth"] = kwargs["optical_depth"]
            if "hi_intensity" in kwarg_keys:
                self.files["hi_intensity"] = kwargs["hi_intensity"]
            if "hi_optical_depth" in kwarg_keys:
                self.files["hi_optical_depth"] = kwargs["hi_optical_depth"]
            if "f_vox" in kwarg_keys:
                self.files["f_vox"] = kwargs["f_vox"]
            if "species_number" in kwarg_keys:
                self.files["species_number"] = kwargs["species_number"]
            if "t_gas" in kwarg_keys:
                self.files["t_gas"] = kwargs["t_gas"]
            if "t_dust" in kwarg_keys:
                self.files["t_dust"] = kwargs["t_dust"]
            if "dust_absorption" in kwarg_keys:
                self.files["dust_absorption"] = kwargs["dust_absorption"]
            if "dust_emissivity" in kwarg_keys:
                self.files["dust_emissivity"] = kwargs["dust_emissivity"]
            if "species_absorption" in kwarg_keys:
                self.files["species_absorption"] = kwargs["species_absorption"]
            if "species_emissivity" in kwarg_keys:
                self.files["species_emissivity"] = kwargs["species_emissivity"]
            if "hi_absorption" in kwarg_keys:
                self.files["hi_absorption"] = kwargs["hi_absorption"]
            if "hi_emissivity" in kwarg_keys:
                self.files["hi_emissivity"] = kwargs["hi_emissivity"]
            if "clump_number" in kwarg_keys:
                self.files["clump_number"] = kwargs["clump_number"]
            if "clump_radius" in kwarg_keys:
                self.files["clump_radius"] = kwargs["clump_radius"]
            if "density" in kwarg_keys:
                self.files["density"] = kwargs["density"]
            if "ensemble_dispersion" in kwarg_keys:
                self.files["ensemble_dispersion"] = kwargs["ensemble_dispersion"]
            if "ensemble_mass" in kwarg_keys:
                self.files["ensemble_mass"] = kwargs["ensemble_mass"]
            if "hi_mass" in kwarg_keys:
                self.files["hi_mass"] = kwargs["hi_mass"]
            if "h2_mass" in kwarg_keys:
                self.files["h2_mass"] = kwargs["h2_mass"]
            if "fuv_absorption" in kwarg_keys:
                self.files["fuv_absorption"] = kwargs["fuv_absorption"]
            if "fuv" in kwarg_keys:
                self.files["fuv"] = kwargs["fuv"]
            if "position" in kwarg_keys:
                self.files["position"] = kwargs["position"]
            if "velocity" in kwarg_keys:
                self.files["velocity"] = kwargs["velocity"]
            if "los_count" in kwarg_keys:
                self.files["los_count"] = kwargs["los_count"]
            if "log" in kwarg_keys:
                self.files["log"] = kwargs["log"]

            return model(self, **kwargs)

        return wrapper

    @change_files
    @close_files
    def load_model(self, directory="", map_units="deg", **kwargs):
        """
        Load all of the data for one model. Any additional information such as
        observing velocities, latitude, and longitude are computed as well.

        **Note** that this can be quite computationally expensive. It is not
        recommended to load multiple models at the same time.

        :param directory: The directory of all of the model. Note that this is
                          appended to `self.base_dir`.
        :param kwargs: optional kwargs to modify the model files used
                       (specified in `self.files`).


        """

        self.dir_path = os.path.join(self.base_dir, directory)
        # Ensure model path exists
        if os.path.exists(self.dir_path):
            # if os.path.exists(self.base_dir + directory + self.files['intensity'] + '.fits'):
            self.files["directory"] = directory
            self.model_files = os.listdir(self.base_dir + directory)
        else:
            raise FileNotFoundError(
                f"Directory {self.base_dir + directory} is not valid!! Check "
                + "to see if `base_dir` was set when initialising or if the "
                + "specified directory was correct."
            )

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
        self.intensity_file = self.__open_file(self.files["intensity"])
        (
            self.map_positions,
            self.intensity_species,
            self.intensity_dust,
        ) = self.__hdul_data(self.intensity_file, idx=[0, 1, 2])
        self.optical_depth_file = self.__open_file(self.files["optical_depth"])
        (self.optical_depth_species, self.optical_depth_dust) = self.__hdul_data(
            self.optical_depth_file, idx=[1, 2]
        )
        self.f_vox_file = self.__open_file(self.files["f_vox"])
        self.f_vox = self.__hdul_data(self.f_vox_file)
        self.hi_intensity_file = self.__open_file(self.files["hi_intensity"])
        (
            self.hi_map_positions,
            self.hi_intensity_species,
            self.hi_intensity_dust,
        ) = self.__hdul_data(self.hi_intensity_file, idx=[0, 1, 2])
        self.hi_optical_depth_file = self.__open_file(self.files["hi_optical_depth"])
        (self.hi_optical_depth_species, self.hi_optical_depth_dust) = self.__hdul_data(
            self.hi_optical_depth_file, idx=[1, 2]
        )
        self.dust_absorption_file = self.__open_file(self.files["dust_absorption"])
        self.dust_absorption = self.__hdul_data(self.dust_absorption_file)
        self.dust_emissivity_file = self.__open_file(self.files["dust_emissivity"])
        self.dust_emissivity = self.__hdul_data(self.dust_emissivity_file)
        self.species_absorption_file = self.__open_file(
            self.files["species_absorption"]
        )
        self.species_absorption = self.__hdul_data(self.species_absorption_file)
        self.species_emissivity_file = self.__open_file(
            self.files["species_emissivity"]
        )
        self.species_emissivity = self.__hdul_data(self.species_emissivity_file)
        self.hi_absorption_file = self.__open_file(self.files["hi_absorption"])
        self.hi_absorption = self.__hdul_data(self.hi_absorption_file)
        self.hi_emissivity_file = self.__open_file(self.files["hi_emissivity"])
        self.hi_emissivity = self.__hdul_data(self.hi_emissivity_file)
        self.species_number_file = self.__open_file(self.files["species_number"])
        self.species_number = self.__hdul_data(self.species_number_file)
        self.t_gas_file = self.__open_file(self.files["t_gas"])
        self.t_gas = self.__hdul_data(self.t_gas_file)
        self.t_dust_file = self.__open_file(self.files["t_dust"])
        self.t_dust = self.__hdul_data(self.t_dust_file)
        self.clump_number_file = self.__open_file(self.files["clump_number"])
        self.clump_number = self.__hdul_data(self.clump_number_file)
        self.clump_radius_file = self.__open_file(self.files["clump_radius"])
        self.clump_radius = self.__hdul_data(self.clump_radius_file)
        self.density_file = self.__open_file(self.files["density"])
        self.density = self.__hdul_data(self.density_file)
        self.ensemble_dispersion_file = self.__open_file(
            self.files["ensemble_dispersion"]
        )
        self.ensemble_dispersion = self.__hdul_data(self.ensemble_dispersion_file)
        self.ensemble_mass_file = self.__open_file(self.files["ensemble_mass"])
        self.ensemble_mass = self.__hdul_data(self.ensemble_mass_file)
        self.hi_mass_file = self.__open_file(self.files["hi_mass"])
        self.hi_mass = self.__hdul_data(self.hi_mass_file)
        self.h2_mass_file = self.__open_file(self.files["h2_mass"])
        self.h2_mass = self.__hdul_data(self.h2_mass_file)
        self.fuv_absorption_file = self.__open_file(self.files["fuv_absorption"])
        self.fuv_absorption = self.__hdul_data(self.fuv_absorption_file)
        self.fuv_file = self.__open_file(self.files["fuv"])
        self.fuv = self.__hdul_data(self.fuv_file)
        self.position_file = self.__open_file(self.files["position"])
        self.position = self.__hdul_data(self.position_file)
        self.velocity_file = self.__open_file(self.files["velocity"])
        self.velocity = self.__hdul_data(self.velocity_file)
        self.los_count = self.__open_file(self.files["los_count"])
        self.log = self.__open_file(self.files["log"])

        # Extract headers and create additional axes
        if isinstance(self.species_absorption_file, fits.HDUList):
            self.info = self.species_absorption_file[0].header["COMMENT"]
            self.ds = float(self.info[1].split()[1])  # position of voxel size
        else:
            self.info = """N/A"""
            self.ds = np.nan
        if isinstance(self.species_number_file, fits.HDUList):
            self.N_species = self.species_number_file[0].header["SPECIES"].split(", ")
        else:
            self.N_species = np.nan
        if isinstance(self.intensity_file, fits.HDUList) and isinstance(
            self.optical_depth_file, fits.HDUList
        ):
            self.species_header = self.intensity_file[1].header
            self.species_header["BUNIT"] = (
                self.intensity_file[1].header["BUNIT"]
                + "/"
                + self.optical_depth_file[1].header["BUNIT"]
            )
            self.dust_header = self.intensity_file[2].header
            self.dust_header["BUNIT"] = (
                self.intensity_file[2].header["BUNIT"]
                + "/"
                + self.optical_depth_file[2].header["BUNIT"]
            )
            self.species = np.array(self.species_header["SPECIES"].split(", "))
            self.dust = np.array(self.dust_header["DUST"].split(", "))
            self.dust_header["BUNIT"] = (
                self.intensity_file[2].header["BUNIT"]
                + "/"
                + self.optical_depth_file[2].header["BUNIT"]
            )
            self.map_lon = np.linspace(
                self.species_header["CRVAL2"]
                - self.species_header["CDELT2"] * (self.species_header["CRPIX2"] - 0.5),
                self.species_header["CRVAL2"]
                + self.species_header["CDELT2"]
                * (self.species_header["NAXIS2"] - self.species_header["CRPIX2"] - 0.5),
                num=self.species_header["NAXIS2"],
            )
            self.map_lat = np.linspace(
                self.species_header["CRVAL3"]
                - self.species_header["CDELT3"] * (self.species_header["CRPIX3"] - 0.5),
                self.species_header["CRVAL3"]
                + self.species_header["CDELT3"]
                * (self.species_header["NAXIS3"] - self.species_header["CRPIX3"] - 0.5),
                num=self.species_header["NAXIS3"],
            )
            self.map_vel = np.linspace(
                self.species_header["CRVAL4"]
                - self.species_header["CDELT4"] * (self.species_header["CRPIX4"]),
                self.species_header["CRVAL4"]
                + self.species_header["CDELT4"]
                * (self.species_header["NAXIS4"] - self.species_header["CRPIX4"] - 1),
                num=self.species_header["NAXIS4"],
            )
        else:
            self.species_header = np.nan
            self.dust_header = np.nan
            self.map_lon = np.nan
            self.map_lat = np.nan
            self.map_vel = np.nan

        # convert from radians to degrees if specified
        # (rounding removes floating point error)
        if map_units == "deg" or map_units == "degrees":
            self.map_lon = np.round(self.map_lon * 180 / np.pi, decimals=8)
            self.map_lat = np.round(self.map_lat * 180 / np.pi, decimals=8)

        return

    def get_dust_wavelengths(self):
        """
        Return dust wavelengths used in model.
        """
        wav = list(map(u.Unit, self.dust))
        return list([w.decompose() for w in wav])

    def get_volume_filling_factor(self):
        """
        Return volume filling factor data.
        """
        if isinstance(self.clump_number, np.ndarray) and isinstance(
            self.clump_radius, np.ndarray
        ):
            return (4 / 3 * np.pi * self.clump_number * self.clump_radius**3).sum(
                1
            ) / self.ds**3
        else:
            raise TypeError

    def get_species_number(
        self, species=None, abun=False, nref=[("H", 1), ("H2", 2)], total=True
    ):
        """
        Return species number data. This is the number of a given species contained in each voxel.
        """
        if species in [None, "all"]:
            species = self.N_species
        elif isinstance(species, str):
            species = [species]

        if abun:
            N_0 = 0
            for sp in nref:
                N_0 += sp[1] * self.species_number[:, :, self.N_species.index(sp[0])]
        else:
            N_0 = 1
        N_species = []
        for sp in species:
            N_species.append(self.species_number[:, :, self.N_species.index(sp)] / N_0)
        if total and len(species) == 1:
            return np.sum(N_species, axis=2)[0]
        elif total:
            return np.sum(N_species, axis=2)
        elif len(species) == 1:
            return N_species[0]
        else:
            return N_species

    def get_abundances(self, *args, **kwargs):
        """
        Return species abundance data.
        """
        return self.get_species_number(*args, **kwargs, abun=True)

    def get_gas_temperature(self, total=True):
        """
        Return gas temperature data.
        """
        if total:
            return (self.ensemble_mass * self.t_gas).sum(1) / self.ensemble_mass.sum(1)
        else:
            return copy(self.t_gas)

    def get_dust_temperature(self, total=True):
        """
        Return dust temperature data.
        """
        if total:
            return (self.ensemble_mass * self.t_dust).sum(1) / self.ensemble_mass.sum(1)
        else:
            return copy(self.t_dust)

    def get_model_species_emissivity(
        self, transition=None, idx=None, include_dust=False
    ):
        """
        Return species emissivity in model (ie. for each voxel).
        """
        if transition is None and idx is None:
            transition = self.species
            idx = range(len(self.species))
        elif isinstance(transition, (tuple, list, np.ndarray)):
            idx = (list(self.species).index(t) for t in transition)
        elif isinstance(idx, (tuple, list, np.ndarray)):
            transition = (self.species[i] for i in idx)
        elif isinstance(transition, str) and transition in self.species:
            transition = (transition,)
            idx = (list(self.species).index(transition[0]),)
        elif isinstance(idx, int):
            transition = (self.species[idx],)
            idx = (idx,)
        else:
            print("Please enter a valid value for transition or idx")
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

    def get_model_species_absorption(
        self, transition=None, idx=None, include_dust=False
    ):
        """
        Return species absorption in model (ie. for each voxel).
        """

        if transition is None and idx is None:
            transition = self.species
            idx = range(len(self.species))
        elif isinstance(transition, (tuple, list, np.ndarray)):
            idx = (list(self.species).index(t) for t in transition)
        elif isinstance(idx, (tuple, list, np.ndarray)):
            transition = (self.species[i] for i in idx)
        elif isinstance(transition, str) and transition in self.species:
            transition = (transition,)
            idx = (list(self.species).index(transition[0]),)
        elif isinstance(idx, int):
            transition = (self.species[idx],)
            idx = (idx,)
        else:
            print("Please enter a valid value for transition or idx")
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

    def get_model_species_intensity(
        self, transition=None, idx=None, include_dust=False, integrated=False
    ):
        """
        Return species intensity in model (ie. for each voxel).
        """

        eps = self.get_model_species_emissivity(
            transition=transition, idx=idx, include_dust=include_dust
        )
        kap = self.get_model_species_absorption(
            transition=transition, idx=idx, include_dust=include_dust
        )
        intensity = np.zeros_like(eps)
        i_nan = kap == 0
        intensity[~i_nan] = (
            eps[~i_nan] / kap[~i_nan] * (1 - np.exp(-kap[~i_nan] * self.ds))
        )
        intensity[i_nan] = eps[i_nan]

        if integrated and (intensity.ndim == 3):
            return np.trapz(intensity, self.map_vel, axis=2)
        elif integrated and (intensity.ndim == 2):
            return np.trapz(intensity, self.map_vel, axis=1)
        else:
            return np.array(intensity)

    def get_model_hi_emissivity(self, include_dust=False):
        """
        Return HI 21cm line emissivity in model (ie. for each voxel).
        """

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
        """
        Return HI 21cm line absorption in model (ie. for each voxel).
        """

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
        """
        Return HI 21cm line emissivity in model (ie. for each voxel).
        """

        eps = self.get_model_hi_emissivity(include_dust=include_dust)
        kap = self.get_model_hi_absorption(include_dust=include_dust)
        intensity = np.zeros_like(eps)
        i_nan = kap == 0
        intensity[~i_nan] = (
            eps[~i_nan] / kap[~i_nan] * (1 - np.exp(-kap[~i_nan] * self.ds))
        )
        intensity[i_nan] = eps[i_nan]

        return intensity

    def get_model_dust_emissivity(self, wavelength=None, idx=None):
        """
        Return dust emissivity in model (ie. for each voxel).
        """

        if wavelength is None and idx is None:
            wavelength = self.dust
            idx = range(len(self.dust))
        elif isinstance(wavelength, (tuple, list, np.ndarray)):
            idx = (self.dust.index(t) for t in wavelength)
        elif isinstance(idx, (tuple, list, np.ndarray)):
            wavelength = (self.dust[i] for i in idx)
        elif isinstance(wavelength, str) and wavelength in self.dust:
            wavelength = (wavelength,)
            idx = (list(self.dust).index(wavelength[0]),)
        elif isinstance(idx, int):
            wavelength = (self.dust[idx],)
            idx = (idx,)
        else:
            print("Please enter a valid value for wavelength or idx")
            return

        emissivity = []

        for i in idx:
            emissivity.append(self.dust_emissivity[:, i])

        if len(emissivity) > 1:
            return np.array(emissivity)
        else:
            return np.array(emissivity[0])

    def get_model_dust_absorption(self, wavelength=None, idx=None):
        """
        Return dust absorption in model (ie. for each voxel).
        """

        if wavelength is None and idx is None:
            wavelength = self.dust
            idx = range(len(self.dust))
        elif isinstance(wavelength, (tuple, list, np.ndarray)):
            idx = (self.dust.index(t) for t in wavelength)
        elif isinstance(idx, (tuple, list, np.ndarray)):
            wavelength = (self.dust[i] for i in idx)
        elif isinstance(wavelength, str) and wavelength in self.dust:
            wavelength = (wavelength,)
            idx = (list(self.dust).index(wavelength[0]),)
        elif isinstance(idx, int):
            wavelength = (self.dust[idx],)
            idx = (idx,)
        else:
            print("Please enter a valid value for wavelength or idx")
            return

        absorption = []

        for i in idx:
            absorption.append(self.dust_absorption[:, i])

        if len(absorption) > 1:
            return np.array(absorption)
        else:
            return np.array(absorption[0])

    def get_model_dust_intensity(self, wavelength=None, idx=None):
        """
        Return dust intensity in model (ie. for each voxel).
        """

        eps = self.get_model_dust_emissivity(wavelength=wavelength, idx=idx)
        kap = self.get_model_dust_absorption(wavelength=wavelength, idx=idx)
        intensity = np.zeros_like(eps)
        i_nan = kap == 0
        intensity[~i_nan] = (
            eps[~i_nan] / kap[~i_nan] * (1 - np.exp(-kap[~i_nan] * self.ds))
        )
        intensity[i_nan] = eps[i_nan]

        return np.array(intensity)

    def get_species_intensity(
        self, transition=None, idx=None, include_dust=False, integrated=False
    ):
        """
        Return species intensity in synthetic observation.
        """

        if transition is None and idx is None:
            transition = self.species
            idx = range(len(self.species))
        elif isinstance(transition, (tuple, list, np.ndarray)):
            idx = (list(self.species).index(t) for t in transition)
        elif isinstance(idx, (tuple, list, np.ndarray)):
            transition = (self.species[i] for i in idx)
        elif isinstance(transition, str) and transition in self.species:
            transition = (transition,)
            idx = (list(self.species).index(transition[0]),)
        elif isinstance(idx, int):
            transition = (self.species[idx],)
            idx = (idx,)
        else:
            print("Please enter a valid value for transition or idx")
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
        """
        Return species optical depth in synthetic observation.
        """

        if transition is None and idx is None:
            transition = self.species
            idx = range(len(self.species))
        elif isinstance(transition, (tuple, list, np.ndarray)):
            idx = (list(self.species).index(t) for t in transition)
        elif isinstance(idx, (tuple, list, np.ndarray)):
            transition = (self.species[i] for i in idx)
        elif isinstance(transition, str) and transition in self.species:
            transition = (transition,)
            idx = (list(self.species).index(transition[0]),)
        elif isinstance(idx, int):
            transition = (self.species[idx],)
            idx = (idx,)
        else:
            print("Please enter a valid value for transition or idx")
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
                optical_depth.append(
                    self.optical_depth_species[:, :, :, i] - optical_depth_dust
                )

        if len(optical_depth) > 1:
            return np.array(optical_depth)
        else:
            return np.array(optical_depth[0])

    def get_hi_intensity(self, include_dust=False, integrated=False):
        """
        Return HI 21cm line intensity in synthetic observation.
        """

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
        """
        Return HI 21cm line optical depth in synthetic observation.
        """

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
            optical_depth = (
                self.hi_optical_depth_species[:, :, :, 0] - optical_depth_dust
            )

        return optical_depth

    def get_dust_intensity(self, wavelength=None, idx=None):
        """
        Return dust intensity in synthetic observation.
        """

        if wavelength is None and idx is None:
            wavelength = self.dust
            idx = range(len(self.dust))
        elif isinstance(wavelength, (tuple, list, np.ndarray)):
            idx = (list(self.dust).index(t) for t in wavelength)
        elif isinstance(idx, (tuple, list, np.ndarray)):
            wavelength = (self.dust[i] for i in idx)
        elif isinstance(wavelength, str) and wavelength in self.dust:
            wavelength = (wavelength,)
            idx = (list(self.dust).index(wavelength[0]),)
        elif isinstance(idx, int):
            wavelength = (self.dust[idx],)
            idx = (idx,)
        else:
            print("Please enter a valid value for wavelength or idx")
            return

        intensity = []

        for i in idx:
            intensity.append(self.intensity_dust[:, :, i])

        if len(intensity) > 1:
            return np.array(intensity)
        else:
            return np.array(intensity[0])

    def get_dust_optical_depth(self, wavelength=None, idx=None):
        """
        Return dust optical depth in synthetic observation.
        """

        if wavelength is None and idx is None:
            wavelength = self.dust
            idx = range(len(self.dust))
        elif isinstance(wavelength, (tuple, list, np.ndarray)):
            idx = (self.dust.index(t) for t in wavelength)
        elif isinstance(idx, (tuple, list, np.ndarray)):
            wavelength = (self.dust[i] for i in idx)
        elif isinstance(wavelength, str) and wavelength in self.dust:
            wavelength = (wavelength,)
            idx = (self.dust.index(wavelength[0]),)
        elif isinstance(idx, int):
            wavelength = (self.dust[idx],)
            idx = (idx,)
        else:
            print("Please enter a valid value for wavelength or idx")
            return

        optical_depth = []

        for i in idx:
            optical_depth.append(self.optical_depth_dust[:, :, i])

        if len(optical_depth) > 1:
            return np.array(optical_depth)
        else:
            return np.array(optical_depth[0])

    def plot_model_quantity(
        self,
        quantity=None,
        transition=None,
        transition2=None,
        ens=None,
        log=False,
        stat="max",
        include_dust=False,
        integrated=False,
        vmin=None,
        vmax=None,
        cmap_kwargs={"cmap": "magma", "marker": "s", "alpha": 0.8, "s": 27},
        cbar_kwargs={},
        label_axes=False,
        verbose=False,
        **kwargs,
    ):
        """
        Return a 3DAxes with the specitied quantity shown in the colour scale.
        """

        mpl.rcParams["text.usetex"] = False
        mpl.rcParams["font.family"] = "Nimbus Roman"
        plt.rcParams["mathtext.fontset"] = "stix"

        if verbose:
            print("cmap")
            print(cmap_kwargs)
            print()
            print("cbar")
            print(cbar_kwargs)
            print()

        if not quantity:
            print("Please specify a property to plot")

        if quantity == "intensity" and transition in self.species:
            # transition2 = transition[1]
            transition = transition
            value = self.get_model_species_intensity(
                transition=transition, include_dust=include_dust, integrated=integrated
            )
            if integrated:
                clabel = r"$\varpi_\nu$ (K km s$^{-1}$)"
            else:
                clabel = r"$I_\nu$ (K)"
        elif quantity == "emissivity" and transition in self.species:
            value = self.species_emissivity[:, :, self.species == transition]
            clabel = r"$\epsilon_\nu$ (K pc$^{-1}$)"
        elif quantity == "absorption" and transition in self.species:
            value = self.species_absorption[:, :, self.species == transition]
            clabel = r"$\kappa_\nu$ (pc$^{-1}$)"
        elif quantity == "intensity" and transition in self.dust:
            # transition2 = transition[1]
            transition = transition
            value = self.get_model_dust_intensity(wavelength=transition)
            clabel = r"$I_\nu$ (K)"
        elif quantity == "emissivity" and transition in self.dust:
            value = self.dust_emissivity[:, self.dust == transition]
            clabel = r"$\epsilon_\nu$ (K pc$^{-1}$)"
        elif quantity == "absorption" and transition in self.dust:
            value = self.dust_absorption[:, self.dust == transition]
            clabel = r"$\kappa_\nu$ (pc$^{-1}$)"
        elif quantity == "intensity" and transition == "HI" and self.hi_model:
            # transition2 = transition[1]
            transition = transition
            value = self.get_model_hi_intensity(
                include_dust=include_dust, integrated=integrated
            )
            if integrated:
                clabel = r"$\varpi_\nu$ (K km s$^{-1}$)"
            else:
                clabel = r"$I_\nu$ (K)"
        elif quantity == "emissivity" and transition == "HI" and self.hi_model:
            value = self.hi_emissivity[:, :, 0]
            clabel = r"$\epsilon_\nu$ (K pc$^{-1}$)"
        elif quantity == "absorption" and transition == "HI" and self.hi_model:
            value = self.hi_absorption[:, :, 0]
            clabel = r"$\kappa_\nu$ (pc$^{-1}$)"
        elif quantity == "FUV" or quantity == "fuv":
            if isinstance(ens, int):
                value = self.fuv[:, ens]
            else:
                value = self.fuv[:, 0]
            clabel = r"$\chi$ ($\chi_\mathrm{D}$)"
        elif quantity in ["density", "n"]:
            if isinstance(ens, int):
                value = self.density[:, ens]
            else:
                value = self.density[:, 0]
            clabel = r"$n_\mathrm{ens}$ (cm$^{-3}$)"
        elif quantity == "dispersion" or quantity == "sigma":
            if isinstance(ens, int):
                value = self.ensemble_dispersion[:, ens]
            else:
                value = self.ensemble_dispersion[:, 0]
            clabel = r"$\sigma_\mathrm{ens}$ (km s$^{-1}$)"
        elif quantity in ["t_gas", "gas temperature"]:
            if isinstance(ens, int):
                value = self.get_gas_temperature(total=False)[:, ens]
            else:
                value = self.get_gas_temperature(total=True)
            clabel = r"$T_\mathrm{gas}$ (K)"
        elif quantity in ["t_dust", "dust temperature"]:
            if isinstance(ens, int):
                value = self.get_dust_temperature(total=False)[:, ens]
            else:
                value = self.get_dust_temperature(total=True)
            clabel = r"$T_\mathrm{dust}$ (K)"
        elif quantity == "f_vox" or quantity == "voxel-filling factor":
            value = self.f_vox[:, 0]
            clabel = r"$f_\mathrm{vox}$"
        elif quantity == "f_vol" or quantity == "volume-filling factor":
            value = self.get_volume_filling_factor()
            clabel = r"$f_\mathrm{voxel}$"
        elif quantity == "m_h" or quantity == "atomic mass":
            if isinstance(ens, int):
                value = self.hi_mass[:, ens]
            else:
                value = self.hi_mass.sum(1)
            clabel = r"$M_\mathrm{H^0}$ ($M_\odot$)"
        elif quantity == "m_h2" or quantity == "molecular mass":
            if isinstance(ens, int):
                value = self.h2_mass[:, ens]
            else:
                value = self.h2_mass.sum(1)
            clabel = r"$M_\mathrm{H^2}$ ($M_\odot$)"
        else:
            print("Quantity not available.")
            return

        if quantity == "intensity" and transition2 in self.species:
            # transition2 = transition[1]
            transition = transition
            value2 = self.get_model_species_intensity(
                transition=transition2, include_dust=include_dust, integrated=integrated
            )
            if integrated:
                clabel = r"$\varpi_\nu$ ratio (K km s$^{-1}$)"
            else:
                clabel = r"$I_\nu$ ratio (K)"
        elif quantity == "emissivity" and transition2 in self.species:
            value2 = self.species_emissivity[:, :, self.species == transition2]
            clabel = r"$\epsilon_\nu$ ratio (K pc$^{-1}$)"
        elif quantity == "absorption" and transition2 in self.species:
            value2 = self.species_absorption[:, :, self.species == transition2]
            clabel = r"$\kappa_\nu$ ratio (pc$^{-1}$)"
        elif quantity == "intensity" and transition2 in self.dust:
            # transition2 = transition[1]
            transition = transition
            value2 = self.get_model_dust_intensity(wavelength=transition2)
            clabel = r"$I_\nu$ ratio (K)"
        elif quantity == "emissivity" and transition2 in self.dust:
            value2 = self.dust_emissivity[:, self.dust == transition2]
            clabel = r"$\epsilon_\nu$ ratio (K pc$^{-1}$)"
        elif quantity == "absorption" and transition2 in self.dust:
            value2 = self.dust_absorption[:, self.dust == transition2]
            clabel = r"$\kappa_\nu$ ratio (pc$^{-1}$)"
        elif quantity == "intensity" and transition2 == "HI" and self.hi_model:
            # transition2 = transition[1]
            transition = transition
            value2 = self.get_model_hi_intensity(
                include_dust=include_dust, integrated=integrated
            )
            if integrated:
                clabel = r"$\varpi_\nu$ ratio (K km s$^{-1}$)"
            else:
                clabel = r"$I_\nu$ ratio (K)"
        elif quantity == "emissivity" and transition2 == "HI" and self.hi_model:
            value2 = self.hi_emissivity[:, :, 0]
            clabel = r"$\epsilon_\nu$ ratio (K pc$^{-1}$)"
        elif quantity == "absorption" and transition2 == "HI" and self.hi_model:
            value2 = self.hi_absorption[:, :, 0]
            clabel = r"$\kappa_\nu$ ratio (pc$^{-1}$)"

        if (
            (quantity in ["emissivity", "absorption"])
            or (quantity == "intensity" and integrated is False)
        ) and not (transition in self.dust):
            if stat == "std" or stat == "sigma":
                value = value.std(1)
                if transition2:
                    value = value / value2.std(1)
            elif stat == "mean":
                value = value.mean(1)
                if transition2:
                    value = value / value2.mean(1)
            elif stat == "max":
                value = value.max(1)
                if transition2:
                    value = value / value2.max(1)
            elif stat == "min":
                value = value.min(1)
                if transition2:
                    value = value / value2.min(1)
            else:
                print(f"{stat} not available")
                return

        if (
            (quantity == "intensity" and integrated is True)
            or (transition in self.dust)
        ) and not transition2 == None:
            value = value / value2

        if log:
            value = np.log10(value)
            clabel = r"log$_{10}$ " + clabel

        X, Y, Z = self.position.T
        lims = (X.min(), X.max())

        fig = plt.figure(figsize=(10, 8))
        ax = fig.add_subplot(111, projection="3d")

        cm = ax.scatter(X, Y, Z, c=value, vmin=vmin, vmax=vmax, **cmap_kwargs)
        cb = fig.colorbar(cm, **cbar_kwargs)

        if label_axes:
            ax.set_xlabel("X (kpc)", fontsize=32)
            ax.set_ylabel("Y (kpc)", fontsize=32)
            ax.set_zlabel("Z (kpc)", fontsize=32)
        else:
            ax.set_axis_off()
        cb.ax.set_ylabel(clabel, fontsize=32)
        ax.tick_params(labelsize=16)
        ax.set_xlim(lims)
        ax.set_ylim(lims)
        ax.set_zlim(lims)
        ax.view_init(elev=90, azim=270)

        return ax

    def radial_plot(
        self,
        quantity="intensity",
        transition=["CO 1"],
        transition2=[],
        idx=0,
        lat=0,
        include_dust=False,
        integrated=False,
        log=False,
        scale=False,
        normalized=False,
        ls="-",
        lw=2,
        color="xkcd:maroon",
        label="",
        fontsize=42,
        labelsize=36,
        legendsize=36,
        legendloc=0,
        bins=36,
        bin_lim=(0, 18000),
        stat="mean",
        voxel_size=None,
        ax=None,
    ):
        """
        Plot a given quantity as a function of galactocentric radius.
        This is only valid for galaxy models of type 'disk'.
        """

        mpl.rcParams["text.usetex"] = False
        mpl.rcParams["font.family"] = "Nimbus Roman"
        plt.rcParams["mathtext.fontset"] = "stix"

        species = self.species
        positions = self.position
        rgal = np.sqrt(
            positions[:, 0] ** 2 + positions[:, 1] ** 2 + positions[:, 2] ** 2
        )

        if isinstance(ax, mpl.axes._axes.Axes):
            fig = ax.get_figure()
        else:
            fig, ax = plt.subplots(figsize=(15, 10))

        i_lat = np.where(self.map_lat == lat)[0][0]
        bins = np.linspace(*bin_lim, num=bins)
        bins_mid = bins[:-1] + self.ds / 2

        if label == "":
            label = quantity

        if scale:
            f_vox = self.f_vox
        else:
            f_vox = 1

        if quantity in ["intensity", "optical depth", "emissivity", "absorption"]:
            for i, t in enumerate(transition):
                if t == "HI":
                    eps = (
                        self.get_model_hi_emissivity(include_dust=include_dust) / f_vox
                    )
                    kap = (
                        self.get_model_hi_absorption(include_dust=include_dust) / f_vox
                    )
                    tau = kap * self.ds
                    intensity = (
                        self.get_model_hi_intensity(include_dust=include_dust) / f_vox
                    )
                else:
                    eps = (
                        self.get_model_species_emissivity(
                            transition=t, include_dust=include_dust
                        )
                        / f_vox
                    )
                    kap = (
                        self.get_model_species_absorption(
                            transition=t, include_dust=include_dust
                        )
                        / f_vox
                    )
                    tau = kap * self.ds
                    intensity = (
                        self.get_model_species_intensity(
                            transition=t, include_dust=include_dust
                        )
                        / f_vox
                    )

                if len(transition2) == 1:
                    if transition2[0] == "HI":
                        eps2 = (
                            self.get_model_hi_emissivity(include_dust=include_dust)
                            / f_vox
                        )
                        kap2 = (
                            self.get_model_hi_absorption(include_dust=include_dust)
                            / f_vox
                        )
                        tau2 = kap2 * self.ds
                        intensity2 = (
                            self.get_model_hi_intensity(include_dust=include_dust)
                            / f_vox
                        )
                    else:
                        eps2 = (
                            self.get_model_species_emissivity(
                                transition=transition2, include_dust=include_dust
                            )
                            / f_vox
                        )
                        kap2 = (
                            self.get_model_species_absorption(
                                transition=transition2, include_dust=include_dust
                            )
                            / f_vox
                        )
                        tau2 = kap2 * self.ds
                        intensity2 = (
                            self.get_model_species_intensity(
                                transition=transition2[0], include_dust=include_dust
                            )
                            / f_vox
                        )
                else:
                    eps2 = None
                    kap2 = None
                    intensity2 = None

                # intensity = eps/kap * (1-np.exp(-kap*self.ds))

                if len(transition2) == 0:
                    if quantity == "emissivity":
                        value = eps.max(1)
                        ylabel = r"$\epsilon$ (K pc$^{-1}$ kpc$^{-1}$)"
                    elif quantity == "absorption":
                        value = kap.max(1)
                        ylabel = r"$\kappa$ (pc$^{-1}$ kpc$^{-1}$)"
                    elif quantity == "optical depth":
                        value = tau.max(1)
                        ylabel = r"$\tau$ (kpc$^{-1}$)"
                    elif integrated:
                        value = np.trapz(intensity, self.map_vel, axis=1)
                        ylabel = r"$W$ (K km s$^{-1}$ kpc$^{-1}$)"
                    else:
                        value = intensity.max(1)
                        ylabel = r"$I$ (K kpc$^{-1}$)"
                else:
                    if quantity == "emissivity":
                        value = eps.max(1) / eps2.max(1)
                        ylabel = r"$R_\epsilon$ (K pc$^{-1}$ kpc$^{-1}$)"
                    elif quantity == "absorption":
                        value = kap.max(1) / kap2.max(1)
                        ylabel = r"$R_\kappa$ (pc$^{-1}$ kpc$^{-1}$)"
                    elif quantity == "optical depth":
                        value = tau.max(1) / tau2.max(1)
                        ylabel = r"$R_\tau$ (kpc$^{-1}$)"
                    elif integrated:
                        value = np.trapz(intensity, self.map_vel, axis=1) / np.trapz(
                            intensity2, self.map_vel, axis=1
                        )
                        ylabel = r"$R_W$ (K km s$^{-1}$ kpc$^{-1}$)"
                    else:
                        value = intensity.max(1) / intensity2.max(1)
                        ylabel = r"$R_I$ (K kpc$^{-1}$)"

                if voxel_size == 1:
                    ylabel = ylabel.replace(r" kpc$^{-1}$", "")
                elif voxel_size is None:
                    voxel_size = self.ds / 1e3

                value_stat, _, _ = binned_statistic(
                    rgal, value, statistic=stat, bins=bins
                )
                ax.plot(bins_mid, value_stat / voxel_size, lw=lw, ls=ls, label=t)

            ax.legend(fontsize=legendsize, loc=legendloc)
            ax.tick_params(labelsize=labelsize)
            ax.set_xlabel(r"$R_\mathrm{gal}$ (kpc)", fontsize=fontsize)
            ax.set_ylabel(ylabel, fontsize=fontsize)
            return ax

        elif quantity == "mass":
            ylabel = r"$M_\mathrm{ens}$ (M$_\odot$ kpc$^{-1}$)"
            # for idx in range(self.ensemble_mass.shape[1]):
            value = f_vox * self.ensemble_mass[:, 0]
            value_stat, _, _ = binned_statistic(rgal, value, statistic=stat, bins=bins)
            label = f"ens 0 mass"
            ax.semilogy(bins_mid, value_stat, ls="--", lw=3, label=label)
            value = f_vox * self.ensemble_mass[:, 1]
            value_stat, _, _ = binned_statistic(rgal, value, statistic=stat, bins=bins)
            label = f"ens 1 mass"
            ax.semilogy(bins_mid, value_stat, ls="--", lw=3, label=label)
            value = f_vox * self.hi_mass.sum(1)
            value_stat, _, _ = binned_statistic(rgal, value, statistic=stat, bins=bins)
            label = f"H mass"
            ax.semilogy(bins_mid, value_stat, lw=1, label=label)
            value = f_vox * self.h2_mass.sum(1)
            value_stat, _, _ = binned_statistic(rgal, value, statistic=stat, bins=bins)
            label = f"H2 mass"
            ax.semilogy(bins_mid, value_stat, lw=1, label=label)
            ax.legend(fontsize=16)
            ax.tick_params(labelsize=16)
            ax.set_xlabel(r"$R_\mathrm{gal}$ (kpc)", fontsize=24)
            ax.set_ylabel(ylabel, fontsize=24)
            return ax

        elif quantity == "ensemble mass":
            value = []
            label_suffix = [" total", " clump", " interclump"]
            value.append(f_vox * self.ensemble_mass.sum(1))
            value.append(f_vox * self.ensemble_mass[:, 0])
            value.append(f_vox * self.ensemble_mass[:, 1])
            ylabel = r"$M_\mathrm{ens}$ (M$_\odot$ kpc$^{-1}$)"

        elif quantity == "hi mass":
            value = []
            label_suffix = [" total", " clump", " interclump"]
            mass_cl = self.f_vox[:, 0] * self.hi_mass[:, 0]
            if self.hi_mass.shape[1] == 3:
                f_vox_int = self.f_vox[:, 1] + self.f_vox[:, 2]
                mass_int = f_vox_int * (self.hi_mass[:, 1] + self.hi_mass[:, 2])
            else:
                f_vox_int = self.f_vox[:, 1]
                mass_int = f_vox_int * self.hi_mass[:, 1]
            value.append(mass_cl + mass_int)
            value.append(mass_cl)
            value.append(mass_int)
            ylabel = r"$M_\mathrm{H^0}$ (M$_\odot$ kpc$^{-1}$)"

        elif quantity == "h2 mass":
            value = []
            label_suffix = [" total", " clump", " interclump"]
            mass_cl = self.f_vox[:, 0] * self.h2_mass[:, 0]
            if self.hi_mass.shape[1] == 3:
                f_vox_int = self.f_vox[:, 1] + self.f_vox[:, 2]
                mass_int = f_vox_int * (self.h2_mass[:, 1] + self.h2_mass[:, 2])
            else:
                f_vox_int = self.f_vox[:, 1]
                mass_int = f_vox_int * self.h2_mass[:, 1]
            value.append(mass_cl + mass_int)
            value.append(mass_cl)
            value.append(mass_int)
            ylabel = r"$M_\mathrm{H_2}$ (M$_\odot$ kpc$^{-1}$)"

        elif quantity in ["X_CO", "Xco", "xco"]:
            value = (
                self.get_species_number(species="H2", total=True)
                / (self.ds * constants.pc * 100) ** 2
                / self.get_model_species_intensity(include_dust=False, integrated=True)[
                    list(self.species).index(transition[0]), :
                ]
                / 2e20
            )
            ylabel = r"$X_\mathrm{CO}$ / $X_\mathrm{CO, MW}$"

        elif quantity == "ensemble density":
            value = self.density[:, idx]
            ylabel = r"$n_\mathrm{ens, " + f"{idx}" + "}$ (cm$^{-3}$ kpc$^{-1}$)"

        elif quantity in ["R_sa", "r_sa"]:
            value = np.nanmax(
                (
                    self.get_model_species_emissivity(
                        transition=transition[0], include_dust=include_dust
                    )
                    / self.get_model_species_intensity(
                        transition=transition[0], include_dust=include_dust
                    )
                ),
                1,
            )
            ylabel = r"$R_\mathrm{" + f"{transition[0]}" + r"}$"

        elif quantity == "ensemble FUV":
            value = self.fuv[:, idx]
            ylabel = r"$\chi$ ($\chi_\mathrm{D}$)"

        elif quantity in ["fvox", "f_vox"]:
            ylabel = r"$f_\mathrm{vox}$"

            value = self.f_vox[:, 0]
            value_stat, _, _ = binned_statistic(rgal, value, statistic=stat, bins=bins)
            label = r"H$_\mathrm{2}$"
            ax.plot(
                bins_mid / 1e3,
                value_stat,
                ls="--",
                lw=3,
                color="xkcd:sapphire",
                label=label,
            )

            if self.f_vox.shape[1] == 3:
                value = self.f_vox[:, 1] + self.f_vox[:, 2]
            else:
                value = self.f_vox[:, 1]
            value_stat, _, _ = binned_statistic(rgal, value, statistic=stat, bins=bins)
            label = r"H$^\mathrm{0}$"
            ax.plot(
                bins_mid / 1e3,
                value_stat,
                ls="--",
                lw=3,
                color="xkcd:maroon",
                label=label,
            )

            ax.legend(fontsize=36)
            ax.tick_params(labelsize=36)
            ax.set_xlabel(r"$R_\mathrm{gal}$ (kpc)", fontsize=42)
            ax.set_ylabel(ylabel, fontsize=42)
            return ax

        else:
            print(f"quantity {quantity} not available...")
            exit()

        if log:
            value = np.log10(value)
            ylabel = r"$\mathrm{log}_{10}$ " + ylabel

        if normalized:
            factor = self.ds / 1e3
        else:
            factor = 1
        if "mass" in quantity:
            for i, v in enumerate(value):
                value_stat, _, _ = binned_statistic(rgal, v, statistic=stat, bins=bins)
                ax.plot(
                    bins_mid / 1e3,
                    value_stat / factor,
                    lw=lw[i],
                    ls=ls[i],
                    color=color[i],
                    label=label + label_suffix[i],
                )
        else:
            value_stat, _, _ = binned_statistic(rgal, value, statistic=stat, bins=bins)
            ax.plot(
                bins_mid / 1000,
                value_stat / factor,
                lw=lw,
                ls=ls,
                color=color,
                label=label,
            )
        ax.legend(fontsize=legendsize, loc=legendloc)
        ax.tick_params(labelsize=labelsize)
        ax.set_xlabel(r"$R_\mathrm{gal}$ (kpc)", fontsize=fontsize)
        ax.set_ylabel(ylabel, fontsize=fontsize)

        return ax
