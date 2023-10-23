"""
A module containing the :code:`VoxelGrid` class.
"""

import os
import importlib as il
import gc
from multiprocessing import Pool

import numpy as np
from astropy.io import fits
import matplotlib.pyplot as plt
from tqdm import tqdm

# import progressbar as pb
from numba import jit
from logging import getLogger
from copy import copy
from pprint import pprint
from time import time

from .voxel import *
from kosmatau3d.models import constants, species, interpolations


class VoxelGrid(object):
    """
    This is a class to handle all of the voxels in KOSMA-tau^3. It contains a
    specified arrangement of voxels, and must coordinate with the Dimensions class
    to make the Shape class functional.
    """

    # PRIVATE

    def __init__(self, shape, suggested_calc=True, verbose=False, debug=False):
        """
        Initialise all attributes.
        """

        self.__shape = shape

        self.__voxel_number = self.__shape.voxelNumber()

        self.__voxels = []

        self.__map = (
            {}
        )  # dictionary object to map the voxel indeces to the correct location

        # self.__species = None

        self.__voxel_species_emissivity = []
        self.__voxel_species_absorption = []
        self.__voxel_dust_emissivity = []
        self.__voxel_dust_absorption = []

        self.__voxel_mass = []
        self.__voxel_density = []
        self.__voxel_crir = []
        self.__voxel_fuv = []
        self.__voxel_taufuv = []

        self.__voxel_velocities = []
        self.__voxel_dispersion = []

        self.__x = []
        self.__y = []
        self.__z = []

        # constants.history = 'r{}_cm{}_d{}_uv{}/'.format(int(constants.voxel_size),
        #                                                 '-'.join(str(f) for f in constants.clumpMassFactor),
        #                                                 constants.densityFactor,
        #                                                 constants.globalUV)

        self.suggested_calc = suggested_calc

        self.__verbose = verbose
        self.__debug = debug

        self.__logger = getLogger(__name__)
        if self.__debug:
            self.__logger.setLevel("DEBUG")
        elif self.__verbose:
            self.__logger.setLevel("INFO")
        else:
            self.__logger.setLevel("WARNING")

        return

    def __initialiseGrid(self):
        """
        Initialise all `Voxel()` instances used in model.
        """
        self.__voxels = []
        for i in range(self.__voxel_number):
            self.__voxels.append(Voxel(i))
        return

    def __str__(self):
        """
        Print `VoxelGrid()` instance information.
        """
        return "VoxelGrid\n  ->{} voxels\n".format(self.__voxel_number)

    def __calculateProperties(self, X, Y, Z, average=False):
        """
        This is a method to calculate the dict to unpack into the argument for `Voxel.setProperties()`.
        """

        if average:
            x, y = np.meshgrid(
                np.linspace(
                    X - 0.5 * constants.voxel_size,
                    X + 0.5 * constants.voxel_size,
                    average,
                ),
                np.linspace(
                    Y - 0.5 * constants.voxel_size,
                    Y + 0.5 * constants.voxel_size,
                    average,
                ),
            )
        else:
            x = np.array([X])
            y = np.array([Y])
        rPol = np.array([x.flatten(), y.flatten()]).T
        rPol = np.linalg.norm(rPol, axis=1)
        # print(f'radius{rPol}')

        # Velocity
        velocity = interpolations.interpolate_galaxy_rotation(rPol)

        if constants.from_earth:
            # Calculate the correction to the voxel velocity vectors
            relativeRpol = np.sqrt(
                (x.flatten() - constants.rgal_earth) ** 2 + y.flatten() ** 2
            )
            relativePhi = np.arctan2(y.flatten(), constants.rgal_earth - x.flatten())
            relativeSigma = np.arccos(
                (rPol**2 + relativeRpol**2 - constants.rgal_earth**2)
                / (2 * rPol * relativeRpol)
            )
            relativeTheta = np.arctan(Z / relativeRpol)

            # Correct the relative velocity of the voxel
            velocityEarth = interpolations.interpolate_galaxy_rotation(
                constants.rgal_earth
            )
            velocityCirc = velocity - velocityEarth * rPol / constants.rgal_earth

            velocity = (
                np.sign(relativePhi)
                * velocityCirc
                * np.sin(relativeSigma)
                * np.cos(relativeTheta)
            )
            # velocity = (np.sign(relativePhi) * velocityCirc * np.sin(relativeSigma))
            # velocity = np.sign(np.arctan2(Y,X))*velocity*np.sin(relativeSigma) - velocityEarth*np.sin(relativePhi)

            if (rPol == 0).any():
                velocity[rPol == 0] = 0
            # self.__velocity = (velocity.mean()) * np.sin(self.__phi)

            # this is leftover from Silke's version
            if constants.disp_gc and (
                np.linalg.norm([x.mean(), y.mean()]) <= constants.disp_r_gc
            ):
                ensembleDispersion = np.float64(constants.disp_gc)
            # elif constants.disp_gmc>0.1:
            #    ensembleDispersion = np.linalg.norm((velocityCirc.std(), constants.clump_dispersion, constants.disp_gmc))
            else:
                ensembleDispersion = np.linalg.norm(
                    (velocity.std(), constants.clump_dispersion, constants.disp_gmc)
                )
            velocity = velocity.mean()

        else:
            velocity = velocity.mean()
            ensembleDispersion = constants.ensemble_dispersion

        # Use this to check the evaluation of the velocity field. It is still not working correctly...
        # print(self.__velocity)

        # ensembleDispersion = interpolations.interpolateVelocityDispersion(rPol)
        #
        # ensembleDispersion = ensembleDispersion.mean()

        # Voxel-filling factor
        # These lines were hard-coded
        # f_vox = [interpolations.interpolate_h2_voxel_filling_factor(rPol, Z),
        #          interpolations.interpolate_hi_voxel_filling_factor(rPol, Z)]
        f_vox = []
        for ens in range(constants.ensembles):
            if constants.interclump_idx[ens]:
                if constants.interclump_wnm_idx[ens]:
                    f_vox.append(
                        constants.interclump_wnm_ratio
                        * interpolations.interpolate_hi_voxel_filling_factor(rPol, Z)
                    )
                else:
                    f_vox.append(
                        (1 - constants.interclump_wnm_ratio)
                        * interpolations.interpolate_hi_voxel_filling_factor(rPol, Z)
                    )
            else:
                f_vox.append(
                    interpolations.interpolate_h2_voxel_filling_factor(rPol, Z)
                )

        # Mass (clump mass MH2+(1-f)*MHI, interclump mass f*MHI)
        # These lines were hard-coded
        # ensembleMass = [(constants.h2_mass_factor
        #                  * interpolations.interpolate_h2_mass(rPol, Z)
        #                  + (1-constants.interclump_hi_ratio)
        #                  * constants.hi_mass_factor
        #                  * interpolations.interpolate_hi_mass(rPol, Z)),
        #                 (constants.interclump_hi_ratio
        #                  * constants.hi_mass_factor
        #                  * interpolations.interpolate_hi_mass(rPol, Z))]
        # ensembleMass = [constants.ensemble_mass_factor[ens]*np.asarray(ensembleMass).mean(1)[ens]
        #                 for ens in range(constants.ensembles)]
        ensemble_mass = []
        for ens in range(constants.ensembles):
            if constants.interclump_idx[ens]:
                if constants.interclump_wnm_idx[ens]:
                    ensemble_mass.append(
                        constants.ensemble_mass_factor[ens]
                        * (
                            constants.interclump_wnm_ratio
                            * constants.interclump_hi_ratio
                            * constants.hi_mass_factor
                            * interpolations.interpolate_hi_mass(rPol, Z)
                        ).mean()
                    )
                else:
                    ensemble_mass.append(
                        constants.ensemble_mass_factor[ens]
                        * (
                            (1 - constants.interclump_wnm_ratio)
                            * constants.interclump_hi_ratio
                            * constants.hi_mass_factor
                            * interpolations.interpolate_hi_mass(rPol, Z)
                        ).mean()
                    )
            else:
                ensemble_mass.append(
                    constants.ensemble_mass_factor[ens]
                    * (
                        constants.h2_mass_factor
                        * interpolations.interpolate_h2_mass(rPol, Z)
                        + constants.ensemble_mass_factor[ens]
                        * (1 - constants.interclump_hi_ratio)
                        * constants.hi_mass_factor
                        * interpolations.interpolate_hi_mass(rPol, Z)
                    ).mean()
                )

        # Ensemble density
        # These lines were hard-coded
        # ensembleDensity = interpolations.interpolate_number_density(rPol)
        # ensembleDensity = constants.density_factor*ensembleDensity.mean()
        # ensembleDensity = [ensembleDensity, np.max([0.01*ensembleDensity, constants.interclump_density])]
        # ensembleDensity = [ensembleDensity, np.max([10, constants.interclump_density])]
        ensemble_density = []
        for ens in range(constants.ensembles):
            if constants.interclump_idx[ens]:
                if constants.interclump_wnm_idx[ens]:
                    # Change this if you want a different WNM density (such as 0.1 cm^-3)
                    ensemble_density.append(np.max([10, constants.interclump_density]))
                else:
                    ensemble_density.append(np.max([10, constants.interclump_density]))
            else:
                ensemble_density.append(
                    constants.density_factor
                    * interpolations.interpolate_number_density(rPol).mean()
                )

        # FUV
        fuv = constants.fuv_factor * interpolations.interpolate_fuv(rPol, Z)
        if not constants.clump_log_fuv is None:
            clump_fuv = 10**constants.clump_log_fuv
        else:
            # clump_fuv = np.clip(fuv, np.max((1, 10**constants.interclumpLogFUV)), None).mean()
            clump_fuv = np.clip(fuv, constants.fuv_ism, None).mean()
        if not constants.interclump_log_fuv is None:
            interclump_fuv = 10**constants.interclump_log_fuv
        else:
            # interclump_fuv = np.clip(fuv, np.max((1, 10**constants.clumpLogFUV)), None).mean()
            interclump_fuv = np.clip(fuv, constants.fuv_ism, None).mean()
        if not constants.interclump_wnm_log_fuv is None:
            interclump_wnm_fuv = 10**constants.interclump_wnm_log_fuv
        else:
            # interclump_fuv = np.clip(fuv, np.max((1, 10**constants.clumpLogFUV)), None).mean()
            interclump_wnm_fuv = np.clip(
                constants.interclump_f_fuv_wnm * fuv, constants.fuv_ism, None
            ).mean()
        # This line was hard-coded
        # FUV = [copy(clump_fuv), copy(interclump_fuv)]
        FUV = []
        for ens in range(constants.ensembles):
            if constants.interclump_idx[ens]:
                if constants.interclump_wnm_idx[ens]:
                    FUV.append(copy(interclump_wnm_fuv))
                else:
                    FUV.append(copy(interclump_fuv))
            else:
                FUV.append(copy(clump_fuv))

        # CRIR
        if rPol.mean() >= constants.r_cmz:
            crir = constants.zeta_sol
        else:
            crir = constants.zeta_cmz

        # Save the properties in private lists
        self.__voxel_velocities.append(velocity)
        self.__voxel_dispersion.append(ensembleDispersion)
        self.__voxel_mass.append(ensemble_mass)
        self.__voxel_density.append(ensemble_density)
        self.__voxel_fuv.append(FUV)
        self.__voxel_crir.append(crir)

        self.__properties = {
            # Model parameters
            "from_grid": True,
            "voxel_size": constants.voxel_size,
            "clump_mass_range": constants.clump_log_mass_range,
            "clump_mass_number": constants.clump_mass_number,
            "clump_n_max": constants.clump_n_max,
            "velocity_range": constants.velocity_bin,
            "velocity_number": constants.velocity_number,
            # Voxel properties
            "velocity": velocity,
            "ensemble_dispersion": ensembleDispersion,
            "voxel_factor": f_vox,
            "ensemble_mass": ensemble_mass,
            "ensemble_density": ensemble_density,
            "fuv": FUV,
            "crir": crir,
            # Calculation
            "suggested_calc": self.suggested_calc,
            "velocity_resolution": 1,
        }

        return

    # PUBLIC

    def getDimensions(self):
        """
        Return dimensions of `Shape` object.
        """
        return self.__shape.getDimensions()

    def calculateEmission(
        self,
        index=0,
        dilled=False,
        kind="linear",
        timed=False,
        verbose=False,
        debug=False,
        multiprocessing=0,
    ):
        """
        This will initialise the grid of voxels and calculate their emission. This has to be
        done in succession for each voxel since the calculations are modular (the temporary
        voxel terms are changed when calculating different voxels). This can be rerun and it
        will reinitialise the grid.
        """

        if timed:
            t0 = time()

        # print(observations.tauCenterline)
        # interpolations.initialise_grid(dilled=dilled)
        # interpolations.initialise_model()
        self.__initialiseGrid()

        if timed:
            t1 = time()
            self.__logger.info("Grid initialised: {:.4f} s".format(t1 - t0))

        print("\nCalculating Grid Emission...")

        x, y, z = self.__shape.voxelCartesianPositions()
        r, phi = self.__shape.voxelPolarPositions()
        print()
        # self.__unusedVoxels = []

        # Setup fits files to stream the voxel emissivity and absorption.
        # number
        dim = [len(constants.abundances), constants.ensembles, self.__voxel_number]
        if verbose:
            print("species number dimension:", dim)
        shdu_ensemble_species_number = self.shdu_header(
            name="Clump species number",
            units="#",
            abundance=True,
            filename="species_number",
            velocity=True,
            dim=dim,
        )
        # species
        dim = [
            len(species.clump_transition_wavelengths),
            constants.velocity_range.size,
            self.__voxel_number,
        ]
        if verbose:
            print("species emission dimension:", dim)
        shdu_voxel_emissivity_species = self.shdu_header(
            name="Clump species emissivity",
            units="K/pc",
            molecules=True,
            filename="species_emissivity",
            velocity=True,
            dim=dim,
        )
        shdu_voxel_absorption_species = self.shdu_header(
            name="Clump species absorption",
            units="1/pc",
            molecules=True,
            filename="species_absorption",
            velocity=True,
            dim=dim,
        )
        # HI
        dim = [1, constants.velocity_range.size, self.__voxel_number]
        if verbose:
            print("species emission dimension:", dim)
        shdu_voxel_emissivity_hi = self.shdu_header(
            name="Clump HI emissivity",
            units="K/pc",
            hi=True,
            filename="hi_emissivity",
            velocity=True,
            dim=dim,
        )
        shdu_voxel_absorption_hi = self.shdu_header(
            name="Clump HI absorption",
            units="1/pc",
            hi=True,
            filename="hi_absorption",
            velocity=True,
            dim=dim,
        )
        # dust
        nDust = constants.wavelengths[constants.n_dust].size
        dim = [nDust, self.__voxel_number]
        if verbose:
            print("dust emission dimension:", dim)
        shdu_voxel_emissivity_dust = self.shdu_header(
            name="Clump dust emissivity",
            units="K/pc",
            dust=True,
            filename="dust_emissivity",
            dim=dim,
        )
        shdu_voxel_absorption_dust = self.shdu_header(
            name="Clump dust absorption",
            units="1/pc",
            dust=True,
            filename="dust_absorption",
            dim=dim,
        )

        dim = [np.sum(constants.clump_mass_number), self.__voxel_number]
        shdu_clumpNj = self.shdu_header(
            name="N_j", units="N/A", filename="clump_number", dim=dim
        )
        shdu_clumpRj = self.shdu_header(
            name="Radius", units="pc", filename="clump_radius", dim=dim
        )

        dim = [len(constants.clump_mass_number), self.__voxel_number]
        shdu_f_vox = self.shdu_header(
            name="f_vox", units="N/A", filename="voxel-filling_factor", dim=dim
        )
        shdu_ensemble_mass = self.shdu_header(
            name="Ensemble mass", units="Msol", filename="voxel_ensemble_mass", dim=dim
        )
        shdu_hi_mass = self.shdu_header(
            name="HI mass", units="Msol", filename="voxel_hi_mass", dim=dim
        )
        shdu_h2_mass = self.shdu_header(
            name="H2 mass", units="Msol", filename="voxel_h2_mass", dim=dim
        )
        shdu_t_gas = self.shdu_header(
            name="T_gas", units="K", filename="voxel_t_gas", dim=dim
        )
        shdu_t_dust = self.shdu_header(
            name="T_dust", units="K", filename="voxel_t_dust", dim=dim
        )
        shdu_ensemble_density = self.shdu_header(
            name="Density", units="cm^-3", filename="voxel_density", dim=dim
        )
        dim = [3, self.__voxel_number]
        shdu_voxel_position = self.shdu_header(
            name="Position", units="pc", filename="voxel_position", dim=dim
        )
        dim = [1, self.__voxel_number]
        shdu_voxel_velocity = self.shdu_header(
            name="Velocity", units="km/s", filename="voxel_velocity", dim=dim
        )
        shdu_voxel_dispersion = self.shdu_header(
            name="Velocity dispersion",
            units="km/s",
            filename="voxel_ensemble_dispersion",
            dim=dim,
        )
        dim = [len(constants.clump_mass_number), self.__voxel_number]
        shdu_fuv = self.shdu_header(
            name="FUV", units="Draine", filename="voxel_fuv", dim=dim
        )
        shdu_taufuv = self.shdu_header(
            name="FUV optical depth",
            units="mag",
            filename="voxel_FUVabsorption",
            dim=dim,
        )

        # test of multiprocessing
        def getProperties(ivoxel):
            i, voxel = ivoxel

            self.__x.append(x[i])
            self.__y.append(y[i])
            self.__z.append(z[i])

            voxel.set_index(i)  # -len(self.__unusedVoxels))
            voxel.set_position(x[i], y[i], z[i], r[i], phi[i])
            self.__calculateProperties(x[i], y[i], z[i], average=3)
            # t0 = time()
            # pprint(self.__properties)
            voxel.set_properties(**self.__properties)
            # print('Initialised:', (time()-t0)/60)
            voxel.calculate_emission()
            # print('Calculated:', (time()-t0)/60)
            # print()

            return voxel

        if multiprocessing:
            pool = Pool(multiprocessing)
            chunksize = int(len(self.__voxels) / 1000 / multiprocessing)
            voxels = pool.imap(getProperties, list(enumerate(self.__voxels)), chunksize)
        else:
            voxels = self.__voxels

        with tqdm(
            total=len(self.__voxels),
            desc="Voxels initialised",
            miniters=1,
            dynamic_ncols=True,
        ) as progress:
            # for i,voxel in enumerate(self.__voxels):
            for i, voxel in enumerate(voxels):
                if timed:
                    t2 = time()

                self.__logger.info("Max X, Radius: {} {}".format(max(x), r[i]))

                if not multiprocessing:
                    voxel = getProperties((i, copy(voxel)))

                X, Y, Z = voxel.get_position()
                self.__x.append(X)
                self.__y.append(Y)
                self.__z.append(Z)

                shdu_clumpNj.write(np.hstack(ensemble.clumpNj).flatten())
                shdu_clumpRj.write(np.hstack(masspoints.clump_radius).flatten())

                if timed:
                    self.__logger.info(
                        "\nVoxel initialised: {:.4f} s".format(time() - t2)
                    )

                # this is likely unneeded now that I am directly writing to a fits file
                self.__voxel_taufuv.append(voxel.get_taufuv())

                shdu_voxel_position.write(voxel.get_position())
                velocity, dispersion = voxel.get_velocity()
                if isinstance(velocity, float):
                    shdu_voxel_velocity.write(np.array([velocity]))
                else:
                    shdu_voxel_velocity.write(np.linalg.norm(velocity))
                shdu_voxel_dispersion.write(np.array([dispersion[0]]))
                shdu_f_vox.write(np.asarray(voxel.get_voxel_filling_factor()))
                shdu_ensemble_mass.write(np.asarray(voxel.get_ensemble_mass()))
                shdu_hi_mass.write(np.asarray(voxel.get_hi_mass()))
                shdu_h2_mass.write(np.asarray(voxel.get_h2_mass()))
                shdu_ensemble_species_number.write(
                    np.asarray(voxel.get_species_number(total=False))
                )
                shdu_t_gas.write(np.asarray(voxel.get_gas_temperature(total=False)))
                shdu_t_dust.write(np.asarray(voxel.get_dust_temperature(total=False)))
                shdu_ensemble_density.write(np.asarray(voxel.get_density()))
                shdu_fuv.write(np.array([voxel.get_fuv()]))
                shdu_taufuv.write(np.asarray([voxel.get_taufuv()]))

                # Save emissivity and absorption
                shdu_voxel_emissivity_species.write(
                    voxel.get_species_emissivity(kind=kind, include_dust=True)
                )
                shdu_voxel_absorption_species.write(
                    voxel.get_species_absorption(kind=kind, include_dust=True)
                )
                shdu_voxel_emissivity_hi.write(
                    voxel.get_species_emissivity(hi=True, kind=kind, include_dust=True)
                )
                shdu_voxel_absorption_hi.write(
                    voxel.get_species_absorption(hi=True, kind=kind, include_dust=True)
                )
                shdu_voxel_emissivity_dust.write(
                    voxel.get_dust_emissivity(minimal=True)
                )
                shdu_voxel_absorption_dust.write(
                    voxel.get_dust_absorption(minimal=True)
                )

                voxels[i] = None
                # self.__voxels[i] = None
                # print(len(self.__voxels))

                if timed:
                    self.__logger.info(
                        "Voxel calculated: {:.4f} s / {:.4f} s".format(
                            time() - t2, time() - t1
                        )
                    )

                progress.update()

            shdu_clumpNj.close()
            shdu_clumpRj.close()
            shdu_voxel_position.close()
            shdu_voxel_velocity.close()
            shdu_voxel_dispersion.close()
            shdu_f_vox.close()
            shdu_ensemble_mass.close()
            shdu_ensemble_species_number.close()
            shdu_t_gas.close()
            shdu_t_dust.close()
            shdu_ensemble_density.close()
            shdu_fuv.close()
            shdu_taufuv.close()
            shdu_voxel_emissivity_species.close()
            shdu_voxel_absorption_species.close()
            shdu_voxel_emissivity_hi.close()
            shdu_voxel_absorption_hi.close()
            shdu_voxel_emissivity_dust.close()
            shdu_voxel_absorption_dust.close()

            progress.close()

            print("Emission calculated successfully.")

        return

    #     def writeDataAlt(self):
    #         # Not yet working.
    #
    #         voxel_species_emissivity = fits.ImageHDU(self.__voxelSpeciesEmissivity, header_species)
    #         voxel_species_absorption = fits.ImageHDU(self.__voxelSpeciesAbsorption, header_species)
    #         voxel_dust_emissivity = fits.ImageHDU(self.__voxelDustEmissivity, header_dust)
    #         voxel_dust_absorption = fits.ImageHDU(self.__voxelDustAbsorption, header_dust)
    #
    #         return

    def writeEmission(self, verbose=False, debug=False):
        """
        NO LONGER USED: writing of data is coupled to calculating the emission!
        This will stream the model parameters and emission to FITS files in the corresponding
        directory for the model.
        """

        print("\nStreaming to fits files...\n")

        if debug:
            # dim = [1, self.__voxelNumber]
            # shdu_mass = self.shdu_header(name='Mass', units='Msol', filename='voxel_mass', dim=dim)

            dim = [len(constants.clump_mass_number), self.__voxelNumber]
            shdu_clump_mass = self.shdu_header(
                name="Ensemble Mass",
                units="Msol",
                filename="voxel_ensemble_mass",
                dim=dim,
            )

            # dim = [1, self.__voxelNumber]
            # shdu_interclump_mass = self.shdu_header(name='Interclump Mass', units='Msol',
            #                                         filename='voxel_interclump_mass', dim=dim)

            dim = [len(constants.clump_mass_number), self.__voxelNumber]
            shdu_density = self.shdu_header(
                name="Density", units="cm^-3", filename="voxel_density", dim=dim
            )

        dim = [3, self.__voxelNumber]
        shdu_position = self.shdu_header(
            name="Position", units="pc", filename="voxel_position", dim=dim
        )

        dim = [1, self.__voxelNumber]
        shdu_velocity = self.shdu_header(
            name="Velocity", units="km/s", filename="voxel_velocity", dim=dim
        )

        # dim = [constants.clumpMaxIndeces[0], self.__voxelNumber]
        # shdu_clump_velocity = self.shdu_header(name='Clump velocity', units='km/s',
        #                                        filename='voxel_clump_velocity', dim=dim)

        # dim = [constants.clumpMaxIndeces[1], self.__voxelNumber]
        # shdu_interclump_velocity = self.shdu_header(name='Interclump Velocity', units='km/s',
        #                                             filename='voxel_interclump_velocity', dim=dim)

        dim = [np.hstack(ensemble.clumpNj).size, self.__voxelNumber]
        shdu_clumpNj = self.shdu_header(
            name="N_j", units="N/A", filename="clump_number", dim=dim
        )
        shdu_clumpRj = self.shdu_header(
            name="Radius", units="pc", filename="clump_radius", dim=dim
        )

        dim = [len(constants.clump_mass_number), self.__voxelNumber]
        shdu_FUV = self.shdu_header(
            name="FUV", units="Draine", filename="voxel_fuv", dim=dim
        )
        shdu_FUVabsorption = self.shdu_header(
            name="tau_FUV", units="mag", filename="voxel_FUVabsorption", dim=dim
        )

        dim = [
            len(constants.clump_mass_number),
            constants.abundances,
            self.__voxelNumber,
        ]
        shdu_species_number = self.shdu_header(
            name="N_species",
            units="#",
            filename="voxel_N_species",
            abundance=True,
            dim=dim,
        )

        # This is for a test of the voxel Emissions before streaming
        wav = np.append(
            constants.wavelengths[constants.n_dust], species.molecule_wavelengths
        )
        nDust = constants.wavelengths[constants.n_dust].size

        # The dimensions of the emissions are the expected velocity range (5 for clumps, 9 for interclumps),
        # the number of wavelengths at which the emission is calculated (#molecules + 333 dust wavelengths),
        # and finally the number of voxels.
        # dim = [len(species.moleculeWavelengths)+nDust, constants.clumpMaxIndeces[0], self.__voxelNumber]
        dim = [
            len(species.molecule_wavelengths),
            constants.velocity_range.size,
            self.__voxelNumber,
        ]
        shdu_clump_emissivity_species = self.shdu_header(
            name="Clump species emissivity",
            units="K/pc",
            molecules=True,
            filename="species_emissivity_clump",
            dim=dim,
        )
        shdu_clump_absorption_species = self.shdu_header(
            name="Clump species absorption",
            units="1/pc",
            molecules=True,
            filename="species_absorption_clump",
            dim=dim,
        )

        dim = [nDust, self.__voxelNumber]
        shdu_clump_emissivity_dust = self.shdu_header(
            name="Clump dust emissivity",
            units="K/pc",
            dust=True,
            filename="dust_emissivity_clump",
            dim=dim,
        )
        shdu_clump_absorption_dust = self.shdu_header(
            name="Clump dust absorption",
            units="1/pc",
            dust=True,
            filename="dust_absorption_clump",
            dim=dim,
        )

        # dim = [len(species.moleculeWavelengths)+nDust, constants.clumpMaxIndeces[1], self.__voxelNumber]
        # shdu_interclump_intensity = self.shdu_header(name='Clump intensity', units='K',
        #                                              filename='intensity_interclump', dim=dim)
        # shdu_interclump_tau = self.shdu_header(name='Clump optical depth', units='1/cm',
        #                                        filename='opticalDepth_interclump', dim=dim)

        with tqdm(
            total=len(self.__voxels),
            desc="Voxel emissions",
            miniters=1,
            dynamic_ncols=True,
        ) as progress:
            for i, voxel in enumerate(self.__voxels):
                gc.collect()

                epsilon_species = voxel.getSpeciesEmissivity(
                    include_dust=True, total=True
                )
                kappa_species = voxel.getSpeciesAbsorption(
                    include_dust=True, total=True
                )
                epsilon_dust = voxel.getDustEmissivity(minimal=True)
                kappa_dust = voxel.getDustAbsorption(minimal=True)

                # Optain the voxel emission data
                # clumpIntensity = intensity[0]
                # clumpTau = opticalDepth[0]
                # clumpVelocity = voxel.getClumpVelocity()[0]
                # # print(clumpVelocity.size)
                # interclumpIntensity = intensity[1]
                # interclumpTau = opticalDepth[1]
                # interclumpVelocity = voxel.getClumpVelocity()[1]

                # print()
                # print(voxel.getPosition())
                # print()
                # plt.loglog(wav[nDust:], clumpIntensity[0,nDust:], marker='x', ms=2, ls='')
                # plt.loglog(wav[:nDust], clumpIntensity[0,:nDust], marker='', ls='-', lw=1)
                # plt.show()

                # Transform the voxel emission data to the maximum size in the model (to deal with numpy nd arrays)
                # while clumpIntensity[:,0].size<constants.clumpMaxIndeces[0]:
                #   clumpIntensity = np.append(clumpIntensity, np.zeros((1,clumpIntensity[0,:].size)), axis=0)
                #   clumpTau = np.append(clumpTau, np.zeros((1,clumpTau[0,:].size)), axis=0)
                #   clumpVelocity = np.append(clumpVelocity, [np.nan], axis=0)

                # while interclumpIntensity[:,0].size<constants.clumpMaxIndeces[1]:
                #   interclumpIntensity = np.append(interclumpIntensity, np.zeros((1,interclumpIntensity[0,:].size)),
                #                                   axis=0)
                #   interclumpTau = np.append(interclumpTau, np.zeros((1,interclumpTau[0,:].size)), axis=0)
                #   interclumpVelocity = np.append(interclumpVelocity, [np.nan], axis=0)

                shdu_position.write(voxel.getPosition())
                shdu_clumpNj.write(np.hstack(ensemble.clumpNj).flatten())
                shdu_clumpRj.write(np.hstack(masspoints.clump_radius).flatten())
                velocity = voxel.getVelocity()[0]
                if isinstance(velocity, float):
                    shdu_velocity.write(np.array([velocity]))
                else:
                    shdu_velocity.write(np.linalg.norm(velocity))
                # shdu_clump_velocity.write(clumpVelocity)
                # shdu_interclump_velocity.write(interclumpVelocity)
                shdu_FUV.write(np.array([voxel.getFUV()]))
                shdu_FUVabsorption.write(np.asarray([voxel.getFUVabsorption()]))

                if debug:
                    # shdu_mass.write(np.array([voxel.getMass()]))
                    shdu_clump_mass.write(np.array([voxel.getEnsembleMass()]))
                    # shdu_interclump_mass.write(np.array([voxel.getInterclumpMass()]))
                    shdu_density.write(np.array([voxel.getDensity()]))

                shdu_clump_emissivity_species.write(np.array(epsilon_species))
                shdu_clump_absorption_species.write(np.array(kappa_species))
                shdu_clump_emissivity_dust.write(np.array(epsilon_dust))
                shdu_clump_absorption_dust.write(np.array(kappa_dust))

                # shdu_interclump_intensity.write(interclumpIntensity)
                # shdu_interclump_tau.write(interclumpTau)

                progress.update()

            progress.close()

        print("\nData files have been written successfully.\n")

        return

    def allVoxels(self):
        """Just in case all of the `Voxel()` instances need to be retrieved"""
        return self.__voxels

    #     # These methods may be used if and when we decide to work with models as objects
    #     # rather than separating the data into files to post-process.
    #
    #     def getVoxelNumber(self):
    #         return self.__voxelNumber
    #
    #     def getVoxelPositions(self):
    #         return np.array([self.__x, self.__y, self.__z])
    #
    #     def totalEmission(self):
    #         # Return the emission from all of the voxels, separated by observed velocity
    #         return np.array([self.__voxelIntensity.sum(1), self.__voxelOpticalDepth.sum(1)])
    #
    #     def clumpEmission(self):
    #         # Return the emission from all of the voxels, separated by observed velocity
    #         return np.array([self.__voxelIntensity[:, 0], self.__voxelOpticalDepth[:, 0]])
    #
    #     def interclumpEmission(self):
    #         # Return the emission from all of the voxels, separated by observed velocity
    #         return np.array([self.__voxelIntensity[:, 1], self.__voxelOpticalDepth[:, 1]])
    #
    #     def getFUV(self):
    #         return self.__voxelFUV
    #
    #     def getFUVabsorption(self):
    #         return self.__voxelFUVabsorption

    def printVoxels(self):
        """Print each voxel."""
        for voxel in self.__voxels:
            voxel.printVoxel()
        return

    def voxelproperties(self, x, y, z):
        """Print the parameters passed to a particular voxel instance depending on x, y, z."""
        self.__calculateProperties(x, y, z)
        print(self.__properties)
        return

    def shdu_header(
        self,
        name="",
        units="",
        molecules=False,
        hi=False,
        dust=False,
        abundance=False,
        velocity=False,
        dim=None,
        kw=25,
        cw=50,
        filename=None,
    ):
        """
        method to setup the header used for a `StreamingHDU` instance.
        """

        if filename == None or dim == None:
            return

        header = fits.Header()

        header["SIMPLE"] = (True, "conforms to FITS standard")
        header["BITPIX"] = (-64, "element size")
        header["NAXIS"] = (len(dim), "number of axes")
        header["EXTEND"] = True
        if molecules:
            header["SPECIES"] = ", ".join(species.clump_transitions)
        elif hi:
            header["SPECIES"] = "HI"
        elif abundance:
            header["SPECIES"] = ", ".join(constants.abundances)
        elif dust:
            header["DUST"] = ", ".join(constants.dust_names[constants.n_dust])
        else:
            pass

        header["NAME"] = name
        header["UNITS"] = units

        for i in range(len(dim)):
            header["NAXIS{}".format(i + 1)] = dim[i]
            if velocity & i == 0:
                header["CTYPE1"] = "Transition/wavelength"
                header["CUNIT1"] = "m"
                header["CRPIX1"] = "N/A"
                header["CRVAL1"] = "N/A"
                header["CDELT1"] = "N/A"
            elif velocity & i == 1:
                header["CTYPE2"] = "Velocity"
                header["CUNIT2"] = "km/s"
                header["CRPIX2"] = (header["NAXIS2"] - 1) / 2
                header["CRVAL2"] = 0
                header["CDELT2"] = (
                    constants.velocity_range[-1] - constants.velocity_range[0]
                ) / (header["NAXIS2"] - 1)
            elif velocity & i == 2:
                header["CTYPE3"] = "Voxel"
                header["CUNIT1"] = "N/A"
                header["CRPIX1"] = "N/A"
                header["CRVAL1"] = "N/A"
                header["CDELT1"] = "N/A"

        # header['NAXIS1'] = self.__constants.velocityBins.size#emission[0].data[0,0,0,:].size
        # header['NAXIS2'] = len(self.__species[0].getMolecules()) + \
        #                    len(self.__species[1].getDust())#emission[0].data[0,0,:,0].size
        # header['NAXIS3'] = 2#emission[0].data[:,0,0,0].size
        # header['NAXIS4'] = self.__voxelNumber#emission[0].data[0,:,0,0].size

        header["COMMENT"] = "".ljust(cw)
        header["COMMENT"] = "voxel_size: {} pc".format(constants.voxel_size).ljust(cw)
        header["COMMENT"] = "".ljust(cw)
        header["COMMENT"] = "".ljust(cw, "-")
        header["COMMENT"] = "Model Information".center(cw)
        header["COMMENT"] = "".ljust(cw, "-")
        header["COMMENT"] = "It is necessary to have a record of the".ljust(cw)
        header["COMMENT"] = "relevant model information in each file,".ljust(cw)
        header["COMMENT"] = "so below we record the relevant information".ljust(cw)
        header["COMMENT"] = "defining the model parameters.".ljust(cw)
        header["COMMENT"] = "".ljust(cw)
        header["COMMENT"] = constants.history.ljust(cw)
        header["COMMENT"] = "".ljust(cw)
        header["COMMENT"] = "".ljust(cw)
        header["COMMENT"] = "Grid files".ljust(cw)
        header["COMMENT"] = "".ljust(cw)
        header["COMMENT"] = constants.tb_grid_file.ljust(cw)
        header["COMMENT"] = constants.tau_grid_file.ljust(cw)
        header["COMMENT"] = constants.tau_fuv_grid_file.ljust(cw)
        header["COMMENT"] = "".ljust(cw)
        header["COMMENT"] = "Input files".ljust(cw)
        header["COMMENT"] = "".ljust(cw)
        header["COMMENT"] = "galactic center: {} pc".format(constants.r_gc).ljust(cw)
        header["COMMENT"] = "full voxels: {}".format(constants.all_full).ljust(cw)
        header["COMMENT"] = constants.h2_surface_density_file.ljust(cw)
        header["COMMENT"] = "  scale galactic center: {}".format(
            constants.mh2_scale_gc
        ).ljust(cw)
        header["COMMENT"] = constants.hi_surface_density_file.ljust(cw)
        header["COMMENT"] = "  scale galactic center: {}".format(
            constants.mhi_scale_gc
        ).ljust(cw)
        header["COMMENT"] = "  same as h2: {}".format(constants.hi_mass_file).ljust(cw)
        header["COMMENT"] = constants.h2_scale_height_file.ljust(cw)
        header["COMMENT"] = constants.hi_scale_height_file.ljust(cw)
        header["COMMENT"] = constants.h_number_density_file.ljust(cw)
        header["COMMENT"] = constants.fuv_file.ljust(cw)
        header["COMMENT"] = "  scale galactic center: {}".format(
            constants.fuv_scale_gc
        ).ljust(cw)
        header["COMMENT"] = "  average spectrum: {}".format(
            constants.average_fuv
        ).ljust(cw)
        header["COMMENT"] = "  integration range: {}".format(constants.l_range).ljust(
            cw
        )
        header["COMMENT"] = constants.velocity_file.ljust(cw)
        header["COMMENT"] = "".ljust(cw)
        header["COMMENT"] = "improved calculation: {}".format(
            self.suggested_calc
        ).ljust(cw)
        header["COMMENT"] = "".ljust(cw)
        header["COMMENT"] = "".ljust(cw)
        header["COMMENT"] = "Model".center(cw)
        header["COMMENT"] = "".ljust(cw)
        header["COMMENT"] = "Model shape: {}".rjust(kw).format(constants.type).ljust(cw)
        x, y, z = self.__shape.getShape()
        header["COMMENT"] = "voxels: {}".rjust(kw).format(self.__voxel_number).ljust(cw)
        header["COMMENT"] = "x (pc): {}".rjust(kw).format(x).ljust(cw)
        header["COMMENT"] = "y (pc): {}".rjust(kw).format(y).ljust(cw)
        header["COMMENT"] = "z (pc): {}".rjust(kw).format(z).ljust(cw)
        header["COMMENT"] = "radius (pc): {}".rjust(kw).format(constants.rgal).ljust(cw)
        header["COMMENT"] = (
            "R_gal,Sun (pc): {}".rjust(kw).format(constants.rgal_earth).ljust(cw)
        )
        # header['COMMENT'] = '  thickness:'.ljust(50)
        header["COMMENT"] = "".ljust(cw)
        header["COMMENT"] = "".ljust(cw)
        header["COMMENT"] = "Ensemble properties".center(cw)
        header["COMMENT"] = "".ljust(cw)
        header["COMMENT"] = (
            "mass range (M_Sun): {}".rjust(kw)
            .format(constants.clump_log_mass_range)
            .ljust(cw)
        )
        header["COMMENT"] = (
            "# clumps: {}".rjust(kw).format(constants.clump_mass_number).ljust(cw)
        )
        header["COMMENT"] = (
            "log m_cl (M_Sun): {}".rjust(kw).format(constants.clump_log_mass).ljust(cw)
        )
        header["COMMENT"] = (
            "disp_gmc (km/s): {}".rjust(kw).format(constants.disp_gmc).ljust(cw)
        )
        header["COMMENT"] = (
            "N_m_cl: {}".rjust(kw).format(constants.clump_n_max).ljust(cw)
        )
        header["COMMENT"] = (
            "v_obs_range (km/s): {}".rjust(kw).format(constants.velocity_bin).ljust(cw)
        )
        header["COMMENT"] = (
            "N_v_obs: {}".rjust(kw).format(constants.velocity_number).ljust(cw)
        )
        header["COMMENT"] = "".ljust(cw)
        header["COMMENT"] = "".ljust(cw)
        header["COMMENT"] = "Parameters".center(cw)
        header["COMMENT"] = "".ljust(cw)
        header["COMMENT"] = (
            "f_Mcl: {}".rjust(kw).format(constants.h2_mass_factor).ljust(cw)
        )
        header["COMMENT"] = (
            "f_Micl: {}".rjust(kw).format(constants.hi_mass_factor).ljust(cw)
        )
        header["COMMENT"] = (
            "f_Mens: {}".rjust(kw).format(constants.ensemble_mass_factor).ljust(cw)
        )
        header["COMMENT"] = (
            "f_n: {}".rjust(kw).format(constants.density_factor).ljust(cw)
        )
        header["COMMENT"] = "f_FUV: {}".rjust(kw).format(constants.fuv_factor).ljust(cw)
        header["COMMENT"] = (
            "f_FUV_WNM: {}".rjust(kw).format(constants.interclump_f_fuv_wnm).ljust(cw)
        )
        header["COMMENT"] = (
            "f_HI: {}".rjust(kw).format(constants.interclump_hi_ratio).ljust(cw)
        )
        header["COMMENT"] = (
            "f_WNM: {}".rjust(kw).format(constants.interclump_wnm_ratio).ljust(cw)
        )
        header["COMMENT"] = (
            "log chi_cl: {}".rjust(kw).format(constants.clump_log_fuv).ljust(cw)
        )
        header["COMMENT"] = (
            "log chi_icl: {}".rjust(kw).format(constants.interclump_log_fuv).ljust(cw)
        )
        header["COMMENT"] = (
            "log chi_wnm: {}".rjust(kw)
            .format(constants.interclump_wnm_log_fuv)
            .ljust(cw)
        )
        header["COMMENT"] = "R_gal,CMZ: {}".rjust(kw).format(constants.r_cmz).ljust(cw)
        header["COMMENT"] = (
            "zeta_CMZ: {}".rjust(kw).format(constants.zeta_cmz).ljust(cw)
        )
        header["COMMENT"] = (
            "zeta_Sun: {}".rjust(kw).format(constants.zeta_sol).ljust(cw)
        )
        header["COMMENT"] = "".ljust(cw)
        header["COMMENT"] = "".ljust(cw)
        header["COMMENT"] = "".ljust(cw)

        directory = constants.HISTORYPATH + constants.directory + constants.history

        if not os.path.exists(directory):
            os.makedirs(directory)

        if ".fits" not in filename:
            filename = directory + filename + ".fits"
        else:
            filename = directory + filename

        if os.path.exists(filename):
            os.remove(filename)

        shdu = fits.StreamingHDU(filename, header)

        return shdu
