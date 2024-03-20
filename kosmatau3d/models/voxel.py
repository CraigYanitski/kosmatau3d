"""
A module containing the :code:`Voxel` class, which is used to model the
synthetic emission from one or more *ensembles*.
The :code:`Voxel` class wraps together all of the computations of
:ref:`models.masspoints.masspoint`, :ref:`models.combinations.combination`,
and :ref:`models.ensemble.ensemble`.
"""

import importlib as il
from copy import copy, deepcopy
from time import time

import numpy as np
from numba import jit
from scipy.interpolate import interp1d
from scipy.optimize import curve_fit
from logging import getLogger, DEBUG, INFO, WARNING

from kosmatau3d.models import (
    constants,
    interpolations,
    observations,
    species,
    ensemble,
    combinations,
    masspoints,
)


class Voxel(object):
    """
    .. Voxel:

    A class to compute the emission from a voxel in :code:`kosmatau3d`.
    It contains ensembles of spherical PDR simulations (to mimic the fractal
    structure of PDRs), a FUV field element, and element absorption, separated
    for dust and transitions.
    It has been tested to contain two ensembles: one for the dense clumpy
    medium and one for the diffuse interclump medium.
    This is to account for dense clumps embedded in a diffuse environment.
    """

    # PRIVATE

    def __init__(self, index=0, debugging=False):
        """ """

        self.__index = index  # index of voxel in VoxelGrid, sort of like its ID in the overall model

        self.__debugging = debugging

        self.__mass = 0
        self.__ensemble_mass = 0
        self.__velocity = 0  # velocity of mass at voxel point
        self.__ensemble_dispersion = 0  # dispersion of velocity at voxel point

        self.__hi_mass = 0
        self.__h2_mass = 0

        self.__clump_velocity_indeces = []
        self.__model_mass = []
        self.__volume_factor = []

        self.__intensity_species = []  # intensity of species transitions in voxel
        self.__optical_depth_species = (
            []
        )  # optical depth of species transitions in voxel
        self.__intensity_dust = []  # intensity of dust continuum in voxel
        self.__optical_depth_dust = []  # optical depth of dust continuum in voxel

        self.__fuv = 0.0
        self.__taufuv = 0.0

        self.__x = 0
        self.__y = 0
        self.__z = 0

        self.__r = 0
        self.__phi = 0

        # testing flags
        self.test_calc = False
        self.test_opacity = False
        self.test_pexp = False
        self.test_fv = False
        self.suggested_calc = True

        return

    def __set_mass(self):
        """Set the total mass."""
        self.__mass = self.__clump_mass.sum()

    def __set_clump_mass(self, r, z):
        """Set the clump mass."""
        mass = [
            interpolations.interpolate_h2_mass(r, z),
            interpolations.interpolate_hi_mass(r, z),
        ]
        self.__clump_mass = constants.clump_mass_factor * np.asarray(mass).mean(1)
        return

    def __setVelocity(self, r):
        """Set the average velocity."""
        velocity = interpolations.interpolate_galaxy_rotation(r)

        if constants.from_earth:
            # Calculate the correction to the voxel velocity vectors
            relativeRpol = np.sqrt(
                (self.__x - constants.r_gal_earth) ** 2 + self.__y**2
            )
            relativePhi = np.arctan2(self.__y, self.__x - constants.r_gal_earth)
            relativeSigma = np.arccos(
                (self.__r**2 + relativeRpol**2 - constants.r_gal_earth**2)
                / (2 * self.__r * relativeRpol)
            )
            sigma = np.arctan2(self.__z, abs(self.__x - constants.r_gal_earth))

            # Correct the relative velocity of the voxel
            velocityEarth = interpolations.interpolate_galaxy_rotation(
                constants.r_gal_Earth
            )
            velocityCirc = (
                velocity.mean() - velocityEarth * self.__r / constants.r_gal_earth
            )

            self.__velocity = (
                np.sign(relativePhi)
                * velocityCirc
                * np.sin(relativeSigma)
                * np.cos(sigma)
            )

            if self.__r == 0:
                self.__velocity = 0
            # self.__velocity = (velocity.mean()) * np.sin(self.__phi)

        else:
            self.__velocity = np.array(velocity)

        # Use this to check the evaluation of the velocity field. It is still
        # not working correctly...
        # print(self.__velocity)

        ensembleDispersion = interpolations.interpolate_velocity_dispersion(r)

        self.__ensemble_dispersion = ensembleDispersion.mean()

        return

    def __set_ensemble_density(self, r):
        """Set the ensemble-averaged density."""
        density = interpolations.interpolate_number_density(r)
        self.__ensemble_density = constants.density_factor * density.mean()
        return

    def __set_taufuv(self):
        """Set the voxel-averaged FUV optical depth."""
        self.__taufuv = interpolations.interpolate_taufuv(
            self.__ensembleDensity, self.__clumpMass + self.__interclumpMass
        )
        return

    def __set_fuv(self, r, z):
        """
        This is in units of the Draine field.
        """
        fuv = interpolations.interpolate_fuv(r, z) / constants.u_draine0
        self.__fuv = np.clip(fuv, 1, None)
        # self.__fuv = FUVfield(fuv)
        return

    def __str__(self):
        """ """
        return "Voxel {}".format(self.__index)

    # PUBLIC

    # def reloadModules(self):
    #   il.reload(Ensemble)
    #   il.reload(FUVfield)
    #   self.__clump.reloadModules()
    #   self.__interclump.reloadModules()
    #   return

    def set_index(self, index):
        """Change the index."""
        self.__index = index
        return

    def get_index(self):
        """Return the voxel index."""
        return self.__index

    def set_position(self, x, y, z, r, phi):
        """Change the voxel location."""
        self.__x = x
        self.__y = y
        self.__z = z
        self.__r = r
        self.__phi = phi
        return

    def get_fuv(self):
        """Return the voxel-averaged far-UV density."""
        return self.__fuv

    def set_properties(
        self,
        voxel_size=1,
        transitions="all",
        abundances=["ELECTR", "H", "H2", "H+"],
        dust="PAH",
        alpha=1.84,
        gamma=2.31,
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
        clump_mass_number=(3, 1),
        clump_mass_range=([0, 2], [-3]),
        clump_n_max=(1, 100),
        interclump_idx=[False, True],
        interclump_wnm_idx=[False, False],
        velocity_range=(-10, 10),
        velocity_number=201,
        velocity=0.0,
        velocity_resolution=10,
        ensemble_mass=(100, 40),
        ensemble_dispersion=2,
        ensemble_density=(1e5, 1911),
        volume_factor=None,
        voxel_factor=1,
        fuv=200,
        crir=2e-16,
        suggested_calc=True,
        from_grid=False,
        new_grid=False,
        change_interpolation=False,
        dilled=False,
        timed=False,
        verbose=False,
        debug=False,
    ):
        """
        This method initialises all of the voxel parameters and, if not already
        evaluated, also initialises the KOSMA-:math:`\\tau` grid and prepares
        interpolation routines.

        There are quite a few kwargs that can be used to initialise a
        :class:`Voxel` instance:

        :param voxel_size: The voxel size in parsecs.
            The default is :code:`1` pc.
        :param transitions: The transitions included in the model.
            This is a list of strings, where each string has the species name
            (as in the grid) followed by the transition number (which is taken
            from KOSMA-:math:`\\tau`) or :code:`"all"` for all transitions
            in the grid.
            It is set to :code:`"all"` by default.
        :param abundances: The species included in the model.
            Independent of the species specified in :code:`transitions`,
            it sets the species for which the fractional abundances are
            calculated for each ensemble.
            The default value is :code:`["ELECTR", "H", "H2", "H+"]`.
        :param dust: This is a string to limit how much of the dust continuum
            is included.
            The available continuum ranges from 1 nm to 3.1 mm.
            This can be limited to the range
            of molecular lines by setting :code:`dust="molecular"` (22 dust
            wavelengths), or to the range of PAH dust features by setting
            :code:`dust="PAH"` (201 dust wavelengths). All dust wavelengths
            (333) will be included in the calculation if :code:`dust=""`.
        :param alpha: The power law distribution required for initial mass
            function :math:`\\mathrm{d}N(M) / \\mathrm{d}M = M^-\\alpha`.
            The default is :code:`1.84` as determined in Heithausen et al.
            (1998).
        :param gamma: The power law mass-size relation satisfying
            :math:`M = C R^\\gamma`.
            The default is :code:`2.31` as determined in Heithausen et al.
            (1998).
        :param velocity_range: The range (list) of the observed radial
            velocities.
            By default it is :code:`[-10, 10]`.
        :param velocity_number: The number of observed velocities in the
            specified range (including endpoints).
            The default is :code:`51`.
        :param clump_mass_number: The number of clump masses for each set.
            This should be a list of integers with the same length as
            clump_mass_range.
            The default is :code:`[3, 1]`.
        :param clump_mass_range: The range (list) of clump masses for each set.
            This is a list of ranges with the same length as clump_mass_number.
            The ranges can have length 1, in which case it only evaluates that
            clump mass.
            The masses are in units of :math:`\mathrm{dex}(M_\odot)`.
            By default it is :code:`[[0, 2], [-2]]`.
        :param clump_density: This changes the default clump density in the
            constants module.
            This overrides the voxel density kwarg, and is useful for creating
            interclumps medium that is not scaled according to the clumps.
            The default is :code:`[None, 1911]`, so the first clump set takes
            the voxel density.
            This is in units of :math:`\mathrm{cm}^{-3}`.
        :param clump_fuv: This changes the default clump FUV in the constants
            module.
            This overrides the voxel FUV kwarg.
            It is a list of length of the number of clump sets.
            Use a value of :code:`None` to use the voxel FUV field.
            The default is :code:`[None, 1]`, in units of the Draine field.
        :param clump_n_max: The maximum number of the largest clump in the
            ensemble.
            This helps to limit the calculation of a large combination of
            clumps by rescaling the voxel (thus eliminating any change to the
            brightness temperature).
            Larger combinations are used if the rescaled voxel is smaller than
            the largest clump.
            The default is :code:`[1, 100]`.

        :param velocity: The average radial velocity of the observed Voxel
            instance.
            This is used with the model velocities to determine the number of
            clumps seen at a given velocity.
            The default is 0 :math:`\mathrm{km\, s^{-1}}`.
        :param ensemble_mass: This can be either a float or list of floats,
            depending on how many ensembles you have.
            The default is 100 :math:`M_\odot`.
        :param ensemble_dispersion: The velocity dispersion of the ensemble.
            This is used with the model velocities to determine the number
            of clumps at a given velocity.
            The default is 1 :math:`\mathrm{km, s^{-1}}`.
        :param ensemble_density: The observed hydrogen density in the voxel.
            This is different from and overridden by any default density
            specified in the constants module.
            The default is 15000 :math:`\mathrm{cm^{-3}}`.
        :param volume_factor:
        :param voxel_factor:
        :param fuv: The far-UV radiation of the voxel. All clumps are
            isotropically radiated with this radiation.
            It is different from and overridden by any default far-UV radiation
            specified in the constants sub-module.
            The default is 20000 :math:`\\chi_\mathrm{D}`, where
            :math:`\\chi_\mathrm{D}` refers to the radiation in the Draine
            field.
        :param crir: The primary cosmic ray ionisation rate with respect to
            molecular hydrogen (:math:`\\zeta_{H_2}`). By default,
            the local rate is used (2e-16 :math:`s^{-1}`).

        :param suggested_calc: This flag is used to specify the corrected
            version of the calculation should be used. It is True by default.
        :param from_file: This flag can be set to retrieve the voxel properties
            from a file. It is False by default.
        :param verbose: This is mainly used for debugging purposes. It prints
            various statements as the code is evaluated.
        :param debug:
        :param timed:

        """
        # print('Voxel instance initialised')

        self.__logger = getLogger(__name__)
        if debug:
            self.__logger.setLevel(DEBUG)
        elif verbose:
            self.__logger.setLevel(INFO)
        else:
            self.__logger.setLevel(WARNING)

        if timed:
            t0 = time()

        self.suggested_calc = suggested_calc

        x, y = np.meshgrid(
            np.linspace(
                self.__x - 0.5 * constants.voxel_size,
                self.__x + 0.5 * constants.voxel_size,
                3,
            ),
            np.linspace(
                self.__y - 0.5 * constants.voxel_size,
                self.__y + 0.5 * constants.voxel_size,
                3,
            ),
        )
        r = np.array([x.flatten(), y.flatten()]).T
        r = np.linalg.norm(r, axis=1)

        if not from_grid:
            if constants.voxel_size != voxel_size:
                constants.voxel_size = voxel_size

            if constants.alpha != alpha or constants.gamma != gamma:
                constants.change_mass_function_parameters(alpha=alpha, gamma=gamma)

            if (
                constants.velocity_bin != velocity_range
                or constants.velocity_number != velocity_number
            ):
                constants.change_velocity_range(velocity_range)
                constants.change_velocity_number(velocity_number)

            if constants.velocity_resolution != velocity_resolution:
                constants.velocity_resolution = velocity_resolution

            if (
                constants.clump_log_mass_range != clump_mass_range
                or constants.clump_mass_number != clump_mass_number
                or constants.clump_n_max != clump_n_max
            ):
                constants.add_clumps(
                    mass_range=clump_mass_range,
                    num=clump_mass_number,
                    n_max=clump_n_max,
                    reset=True,
                )

            if constants.interclump_idx != interclump_idx:
                constants.set_interclump_ensemble(interclump_idx)
                change_interpolation = True

            if constants.interclump_wnm_idx != interclump_idx:
                constants.set_interclump_wnm(interclump_wnm_idx)

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

            if (
                new_grid
                or change_interpolation
                or not interpolations.initialised
                or not observations.grid_initialised
                or dust != constants.dust
            ):
                constants.change_dust_wavelengths(dust)
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
                species.add_transitions(transitions)
                species.add_transitions(transitions, interclump=True)
                interpolations.initialise_grid(dilled=dilled)

        # print('interpolation initialised:', (time()-t0)/60)

        if timed:
            t1 = time()
            self.__logger.info("Model setup: {}".format(t1 - t0))

        masspoints.reinitialise()
        combinations.reinitialise()
        ensemble.reinitialise()

        self.__velocity = velocity

        if isinstance(ensemble_mass, (list, tuple, np.ndarray)):
            self.__ensemble_mass = ensemble_mass
        else:
            self.__ensemble_mass = [ensemble_mass] * len(constants.clump_mass_number)

        if isinstance(ensemble_dispersion, (list, tuple, np.ndarray)):
            self.__ensemble_dispersion = ensemble_dispersion
        else:
            self.__ensemble_dispersion = [ensemble_dispersion] * len(
                constants.clump_mass_number
            )

        if isinstance(voxel_factor, (int, float)):
            self.__voxel_filling_factor = [voxel_factor] * len(
                constants.clump_mass_number
            )
        else:
            self.__voxel_filling_factor = voxel_factor

        if isinstance(volume_factor, (float, int)):
            volume_factor = [volume_factor] * len(constants.clump_mass_number)
        if volume_factor:
            ensemble_density = [
                ensemble_mass[ens]
                * constants.mass_solar
                / constants.mass_h
                / volume_factor[ens]
                / constants.voxel_size**3
                / constants.pc**3
                / 100**3
                for ens in range(len(constants.clump_mass_number))
            ]
        if isinstance(ensemble_density, (list, tuple, np.ndarray)):
            self.__ensemble_density = ensemble_density
        else:
            self.__ensemble_density = [ensemble_density] * len(
                constants.clump_mass_number
            )

        if isinstance(fuv, (list, tuple, np.ndarray)):
            self.__fuv = fuv
        else:
            self.__fuv = [fuv] * len(clump_mass_number)

        self.__crir = crir

        velocity = self.__velocity

        # This will allow the code to reuse the standard clump density constants for voxel sets (ie. do not alter the model constants)
        # density = copy(constants.clumpDensity)
        # fuv = copy(constants.clumpFUV)
        # for i in range(len(density)):
        # if density[i]=='auto': density[i] = self.__ensembleDensity
        # if fuv[i]==None: fuv[i] = self.__FUV

        # print('properties set:', (time()-t0)/60)

        masspoints.set_masspoint_data(
            density=self.__ensemble_density, fuv=self.__fuv, crir=self.__crir
        )
        # print('  masspoints:', (time()-t0)/60)
        ensemble.initialise(
            ensemble_dispersion=self.__ensemble_dispersion,
            ensemble_mass=self.__ensemble_mass,
            suggested_calc=self.suggested_calc,
        )
        # print('  ensemble:', (time()-t0)/60)
        combinations.initialise(
            clump_combination=[
                ensemble.clumpCombinations[ens][ensemble.clumpLargestIndex[ens]]
                for ens in range(len(constants.clump_mass_number))
            ],
            total_combination=[
                ensemble.CLmaxCombinations[ens][0]
                for ens in range(len(constants.clump_mass_number))
            ],
        )

        # self.__hi_mass =

        # print('probabilities calculated:', (time()-t0)/60)
        # print('ensemble velocities:', ensemble.clumpVelocities)

        if timed:
            t2 = time()
            self.__logger.info("Modules initialised:".format(t2 - t1))

        self.__hi_mass = [0.0] * constants.ensembles
        self.__h2_mass = [0.0] * constants.ensembles

        for ens in range(len(constants.clump_mass_number)):
            if self.suggested_calc:
                if self.__ensemble_dispersion[ens] > constants.clump_dispersion:
                    dispersion = np.sqrt(
                        self.__ensemble_dispersion[ens] ** 2
                        - constants.clump_dispersion**2
                    )
                    self.__model_mass.append(
                        (
                            ensemble.clumpDeltaNji[ens].sum(1)
                            / np.sqrt(2 * np.pi)
                            / dispersion
                            * self.__ensemble_dispersion[ens]
                            / constants.clump_dispersion
                            * 10 ** constants.clump_log_mass[ens]
                        ).sum()
                    )
                else:
                    self.__model_mass.append(
                        (
                            ensemble.clumpDeltaNji[ens].sum(1)
                            / np.sqrt(2 * np.pi)
                            / constants.clump_dispersion
                            * 10 ** constants.clump_log_mass[ens]
                        ).sum()
                    )
            else:
                self.__model_mass.append(
                    (
                        ensemble.clumpDeltaNji[ens].sum(1)
                        * 10 ** constants.clump_log_mass[ens]
                    ).sum()
                )
            self.__volume_factor.append(
                (
                    ensemble.clumpNj[ens]
                    * (masspoints.clump_radius[ens] ** 3 * 4 * np.pi / 3)
                ).sum()
                / constants.voxel_size**3
            )
            if (
                abs(self.__model_mass[ens] - self.__ensemble_mass[ens])
                > 0.1 * self.__ensemble_mass[ens]
                and not suggested_calc
            ):
                self.__logger.error(
                    "ERROR: Voxel {} mass difference for clump set {} greater than 10%".format(
                        self.__index, ens + 1
                    )
                )

        # clumpAfuv,interclumpAfuv = combinations.getAfuv()
        # clumpAfuv *= ensemble.CLmaxProbability
        # interclumpAfuv *= ensemble.ICmaxProbability

        self.__clump_velocities = copy(ensemble.clumpVelocities)
        self.__clump_velocity_indeces = copy(
            ensemble.clumpIndeces
        )  # all list entries should be the same

        if suggested_calc:
            self.__intensity_species = [
                np.zeros(
                    (len(ensemble.clumpIndeces[_]), len(species.clump_transitions)),
                    dtype=np.float64,
                )
                for _ in range(constants.ensembles)
            ]
            self.__optical_depth_species = [
                np.zeros(
                    (len(ensemble.clumpIndeces[_]), len(species.clump_transitions)),
                    dtype=np.float64,
                )
                for _ in range(constants.ensembles)
            ]
            self.__intensity_dust = [
                np.zeros(
                    (
                        len(ensemble.clumpIndeces[_]),
                        constants.wavelengths[constants.n_dust].size,
                    ),
                    dtype=np.float64,
                )
                for _ in range(constants.ensembles)
            ]
            self.__optical_depth_dust = [
                np.zeros(
                    (
                        len(ensemble.clumpIndeces[_]),
                        constants.wavelengths[constants.n_dust].size,
                    ),
                    dtype=np.float64,
                )
                for _ in range(constants.ensembles)
            ]
            self.__intensity_hi = [
                np.zeros((len(ensemble.clumpIndeces[_]), 1), dtype=np.float64)
                for _ in range(constants.ensembles)
            ]
            self.__optical_depth_hi = [
                np.zeros((len(ensemble.clumpIndeces[_]), 1), dtype=np.float64)
                for _ in range(constants.ensembles)
            ]
        else:
            self.__intensity_species = [
                np.zeros(
                    (constants.velocity_range.size, len(species.clump_transitions)),
                    dtype=np.float64,
                )
                for _ in range(constants.ensembles)
            ]
            self.__optical_depth_species = [
                np.zeros(
                    (constants.velocity_range.size, len(species.clump_transitions)),
                    dtype=np.float64,
                )
                for _ in range(constants.ensembles)
            ]
            self.__intensity_dust = [
                np.zeros(
                    (
                        constants.velocity_range.size,
                        constants.wavelengths[constants.n_dust].size,
                    ),
                    dtype=np.float64,
                )
                for _ in range(constants.ensembles)
            ]
            self.__optical_depth_dust = [
                np.zeros(
                    (
                        constants.velocity_range.size,
                        constants.wavelengths[constants.n_dust].size,
                    ),
                    dtype=np.float64,
                )
                for _ in range(constants.ensembles)
            ]

        self.__emissivity_species = deepcopy(self.__intensity_species)
        self.__absorption_species = deepcopy(self.__optical_depth_species)
        self.__emissivity_dust = deepcopy(self.__intensity_dust)
        self.__absorption_dust = deepcopy(self.__optical_depth_dust)
        self.__emissivity_hi = deepcopy(self.__intensity_hi)
        self.__absorption_hi = deepcopy(self.__optical_depth_hi)

        # This gives an error if there are too many clumps in a line-of-sight; tau_FUV is too large for this equation...
        Afuv = combinations.get_fuv_extinction()
        self.__taufuv = [
            -np.log((ensemble.CLmaxProbability[ens].prod(1) * Afuv[ens]).sum())
            for ens in range(len(constants.clump_mass_number))
        ]

        if timed:
            self.__logger.info(
                "setProperties() time of execution: {}".format(time() - t0)
            )

        # print('voxel initialised:', (time()-t0)/60)
        self.__logger.info("voxel initialised")

        return

    def get_position(self):
        """Return voxel position."""
        return np.array([self.__x, self.__y, self.__z])

    def get_density(self):
        """Return the ensemble-averaged density."""
        return self.__ensemble_density

    def get_ensemble_mass(self, total=False):
        """Return the ensemble mass."""
        if total:
            return np.sum(self.__ensemble_mass)
        else:
            return self.__ensemble_mass

    def get_model_mass(self, total=False):
        """Return the mass represented by the model."""
        if total:
            return np.sum(self.__model_mass)
        else:
            return self.__model_mass

    def get_hi_mass(self, total=False):
        """Return the :math:`\mathrm{HI}` mass represented in the model."""
        if total:
            return np.sum(self.__hi_mass)
        else:
            return self.__hi_mass

    def get_h2_mass(self, total=False):
        """Return the :math:`\mathrm{H}_2` mass represented in the model."""
        if total:
            return np.sum(self.__h2_mass)
        else:
            return self.__h2_mass

    def get_volume_filling_factor(self):
        """Return the volume-filling factor."""
        return self.__volume_factor

    def get_voxel_filling_factor(self):
        """
        Return the voxel-filling factor (the fraction of the voxel filled
        by gas).
        """
        return self.__voxel_filling_factor

    def get_velocity(self):
        """
        Return a tuple of the voxel velocity and ensemble-averaged dispersion.
        """
        return self.__velocity, self.__ensemble_dispersion

    def get_clump_velocity(self):
        """
        Return a tuple of the velocity bins used in the calculation and the
        indeces.
        """
        return self.__clump_velocities, self.__clump_velocity_indeces

    def get_taufuv(self):
        """return the voxel-averaged far-UV optical depth."""
        return self.__taufuv

    def get_species_number(
        self, species=None, abun=False, nref=[("H", 1), ("H2", 2)], total=True
    ):
        """
        Return the either the total number of species or fractional abundance
        in each ensemble. The default functionality is to express the fractional
        in terms of total hydrogen (:math:`N_\mathrm{H} = N_\mathrm{H}^0 +
        2N_\mathrm{H}_2`).
        """
        if species in [None, "all"]:
            species = constants.abundances
        elif isinstance(species, str):
            species = [species]

        if nref == None:
            nref = ["H", "H2"]
        elif isinstance(nref, str):
            nref = [nref]

        N_species = []

        for ens in range(constants.ensembles):
            if abun:
                N_0 = 0
                for sp in nref:
                    N_0 += (
                        sp[1]
                        * self.clump_N_species[ens][constants.abundances.index(sp[0])]
                    )
            else:
                N_0 = 1

            N_sp = []
            for sp in species:
                N_sp.append(
                    self.clump_N_species[ens][constants.abundances.index(sp)] / N_0
                )

            N_species.append(copy(N_sp))

        if total:
            return np.sum(N_species, axis=0)
        else:
            return N_species
    def get_column_density(self, *args, **kwargs):

        """Return column densities."""
        area = (constants.voxel_size * 3.086* (10**(18)))**2
        cd = ((self.get_species_number(*args, **kwargs))/area)
        return cd

    def get_abundances(self, *args, **kwargs):
        """Return species abundances."""
        return self.get_species_number(*args, **kwargs, abun=True)

    def get_dust_temperature(self, total=True):
        """Return ensemble-averaged dust temperature."""
        if total:
            return (
                np.asarray(self.__ensemble_mass) * np.asarray(self.t_dust)
            ).sum() / np.asarray(self.__ensemble_mass).sum()
        else:
            return self.t_dust

    def get_gas_temperature(self, total=True):
        """Return ensemble-averaged gas temperature."""
        if total:
            return (
                np.asarray(self.__ensemble_mass) * np.asarray(self.t_gas)
            ).sum() / np.asarray(self.__ensemble_mass).sum()
        else:
            return self.t_gas

    # @jit
    def calculate_emission(
        self,
        taufuv=0,
        test_calc=False,
        test_opacity=False,
        test_pexp=False,
        test_fv=False,
        verbose=False,
        timed=False,
    ):
        """Calculate voxel-averaged absorption and emmisivity."""
        self.test_calc = test_calc
        self.test_opacity = test_opacity
        self.test_pexp = test_pexp
        self.test_fv = test_fv

        # if timed:
        t0 = time()

        # reinitialise properties
        self.__intensity_species = [
            np.zeros(
                (len(ensemble.clumpIndeces[_]), len(species.clump_transitions)),
                dtype=np.float64,
            )
            for _ in range(constants.ensembles)
        ]
        self.__optical_depth_species = [
            np.zeros(
                (len(ensemble.clumpIndeces[_]), len(species.clump_transitions)),
                dtype=np.float64,
            )
            for _ in range(constants.ensembles)
        ]
        self.__intensity_dust = [
            np.zeros(
                (
                    len(ensemble.clumpIndeces[_]),
                    constants.wavelengths[constants.n_dust].size,
                ),
                dtype=np.float64,
            )
            for _ in range(constants.ensembles)
        ]
        self.__optical_depth_dust = [
            np.zeros(
                (
                    len(ensemble.clumpIndeces[_]),
                    constants.wavelengths[constants.n_dust].size,
                ),
                dtype=np.float64,
            )
            for _ in range(constants.ensembles)
        ]
        self.__intensity_hi = [
            np.zeros((len(ensemble.clumpIndeces[_]), 1), dtype=np.float64)
            for _ in range(constants.ensembles)
        ]
        self.__optical_depth_hi = [
            np.zeros((len(ensemble.clumpIndeces[_]), 1), dtype=np.float64)
            for _ in range(constants.ensembles)
        ]

        masspoints.calculate_emission(taufuv=taufuv, timed=timed)

        if timed:
            t1 = time()
            self.__logger.info("Masspoint emission calculated: {}".format(t1 - t0))
        # print('masspoint emission calculated:', (time()-t0)/60)

        combinations.calculate_emission(
            test_calc=test_calc,
            test_opacity=test_opacity,
            test_fv=test_fv,
            f_v=self.__volume_factor,
            suggested_calc=self.suggested_calc,
        )

        if timed:
            t2 = time()
            self.__logger.info("Combination emission calculated: {}".format(t2 - t1))
        # print('combination emission calculated:', (time()-t0)/60)

        clumpIntensity = [[] for _ in range(len(constants.clump_mass_number))]
        clumpOpticalDepth = [[] for _ in range(len(constants.clump_mass_number))]
        clumpIntensityDust = [[] for _ in range(len(constants.clump_mass_number))]
        clumpOpticalDepthDust = [[] for _ in range(len(constants.clump_mass_number))]

        clumpIntensityHI = [[] for _ in range(len(constants.clump_mass_number))]
        clumpOpticalDepthHI = [[] for _ in range(len(constants.clump_mass_number))]

        # interclumpIntensity = []
        # interclumpOpticalDepth = []

        iDust = constants.wavelengths[constants.n_dust].size
        f_ds = [
            np.maximum(1, self.__volume_factor[ens])
            for ens in range(len(constants.clump_mass_number))
        ]

        # Clump
        for ens in range(constants.ensembles):
            self.__hi_mass[ens] = (
                masspoints.clump_hi_mass[ens] * ensemble.clumpNj[ens]
            ).sum()
            self.__h2_mass[ens] = (
                masspoints.clump_h2_mass[ens] * ensemble.clumpNj[ens]
            ).sum()

            vel = constants.velocity_range  # v_obs
            clumpVel = ensemble.clumpVelocities[ens]  # v_sys
            # print(vel.shape, clumpVel.shape)

            if self.suggested_calc:
                # print('new calculation:', (time()-t0)/60)
                # dust calculation
                optical_depth_comb_dust = copy(
                    combinations.clump_dust_optical_depth[ens]
                )
                intensity_comb_dust = copy(combinations.clump_dust_intensity[ens])
                i_nan = optical_depth_comb_dust < 1e-16
                intensity_comb_dust[~i_nan] = (
                    intensity_comb_dust[~i_nan]
                    / optical_depth_comb_dust[~i_nan]
                    # *constants.voxel_size
                    * (1 - np.exp(-optical_depth_comb_dust[~i_nan]))
                )
                intensity_comb_dust[i_nan] = (
                    intensity_comb_dust[i_nan] / constants.voxel_size
                )
                p_tot = ensemble.CLmaxProbability[ens].prod(
                    1
                ) / ensemble.CLmaxProbability[ens].prod(1).sum(0)
                clumpIntensityDust[ens].append((p_tot * intensity_comb_dust.T).sum(1))
                optical_depth_exp_dust = np.exp(-optical_depth_comb_dust)
                optical_depth_exp_dust[i_nan] = 1
                optical_depth_prob_dust = (p_tot * optical_depth_exp_dust.T).sum(1)
                i_zero = (optical_depth_exp_dust % 1 == 0).all(0)
                clumpOpticalDepthDust[ens].append(-np.log(optical_depth_prob_dust))
                clumpOpticalDepthDust[ens][-1][i_zero] = 0
                # clumpOpticalDepthDust[ens].append(-np.log((p_tot * np.exp(-optical_depth_comb_dust.T)).sum(1)))
                # transition calculation
                optical_depth_comb = copy(combinations.clump_species_optical_depth[ens])
                intensity_comb = copy(combinations.clump_species_intensity[ens])
                optical_depth_hi_comb = copy(combinations.clump_hi_tau[ens])
                intensity_hi_comb = copy(combinations.clump_hi_tb[ens])
                i_nan = optical_depth_comb < 1e-16
                i_nan_hi = optical_depth_hi_comb < 1e-16
                optical_depth_exp = np.exp(-optical_depth_comb)
                optical_depth_exp[i_nan] = 1
                optical_depth_hi_exp = np.exp(-optical_depth_hi_comb)
                optical_depth_hi_exp[i_nan_hi] = 1
                intensity_comb[~i_nan] = (
                    intensity_comb[~i_nan]
                    / optical_depth_comb[~i_nan]  # *constants.voxel_size
                    * (1 - np.exp(-optical_depth_comb[~i_nan]))
                )
                intensity_comb[i_nan] = intensity_comb[i_nan] / constants.voxel_size
                intensity_hi_comb[~i_nan_hi] = (
                    intensity_hi_comb[~i_nan_hi]
                    / optical_depth_hi_comb[~i_nan_hi]  # *constants.voxel_size
                    * (1 - np.exp(-optical_depth_hi_comb[~i_nan_hi]))
                )
                intensity_hi_comb[i_nan_hi] = (
                    intensity_hi_comb[i_nan_hi] / constants.voxel_size
                )
                for i, probability in enumerate(ensemble.clumpProbability[ens]):
                    p_i = probability.prod(1) / probability.prod(1).sum(0)
                    # clumpIntensity[ens].append((probability.prod(1)#/probability.prod(1).sum(0)
                    #                            intensity_comb.T))
                    clumpIntensity[ens].append((p_i * intensity_comb.T).sum(1))
                    optical_depth_prob = (p_i * optical_depth_exp.T).sum(1)
                    clumpIntensityHI[ens].append((p_i * intensity_hi_comb.T).sum(1))
                    optical_depth_hi_prob = (p_i * optical_depth_hi_exp.T).sum(1)
                    i_zero = (optical_depth_comb % 1 == 0).all(0)
                    i_zero_hi = (optical_depth_hi_comb % 1 == 0).all(0)
                    clumpOpticalDepth[ens].append(-np.log(optical_depth_prob))
                    clumpOpticalDepth[ens][-1][i_zero] = 0
                    clumpOpticalDepthHI[ens].append(-np.log(optical_depth_hi_prob))
                    clumpOpticalDepthHI[ens][-1][i_zero_hi] = 0
                    # clumpOpticalDepth[ens].append(-np.log((probability.prod(1)#/probability.prod(1).sum(0)
                    #                                       * np.exp(-optical_depth_comb.T)).sum(1)))

                # Transitions
                self.__intensity_species[ens] += np.asarray(clumpIntensity[ens])
                self.__optical_depth_species[ens] += np.asarray(clumpOpticalDepth[ens])
                self.__absorption_species[ens] = (
                    self.__optical_depth_species[ens] / constants.voxel_size / f_ds[ens]
                )
                i_undef = self.__optical_depth_species[ens] < 1e-16
                self.__emissivity_species[ens][~i_undef] = (
                    self.__intensity_species[ens][~i_undef]
                    * self.__absorption_species[ens][~i_undef]
                    / (1 - np.exp(-self.__optical_depth_species[ens][~i_undef]))
                )
                self.__emissivity_species[ens][i_undef] = self.__intensity_species[ens][
                    i_undef
                ]

                # Dust
                self.__intensity_dust[ens] += np.asarray(clumpIntensityDust[ens])
                self.__optical_depth_dust[ens] += np.asarray(clumpOpticalDepthDust[ens])
                self.__absorption_dust[ens] = (
                    self.__optical_depth_dust[ens] / constants.voxel_size / f_ds[ens]
                )
                i_undef = self.__absorption_dust[ens] == 0
                self.__emissivity_dust[ens][~i_undef] = (
                    self.__intensity_dust[ens][~i_undef]
                    * self.__absorption_dust[ens][~i_undef]
                    / (1 - np.exp(-self.__optical_depth_dust[ens][~i_undef]))
                )
                self.__emissivity_dust[ens][i_undef] = self.__intensity_dust[ens][
                    i_undef
                ]

                # HI
                self.__intensity_hi[ens] += np.asarray(clumpIntensityHI[ens])
                self.__optical_depth_hi[ens] += np.asarray(clumpOpticalDepthHI[ens])
                self.__absorption_hi[ens] = (
                    self.__optical_depth_hi[ens] / constants.voxel_size / f_ds[ens]
                )
                i_undef = self.__optical_depth_hi[ens] < 1e-16
                self.__emissivity_hi[ens][~i_undef] = (
                    self.__intensity_hi[ens][~i_undef]
                    * self.__absorption_hi[ens][~i_undef]
                    / (1 - np.exp(-self.__optical_depth_hi[ens][~i_undef]))
                )
                self.__emissivity_hi[ens][i_undef] = self.__intensity_hi[ens][i_undef]

            # The old computation that will soon be removed
            else:
                nv = np.abs(vel - self.__velocity) <= 4 * np.maximum(
                    self.__ensemble_dispersion[ens], constants.clump_dispersion
                )
                factor = (
                    1
                    / np.sqrt(2 * np.pi * constants.clump_dispersion**2)
                    * np.exp(
                        -(
                            (
                                vel[nv].reshape(1, -1)
                                - clumpVel.reshape(-1, 1)
                                - self.__velocity
                            )
                            ** 2
                        )
                        / 2
                        / constants.clump_dispersion**2
                    )
                )

                # clumpIntensity[ens] = np.zeros((ensemble.clumpVelocities[ens].size,
                #                                 constants.velocityRange.size,
                #                                 combinations.clumpIntensity[ens].shape[0],
                #                                 combinations.clumpIntensity[ens].shape[1]))
                # clumpIntensityDust[ens] = np.zeros((constants.velocityRange.size,
                #                                     combinations.clumpIntensity[ens].shape[0],
                #                                     combinations.clumpIntensity[ens].shape[1]))
                # clumpOpticalDepth[ens] = np.zeros((ensemble.clumpVelocities[ens].size,
                #                                    constants.velocityRange.size,
                #                                    combinations.clumpIntensity[ens].shape[0],
                #                                    combinations.clumpIntensity[ens].shape[1]))
                # clumpOpticalDepthDust[ens] = np.zeros((constants.velocityRange.size,
                #                                        combinations.clumpIntensity[ens].shape[0],
                #                                        combinations.clumpIntensity[ens].shape[1]))

                p_tot = ensemble.CLmaxProbability[ens].prod(
                    1
                ) / ensemble.CLmaxProbability[ens].prod(1).sum(0)

                if nv.any():
                    for i, probability in enumerate(ensemble.clumpProbability[ens]):
                        if timed:
                            t3 = time()
                            self.__logger.info(
                                "Start I_xi calculation:".format(t3 - t2)
                            )

                        p_i = probability.prod(1) / probability.prod(1).sum(0)

                        intensity = copy(combinations.clump_species_intensity[ens])
                        # shape (v_obs, combination, wavelength)
                        intensity = np.array(
                            [
                                intensity * factor[ensemble.clumpIndeces[ens][i], j]
                                for j in range(factor.shape[1])
                            ]
                        )
                        # shape (v_sys, v_obs, combination, wavelength)
                        clumpIntensity[ens].append(
                            np.array(
                                [
                                    (p_i * intensity[j].T).T
                                    for j in range(factor.shape[1])
                                ]
                            )
                        )
                        # shape (combination, wavelength)
                        intensityDust = copy(combinations.clump_dust_intensity[ens])
                        # shape (v_obs, combination, wavelength)
                        clumpIntensityDust[ens].append((p_tot * intensityDust.T).T)

                        opticalDepth = copy(
                            combinations.clump_species_optical_depth[ens]
                        )
                        opticalDepth = np.array(
                            [
                                opticalDepth * factor[ensemble.clumpIndeces[ens][i], j]
                                for j in range(factor.shape[1])
                            ]
                        )
                        opticalDepthDust = copy(
                            combinations.clump_dust_optical_depth[ens]
                        )
                        if self.test_pexp:
                            clumpOpticalDepthDust[ens].append(
                                (p_tot * opticalDepthDust.T).T
                            )
                            clumpOpticalDepth[ens].append(
                                np.array(
                                    [
                                        (p_i * opticalDepth[j].T).T
                                        for j in range(factor.shape[1])
                                    ]
                                )
                            )
                        elif self.test_opacity:
                            clumpOpticalDepthDust[ens].append(
                                (
                                    p_tot
                                    * np.exp(
                                        -opticalDepthDust.T
                                        * constants.voxel_size
                                        * f_ds[ens]
                                    )
                                ).T
                            )
                            clumpOpticalDepth[ens].append(
                                np.array(
                                    [
                                        (
                                            p_i
                                            * np.exp(
                                                -opticalDepth[j].T
                                                * constants.voxel_size
                                                * f_ds[ens]
                                            )
                                        ).T
                                        for j in range(factor.shape[1])
                                    ]
                                )
                            )
                        else:
                            clumpOpticalDepthDust[ens].append(
                                (p_tot * np.exp(-opticalDepthDust.T)).T
                            )
                            clumpOpticalDepth[ens].append(
                                np.array(
                                    [
                                        (p_i * np.exp(-opticalDepth[j].T)).T
                                        for j in range(factor.shape[1])
                                    ]
                                )
                            )

                        if timed:
                            t4 = time()
                            self.__logger.info("End I_xi calculation:".format(t4 - t3))

                    # All of these have shape (velocity, wavelength)
                    self.__intensity_species[ens][nv, :] = self.__intensity_species[
                        ens
                    ][nv, :] + (np.array(clumpIntensity[ens]).sum(2)).sum(0).astype(
                        constants.dtype
                    )
                    # self.__intensityDust[ens][:,:] = self.__intensityDust[ens][:,:] +
                    #                                  np.array([np.array(clumpIntensityDust[ens]).sum(1).sum(0)
                    #                                  for _ in range(factor.shape[1])]).astype(constants.dtype)
                    # self.__opticalDepthDust[ens][:,:] = self.__opticalDepthDust[ens][:,:] +
                    #                                     np.array([-np.log(np.array(clumpOpticalDepthDust[ens]).sum(1)).sum(0)
                    #                                     for _ in range(factor.shape[1])]).astype(constants.dtype)
                    #                                     (factor.shape[1])]).astype(constants.dtype)
                    self.__intensity_dust[ens][:, :] = self.__intensity_dust[ens][
                        :, :
                    ] + np.array(
                        [
                            np.array(clumpIntensityDust[ens]).sum(1).max(0)
                            for _ in range(self.__intensity_dust[ens].shape[0])
                        ]
                    ).astype(
                        constants.dtype
                    )
                    if self.test_pexp:
                        self.__optical_depth_species[ens][
                            nv, :
                        ] = self.__optical_depth_species[ens][nv, :] + (
                            np.array(clumpOpticalDepth[ens]).sum(2)
                        ).sum(
                            0
                        ).astype(
                            constants.dtype
                        )
                        self.__optical_depth_dust[ens][
                            :, :
                        ] = self.__optical_depth_dust[ens][:, :] + np.array(
                            [
                                np.array(clumpOpticalDepthDust[ens]).sum(1).max(0)
                                for _ in range(self.__optical_depth_dust[ens].shape[0])
                            ]
                        ).astype(
                            constants.dtype
                        )
                    else:
                        self.__optical_depth_species[ens][
                            nv, :
                        ] = self.__optical_depth_species[ens][nv, :] + (
                            -np.log(np.array(clumpOpticalDepth[ens]).sum(2))
                        ).sum(
                            0
                        ).astype(
                            constants.dtype
                        )
                        self.__optical_depth_dust[ens][
                            :, :
                        ] = self.__optical_depth_dust[ens][:, :] + np.array(
                            [
                                -np.log(
                                    np.array(clumpOpticalDepthDust[ens]).sum(1)
                                ).max(0)
                                for _ in range(self.__optical_depth_dust[ens].shape[0])
                            ]
                        ).astype(
                            constants.dtype
                        )

                    if iDust > 10:
                        intensityDustInterp = interp1d(
                            constants.wavelengths[constants.n_dust],
                            self.__intensity_dust[ens].max(0),
                            fill_value="extrapolate",
                        )
                        opticalDepthDustInterp = interp1d(
                            constants.wavelengths[constants.n_dust],
                            self.__optical_depth_dust[ens].max(0),
                            fill_value="extrapolate",
                        )

                    for i, transition in enumerate(
                        species.clump_transition_wavelengths
                    ):
                        if iDust > 10:
                            self.__intensity_species[ens][:, i] += intensityDustInterp(
                                transition
                            )
                            self.__optical_depth_species[ens][
                                :, i
                            ] += opticalDepthDustInterp(transition)
                        else:
                            self.__intensity_species[ens][
                                :, i
                            ] += self.__intensity_dust[ens].max()
                            self.__optical_depth_species[ens][
                                :, i
                            ] += self.__optical_depth_dust[ens].max()

                else:
                    self.__logger.warning(
                        "Voxel with velocity {} not within given observing velocity range.".format(
                            self.__velocity
                        )
                    )

        if self.suggested_calc:
            # Handled above
            pass
        elif self.test_calc:
            self.__emissivity_species = deepcopy(self.__intensity_species)
            self.__emissivity_dust = deepcopy(self.__intensity_dust)
        elif self.test_fv:
            f_ds = [
                np.maximum(1, self.__volume_factor[ens])
                for ens in range(len(constants.clump_mass_number))
            ]
            self.__emissivity_species = [
                self.__intensity_species[ens] / constants.voxel_size / f_ds[ens]
                for ens in range(len(constants.clump_mass_number))
            ]
            self.__emissivity_dust = [
                self.__intensity_dust[ens] / constants.voxel_size / f_ds[ens]
                for ens in range(len(constants.clump_mass_number))
            ]
        else:
            self.__emissivity_species = [
                self.__intensity_species[ens] / constants.voxel_size
                for ens in range(len(constants.clump_mass_number))
            ]
            self.__emissivity_dust = [
                self.__intensity_dust[ens] / constants.voxel_size
                for ens in range(len(constants.clump_mass_Number))
            ]
        if self.suggested_calc:
            # Handled above
            pass
        elif self.test_opacity:
            self.__absorption_species = deepcopy(self.__optical_depth_species)
            self.__absorption_dust = deepcopy(self.__optical_depth_dust)
        elif self.test_fv:
            f_ds = [
                np.maximum(1, self.__volume_factor[ens])
                for ens in range(len(constants.clump_mass_number))
            ]
            self.__absorption_species = [
                self.__optical_depth_species[ens] / constants.voxel_size / f_ds[ens]
                for ens in range(len(constants.clump_mass_number))
            ]
            self.__absorption_dust = [
                self.__optical_depth_dust[ens] / constants.voxel_size / f_ds[ens]
                for ens in range(len(constants.clump_mass_number))
            ]
        else:
            self.__absorption_species = [
                self.__optical_depth_species[ens] / constants.voxel_size
                for ens in range(len(constants.clump_mass_number))
            ]
            self.__absorption_dust = [
                self.__optical_depth_dust[ens] / constants.voxel_size
                for ens in range(len(constants.clump_mass_number))
            ]

        self.t_gas = [
            (
                ensemble.clumpNj[_]
                * 10 ** constants.clump_log_mass[_]
                * masspoints.clump_t_gas[_]
            ).sum()
            / (ensemble.clumpNj[_] * 10 ** constants.clump_log_mass[_]).sum()
            for _ in range(constants.ensembles)
        ]
        self.t_dust = [
            (
                ensemble.clumpNj[_]
                * 10 ** constants.clump_log_mass[_]
                * masspoints.clump_t_dust[_]
            ).sum()
            / (ensemble.clumpNj[_] * 10 ** constants.clump_log_mass[_]).sum()
            for _ in range(constants.ensembles)
        ]

        self.clump_N_species = [
            (ensemble.clumpNj[_].T * masspoints.clump_N_species[_]).sum(0)
            for _ in range(constants.ensembles)
        ]

        self.__logger.info(
            "NaN in species intensity: {}".format(
                np.isnan(self.__intensity_species).any()
            )
        )
        self.__logger.info(
            "NaN in species optical depth: {}".format(
                np.isnan(self.__optical_depth_species).any()
            )
        )
        self.__logger.info(
            "NaN in dust intensity: {}".format(np.isnan(self.__intensity_dust).any())
        )
        self.__logger.info(
            "NaN in dust optical depth: {}".format(
                np.isnan(self.__optical_depth_dust).any()
            )
        )

        self.__logger.info("Voxel emission calculated.")

        if timed:
            self.__logger.info(
                "calculateEmission() time of execution: {}".format(time() - t0)
            )

        return

    def gaussian(self, x, a, sigma):
        """Return Gaussian."""
        return a * np.exp(-(x**2) / (2 * sigma**2))

    def two_gaussians(self, x, a1, a2, sigma):
        """Return two Gaussians stacked horizontally."""
        g1 = self.gaussian(x, a1, sigma)
        g2 = self.gaussian(x, a2, sigma)
        return np.hstack((g1, g2))

    def get_species_emissivity(
        self, kind="linear", include_dust=False, fit=True, total=True, hi=False
    ):
        """Return line emission emissivity."""
        if self.suggested_calc:
            epsilon_species = []
            if include_dust:
                if False:
                    epsilon_dust = [self.get_dust_emissivity(total=total)]
                else:
                    epsilon_dust = self.get_dust_emissivity(total=False)
            else:
                pass
            for ens in range(len(constants.clump_mass_number)):
                if fit:
                    clump_vel = (
                        self.__clump_velocities[ens][self.__clump_velocity_indeces[ens]]
                        + self.__velocity
                    )
                    if kind == "gaussian":
                        eps = []
                        if hi:
                            diff_eps = self.__emissivity_hi[ens].var(0)
                            diff_kap = self.__absorption_hi[ens].var(0)
                        else:
                            diff_eps = self.__emissivity_species[ens].var(0)
                            diff_kap = self.__absorption_species[ens].var(0)
                        for i in range(self.__emissivity_species[ens].shape[1]):
                            if diff_eps[i] or diff_kap[i]:
                                a_eps, a_kap, sigma = curve_fit(
                                    self.two_gaussians,
                                    clump_vel,
                                    np.hstack(
                                        (
                                            self.__emissivity_species[ens][:, i],
                                            self.__absorption_species[ens][:, i],
                                        )
                                    ),
                                    p0=(1, 1, self.__ensembleDispersion[ens]),
                                )[0]
                                eps.append(
                                    self.gaussian(
                                        constants.velocity_range, a_eps, sigma
                                    )
                                )
                            else:
                                eps.append(np.zeros_like(constants.velocity_range))
                        epsilon_species.append(np.asarray(eps).T)
                    else:
                        # eps = []
                        # for i in range(self.__emissivity_species[ens].shape[1]):
                        #     fn = interp1d(clump_vel, self.__emissivity_species[ens][:, i],
                        #                   kind=kind, fill_value=0, bounds_error=False)
                        #     eps.append(fn(constants.velocityRange))
                        # epsilon_species.append(np.asarray(eps).T)
                        if hi:
                            fn = interp1d(
                                clump_vel,
                                self.__emissivity_hi[ens],
                                kind=kind,
                                fill_value=0,
                                bounds_error=False,
                                axis=0,
                            )
                        else:
                            fn = interp1d(
                                clump_vel,
                                self.__emissivity_species[ens],
                                kind=kind,
                                fill_value=0,
                                bounds_error=False,
                                axis=0,
                            )
                        epsilon_species.append(fn(constants.velocity_range))
                    if include_dust:
                        if np.where(constants.n_dust)[0].size > 10:
                            epsilon_dust_interp = interp1d(
                                constants.wavelengths[constants.n_dust],
                                epsilon_dust[ens].max(0),
                                fill_value="extrapolate",
                            )
                            if hi:
                                epsilon_species[ens][:, 0] += epsilon_dust_interp(
                                    0.21106
                                )
                            else:
                                for i, transition in enumerate(
                                    species.clump_transition_wavelengths
                                ):
                                    epsilon_species[ens][:, i] += epsilon_dust_interp(
                                        transition
                                    )
                        else:
                            if hi:
                                epsilon_species[ens][:, 0] += epsilon_dust[ens].mean()
                            else:
                                for i in range(len(species.clump_transitions)):
                                    epsilon_species[ens][:, i] += epsilon_dust[
                                        ens
                                    ].mean()
                else:
                    if hi:
                        epsilon_species.append(copy(self.__emissivity_hi[ens]))
                    else:
                        epsilon_species.append(copy(self.__emissivity_species[ens]))
        else:
            epsilon_species = deepcopy(self.__emissivity_species)
            if not include_dust:
                epsilon_dust = self.__emissivity_dust
                for ens in range(len(constants.clump_mass_number)):
                    if np.where(constants.n_dust)[0].size > 1:
                        epsilon_dust_interp = interp1d(
                            constants.wavelengths[constants.n_dust],
                            epsilon_dust[ens].max(0),
                            fill_value="extrapolate",
                        )
                        for i, transition in enumerate(
                            species.clump_transition_wavelengths
                        ):
                            if epsilon_dust[ens].shape[1] > 1:
                                epsilon_species[ens][:, i] -= epsilon_dust_interp(
                                    transition
                                )
                    else:
                        for i in range(len(species.clump_transitions)):
                            epsilon_species[ens][:, i] -= epsilon_dust[ens].max(0)
        if total:
            return np.asarray(epsilon_species).sum(0)
        else:
            return np.asarray(epsilon_species)

    def get_species_absorption(
        self, kind="linear", include_dust=False, fit=True, total=True, hi=False
    ):
        """Return line emission absorption."""
        if self.suggested_calc:
            kappa_species = []
            if include_dust:
                if False:
                    kappa_dust = [self.get_dust_absorption(total=total)]
                else:
                    kappa_dust = self.get_dust_absorption(total=False)
            else:
                pass
            for ens in range(len(constants.clump_mass_number)):
                if fit:
                    clump_vel = (
                        self.__clump_velocities[ens][self.__clump_velocity_indeces[ens]]
                        + self.__velocity
                    )
                    if kind == "gaussian":
                        kap = []
                        if hi:
                            diff_eps = self.__emissivity_hi[ens].var(0)
                            diff_kap = self.__absorption_hi[ens].var(0)
                        else:
                            diff_eps = self.__emissivity_species[ens].var(0)
                            diff_kap = self.__absorption_species[ens].var(0)
                        for i in range(self.__absorption_species[ens].shape[1]):
                            if diff_eps[i] or diff_kap[i]:
                                if hi:
                                    a_eps, a_kap, sigma = curve_fit(
                                        self.two_gaussians,
                                        clump_vel,
                                        np.hstack(
                                            (
                                                self.__emissivity_hi[ens][:, i],
                                                self.__absorption_hi[ens][:, i],
                                            )
                                        ),
                                        p0=(1, 1, self.__ensembleDispersion[ens]),
                                    )[0]
                                else:
                                    a_eps, a_kap, sigma = curve_fit(
                                        self.two_gaussians,
                                        clump_vel,
                                        np.hstack(
                                            (
                                                self.__emissivity_species[ens][:, i],
                                                self.__absorption_species[ens][:, i],
                                            )
                                        ),
                                        p0=(1, 1, self.__ensembleDispersion[ens]),
                                    )[0]
                                kap.append(
                                    self.gaussian(
                                        constants.velocity_range, a_kap, sigma
                                    )
                                )
                            else:
                                kap.append(np.zeros_like(constants.velocity_range))
                        kappa_species.append(np.asarray(kap).T)
                    else:
                        # kap = []
                        # for i in range(self.__absorption_species[ens].shape[1]):
                        #     fn = interp1d(clump_vel, self.__absorption_species[ens][:, i],
                        #                   kind=kind, fill_value=0, bounds_error=False)
                        #     kap.append(fn(constants.velocityRange))
                        # kappa_species.append(np.asarray(kap).T)
                        if hi:
                            fn = interp1d(
                                clump_vel,
                                self.__absorption_hi[ens],
                                kind=kind,
                                fill_value=0,
                                bounds_error=False,
                                axis=0,
                            )
                        else:
                            fn = interp1d(
                                clump_vel,
                                self.__absorption_species[ens],
                                kind=kind,
                                fill_value=0,
                                bounds_error=False,
                                axis=0,
                            )
                        kappa_species.append(fn(constants.velocity_range))
                    if include_dust:
                        if np.where(constants.n_dust)[0].size > 10:
                            kappa_dust_interp = interp1d(
                                constants.wavelengths[constants.n_dust],
                                kappa_dust[ens].max(0),
                                fill_value="extrapolate",
                            )
                            if hi:
                                kappa_species[ens][:, 0] += kappa_dust_interp(0.21106)
                            else:
                                for i, transition in enumerate(
                                    species.clump_transition_wavelengths
                                ):
                                    kappa_species[ens][:, i] += kappa_dust_interp(
                                        transition
                                    )
                        else:
                            if hi:
                                kappa_species[ens][:, 0] += kappa_dust[ens].mean()
                            else:
                                for i in range(len(species.clump_transitions)):
                                    kappa_species[ens][:, i] += kappa_dust[ens].mean()
                else:
                    if hi:
                        kappa_species.append(copy(self.__absorption_hi[ens]))
                    else:
                        kappa_species.append(copy(self.__absorption_species[ens]))
        else:
            kappa_species = deepcopy(self.__absorption_species)
            if not include_dust:
                kappa_dust = self.__absorption_dust
                for ens in range(len(constants.clump_mass_number)):
                    if np.where(constants.n_dust)[0].size > 1:
                        kappa_dust_interp = interp1d(
                            constants.wavelengths[constants.n_dust],
                            kappa_dust[ens].max(0),
                            fill_value="extrapolate",
                        )
                        for i, transition in enumerate(
                            species.clump_transition_wavelengths
                        ):
                            if kappa_dust[ens].shape[1] > 1:
                                kappa_species[ens][:, i] -= kappa_dust_interp(
                                    transition
                                )
                    else:
                        for i in range(len(species.clump_transitions)):
                            kappa_species[ens][i] -= kappa_dust[ens].max(0)
        if total:
            return np.asarray(kappa_species).sum(0)
        else:
            return np.asarray(kappa_species)

    def get_species_optical_depth(
        self, kind="linear", include_dust=False, fit=True, total=True, hi=False
    ):
        """Return line emission optical depth."""
        if self.suggested_calc:
            tau_species = []
            if include_dust:
                if False:
                    tau_dust = [self.get_dust_optical_depth(total=total)]
                else:
                    tau_dust = self.get_dust_optical_depth(total=False)
            for ens in range(len(constants.clump_mass_number)):
                if fit:
                    clump_vel = (
                        self.__clump_velocities[ens][self.__clump_velocity_indeces[ens]]
                        + self.__velocity
                    )
                    if kind == "gaussian":
                        tau = []
                        if hi:
                            diff_eps = self.__emissivity_hi[ens].var(0)
                            diff_kap = self.__absorption_hi[ens].var(0)
                        else:
                            diff_eps = self.__emissivity_species[ens].var(0)
                            diff_kap = self.__absorption_species[ens].var(0)
                        for i in range(self.__opticalDepthSpecies[ens].shape[1]):
                            if diff_eps[i] or diff_kap[i]:
                                if hi:
                                    a_eps, a_kap, sigma = curve_fit(
                                        self.two_gaussians,
                                        clump_vel,
                                        np.hstack(
                                            (
                                                self.__emissivity_species[ens][:, i],
                                                self.__absorption_species[ens][:, i],
                                            )
                                        ),
                                        p0=(1, 1, self.__ensembleDispersion[ens]),
                                    )[0]
                                else:
                                    a_eps, a_kap, sigma = curve_fit(
                                        self.two_gaussians,
                                        clump_vel,
                                        np.hstack(
                                            (
                                                self.__emissivity_species[ens][:, i],
                                                self.__absorption_species[ens][:, i],
                                            )
                                        ),
                                        p0=(1, 1, self.__ensembleDispersion[ens]),
                                    )[0]
                                tau.append(
                                    self.gaussian(
                                        constants.velocity_range, a_kap, sigma
                                    )
                                    * constants.voxel_size
                                )
                            else:
                                tau.append(np.zeros_like(constants.velocity_range))
                        tau_species.append(np.asarray(tau).T)
                    else:
                        # tau = []
                        # for i in range(self.__opticalDepthSpecies[ens].shape[1]):
                        #     fn = interp1d(clump_vel, self.__opticalDepthSpecies[ens][:, i],
                        #                   kind=kind, fill_value=0, bounds_error=False)
                        #     tau.append(fn(constants.velocityRange))
                        # tau_species.append(np.asarray(tau).T)
                        if hi:
                            fn = interp1d(
                                clump_vel,
                                self.__optical_depth_hi[ens],
                                kind=kind,
                                fill_value=0,
                                bounds_error=False,
                                axis=0,
                            )
                        else:
                            fn = interp1d(
                                clump_vel,
                                self.__optical_depth_species[ens],
                                kind=kind,
                                fill_value=0,
                                bounds_error=False,
                                axis=0,
                            )
                        tau_species.append(fn(constants.velocity_range))
                    if include_dust:
                        if np.where(constants.n_dust)[0].size > 10:
                            tau_dust_interp = interp1d(
                                constants.wavelengths[constants.n_dust],
                                tau_dust[ens].max(0),
                                fill_value="extrapolate",
                            )
                            if hi:
                                tau_species[ens][:, 0] += tau_dust_interp(0.21106)
                            else:
                                for i, transition in enumerate(
                                    species.clump_transition_wavelengths
                                ):
                                    tau_species[ens][:, i] += tau_dust_interp(
                                        transition
                                    )
                        else:
                            if hi:
                                tau_species[ens][:, 0] += tau_dust[ens].mean()
                            else:
                                for i in range(len(species.clump_transitions)):
                                    tau_species[ens][:, i] += tau_dust[ens].mean()
                else:
                    if hi:
                        tau_species.append(copy(self.__optical_depth_hi[ens]))
                    else:
                        tau_species.append(copy(self.__optical_depth_species[ens]))
        # else:
        #    return
        if total:
            return np.asarray(tau_species).sum(0)
        else:
            return np.asarray(tau_species)

    def get_species_intensity(
        self,
        integrated=False,
        kind="linear",
        include_dust=False,
        fit=True,
        total=True,
        hi=False,
    ):
        """
        Return line emission intensity. This integrated the radiative transfer
        equation over the length of the voxel.
        """
        if self.suggested_calc:
            if total and not integrated and not fit:
                print(
                    "Cannot return total intensity without first fitting. Returning fitted intensity instead."
                )
                fit = True
            intensity_final = []
            if False:  # total:
                epsilon_species = [
                    self.get_species_emissivity(
                        kind=kind,
                        include_dust=include_dust,
                        fit=fit,
                        total=total,
                        hi=hi,
                    )
                ]
                kappa_species = [
                    self.get_species_absorption(
                        kind=kind,
                        include_dust=include_dust,
                        fit=fit,
                        total=total,
                        hi=hi,
                    )
                ]
                tau_species = [
                    self.get_species_optical_depth(
                        kind=kind,
                        include_dust=include_dust,
                        fit=fit,
                        total=total,
                        hi=hi,
                    )
                ]
                ensembles = 1
            else:
                epsilon_species = self.get_species_emissivity(
                    kind=kind, include_dust=True, fit=fit, total=False, hi=hi
                )
                kappa_species = self.get_species_absorption(
                    kind=kind, include_dust=True, fit=fit, total=False, hi=hi
                )
                tau_species = self.get_species_optical_depth(
                    kind=kind, include_dust=True, fit=fit, total=False, hi=hi
                )
                intensity_dust = self.get_dust_intensity(fit=fit, total=False)
                ensembles = constants.ensembles
            for ens in range(ensembles):
                clump_vel = (
                    self.__clump_velocities[ens][self.__clump_velocity_indeces[ens]]
                    + self.__velocity
                )
                if fit:
                    vel = constants.velocity_range
                    # if kind == 'gaussian':
                    i_nan = kappa_species[ens] < 1e-16
                    intensity = copy(epsilon_species[ens])
                    intensity[~i_nan] = (
                        intensity[~i_nan]
                        / kappa_species[ens][~i_nan]
                        * (1 - np.exp(-tau_species[ens][~i_nan]))
                    )
                    # else:
                    #     intensity_species = self.__intensity_species[ens]
                    #     intensity_interp = interp1d(clump_vel, intensity_species,
                    #                                 kind=kind, fill_value=0, bounds_error=False, axis=0)
                    #     intensity = intensity_interp(vel)
                else:
                    vel = clump_vel
                    intensity = copy(self.__intensity_species[ens])
                if integrated:
                    if np.where(constants.n_dust)[0].size > 10:
                        intensity_dust_interp = interp1d(
                            constants.wavelengths[constants.n_dust],
                            intensity_dust[ens].max(0),
                            fill_value="extrapolate",
                        )
                        for i, transition in enumerate(
                            species.clump_transition_wavelengths
                        ):
                            intensity[:, i] -= intensity_dust_interp(transition)
                    else:
                        for i in range(len(species.clump_transitions)):
                            intensity[:, i] -= intensity_dust[ens].mean()
                    intensity_final.append(np.trapz(intensity, vel, axis=0))
                    if include_dust:
                        if np.where(constants.n_dust)[0].size > 10:
                            intensity_dust_interp = interp1d(
                                constants.wavelengths[constants.n_dust],
                                intensity_dust[ens].max(0),
                                fill_value="extrapolate",
                            )
                            if hi:
                                intensity_final[-1][0] += intensity_dust_interp(0.21106)
                            else:
                                for i, transition in enumerate(
                                    species.clump_transition_wavelengths
                                ):
                                    intensity_final[-1][i] += intensity_dust_interp(
                                        transition
                                    )
                        else:
                            if hi:
                                intensity_final[-1][0] += intensity_dust[ens].mean()
                            else:
                                for i in range(len(species.clump_transitions)):
                                    intensity_final[-1][i] += intensity_dust[ens].mean()
                else:
                    if include_dust:
                        intensity_final.append(copy(intensity))
                    else:
                        intensity_dust = self.get_dust_intensity(fit=fit, total=False)
                        intensity_final.append(copy(intensity))
                        if np.where(constants.n_dust)[0].size > 10:
                            intensity_dust_interp = interp1d(
                                constants.wavelengths[constants.n_dust],
                                intensity_dust[ens].max(0),
                                fill_value="extrapolate",
                            )
                            if hi:
                                intensity_final[-1][:, 0] -= intensity_dust_interp(
                                    0.21106
                                )
                            else:
                                for i, transition in enumerate(
                                    species.clump_transition_wavelengths
                                ):
                                    intensity_final[-1][:, i] -= intensity_dust_interp(
                                        transition
                                    )
                        else:
                            if hi:
                                intensity_final[-1][:, 0] -= intensity_dust[ens].mean()
                            else:
                                for i in range(len(species.clump_transitions)):
                                    intensity_final[-1][:, i] -= intensity_dust[
                                        ens
                                    ].mean()
        else:
            if total:
                epsilon = [self.get_species_emissivity(include_dust=True, total=total)]
                kappa = [self.get_species_absorption(include_dust=True, total=total)]
                intensity_dust = [self.get_dust_intensity(total=total)]
            else:
                epsilon = self.get_species_emissivity(include_dust=True, total=total)
                kappa = self.get_species_absorption(include_dust=True, total=total)
                intensity_dust = self.get_dust_intensity(total=total)
            intensity = np.zeros_like(epsilon)

            for ens in range(constants.ensembles):
                if self.__volume_factor[ens] > 1 and self.test_fv:
                    ds = self.__volume_factor[ens] * constants.voxel_size
                else:
                    ds = constants.voxel_size
                intensity[ens] = (
                    epsilon[ens] / kappa[ens] * (1 - np.exp(-kappa[ens] * ds))
                )
                i_nan = np.isnan(intensity[ens])
                intensity[ens][i_nan] = epsilon[ens][i_nan]

                if not include_dust or (include_dust and integrated):
                    if np.where(constants.n_dust)[0].size > 1:
                        intensity_dust_interp = interp1d(
                            constants.wavelengths[constants.n_dust],
                            intensity_dust[ens].max(0),
                            fill_value="extrapolate",
                        )
                        for i, transition in enumerate(
                            species.clump_transition_wavelengths
                        ):
                            if intensity_dust[ens].shape[1] > 1:
                                intensity[ens][:, i] -= intensity_dust_interp(
                                    transition
                                )
                    else:
                        intensity[ens] -= intensity_dust[ens].max(0)

            if integrated:
                intensity_final = np.zeros((intensity.shape[0], intensity.shape[2]))
                for ens in range(len(constants.clump_mass_number)):
                    intensity_final[ens] = np.trapz(
                        intensity[ens], constants.velocity_range, axis=0
                    )
                    if include_dust:
                        if np.where(constants.n_dust)[0].size > 1:
                            intensity_dust_interp = interp1d(
                                constants.wavelengths[constants.n_dust],
                                intensity_dust[ens].max(0),
                                fill_value="extrapolate",
                            )
                            for i, transition in enumerate(
                                species.clump_transition_wavelengths
                            ):
                                intensity_final[ens][i] += intensity_dust_interp(
                                    transition
                                )
                        else:
                            intensity_final[ens] += intensity_dust[ens].max(0)
            else:
                intensity_final = intensity
        if total:
            return np.asarray(intensity_final).sum(0)
        else:
            return np.asarray(intensity_final)

    def get_dust_emissivity(self, fit=True, total=True, minimal=False):
        """Return the dust continuum emissivity."""
        if self.suggested_calc:
            if fit:
                emissivity = [
                    (
                        [self.__emissivity_dust[ens].max(0)]
                        * constants.velocity_range.size
                    )
                    for ens in range(constants.ensembles)
                ]
            else:
                emissivity = copy(self.__emissivity_dust)
        else:
            emissivity = self.__emissivity_dust
        if minimal:
            return np.asarray(emissivity).sum(0).max(0)
        elif total:
            return np.asarray(emissivity).sum(0)
        else:
            return np.asarray(emissivity)

    def get_dust_absorption(self, fit=True, total=True, minimal=False):
        """Return the dust continuum absorption."""
        if self.suggested_calc:
            if fit:
                absorption = [
                    (
                        [self.__absorption_dust[ens].max(0)]
                        * constants.velocity_range.size
                    )
                    for ens in range(constants.ensembles)
                ]
            else:
                absorption = copy(self.__absorption_dust)
        else:
            absorption = self.__absorption_dust
        if minimal:
            return np.asarray(absorption).sum(0).max(0)
        elif total:
            return np.asarray(absorption).sum(0)
        else:
            return np.asarray(absorption)

    def get_dust_optical_depth(self, fit=True, total=True, minimal=False):
        """Return the dust continuum optical depth."""
        if self.suggested_calc:
            if fit:
                optical_depth = [
                    (
                        [self.__optical_depth_dust[ens].max(0)]
                        * constants.velocity_range.size
                    )
                    for ens in range(constants.ensembles)
                ]
            else:
                optical_depth = copy(self.__optical_depth_dust)
        else:
            optical_depth = self.__optical_depth_dust
        if minimal:
            return np.asarray(optical_depth).sum(0).max(0)
        elif total:
            return np.asarray(optical_depth).sum(0)
        else:
            return np.asarray(optical_depth)

    def get_dust_intensity(self, fit=True, total=True, minimal=False):
        """
        Return the dust continuum intensity. This integrates the radiative
        transfer equation over the length of the voxel.
        """
        if self.suggested_calc:
            if fit:
                if total:
                    intensity = [
                        np.sum(
                            [
                                (
                                    [self.__intensity_dust[ens].max(0)]
                                    * constants.velocity_range.size
                                )
                                for ens in range(constants.ensembles)
                            ],
                            axis=0,
                        )
                    ]
                else:
                    intensity = np.array(
                        [
                            (
                                [self.__intensity_dust[ens].max(0)]
                                * constants.velocity_range.size
                            )
                            for ens in range(constants.ensembles)
                        ]
                    )
            else:
                if total:
                    intensity = np.sum(self.__intensity_dust, axis=0)
                else:
                    intensity = copy(self.__intensity_dust)
        else:
            if False:  # total:
                epsilon = [self.get_dust_emissivity(total=total)]
                kappa = [self.get_dust_absorption(total=total)]
            else:
                epsilon = self.get_dust_emissivity(total=total)
                kappa = self.get_dust_absorption(total=total)
            intensity = np.zeros_like(epsilon)
            for ens in range(len(constants.clump_mass_number)):
                if self.__volume_factor[ens] > 1 and self.test_fv:
                    ds = self.__volume_factor[ens] * constants.voxel_size
                else:
                    ds = constants.voxel_size
                intensity[ens] = (
                    epsilon[ens] / kappa[ens] * (1 - np.exp(-kappa[ens] * ds))
                )
                i_nan = np.isnan(intensity[ens])
                intensity[ens][i_nan] = epsilon[ens][i_nan]
        if minimal:
            return np.asarray(intensity).sum(0).max(0)
        elif total:
            return np.asarray(intensity).sum(0)
        else:
            return np.asarray(intensity)

    def plot_clump_number(self, effective=False, mass=False):
        """
        Plot the number of clumpy used as a function of velocity bin.

        This requires further testing...
        """
        import matplotlib.pyplot as plt

        xlabel = r"$v_i \ \left( \mathrm{\frac{km}{s}} \right)$"
        if effective or not self.suggested_calc:
            ylabel = r"$\mathcal{N}_mathrm{eff}$"
            title = "Clump ensemble {} effective number"
        else:
            ylabel = r"$\mathcal{N}$"
            title = "Clump ensemble {} number"

        for ens in range(constants.ensembles):
            dNj_dvi = combinations.clump_combination[ens].sum(1)  # Sum over clump types
            if effective or not self.suggested_calc:
                pass
            else:
                dNj_dvi = dNj_dvi / np.sqrt(2 * np.pi * constants.clump_dispersion**2)
            pji = copy(ensemble.clumpProbability[ens])
            vel = copy(ensemble.clumpVelocities[ens][ensemble.clumpIndeces[ens]])
            N = []
            for i, probability in enumerate(pji):
                pi = probability.prod(1) / probability.prod(1).sum(
                    0
                )  # Product over clump types and normalise
                N.append((pi * dNj_dvi).sum())  # Sum over combinations

            print("Total number: {:.4f}".format(np.trapz(N, vel)))

            fig, ax = plt.subplots(1, 1, figsize=(10, 5))
            ax.plot(vel, N, marker="o", ms=2, ls="-", lw=1)
            ax.set_xlabel(xlabel, fontsize=20)
            ax.set_ylabel(ylabel, fontsize=20)
            ax.set_title(title.format(ens), fontsize=20)
            plt.setp(ax.get_xticklabels(), fontsize=14)
            plt.setp(ax.get_yticklabels(), fontsize=14)
            fig.tight_layout()
            plt.show()

        return ax

    def plot_transition(
        self,
        transition="",
        quantity="intensity",
        kind="linear",
        include_dust=False,
        total=False,
        label="",
        title="",
        logscale=False,
    ):
        """
        Plot the transition emission with respect to the voxel velocity
        structure.
        """

        import matplotlib.pyplot as plt

        if transition in species.clump_transitions and isinstance(transition, str):
            transition = [transition]

        elif isinstance(transition, list):
            pass

        else:
            transition = species.clump_transitions

        if quantity == "emissivity":
            value = self.get_species_emissivity(
                kind=kind, include_dust=include_dust, total=total
            )
            ylabel = r"$\epsilon_{\lambda} \ \left( \mathrm{\frac{K}{pc}} \right)$"
        elif quantity == "absorption":
            value = self.get_species_absorption(
                kind=kind, include_dust=include_dust, total=total
            )
            ylabel = r"$\kappa_{\lambda} \ \left( \mathrm{\frac{1}{pc}} \right)$"
        elif quantity == "optical depth":
            value = self.get_species_optical_depth(
                kind=kind, include_dust=include_dust, total=total
            )
            ylabel = r"$\tau_{\lambda} \ \left( \mathrm{N/A} \right)$"
        elif quantity == "intensity":
            value = self.get_species_intensity(
                kind=kind, include_dust=include_dust, total=total
            )
            ylabel = r"$I_{\lambda} \ \left( \mathrm{K} \right)$"
        else:
            self.__logger.warning("Quantity {} not available.".format(quantity))
            return

        if total:
            value = [value]

        fig, axes = plt.subplots(len(value), figsize=(10, 5 * len(value)))

        for ens in range(len(value)):
            if isinstance(axes, np.ndarray):
                ax = axes[ens]

            else:
                ax = axes

            vel = constants.velocity_range  # [self.__clumpVelocityIndeces[ens]]
            labels = []

            for n, trans in enumerate(transition):
                if trans not in species.clump_transitions:
                    self.__logger.warning("Species {} not in model.".format(trans))
                    continue

                if constants.interclump_idx[ens]:
                    i = np.where(trans == np.asarray(species.interclump_transitions))[
                        0
                    ][0]
                else:
                    i = np.where(trans == np.asarray(species.clump_transitions))[0][0]

                if logscale:
                    ax.semilogy(vel, value[ens][:, i], marker="o", ms=2, ls="-", lw=1)
                else:
                    ax.plot(vel, value[ens][:, i], marker="o", ms=2, ls="-", lw=1)

                if isinstance(label, list):
                    labels.append(label[n])
                elif label:
                    labels.append(label)
                else:
                    labels.append(trans)

            ax.legend(
                labels=labels, fontsize=16, bbox_to_anchor=(1, 1), loc="upper left"
            )
            ax.set_xlabel(r"$V_r \ ( \mathrm{\frac{km}{s}} )$", fontsize=24)
            ax.set_ylabel(ylabel, fontsize=24)
            plt.setp(ax.get_xticklabels(), fontsize=16)
            plt.setp(ax.get_yticklabels(), fontsize=16)

            if title:
                ax.set_title(title, fontsize=20)
            else:
                ax.set_title("Clump set {} {}".format(ens + 1, quantity), fontsize=24)

        fig.tight_layout()

        # plt.show()

        return axes

    def plot_spectrum(
        self,
        quantity="intensity",
        kind="linear",
        integrated=False,
        total=False,
        vel=None,
        title="",
    ):
        """
        Plot the either the intesity or optical depth spectrum at the voxel
        velocity (the velocity with the largest ensemble).
        """

        import matplotlib as mpl
        import matplotlib.pyplot as plt
        from matplotlib.lines import Line2D

        # mpl.rcParams['text.usetex'] = True

        species_colour = {
            "C+": "xkcd:orange",
            "C": "xkcd:red",
            "O": "xkcd:violet",
            "CO": "xkcd:maroon",
            "13C+": "xkcd:light blue",
            "13C": "xkcd:blue",
            "13CO": "xkcd:sapphire",
            "HCO+": "xkcd:neon green",
            "H13CO+": "xkcd:forrest green",
            "H3O+": "xkcd:vomit",
        }

        nDust = constants.wavelengths[constants.n_dust].size

        if quantity == "emissivity":
            valueDust = self.get_dust_emissivity(total=False)
            valueSpecies = self.get_species_emissivity(
                kind=kind, include_dust=True, total=False
            )
            ylabel = r"$\epsilon_{\lambda} \ \left( \mathrm{\frac{K}{pc}} \right)$"

        elif quantity == "absorption":
            valueDust = self.get_dust_absorption(total=False)
            valueSpecies = self.get_species_absorption(
                kind=kind, include_dust=True, total=False
            )
            ylabel = r"$\kappa_{\lambda} \ \left( \mathrm{\frac{1}{pc}} \right)$"

        elif quantity == "intensity":
            valueDust = self.get_dust_intensity(total=False)
            if integrated:
                valueSpecies = self.get_species_intensity(
                    kind=kind, include_dust=True, integrated=integrated, total=False
                )
                ylabel = r"$\varpi_{\lambda} \ \left( \mathrm{K \frac{km}{s}} \right)$"
            else:
                valueSpecies = self.get_species_intensity(
                    kind=kind, include_dust=True, total=False
                )
                ylabel = r"$I_{\lambda} \ \left( \mathrm{K} \right)$"

        else:
            self.__logger.warning("Quantity {} not available.".format(quantity))
            return

        fig, axes = plt.subplots(
            constants.ensembles, figsize=(10, 5 * constants.ensembles)
        )

        for ens in range(constants.ensembles):
            if isinstance(axes, np.ndarray):
                ax = axes[ens]

            else:
                ax = axes

            # if constants.interclump_idx[ens]:
            #     constants.resort_wavelengths(interclump=True)
            # else:
            #     constants.resort_wavelengths()

            if integrated and quantity == "intensity":
                vel = 0
                valueSpecies = np.asarray([valueSpecies])
            elif vel == None:
                vel = int(valueSpecies[ens].shape[0] / 2)

            if valueDust[ens].shape[1] > 1:
                dustInterp = interp1d(
                    constants.wavelengths[constants.n_dust],
                    valueDust[ens][vel, :],
                    fill_value="extrapolate",
                )
                ax.loglog(
                    1e6 * constants.wavelengths[constants.n_dust],
                    valueDust[ens][vel, :],
                    ls="--",
                    lw=0.5,
                    c="xkcd:black",
                )
            else:
                dust_spectrum = np.full(2, valueDust[ens][vel, 0])
                ax.loglog(
                    1e6 * constants.wavelengths[[0, 20]],
                    dust_spectrum,
                    ls="--",
                    lw=0.5,
                    c="xkcd:black",
                )

            transition_species = []
            colour = []

            for i, transition in enumerate(species.clump_transitions):
                if valueDust[ens].shape[1] > 1:
                    valTemp = dustInterp(species.clump_transition_wavelengths[i])
                else:
                    valTemp = valueDust[ens][vel, 0]
                if not transition.split()[0] in transition_species:
                    colour.append(species_colour[transition.split()[0]])
                    transition_species.append(transition.split()[0])
                ax.loglog(
                    [1e6 * species.clump_transition_wavelengths[i]] * 2,
                    [valTemp, valueSpecies[ens][vel, i]],
                    ls="-",
                    lw=1,
                    c=colour[-1],
                )

            lines = [
                Line2D([0], [0], color=species_colour[transition], lw=1)
                for transition in transition_species
            ]
            leg = ax.legend(
                labels=transition_species,
                handles=lines,
                fontsize=16,
                bbox_to_anchor=(1, 1),
                loc="upper left",
            )

            ax.set_xlabel(r"$\lambda \ (\mu \mathrm{m})$", fontsize=24)
            ax.set_ylabel(ylabel, fontsize=24)
            ax.set_ylim(
                [
                    valueDust[ens][vel, :].min() * 10e-3,
                    valueDust[ens][vel, :].max() * 10e3,
                ]
            )
            ax.tick_params(labelsize=16)
            # plt.setp(ax.get_xticklabels(), fontsize=14)
            # plt.setp(ax.get_yticklabels(), fontsize=14)

            if title:
                ax.set_title(title, fontsize=24)
            else:
                ax.set_title(
                    "Clump Set {} {} spectrum".format(ens + 1, quantity), fontsize=24
                )

        # fig.tight_layout()

        # plt.show()

        return axes

    def printVoxel(self):
        """
        Coming soon...
        """
        return
