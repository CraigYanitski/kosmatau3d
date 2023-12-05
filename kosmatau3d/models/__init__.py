"""
The :code:`models` subpackage contains all of the code relevant for creating
voxels.
There are currently two ways to model PDRs:

- single-voxel model (:code:`class Voxel()`): this is useful to a homogeneous
  voxel to compare to pixels in an observation.
  Note that the voxel does not necessarily need to be cubic; PDR clumps can 
  *overfill* a voxel resulting in a homogeneous *column*. 
  This functionality will require more attention if you plan to integrate 
  the radiative transfer equation.
- full 3-D model (:code:`class Model()`): The full three-dimensional model used 
  in the first applications of this code (Andree-Labsch et al. 2017, Yanitski 
  et al. 2023). 
  This requires knowledge of the three-dimensional structure and distribution 
  of interstellar gas in the object being modelled. Currently this is setup to 
  model a Milky-Way-type galaxy.

Most of the calculations are split into subsubpackages, although there is an
object containing all of the evaluated voxels.
This structure to the code was implemented with the intention of creating a
model in memory to adjust the orientation of the model and thus the synthetic
observation of it.
While the Milky Way models are too large to utilise this feature, it is left
in the source code for future  development of the code.
"""

import dill
import os
import sys

from copy import copy

from .model import *
from .voxel import *

from kosmatau3d.models import constants
from kosmatau3d.models import species
from kosmatau3d.models import shape
from kosmatau3d.models import observations
from kosmatau3d.models import interpolations
from kosmatau3d.models import ensemble
from kosmatau3d.models import combinations
from kosmatau3d.models import masspoints
from kosmatau3d.models import radiativeTransfer
from kosmatau3d.models import plotting


def dill_grid():
    """
    A function that can be used to serialise and save the interpolation 
    functions to make debugging easier.
    Use this function to load the grid, initialise all of the interpolation
    functions, and save them in the package folder.
    This helps to reduce computation time when first running a voxel.
    After the dilled files are created, load the grid from them by passing
    :code:`dill=True` as a kwarg when running :code:`Voxel.set_properties()`.
    """
    model_directory = copy(constants.directory)
    constants.directory = ""
    observations.methods.initialise()
    constants.changeDustWavelengths("all")
    species.addMolecules("all")
    interpolations.initialise_grid()
    if not os.path.exists(constants.GRIDPATH + "dilled/"):
        os.mkdir(constants.GRIDPATH + "dilled/")
    with open(constants.GRIDPATH + "dilled/intensity_interpolation", "wb") as file:
        dill.dump(interpolations.intensityInterpolation, file)
    with open(constants.GRIDPATH + "dilled/tau_interpolation", "wb") as file:
        dill.dump(interpolations.tauInterpolation, file)
    with open(constants.GRIDPATH + "dilled/dust_intensity_interpolation", "wb") as file:
        dill.dump(interpolations.dustIntensityInterpolation, file)
    with open(constants.GRIDPATH + "dilled/dust_tau_interpolation", "wb") as file:
        dill.dump(interpolations.dustTauInterpolation, file)
    with open(constants.GRIDPATH + "dilled/taufuv_interpolation", "wb") as file:
        dill.dump(interpolations.taufuv_interpolation, file)
    with open(
        constants.GRIDPATH + "dilled/hi_column_density_interpolation", "wb"
    ) as file:
        dill.dump(interpolations.hi_column_density_interpolation, file)
    with open(
        constants.GRIDPATH + "dilled/h2_column_density_interpolation", "wb"
    ) as file:
        dill.dump(interpolations.h2_column_density_interpolation, file)
    with open(constants.GRIDPATH + "dilled/tg_interpolation", "wb") as file:
        dill.dump(interpolations.tg_interpolation, file)
    constants.directory = model_directory
    return


def help():
    """A currently defunct function to show information about this submodule."""
    print()
    return
