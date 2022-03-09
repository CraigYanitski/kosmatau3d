import dill
import sys

from copy import copy

from .Model import *
from .Voxel import *

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
    model_directory = copy(constants.directory)
    constants.directory = ''
    observations.methods.initialise()
    constants.changeDustWavelengths('all')
    species.addMolecules('all')
    interpolations.initialise_grid()
    if not os.path.exists(constants.GRIDPATH + 'dilled/'):
        os.mkdir(constants.GRIDPATH + 'dilled/')
    with open(constants.GRIDPATH + 'dilled/intensity_interpolation', 'wb') as file:
        dill.dump(interpolations.intensityInterpolation, file)
    with open(constants.GRIDPATH + 'dilled/tau_interpolation', 'wb') as file:
        dill.dump(interpolations.tauInterpolation, file)
    with open(constants.GRIDPATH + 'dilled/dust_intensity_interpolation', 'wb') as file:
        dill.dump(interpolations.dustIntensityInterpolation, file)
    with open(constants.GRIDPATH + 'dilled/dust_tau_interpolation', 'wb') as file:
        dill.dump(interpolations.dustTauInterpolation, file)
    with open(constants.GRIDPATH + 'dilled/fuv_extinction_interpolation', 'wb') as file:
        dill.dump(interpolations.FUVextinctionInterpolation, file)
    constants.directory = model_directory
    return

def help():

    print()
    return