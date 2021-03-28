import numpy as np

from kosmatau3d.models import constants
from .methods import *


'''
This module will contain the input data needed to properly simulate the PDR. All of the information specific
to the object being simulated should be in their own folder in INPUTPATH. The KOSMA-tau grid data is located
in the folder 'grid'.
'''


# Grid
tbCenterline = None
tauCenterline = None
rhoMassAFUV = None
speciesData = None
eTildeReal = None
eTildeImaginary = None

# Model
clumpMassProfile = None
interclumpMassProfile = None
densityProfile = None
FUVfield = None
rotationProfile = None
