import numpy as np

from .. import constants
from .methods import *
'''
This module will contain the input data needed to properly simulate the PDR. All of the information specific
to the object being simulated should be in their own folder in INPUTPATH. The KOSMA-tau grid data is located
in the folder 'grid'.
'''

# def __str__(self):
#   if not len(speciesData):
#     return 'There are no available transitions yet.'
#   else:
#     printout = 'Available transitions:'
#     i = np.isfinite(speciesData[0])&np.isfinite(speciesData[1])&np.isfinite(speciesData[2])&np.isfinite(speciesData[3])
#     transitions = speciesData[2]
#     elements = speciesData[1]
#     for i in range(len(elements)):
#       printout += '\n  ->{} {}'.format(element[i], transitions[i])
#     return printout

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
densityProfile =  None
FUVfield = None
rotationProfile = None
