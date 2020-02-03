import os
import inspect
import numpy as np
import constants
from methods import *
'''
This module will contain the input data needed to properly simulate the PDR. All of the information specific
to the object being simulated should be in their own folder in INPUTPATH. The KOSMA-tau grid data is located
in the folder 'grid'.
'''
directory = 'MilkyWay'
if directory[-1]!='/': directory = directory + '/'

scale = constants.scale
filename = inspect.getframeinfo(inspect.currentframe()).filename
KOSMAPATH = os.path.abspath(os.path.dirname(filename)+'/../')
INPUTPATH = KOSMAPATH + '/input/'#'/home/craig/Desktop/Köln/kosma-tau^3-develop/kosma-tau-3d/input/'#
GRIDPATH = KOSMAPATH + '/grid/'#'/home/craig/Desktop/Köln/kosma-tau^3-develop/kosma-tau-3d/grid/'#


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

# PUBLIC
clumpMassProfile = clumpMassProfile()
interclumpMassProfile = interclumpMassProfile()
densityProfile =  densityProfile()
FUVfield = FUVfield()
rotationProfile = rotationProfile()
tauCenterline = tauCenterline()
tbCenterline = tbCenterline()
rhoMassAFUV = rhoMassAFUV()
speciesData = speciesData()
eTildeReal = eTildeReal()
eTildeImaginary = eTildeImaginary()
