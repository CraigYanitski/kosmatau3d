# import numpy as np
# from astropy.io import fits
# import matplotlib.pyplot as plt
# from mpl_toolkits.mplot3d import Axes3D

from .orientation import *

from .. import constants

import scipy.interpolate as interpolate
from scipy.stats import norm

positions = []
intensitySpecies = []
intensityDust = []
losVoxels = []
voxelVelocities = []
voxelPositions = []
i_vox = []
sightlines = []
x1LoS = 0
x2LoS = 0
x3LoS = []

# tempClumpVelocity = []
tempSpeciesEmissivity = []
tempSpeciesAbsorption = []
tempDustEmissivity = []
tempDustAbsorption = []
tempPosition = []
# tempInterclumpVelocity = []
# tempInterclumpEmission = []
# tempInterclumpPosition = []

# This is used for the integration of each individual sightline
epsilonSpecies = []
epsilonStepSpecies = 0
kappaSpecies = []
kappaStepSpecies = 0
epsilonDust = []
epsilonStepDust = 0
kappaDust = []
kappaStepDust = 0

# These can be used for debugging the code. It will store the indeces, kappa, and epsilon values for 25 sightlines.
allIndeces = []

allKSpecies = []
allESpecies = []
allKstepSpecies = []
allEstepSpecies = []
allArealSpecies = []
allBrealSpecies = []
allAimagSpecies = []
allBimagSpecies = []
allIntensitySpecies = []

allKstepDust = []
allEstepDust = []
allArealDust = []
allBrealDust = []
allAimagDust = []
allBimagDust = []
allIntensityDust = []

# Here are the values for interpolating the real and imaginary emission
eTildeRealObs = orientation.eTildeReal()
eTildeImaginaryObs = orientation.eTildeImaginary()
eTildeReal = interpolate.interp1d(eTildeRealObs[0], eTildeRealObs[1], kind='linear')
eTildeImaginary = interpolate.interp1d(eTildeImaginaryObs[0], eTildeImaginaryObs[1], kind='linear')