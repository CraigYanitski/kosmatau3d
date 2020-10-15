import numpy as np
from astropy.io import fits
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D

from .orientation import *

from .. import constants

import scipy.interpolate as interpolate
from scipy.stats import norm

positions = []
intensity = []
losVoxels = []
x1LoS = 0
x2LoS = 0
x3LoS = []

tempClumpVelocity = []
tempClumpEmission = []
tempClumpPosition = []
tempInterclumpVelocity = []
tempInterclumpEmission = []
tempInterclumpPosition = []

# This is used for the integration of each individual sightline
epsilon = []
epsilonStep = 0
kappa = []
kappaStep = 0

# These can be used for debugging the code. It will store the indeces, kappa, and epsilon values for 25 sightlines.
allIndeces = []
allK = []
allE = []
allKstep = []
allEstep = []
allAreal = []
allBreal = []
allAimag = []
allBimag = []
allIntensity = []

# Here are the values for interpolating the real and imaginary emission
eTildeRealObs = orientation.eTildeReal()
eTildeImaginaryObs = orientation.eTildeImaginary()
eTildeReal = interpolate.interp1d(eTildeRealObs[0], eTildeRealObs[1], kind='linear')
eTildeImaginary = interpolate.interp1d(eTildeImaginaryObs[0], eTildeImaginaryObs[1], kind='linear')