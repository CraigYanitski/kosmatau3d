# import numpy as np
# from astropy.io import fits
# import matplotlib.pyplot as plt
# from mpl_toolkits.mplot3d import Axes3D
import scipy.interpolate as interpolate
from scipy.stats import norm
from logging import getLogger

from .orientation import *
from kosmatau3d.models import constants


logger = getLogger(__name__)

pencil_beam = True
d_lon = 1
d_lat = 1

positions = []

intensity_species = []
intensity_dust = []

opticaldepth_species = []
opticaldepth_dust = []

los_voxels = []
voxel_velocities = []
voxel_positions = []
i_vox = []

sightlines = []

x1LoS = 0
x2LoS = 0
x3LoS = []

# tempClumpVelocity = []
temp_species_emissivity = []
temp_species_absorption = []

temp_dust_emissivity = []
temp_dust_absorption = []

temp_position = []

# tempInterclumpVelocity = []
# tempInterclumpEmission = []
# tempInterclumpPosition = []

# This is used for the integration of each individual sightline
i0_species = 0
e_species = []
de_species = 0
k_species = []
dk_species = 0

i0_dust = 0
e_dust = []
de_dust = 0
k_dust = []
dk_dust = 0

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
#eTildeRealObs = rt.e_tilde_real()
#eTildeImaginaryObs = rt.e_tilde_imaginary()
#eTildeReal = interpolate.interp1d(eTildeRealObs[0], eTildeRealObs[1], kind='linear')
#eTildeImaginary = interpolate.interp1d(eTildeImaginaryObs[0], eTildeImaginaryObs[1], kind='linear')
