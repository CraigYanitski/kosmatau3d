from numba import jit_module
import numpy as np
import inspect
import os

from numba.errors import NumbaDeprecationWarning, NumbaPendingDeprecationWarning
import warnings
warnings.simplefilter('ignore', category=NumbaDeprecationWarning)
warnings.simplefilter('ignore', category=NumbaPendingDeprecationWarning)

from .change import *
'''
This is a module to contain all of the constants and parameters used throughout the program. There are definitions
to change the model parameters for when this is needed.
'''

# Directory information
directory = ''
history = ''
filename = inspect.getframeinfo(inspect.currentframe()).filename
KOSMAPATH = os.path.abspath(os.path.dirname(filename)+'/../../')
INPUTPATH = KOSMAPATH + '/input/'#'/home/craig/Desktop/Köln/kosma-tau^3-develop/kosma-tau-3d/input/'#
GRIDPATH = KOSMAPATH + '/grid/'#'/home/craig/Desktop/Köln/kosma-tau^3-develop/kosma-tau-3d/grid/'#
HISTORYPATH = KOSMAPATH + '/history/'#'/home/craig/Desktop/Köln/kosma-tau^3-develop/kosma-tau-3d/history/'#

# Interpolation style (it accepts 'linear' or 'cubic'/'radial')
interpolation = 'linear'

# Factors (and constant) for fine-tuning the input data
clumpMassFactor = 1
interclumpMassFactor = 1
densityFactor = 1
FUVFactor = 1

interclumpLogFUV = 1

# Standard constants
massH = 1.007276466812*1.6605*10**-27  #in [kg]
massSolar = 1.98852*10**30  #in [kg]
pc = 3.08567758149*10**16  #in [m]

h = 6.62606957*10**-34 #in [J s]
c = 2.99792458*10**8  #in [m/s]
kB = 1.3806488*10**-23  #in [J/K]

# For observing the Milky Way from Earth
fromEarth = True
rGalEarth = 8178 #from Abuter and the GRAVITY collaboration (2019)
rGal = 18000
hd = 1000

# Integrated map properties
mapShape = (150, 100)   #shape in pixels of the desired map
mapSize = (360, 180)   #size in degrees of the desired map
mapCenter = (0,0)       #center of the desired map in degrees

# Model species
molecules = ['C+', 'C', 'CO', '13CO', '13C+', '13C', 'HCO+', 'H13CO+', 'H13CO', 'H3O+', 'C18O']
dust = True

# Model characteristics
nSigma = 3
resolution = 1000
velocityNumber = 121
velocityBin = [-300, 300]
velocityRange = np.linspace(velocityBin[0], velocityBin[-1], num=velocityNumber)
velocityStep = velocityRange[1] - velocityRange[0]

# Clump characteristics
clumpDispersion = 1.67/2.3548

clumpMassNumber=4
clumpLogMassRange = [-1,2]
clumpLogMass = np.linspace(clumpLogMassRange[0], clumpLogMassRange[-1], num=clumpMassNumber)
clumpLogMass.resize(1,clumpMassNumber)
clumpMaxIndeces = 0
clumpNmax = 1

interclumpDensity = 1911  #as defined in the old version of this code
interclumpMassNumber = 2
interclumpLogMassRange = [-3,-2]
interclumpLogMass = np.linspace(interclumpLogMassRange[0], interclumpLogMassRange[-1], num=interclumpMassNumber)
interclumpLogMass.resize(1,interclumpMassNumber)
interclumpMaxIndeces = 0
interclumpNmax = 100

ensembleDispersion = 0#10./2.3548

# Initial mass function parameters
alpha = 1.84
gamma = 2.31

# Statistics
probability = 'binomial'
nGauss = 1000
pnGauss = 5

# Grid boundaries
densityLimits = [3, 7]
massLimits = [-3, 3]
uvLimits = [0, 6]

# UV adjustment
normUV = 2.89433*10**39
globalUV = 10

# Finally, these are the wavelengths at which the dust emission will be computed
wavelengths = np.array([3100., 2400., 1800., 1300., 1000., 850., 700., 550., 420., 300., 240., 188.4, 177.8, 167.9, 158.5, 149.6, \
                       141.3, 133.4, 125.9, 118.9, 112.2, 105.9, 100., 94.41, 89.13, 84.14, 79.43, 74.99, 70.79, 66.83, 63.1, 59.57, \
                       56.23, 53.09, 50.12, 47.32, 44.67, 44.28, 42.75, 41.33, 39.99, 38.74, 37.57, 36.47, 35.42, 34.44, 33.51, \
                       32.63, 31.79, 31., 30.24, 29.52, 28.83, 28.18, 27.55, 26.95, 26.38, 25.83, 25.3, 24.8, 24.31, 23.84, 23.39, \
                       22.96, 22.54, 22.14, 21.75, 21.38, 21.01, 20.66, 20.32, 20., 19.68, 19.37, 19.07, 18.79, 18.5, 18.23, 17.97, \
                       17.71, 17.46, 17.22, 16.98, 16.75, 16.53, 16.31, 16.1, 15.89, 15.69, 15.5, 15.31, 15.12, 14.94, 14.42, 13.93, \
                       13.78, 13.62, 13.48, 13.33, 13.19, 13.05, 12.91, 12.78, 12.65, 12.52, 12.4, 12.28, 12.15, 12.04, 11.92, 11.81, \
                       11.7, 11.59, 11.48, 11.37, 11.27, 11.17, 11.07, 10.97, 10.88, 10.78, 10.69, 10.6, 10.51, 10.42, 10.33, 10.16, \
                       10.08, 9.998, 9.918, 9.84, 9.762, 9.686, 9.611, 9.537, 9.464, 9.392, 9.322, 9.252, 9.184, 9.116, 9.05, 8.984, \
                       8.919, 8.856, 8.793, 8.731, 8.67, 8.61, 8.55, 8.492, 8.434, 8.377, 8.321, 8.265, 8.211, 8.157, 8.103, 8.051, \
                       7.999, 7.947, 7.897, 7.847, 7.797, 7.749, 7.701, 7.653, 7.606, 7.56, 7.514, 7.424, 7.38, 7.336, 7.293, 7.25, \
                       7.208, 7.166, 7.125, 7.085, 7.044, 7.005, 6.965, 6.926, 6.888, 6.85, 6.812, 6.775, 6.738, 6.702, 6.595, 6.491, \
                       6.391, 6.293, 6.199, 6.107, 5.989, 5.904, 5.793, 5.61, 5.39, 5.209, 4.999, 4.805, 4.592, 4.396, 4.203, 3.999, \
                       3.936, 3.874, 3.815, 3.757, 3.701, 3.647, 3.594, 3.542, 3.492, 3.444, 3.397, 3.342, 3.324, 3.306, 3.297, 3.29, \
                       3.28, 3.246, 3.204, 3.099, 3.002, 2.8, 2.661, 2.512, 2.371, 2.239, 2.113, 1.995, 1.884, 1.778, 1.679, 1.585, \
                       1.496, 1.413, 1.334, 1.259, 1.189, 1.122, 1.059, 1., 0.9441, 0.8913, 0.8414, 0.7943, 0.7499, 0.7079, 0.6683, \
                       0.631, 0.5957, 0.5623, 0.5309, 0.5012, 0.4732, 0.4467, 0.4217, 0.3981, 0.3758, 0.3548, 0.335, 0.3162, 0.2985, \
                       0.2818, 0.2661, 0.2512, 0.2371, 0.2239, 0.2113, 0.1995, 0.1884, 0.1778, 0.1679, 0.1585, 0.1496, 0.1413, 0.1334, \
                       0.1259, 0.1216, 0.1189, 0.1122, 0.1059, 0.1, 0.09441, 0.08913, 0.08414, 0.07943, 0.07499, 0.07079, 0.06683, \
                       0.0631, 0.05957, 0.0584, 0.05623, 0.05309, 0.05012, 0.04732, 0.04467, 0.04217, 0.03981, 0.03758, 0.03548, \
                       0.0335, 0.03162, 0.0304, 0.02985, 0.02818, 0.02661, 0.02512, 0.02371, 0.02239, 0.01995, 0.01778, 0.01585, \
                       0.01413, 0.01259, 0.01122, 0.0106, 0.01, 0.009, 0.008, 0.007, 0.006, 0.005, 0.004, 0.003, 0.002, 0.001]) * 10**-6
# These are the flag and indeces to limit how many dust emission lines are calculated
dustWavelengths = ''
limitMolecular = 10**-4
limitPAH = 5*10**-6
nDust = wavelengths>0
# [DEACTIVATED] These are variables to sort the wavelengths
sortedWavelengths = []
sortedIndeces = wavelengths.argsort()

hclambda = h*c/wavelengths #to aid in calculations

jit_module(nopython=False)
