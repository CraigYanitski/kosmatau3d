import numpy as np
'''
This is a module to contain all of the constants and parameters used throughout the program. There are definitions
to change the model parameters for when this is needed.
'''

# Standard constants
massH = 1.008*1.6605*10**-24  #in [g]
massSolar = 1.98892*10**33  #in [g]
c = 2.998*10**10  #in [cm/s]
kB = 1.3806*10**-16  #in [erg/K]
pc = 3.0856776*10**18  #in [cm]

# Model characteristics
nsigma = 3
pixelWidth = 1000
velocityNumber = 101
velocityBin = [-360, 360]
velocityRange = np.linspace(velocityBin[0], velocityBin[1], num=velocityNumber)

# Clump characteristics
clumpMassNumber=2
clumpMassRange = [1,2]
clumpMass = np.linspace(clumpMassRange[0], clumpMassRange[1], num=clumpMassNumber)
clumpDispersion = 1.67/2.3548

interclumpMassNumber = 1
interclumpMassRange = [-2,-1]
interclumpMass = np.linspace(interclumpMassRange[0], interclumpMassRange[1], num=interclumpMassNumber)
ensembleDispersion = 10./2.3548

# Initial mass function parameters
alpha = 1.84
gamma = 2.31

# Statistics
nGauss = 1000
pnGauss = 5

# Grid boundaries
densityLimits = [3, 7]
massLimits = [-3, 3]
uvLimits = [0, 6]

# UV adjustment
normUV = 2.89433*10**39
globalUV = 10
