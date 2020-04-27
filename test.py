#!/use/bin python3

import os
import sys
import inspect
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D

pathname = inspect.getframeinfo(inspect.currentframe()).filename
sys.path.append(os.path.abspath(os.path.dirname(pathname))+'/classes/')

from Model import *
import constants

complete = False

modelFlag = False
rtFlag = True
cyplotFlag = False

x = 36
y = 36
z = 2

shape = 'disk'

resolution = 1000

focus = 'molecular'

constants.changeDirectory('MilkyWay')

modelFolder = 'r1000_n3015/'

# Factors
constants.clumpMassFactor = 1
constants.interclumpMassFactor = 1
constants.FUVFactor = 1
constants.DensityFactor = 1

# Constant
constants.interclumpLogFUV = 1

# Model masses
constants.clumpLogMassNum = 4
constants.clumpLogMassRange = [-1, 2]
constants.interclumpLogMassNum = 2
constants.interclumpLogMassRange = [-3, -2]

# Limit dust calculation

print('KOSMA-tau^3')

species = ['13CO 10', 'C+ 1', 'CO 1', 'CO 2', 'CO 3', 'CO 4', 'CO 5', 'CO 6', 'CO 7', 'CO 8', 'CO 9', 'CO 10', '13CO 1', '13CO 2', '13CO 3', '13CO 4', '13CO 5', '13CO 6', '13CO 7', '13CO 8', '13CO 9', '13CO 10', 'O 2']
kosma = Model(x, y, z, modelType=shape, resolution=resolution, dustLim=focus)
kosma.addSpecies(species)

if modelFlag:

  kosma.calculateModel()
  kosma.writeEmission()
  
  modelFolder = constants.history

# Calculate integrated intensity maps

if rtFlag:

  import radiativeTransfer

  # radiativeTransfer.plotModel(plot='velocity', directory='r1000.0_n3015/')

  radiativeTransfer.calculateObservation(directory=modelFolder, dim='spherical', plot=True)

# Create cygrid images

if cyplotFlag:

  import cyplot

  images,wcs = cyplot.convertMap(modelFolder, input_coord='spherical')
  print(images.shape)

  np.nan_to_num(images)

  for i in range(images[:,0,0,0].size):
    plt.imshow(images[i,333,:,:])
    plt.show()
    