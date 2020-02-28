import numpy as np
from astropy.io import fits
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D

from radiativeTransfer import orientation

import constants

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

epsilon = []
epsilonStep = 0
kappa = []
kappaStep = 0

eTildeRealObs = orientation.eTildeReal()
eTildeImaginaryObs = orientation.eTildeImaginary()
eTildeReal = interpolate.interp1d(eTildeRealObs[0], eTildeRealObs[1], kind='linear')
eTildeImaginary = interpolate.interp1d(eTildeImaginaryObs[0], eTildeImaginaryObs[1], kind='linear')

def plotModel(plot='total intensity', ce=[], ie=[], directory='/home/craig/projects/pdr/KOSMA-tau^3/history/MilkyWay/resolution500_size36000/', species='13C+ 1', debug=False):
  allSpecies = ['13C 1', '13C 2', '13C 3', '13C+ 1', 'C 1', 'C 2', 'C+ 1', 'CO 1', 'CO 2', 'CO 3', 'CO 4', 'CO 5', 'CO 6', 'CO 7', '13CO 1', '13CO 2', '13CO 3', '13CO 4', '13CO 5', '13CO 6', '13CO 7', 'O 2', 'Dust 1', 'Dust 2', 'Dust 3', 'Dust 4', 'Dust 5', 'Dust 6', 'Dust 7', 'Dust 8', 'Dust 9', 'Dust 10', 'Dust 11', 'Dust 12', 'Dust 13', 'Dust 14', 'Dust 15', 'Dust 16', 'Dust 17', 'Dust 18', 'Dust 19', 'Dust 20', 'Dust 21', 'Dust 22', 'Dust 23', 'Dust 24', 'Dust 25', 'Dust 26', 'Dust 27', 'Dust 28', 'Dust 29', 'Dust 30', 'Dust 31', 'Dust 32', 'Dust 33', 'Dust 34', 'Dust 35', 'Dust 36', 'Dust 37', 'Dust 38', 'Dust 39', 'Dust 40']

  directory = constants.HISTORYPATH + constants.directory

  voxelPositions = fits.open(directory+'voxel_position.fits')[0].data
  fuv = fits.open(directory+'voxel_fuv.fits')[0].data
  Afuv = fits.open(directory+'voxel_Afuv.fits')[0].data
  velocity = fits.open(directory+'voxel_velocity.fits')[0].data
  if len(ce)==0 and len(ie)==0:
    clumpEmission = fits.open(directory+'emission_clump.fits')[0].data
    interclumpEmission = fits.open(directory+'emission_interclump.fits')[0].data

    maxClumpEmission = np.zeros(clumpEmission[:,:,:,0].shape)
    maxInterclumpEmission = np.zeros(interclumpEmission[:,:,:,0].shape)
    for emission in range(2):
      for voxel in range(clumpEmission[:,0,0,0].size):
        for species in range(clumpEmission[0,0,:,0].size):
          maxClumpEmission[voxel,emission,species] = norm.fit(clumpEmission[voxel,emission,species,:])[0]
          maxInterclumpEmission[voxel,emission,species] = norm.fit(interclumpEmission[voxel,emission,species,:])[0]
  else:
    maxClumpEmission = ce
    maxInterclumpEmission = ie

  #positions = self.__grid.getVoxelPositions()
  limits = [voxelPositions.min(), voxelPositions.max()]
  if debug: print(self.__grid.totalEmission().shape)
  if plot=='total intensity':
    weights = (maxClumpEmission[:,0,:]+maxInterclumpEmission[:,0,:]).max(1)
    scale = r'$I \ (\chi)$'
    plotTitle = 'PDR Intensity within the Milky Way'
  elif plot=='total optical depth':
    weights = (maxClumpEmission[:,1,:]+maxInterclumpEmission[:,1,:]).max(1)
    scale = r'$\tau$'
    plotTitle = 'PDR Optical Depth within the Milky Way'
  elif plot=='species total intensity':
    i = allSpecies==species
    weights = maxClumpEmission[:,0,i]+maxInterclumpEmission[:,0,i]
    weights = weights/max(weights)
    scale = r'$I \ (\chi)$'
    plotTitle = species + ' intensity within the Milky Way'
  elif plot=='species total optical depth':
    i = allSpecies==species
    weights = maxClumpEmission[:,1,i]+maxInterclumpEmission[:,1,i]
    weights = weights/weights.max()
    scale = r'$\tau$'
    plotTitle = species + ' optical depth within the Milky Way'
  elif plot=='clump intensity':
    weights = (maxClumpEmission[:,0,:]).sum(1)
    scale = r'$I \ (\chi)$'
    plotTitle = 'Clump intensity within the Milky Way'
  elif plot=='clump optical depth':
    weights = (maxClumpEmission[:,1,:]).sum(1)
    scale = r'$\tau$'
    plotTitle = 'Clump optical depth within the Milky Way'
  elif plot=='interclump intensity':
    weights = (maxInterclumpEmission[:,0,:]).sum(1)
    scale = r'$I \ (\chi)$'
    plotTitle = 'Interclump intensity within the Milky Way'
  elif plot=='interclump optical depth':
    weights = (maxInterclumpEmission[:,1,:]).sum(1)
    scale = r'$\tau$'
    plotTitle = 'Interclump optical depth within the Milky Way'
  elif plot=='FUV':
    weights = (fuv)
    scale = r'$FUV \ (\chi)$'
    plotTitle = 'FUV within the Milky Way'
  elif plot=='Afuv':
    weights = (Afuv)
    scale = r'$\tau_{FUV}$'
    plotTitle = r'$A_{FUV}$ within the Milky Way'
  elif plot=='velocity':
    weights = (velocity)
    scale = r'$v_{rot} \ (\frac{km}{s})$'
    plotTitle = 'Rotational velocity within the Milky Way'
  fig = plt.figure()
  ax = fig.add_subplot(111, projection='3d')
  model = ax.scatter(voxelPositions[:,0], voxelPositions[:,1], voxelPositions[:,2], c=weights, cmap=plt.cm.hot, marker='s', s=27, alpha=0.5, linewidths=0)
  ax.set_xlim(limits)
  ax.set_ylim(limits)
  ax.set_zlim(limits)
  cbar = plt.colorbar(model)
  ax.set_title(plotTitle)
  ax.set_xlabel('X (pc)')
  ax.set_ylabel('Y (pc)')
  ax.set_zlabel('Z (pc)')
  cbar.set_label(scale, rotation=0)
  plt.show()
  return maxClumpEmission,maxInterclumpEmission