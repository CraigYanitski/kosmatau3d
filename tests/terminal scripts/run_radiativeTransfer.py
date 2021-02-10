from kosmatau3d import models

import numpy as np
from astropy.io import fits

# Edit this directory to point to the model you want to use.
# directory = r'C:\Users\cyani\projects\pdr\KT3_history\MilkyWay\r200_cm1-1_d1_uv10'

directory = r'C:\Users\cyani\projects\pdr\KT3_history\MilkyWay\r1000_cm1-1_d1_uv10'

models.constants.velocityRange = np.linspace(-300, 300, 500)

models.radiativeTransfer.voxelPositions = fits.open(directory+'/voxel_position.fits', mode='denywrite')
models.radiativeTransfer.voxelVelocities = fits.open(directory+'/voxel_velocity.fits', mode='denywrite')
models.radiativeTransfer.tempSpeciesEmissivity = fits.open(directory+'/species_emissivity.fits', mode='denywrite')# in K/pc
models.radiativeTransfer.tempSpeciesAbsorption = fits.open(directory+'/species_absorption.fits', mode='denywrite')# in 1/pc
models.radiativeTransfer.tempDustEmissivity = fits.open(directory+'/dust_emissivity.fits', mode='denywrite')# in K/pc
models.radiativeTransfer.tempDustAbsorption = fits.open(directory+'/dust_absorption.fits', mode='denywrite')# in 1/pc

# Calculate integrated intensity maps (adjust sl to change the number of sightlines you want in the map; [longitude, latitude]).
if __name__=='__main__':
  
  print('Main script')
  
  models.radiativeTransfer.calculateObservation(directory=directory, dim='spherical', nsl=[50,25], terminal=True, debug=False, multiprocessing=4)