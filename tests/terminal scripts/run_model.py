#python3.7

import numpy as np
import sys
from astropy.io import fits

from kosmatau3d import models

timed = False
debug = False

# print('KOSMA-tau^3')

# Edit these parameters according to the model you want to produce.
parameters = {
                # Model information
                'history_path': r'\users\cyani\projects\pdr\KT3_history',
                'directory': r'\MilkyWay',
                'x': 36,
                'y': 36,
                'z': 0.5,
                'modelType': 'disk',
                
                # Model parameters
                'resolution': 500,
                # 'molecules': 'all',
                'molecules' : ['C+ 1', 'C 1', 'C 2', 'C 3', 'CO 1', 'CO 2', 'CO 3', 'CO 4', 'CO 5', '13C+ 1', '13C 1', '13C 2', '13CO 1', '13CO 2', '13CO 3', '13CO 4', '13CO 5', 'HCO+ 1'],
                # 'dust' : 'PAH',
                'dust': '3.1mm',
                'clumpMassRange': [[0, 2], [-2]],
                'clumpMassNumber': [3, 1],
                'clumpNmax': [1, 100],
                'velocityRange': [-360, 360],
                'velocityNumber': 101,
                
                # Property factors
                'clumpMassFactor': [1, 1],
                'FUVfactor': 1,
                'densityFactor': 1,
                'globalUV': 10
              }

if __name__=='__main__':
  kosma = models.Model(**parameters)
  
  kosma.calculateModel(timed=timed, debug=debug)

directory = parameters['history_path'] + parameters['directory'] + r'\r{}_cm{}-{}_d{}_uv{}'.format(parameters['resolution'], parameters['clumpMassFactor'][0], parameters['clumpMassFactor'][1], parameters['densityFactor'], parameters['globalUV'])

# Edit the directory in this condition to allow multiprocessing in Windows. Linux does not need this hard encoding.
if ('win' in sys.platform) and (__name__!='__main__'):
  
  models.radiativeTransfer.voxelPositions = fits.open(directory+'/voxel_position.fits', mode='denywrite')
  models.radiativeTransfer.voxelVelocities = fits.open(directory+'/voxel_velocity.fits', mode='denywrite')
  models.radiativeTransfer.tempSpeciesEmissivity = fits.open(directory+'/species_emissivity.fits', mode='denywrite')# in K/pc
  models.radiativeTransfer.tempSpeciesAbsorption = fits.open(directory+'/species_absorption.fits', mode='denywrite')# in 1/pc
  models.radiativeTransfer.tempDustEmissivity = fits.open(directory+'/dust_emissivity.fits', mode='denywrite')# in K/pc
  models.radiativeTransfer.tempDustAbsorption = fits.open(directory+'/dust_absorption.fits', mode='denywrite')# in 1/pc
  
if __name__=='__main__':
  models.radiativeTransfer.calculateObservation(directory=directory, dim='spherical', slRange=[(-np.pi,np.pi), (-np.pi/8,np.pi/8)], nsl=[360,45], terminal=True, debug=False, multiprocessing=4)
