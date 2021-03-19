# python3.7

import numpy as np

from kosmatau3d import models

timed = False
debug = False

print('KOSMA-tau^3')

# Edit these parameters according to the model you want to produce.
parameters = {
              # Model information
              'history_path': r'/media/yanitski/4,0 TB Hard Disk/yanitski/projects/pdr/KT3_history',
              'directory': r'/MilkyWay',
              'x': 36,
              'y': 36,
              'z': 0.5,
              'modelType': 'disk',

              # Model parameters
              'resolution': 200,
              # 'molecules': 'all',
              'molecules': ['C+ 1', 'C 1', 'C 2', 'CO 1', 'CO 2', 'CO 3', 'CO 4', 'CO 5', '13C+ 1', '13C 1', '13C 2', '13CO 1',
                            '13CO 2', '13CO 3', '13CO 4', '13CO 5', 'HCO+ 1'],
              # 'dust': 'PAH',
              'dust': '3.1mm',
              'clumpMassRange': [[0, 2], [-2]],
              'clumpMassNumber': [3, 1],
              'clumpNmax': [1, 100],
              'velocityRange': [-300, 300],
              'velocityNumber': 601,

              # Property factors
              'clumpMassFactor': [1, 1],
              'FUVfactor': 1,
              'densityFactor': 1,
              'globalUV': 10
              }

for f_mass1 in [1]:#[0.5,1.0,2.0,4.0]:
  for f_mass2 in [1]:#[1.0,2.0]:
    for f_density in [1]:#[0.5,1.0,2.0, 4.0]:
      for f_uv in [100]:#[10, 100]:

          parameters['clumpMassFactor'] = [f_mass1, f_mass2]
          parameters['densityFactor'] = f_density
          parameters['globalUV'] = f_uv

          kosma = models.Model(**parameters)
          kosma.calculateModel(timed=timed, debug=debug, multiprocessing=0)

          directory = parameters['history_path'] + parameters['directory'] + '/r{}_cm{}-{}_d{}_uv{}'.format(parameters['resolution'], f_mass1, f_mass2, f_density, f_uv)
          models.radiativeTransfer.calculateObservation(directory=directory, dim='spherical', slRange=[(-np.pi, np.pi), (-np.pi/4, np.pi/4)],
                                                        nsl=[360, 15], terminal=True, debug=False, multiprocessing=6)
