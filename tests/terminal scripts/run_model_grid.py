# python3.7

import numpy as np

from kosmatau3d import models

timed = False
debug = False

print('KOSMA-tau^3')

# Edit these parameters according to the model you want to produce.
parameters = {
              # Model information
              'history_path': r'/mnt/hpc_backup/yanitski/projects/pdr/KT3_history',
              'directory': r'/MilkyWay',
              'folder': '',
              'x': 36,
              'y': 36,
              'z': 0,
              'modelType': 'disk',

              # Model parameters
              'resolution': 400,
              # 'molecules': 'all',
              'molecules': ['C+ 1', 'C 1', 'C 2', 'CO 1', 'CO 2', 'CO 3', 'CO 4', 'CO 5', 'CO 6', 'CO 7', 'CO 8',
                            '13C+ 1', '13C 1', '13C 2', '13CO 1', '13CO 2', '13CO 3', '13CO 4', '13CO 5', 'HCO+ 1'],
              # 'dust': 'PAH',
              'dust': '240um',
              'clumpMassRange': [[0, 2], [-2]],
              'clumpMassNumber': [3, 1],
              'clumpNmax': [1, 100],
              'velocityRange': [-350, 350],
              'velocityNumber': 701,

              # Property factors
              'clumpMassFactor': [1, 1],
              'FUVfactor': 1,
              'densityFactor': 1,
              'globalUV': 10,
              'r_cmz': 0,

              # Logging
              'timed': False,
              'debug': False,
              'verbose': False
              }

f_mass1 = 1.0
f_mass2 = 1.0
f_density = 1.0
f_uv = 10
mass_cl = [0, 2]
mass_icl = [-2]
r_cmz = 0
models.constants.fuv_ism = 2

for resolution in [500]:
        # for f_mass1 in [0.25, 0.5, 1.0, 2.0, 4.0]:
        #     for f_mass2 in [0.25, 0.5, 1.0, 2.0, 4.0]:
        # if (mass_icl == [-2]) and not (f_mass2 == 4.0):
        #     continue
        # for mass_cl in [[0, 2], [0, 3], [-1, 2], [-1, 3]]:
        #     for mass_icl in [[-2], [-3, -2], [-3, -1], [-2, -1]]:
        #         for r_cmz in np.arange(2200, 3001, 200, dtype=np.int):
        #             for f_uv in [10, 50, 100]:
        #                 for f_density in [0.5, 1.0, 2.0, 4.0]:
        #                     for f_uv in [10, 100]:

        parameters['resolution'] = resolution
        parameters['clumpMassFactor'] = [f_mass1, f_mass2]
        parameters['clumpMassRange'] = [mass_cl, mass_icl]
        parameters['clumpMassNumber'] = [len(mass_cl), len(mass_icl)]
        parameters['densityFactor'] = f_density
        parameters['globalUV'] = f_uv
        parameters['r_cmz'] = r_cmz
        parameters['folder'] = 'r{}/'.format(int(resolution))


        kosma = models.Model(**parameters)
        print('    Model {}'.format(models.constants.history))
        kosma.calculateModel(multiprocessing=0)

        directory = parameters['history_path'] + parameters['directory'] + '/' + parameters['folder']
        models.radiativeTransfer.calculateObservation(directory=directory, dim='spherical', multiprocessing=6,
                                                      slRange=[(-np.pi, np.pi), (-2*np.pi/180, 2*np.pi/180)],
                                                      nsl=[361, 5], terminal=True, debug=False)
