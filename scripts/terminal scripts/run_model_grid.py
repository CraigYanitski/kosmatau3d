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
              'resolution': 500,
              # 'molecules': 'all',
              'molecules': ['C+ 1', 'C 1', 'C 2', 'CO 1', 'CO 2', 'CO 3', 'CO 4', 'CO 5', 'CO 6', 'CO 7', 'CO 8',
                            '13C+ 1', '13C 1', '13C 2', '13CO 1', '13CO 2', '13CO 3', '13CO 4', '13CO 5', '13CO 6',
                            '13CO 7', '13CO 8', 'HCO+ 1', 'HCO+ 2', 'HCO+ 3', 'HCO+ 4', 'HCO+ 5'],
              # 'dust': 'PAH',
              'dust': ['240um', '550um'],
              'clumpMassRange': [[0, 2], [-2]],
              'clumpMassNumber': [3, 1],
              'clumpNmax': [1, 100],
              'clumpLogFUV' : None,
              'interclumpLogFUV' : None,
              'velocityRange': [-350, 350],
              'velocityNumber': 701,

              # Property factors
              'clumpMassFactor': [1, 1],
              'FUVfactor': 1,
              'densityFactor': 1,
              'r_cmz': 0,

              # Logging
              'timed': False,
              'debug': False,
              'verbose': False
              }

f_mass1 = 1.0
f_mass2 = 1.0
f_hi = 1
f_density = 1.0
f_uv = 1
mass_cl = [0, 2]
mass_icl = [-2]
r_cmz = 0
models.constants.fuv_ism = 2

for resolution in [500]:
    # for f_hi in [1.0, 0.8, 0.6, 0.4]:
    #     for f_uv in 10**np.linspace(-1, 3, num=5):
    #         for fuv in [1.0, 1.4, 1.8, 2.2]:
                # for f_hi in [1.0, 0.8, 0.6, 0.4]:
                #     for f_uv in 10**np.linspace(-1, 3, num=5):
                #         for fuv in [1.0, 1.4, 1.8, 2.2]:
                # for fuv in [1.4, 1.6, 1.8, 2.0]:
                #     for f_hi in [1.0, 0.8, 0.6]:
                # for f_mass1 in [0.25, 0.5, 1.0, 2.0, 4.0]:
                #     for f_mass2 in [0.25, 0.5, 1.0, 2.0, 4.0]:
                # for f_uv in 10**np.linspace(1.5, 2.5, num=5):
                #     for r_cmz in np.arange(0, 3001, 500, dtype=int):
                # if (mass_icl == [-2]) and not (f_mass2 == 4.0):
                #     continue
                # for mass_cl in [[0, 2], [0, 3], [-1, 2], [-1, 3]]:
                #     for mass_icl in [[-2], [-3, -2], [-3, -1], [-2, -1]]:
                #             for f_uv in [10, 50, 100]:
                #                 for f_density in [0.5, 1.0, 2.0, 4.0]:
                #                     for f_uv in [10, 100]:

                # if f_hi == 1.0 and f_uv < 1000:
                #     continue
                # if fuv < 1.8:
                #     continue

                parameters['resolution'] = resolution
                parameters['clumpMassFactor'] = [f_mass1, f_mass2]
                parameters['interclump_hifactor'] = f_hi
                parameters['clumpMassRange'] = [mass_cl, mass_icl]
                parameters['clumpMassNumber'] = [len(mass_cl), len(mass_icl)]
                parameters['clumpLogFUV'] = None#fuv
                parameters['interclumpLogFUV'] = None 
                parameters['densityFactor'] = f_density
                parameters['FUVfactor'] = f_uv
                parameters['r_cmz'] = r_cmz
                parameters['folder'] = 'r{}_modtest/'.format(int(resolution))


                kosma = models.Model(**parameters)
                print('\n    -> Model {}'.format(models.constants.history))
                print('       ' + '-'*len('Model {}'.format(models.constants.history)))
                kosma.calculateModel(dilled=True, multiprocessing=0)

                directory = parameters['history_path'] + parameters['directory'] + '/' + parameters['folder']
                models.radiativeTransfer.calculateObservation(directory=directory, dim='spherical', multiprocessing=8,
                                                              slRange=[(-np.pi, np.pi), (-2*np.pi/180, 2*np.pi/180)],
                                                              nsl=[361, 5], terminal=True, debug=False)
