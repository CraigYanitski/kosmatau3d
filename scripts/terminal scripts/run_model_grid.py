# python3.8

import numpy as np

from kosmatau3d import models

timed = False
debug = False

print('KOSMA-tau^3')

# Edit these parameters according to the model you want to produce.
parameters = {
              # Data files
              'tau_grid_file': 'clump_tau_LineCenter.dat',
              'tb_grid_file': 'clump_Tmb_LineCenter.dat',
              'tau_fuv_grid_file': 'RhoMassAFUV.dat',
              'h2_surface_density_file': 'h2_surface-density_marasco-bacchini.dat', 
              'hi_surface_density_file': 'hi_surface-density_marasco-bacchini.dat', 
              'h2_scale_height_file': 'h2_scale-height_bacchini.dat', 
              'hi_scale_height_file': 'hi_scale-height_bacchini.dat', 
              'like_clumps': False, 
              'h_number_density_file': 'h_number-density_cubick.dat',
              'fuv_file': 'galactic_FUV_complete.dat',
              'average_fuv': False,
              'l_range': (912, 2066),
              
              # Model information
              'history_path': r'/mnt/hpc_backup/yanitski/projects/pdr/KT3_history',
              'directory': r'/MilkyWay',
              'folder': '',
              'x': 36,
              'y': 36,
              'z': 0,
              'model_type': 'disk',

              # Model parameters
              'resolution': 400,
              # 'molecules': 'all',
              'molecules': ['C+ 1', 
                            'C 1', 'C 2', 
                            'CO 1', 'CO 2', 'CO 3', 'CO 4', 'CO 5', 
                            'CO 6', 'CO 7', 'CO 8', 
                            '13C+ 1', 
                            '13C 1', '13C 2', 
                            '13CO 1', '13CO 2', '13CO 3', '13CO 4', '13CO 5', 
                            '13CO 6', '13CO 7', '13CO 8', 
                            'HCO+ 1', 'HCO+ 2', 'HCO+ 3', 'HCO+ 4', 'HCO+ 5'],
              # 'dust': 'PAH',
              'dust': ['240um', '550um'],
              'clump_mass_range': [[0, 2], [-2]],
              'clump_mass_number': [3, 1],
              'clump_n_max': [1, 100],
              'clump_log_fuv' : None,
              'interclump_log_fuv' : None,
              'velocity_range': [-350, 350],
              'velocity_number': 701,

              # Flags
              'suggested_calc': True,
              'dilled': False,

              # Property factors
              'hi_mass_factor': 1,
              'h2_mass_factor': 1,
              'ensemble_mass_factor': [1, 1],
              'fuv_factor': 1,
              'density_factor': 1,
              'interclump_hi_ratio': 1,
              'r_cmz': 0,

              # Logging
              'timed': False,
              'debug': False,
              'verbose': False
              }

f_mass1 = 1.0
f_mass2 = 1.0
f_m_ens = [1, 1]
f_hi = 1
f_density = 1.0
f_uv = 1
mass_cl = [0, 2]
mass_icl = [-2]
r_cmz = 0
models.constants.fuv_ism = 1
fuv = 1.8

for resolution in [400]:
                    #for f_mass1 in [0.25, 0.5, 1.0, 2.0, 4.0]:
                    #    for mass_cl in [[0, 2], [0, 3], [-1, 2], [-1, 3]]:
                    #                mass_cl_str = '_'.join([str(_) for _ in mass_cl])
                    #for f_mass2 in [0.25, 0.5, 1.0, 2.0, 4.0]:
                    #    for mass_icl in [[0, 2], [0, 3], [-1, 2], [-1, 3]]:
                    #                mass_icl_str = '_'.join([str(_) for _ in mass_icl])
                    #for mass_cl in [[0, 2], [0, 3], [-1, 2], [-1, 3]]:
                    #    for mass_icl in [[-2], [-3, -2], [-3, -1], [-2, -1]]:
                    #                mass_cl_str = '_'.join([str(_) for _ in mass_cl])
                    #                mass_icl_str = '_'.join([str(_) for _ in mass_icl])
                    #for f_uv in 10**np.linspace(0, 2, num=17):
                    #    for r_cmz in np.arange(0, 1001, 50, dtype=int):
                    #for f_hi in [0.7, 0.6, 0.5, 0.4, 0.3]:
                    #    for f_uv in 10**np.linspace(0, 2, num=17):
                    #for f_mass1 in [1.0, 2.0]:
                    #    for f_mass2 in [0.5, 1.0, 2.0, 4.0]:
                    #        for f_density in [0.5, 1.0, 2.0, 4.0]:
                    #            for f_uv in [10, 100]:
                    # for mass_cl in [[0, 2], [0, 3], [-1, 2], [-1, 3]]:
                    #     for mass_icl in [[-2], [-3, -2], [-3, -1], [-2, -1]]:
                    # for f_hi in [1.0, 0.8, 0.6, 0.4]:
                    #     for f_uv in 10**np.linspace(1, 3, num=5):
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
                    # for mass_cl in [[0, 2], [0, 3], [-1, 2], [-1, 3]]:
                    #     for mass_icl in [[-2], [-3, -2], [-3, -1], [-2, -1]]:
                    #             for f_uv in [10, 50, 100]:
                    #                 for f_density in [0.5, 1.0, 2.0, 4.0]:
                    #                     for f_uv in [10, 100]:

                    # if (mass_icl == [-2]) and not (f_mass2 == 4.0):
                    #     continue
                    # if f_hi == 1.0 and f_uv < 1000:
                    #     continue
                    # if fuv < 1.8:
                    #     continue

                    # forces the clump and interclump ensemble masses to be
                    # equivalent to the H2 mass
                    parameters['like_clumps'] = False
                    # uses the average fuv energy density
                    parameters['average_fuv'] = False
                    # set the FUV range to integrate
                    parameters['l_range'] = (912, 2066)

                    parameters['resolution'] = resolution
                    parameters['clump_mass_range'] = [mass_cl, mass_icl]
                    clump_mass_number = []
                    if len(mass_cl) > 1:
                        clump_mass_number.append(int(np.max(mass_cl)-np.min(mass_cl))+1)
                    else:
                        clump_mass_number.append(len(mass_cl))
                    if len(mass_icl) > 1:
                        clump_mass_number.append(int(np.max(mass_icl)-np.min(mass_icl))+1)
                    else:
                        clump_mass_number.append(len(mass_icl))
                    print(clump_mass_number)
                    parameters['clump_mass_number'] = clump_mass_number
                    parameters['clump_log_fuv'] = None
                    parameters['interclump_log_fuv'] = None
                    parameters['folder'] = f'r{int(resolution)}_convergence/'

                    parameters['hi_mass_factor'] = f_mass1
                    parameters['h2_mass_factor'] = f_mass2
                    parameters['ensemble_mass_factor'] = f_m_ens
                    parameters['density_factor'] = f_density
                    parameters['fuv_factor'] = f_uv
                    parameters['interclump_hi_ratio'] = f_hi
                    parameters['r_cmz'] = r_cmz

                    parameters['suggested_calc'] = True
                    parameters['new_grid'] = False


                    kosma = models.Model(**parameters)
                    print('\n    -> Model {}'.format(models.constants.history))
                    print('       ' + '-'*len('Model {}'.format(models.constants.history)))
                    kosma.calculateModel(kind='zero', multiprocessing=0)

                    # print(models.species.molecules)

                    directory = parameters['history_path'] + parameters['directory'] \
                                + '/' + parameters['folder']
                    models.radiativeTransfer.calculateObservation(directory=directory, dim='spherical',
                                                                  multiprocessing=8, vel_pool=False, 
                                                                  slRange=[(-np.pi, np.pi), 
                                                                           (-2*np.pi/180, 2*np.pi/180)],
                                                                  nsl=[721, 9], terminal=True, debug=False)

