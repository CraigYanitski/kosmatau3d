# python3.8

import argparse
import numpy as np

from kosmatau3d import models
from pprint import pprint


# Parse arguments
parser = argparse.ArgumentParser()

parser.add_argument('-f', '--folder', type=str, default='/mnt/hpc_backup/yanitski/projects/pdr/KT3_history', 
                    help='folder containing kosmatau3d models')
parser.add_argument('-g', '--grid', type=str, default='convergence', 
                    choices=['convergence', 'f_cm, cm', 'f_icm, icm', 'cm, icm', 'r_cmz, f_fuv',
                             'f_hi, f_fuv', 'f_cl, f_icl, f_d, f_fuv'], 
                    help='parameters varied in grid')
parser.add_argument('-r', '--resolution', type=int, default=400, 
                    help='voxel size in the model')
parser.add_argument('-m', '--mp', type=int, default=8, 
                    help='number of threads to open during RT \n  set to 0 for debugging')
parser.add_argument('-e', '--eval', type=str, default='all', choices=['all', 'model', 'RT'], 
                    help='which portion of the code to evaluate')
parser.add_argument('-d', '--debug', type=bool, default=False, 
                    help='use debugging mode')
parser.add_argument('-t', '--timed', type=bool, default=False, 
                    help='use timing mode')
parser.add_argument('-v', '--verbose', type=bool, default=False, 
                    help='use verbose mode')

args = parser.parse_args()

if args.grid == 'convergence':
    folder = 'r{}_convergence/'
    params = (400, 200, 100)
    param_keys = ('resolution', )
    param_folders = list(folder.format(*_) for _ in zip(params))
elif args.grid == 'f_cm, cm':
    folder = 'r{}_fcm{}_cm{}/'
    fmass1 = [0.25, 0.5, 1.0, 2.0, 4.0]
    mass_cl = [[0, 2], [0, 3], [-1, 2], [-1, 3]]
    param_keys = ('fmass1', 'clump_mass_range', 'clump_mass_number')
    params = list(_.flatten() for _ in np.meshgrid(fmass1, mass_cl))
    for _ in range(len(params[1])):
        params[1][_] = [params[_], [-2]]
    params.append(list([int(np.max(_)-np.min(_))+1, 1] for _ in params[1]))
    mass_cl_strs = list('_'.join([str(_) for _ in mass]) for mass in params[1])
    param_folders = list(folder.format(*_) for _ in zip((params[0], mass_cl_strs)))
elif args.grid == 'f_icm, icm':
    folder = 'r{}_ficm{}_icm{}/'
    fmass2 = [0.25, 0.5, 1.0, 2.0, 4.0]
    mass_icl = [[-2], [-3, -2], [-3, -1], [-2, -1]]
    params = list(_.flatten() for _ in np.meshgrid(fmass2, mass_icl))
    mass_icl_strs = list('_'.join([str(_) for _ in mass]) for mass in params[1])
    param_folders = list(folder.format(*_) for _ in zip((params[0], mass_icl_strs)))
elif args.grid == 'cm, icm':
    folder = 'r{}_cm{}_icm{}/'
    fmass1 = [0.25, 0.5, 1.0, 2.0, 4.0]
    fmass2 = [0.25, 0.5, 1.0, 2.0, 4.0]
    mass_cl = [[0, 2], [0, 3], [-1, 2], [-1, 3]]
    mass_icl = [[-2], [-3, -2], [-3, -1], [-2, -1]]
    params = list(_.flatten() for _ in np.meshgrid(mass_cl, mass_icl))
    mass_cl_strs = list('_'.join([str(_) for _ in mass]) for mass in params[0])
    mass_icl_strs = list('_'.join([str(_) for _ in mass]) for mass in params[1])
    param_folders = list(folder.format(*_) for _ in zip((mass_cl_strs, mass_icl_strs)))
elif args.grid == 'r_cmz, f_fuv':
    folder = 'r{}_r_cmz{}_f_fuv{:.1f}/'
    r_cmz = np.arange(0, 1001, 50, dtype=int)
    f_fuv = 10**np.linspace(0, 2, num=17)
    params = list(_.flatten() for _ in np.meshgrid(r_cmz, f_fuv))
    param_folders = list(folder.format(*_) for _ in zip(params))
elif args.grid == 'f_hi, f_fuv':
    folder = 'r{}_fhi{}_fuv{:.1f}/'
    f_hi = [0.7, 0.6, 0.5, 0.4, 0.3]
    f_fuv = 10**np.linspace(0, 2, num=17)
    params = list(_.flatten() for _ in np.meshgrid(f_hi, f_fuv))
    param_folders = list(folder.format(*_) for _ in zip(params))
elif args.grid == 'f_cl, f_icl, f_d, f_fuv':
    folder = 'r{}_cm{}-{}_d{}_uv{}/'
    f_cl = [0.5, 1.0, 2.0, 4.0]
    f_icl = [1.0, 2.0]
    f_d = [0.5, 1.0, 2.0, 4.0]
    f_fuv = [1, 10, 100]
    params = list(_.flatten() for _ in np.meshgrid(f_cl, f_icl, f_d, f_fuv))
    param_folders = list(folder.format(*_) for _ in zip(params))
else:
    exit('Invalid grid queried; exiting...')

debug = args.debug

#print('KOSMA-tau^3')

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
              'new_grid': False,
              
              # Model information
              'history_path': args.folder,
              'directory': r'/MilkyWay',
              'folder': '',
              'x': 36,
              'y': 36,
              'z': 0,
              'model_type': 'disk',

              # Model parameters
              'resolution': args.resolution,
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
              'timed': args.timed,
              'debug': args.debug,
              'verbose': args.verbose
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

for i, param in enumerate(zip(params)):
                    for _, p in enumerate(param_keys):
                        parameters[p] = param[_]
                    parameters['folder'] = param_folders[i] 
# for resolution in [400]:
#     for mass_cl in [[0, 2], [0, 3], [-1, 2], [-1, 3]]:
#         for mass_icl in [[-2], [-3, -2], [-3, -1], [-2, -1]]:
#                     mass_cl_str = '_'.join([str(_) for _ in mass_cl])
#                     mass_icl_str = '_'.join([str(_) for _ in mass_icl])
#                     
#                     parameters['folder'] = f'r{int(resolution)}_cm{mass_cl_str}_icm{mass_icl_str}/'
                    

                    #f'r{int(resolution)}_ficm{f_mass1}_icm{mass_cl_str}/'
                    #for f_mass1 in [0.25, 0.5, 1.0, 2.0, 4.0]:
                    #    for mass_cl in [[0, 2], [0, 3], [-1, 2], [-1, 3]]:
                    #                mass_cl_str = '_'.join([str(_) for _ in mass_cl])
                    #f'r{int(resolution)}_ficm{f_mass2}_icm{mass_icl_str}/'
                    #for f_mass2 in [0.25, 0.5, 1.0, 2.0, 4.0]:
                    #    for mass_icl in [[-2], [-3, -2], [-3, -1], [-2, -1]]:
                    #                mass_icl_str = '_'.join([str(_) for _ in mass_icl])
                    #
                    #for mass_cl in [[0, 2], [0, 3], [-1, 2], [-1, 3]]:
                    #    for mass_icl in [[-2], [-3, -2], [-3, -1], [-2, -1]]:
                    #                mass_cl_str = '_'.join([str(_) for _ in mass_cl])
                    #                mass_icl_str = '_'.join([str(_) for _ in mass_icl])
                    #
                    #for f_uv in 10**np.linspace(0, 2, num=17):
                    #    for r_cmz in np.arange(0, 1001, 50, dtype=int):
                    #
                    #for f_hi in [0.7, 0.6, 0.5, 0.4, 0.3]:
                    #    for f_uv in 10**np.linspace(0, 2, num=17):
                    #
                    #for f_mass1 in [1.0, 2.0]:
                    #    for f_mass2 in [0.5, 1.0, 2.0, 4.0]:
                    #        for f_density in [0.5, 1.0, 2.0, 4.0]:
                    #            for f_uv in [10, 100]:
                    #
                    #for mass_cl in [[0, 2], [0, 3], [-1, 2], [-1, 3]]:
                    #    for mass_icl in [[-2], [-3, -2], [-3, -1], [-2, -1]]:
                    #
                    #for f_hi in [1.0, 0.8, 0.6, 0.4]:
                    #    for f_uv in 10**np.linspace(1, 3, num=5):
                    #        for fuv in [1.0, 1.4, 1.8, 2.2]:
                    # 
                    #for f_hi in [1.0, 0.8, 0.6, 0.4]:
                    #    for f_uv in 10**np.linspace(-1, 3, num=5):
                    #        for fuv in [1.0, 1.4, 1.8, 2.2]:
                    #
                    #for fuv in [1.4, 1.6, 1.8, 2.0]:
                    #    for f_hi in [1.0, 0.8, 0.6]:
                    #
                    #for f_mass1 in [0.25, 0.5, 1.0, 2.0, 4.0]:
                    #    for f_mass2 in [0.25, 0.5, 1.0, 2.0, 4.0]:
                    #
                    #for f_uv in 10**np.linspace(1.5, 2.5, num=5):
                    #    for r_cmz in np.arange(0, 3001, 500, dtype=int):
                    #
                    #for mass_cl in [[0, 2], [0, 3], [-1, 2], [-1, 3]]:
                    #    for mass_icl in [[-2], [-3, -2], [-3, -1], [-2, -1]]:
                    #        for f_uv in [10, 50, 100]:
                    #            for f_density in [0.5, 1.0, 2.0, 4.0]:
                    #                for f_uv in [10, 100]:

                    # if (mass_icl == [-2]) and not (f_mass2 == 4.0):
                    #     continue
                    # if f_hi == 1.0 and f_uv < 1000:
                    #     continue
                    # if fuv < 1.8:
                    #     continue

                    # # forces the clump and interclump ensemble masses to be
                    # # equivalent to the H2 mass
                    # parameters['like_clumps'] = False
                    # # uses the average fuv energy density
                    # parameters['average_fuv'] = False
                    # # set the FUV range to integrate
                    # parameters['l_range'] = (912, 2066)

                    # parameters['resolution'] = resolution
                    # parameters['clump_mass_range'] = [mass_cl, mass_icl]
                    # clump_mass_number = []
                    # if len(mass_cl) > 1:
                    #     clump_mass_number.append(int(np.max(mass_cl)-np.min(mass_cl))+1)
                    # else:
                    #     clump_mass_number.append(len(mass_cl))
                    # if len(mass_icl) > 1:
                    #     clump_mass_number.append(int(np.max(mass_icl)-np.min(mass_icl))+1)
                    # else:
                    #     clump_mass_number.append(len(mass_icl))
                    # print(clump_mass_number)
                    # parameters['clump_mass_number'] = clump_mass_number
                    # parameters['clump_log_fuv'] = None
                    # parameters['interclump_log_fuv'] = None

                    # parameters['hi_mass_factor'] = f_mass1
                    # parameters['h2_mass_factor'] = f_mass2
                    # parameters['ensemble_mass_factor'] = f_m_ens
                    # parameters['density_factor'] = f_density
                    # parameters['fuv_factor'] = f_uv
                    # parameters['interclump_hi_ratio'] = f_hi
                    # parameters['r_cmz'] = r_cmz

                    # parameters['suggested_calc'] = True
                    # parameters['new_grid'] = False

                    pprint(parameters)
                    kosma = models.Model(**parameters)
                    print('\n    -> Model {}'.format(models.constants.history))
                    print('       ' + '-'*len('Model {}'.format(models.constants.history)))
                    #kosma.calculateModel(kind='zero', multiprocessing=0)

                    # print(models.species.molecules)

                    directory = parameters['history_path'] + parameters['directory'] \
                                + '/' + parameters['folder']
                    models.radiativeTransfer.calculateObservation(directory=directory, dim='spherical',
                                                                  multiprocessing=8, vel_pool=False, 
                                                                  slRange=[(-np.pi, np.pi), 
                                                                           (-2*np.pi/180, 2*np.pi/180)],
                                                                  nsl=[721, 9], terminal=True, debug=False)

