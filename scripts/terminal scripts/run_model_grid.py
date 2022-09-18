# python3.8

import argparse
import numpy as np

from kosmatau3d import models
from pprint import pprint


# ARGUMENTS
#, fontsize=16 Parse arguments
parser = argparse.ArgumentParser()

parser.add_argument('-f', '--folder', type=str, default='/mnt/yanitski_backup/yanitski/projects/pdr/KT3_history', 
                    help='folder containing kosmatau3d models')
parser.add_argument('-g', '--grid', type=str, default='convergence', 
                    choices=['convergence', 'f_cm-cm', 'f_icm-icm', 'cm-icm', 'f_cm-f_icm', 
                             'r_cmz-f_fuv', 'f_hi-f_fuv', 'f_cl-f_icl-f_d-f_fuv', 
                             'f_fuv_gc', 'disp_gc', 'fuv_gc-disp_gc'], 
                    help='parameters varied in grid')
parser.add_argument('-r', '--resolution', type=int, default=400, 
                    help='voxel size in the model')
parser.add_argument('-m', '--mp', type=int, default=8, 
                    help='number of threads to open during RT \n  set to 0 for debugging')
parser.add_argument('-e', '--exec', type=str, default='all', choices=['all', 'model', 'RT'], 
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
    res = (400, 200, 100)
    grid_flag = (True, )
    param_keys = ('resolution', 'new_grid')
    params = list(_.flatten() for _ in np.meshgrid(res, grid_flag))
    param_folders = list(folder.format(*_) for _ in zip(*params))
elif args.grid == 'f_cm-cm':
    folder = f'r{args.resolution}' + '_fcm{}_cm{}/'
    fmass1 = [0.25, 0.5, 1.0, 2.0, 4.0]
    mass_cl = [[0, 2], [0, 3], [-1, 2], [-1, 3]]
    param_keys = ('h2_mass_factor', 'clump_mass_range', 'clump_mass_number')
    fmass1_mesh, i_cl_mesh = list(_.flatten() for _ in np.meshgrid(fmass1, np.arange(len(mass_cl))))
    mass_bin_mesh = []
    for _ in range(i_cl_mesh.size):
        mass_bin_mesh.append([mass_cl[i_cl_mesh[_]], [-2]])
    mass_num_mesh = list(list(int(np.max(mbin)-np.min(mbin))+1 for mbin in _) for _ in mass_bin_mesh)
    params = (fmass1_mesh, mass_bin_mesh, mass_num_mesh)
    mass_cl_strs = list('_'.join([str(_) for _ in mass[0]]) for mass in params[1])
    param_folders = list(folder.format(*_) for _ in zip(*(params[0], mass_cl_strs)))
elif args.grid == 'f_icm-icm':
    folder = f'r{args.resolution}' + '_ficm{}_icm{}/'
    fmass2 = [0.25, 0.5, 1.0, 2.0, 4.0]
    mass_icl = [[-2], [-3, -2], [-3, -1], [-2, -1]]
    param_keys = ('hi_mass_factor', 'clump_mass_range', 'clump_mass_number')
    fmass2_mesh, i_icl_mesh = list(_.flatten() for _ in np.meshgrid(fmass2, np.arange(len(mass_icl))))
    mass_bin_mesh = []
    for _ in range(i_icl_mesh.size):
        mass_bin_mesh.append([[0, 2], mass_icl[i_icl_mesh[_]]])
    mass_num_mesh = list(list(int(np.max(mbin)-np.min(mbin))+1 for mbin in _) for _ in mass_bin_mesh)
    params = (fmass2_mesh, mass_bin_mesh, mass_num_mesh)
    mass_icl_strs = list('_'.join([str(_) for _ in mass[1]]) for mass in params[1])
    param_folders = list(folder.format(*_) for _ in zip(*(params[0], mass_icl_strs)))
elif args.grid == 'cm-icm':
    folder = f'r{args.resolution}' + '_cm{}_icm{}/'
    mass_cl = [[0, 2], [0, 3], [-1, 2], [-1, 3]]
    mass_icl = [[-2], [-3, -2], [-3, -1], [-2, -1]]
    param_keys = ('clump_mass_range', 'clump_mass_number')
    i_cl_mesh, i_icl_mesh = list(_.flatten() for _ in np.meshgrid(np.arange(len(mass_cl)), 
                                                                  np.arange(len(mass_icl))))
    mass_bin_mesh = []
    for _ in range(i_cl_mesh.size):
        mass_bin_mesh.append([mass_cl[i_cl_mesh[_]], mass_icl[i_icl_mesh[_]]])
    mass_num_mesh = list(list(int(np.max(mbin)-np.min(mbin))+1 for mbin in _) for _ in mass_bin_mesh)
    params = (mass_bin_mesh, mass_num_mesh)
    mass_cl_strs = list('_'.join([str(_) for _ in mass[0]]) for mass in params[0])
    mass_icl_strs = list('_'.join([str(_) for _ in mass[1]]) for mass in params[0])
    param_folders = list(folder.format(*_) for _ in zip(*(mass_cl_strs, mass_icl_strs)))
elif args.grid == 'f_cm-f_icm':
    folder = f'r{args.resolution}' + '_f_cm{}_f_icm{:}/'
    fmass1 = [0.25, 0.5, 1.0, 2.0, 4.0]
    fmass2 = [0.25, 0.5, 1.0, 2.0, 4.0]
    param_keys = ('h2_mass_factor', 'hi_mass_factor')
    params = list(_.flatten() for _ in np.meshgrid(fmass1, fmass2))
    param_folders = list(folder.format(*_) for _ in zip(*params))
elif args.grid == 'r_cmz-f_fuv':
    folder = f'r{args.resolution}' + '_rcmz{}_f_fuv{:.1f}/'
    r_cmz = np.arange(0, 1001, 50, dtype=int)
    f_fuv = 10**np.linspace(0, 1.5, num=13)
    param_keys = ('r_cmz', 'fuv_factor')
    params = list(_.flatten() for _ in np.meshgrid(r_cmz, f_fuv))
    param_folders = list(folder.format(*_) for _ in zip(*params))
elif args.grid == 'f_hi-f_fuv':
    folder = f'r{args.resolution}' + '_fhi{}_fuv{:.1f}/'
    f_hi = [1.0, 0.9, 0.8, 0.7, 0.6, 0.5, 0.4, 0.3]
    f_fuv = 10**np.linspace(0, 1.5, num=13)
    param_keys = ('interclump_hi_ratio', 'fuv_factor')
    params = list(_.flatten() for _ in np.meshgrid(f_hi, f_fuv))
    param_folders = list(folder.format(*_) for _ in zip(*params))
elif args.grid == 'f_cl-f_icl-f_d-f_fuv':
    folder = f'r{args.resolution}' + '_cm{}-{}_d{}_uv{}/'
    f_cl = [0.5, 1.0, 2.0, 4.0]
    f_icl = [1.0, 2.0]
    f_d = [0.5, 1.0, 2.0, 4.0]
    f_fuv = [1, 10, 100]
    param_keys = ('h2_mass_factor', 'hi_mass_factor', 'density_factor', 'fuv_factor')
    params = list(_.flatten() for _ in np.meshgrid(f_cl, f_icl, f_d, f_fuv))
    param_folders = list(folder.format(*_) for _ in zip(*params))
elif args.grid == 'f_fuv_gc':
    folder = f'r{args.resolution}' + '_f_fuv_gc{:.2f}/'
    params = (10**np.array([0, 0.5, 1, 1.5, 2]), )
    param_keys = ('scale_gc', )
    param_folders = list(folder.format(*_) for _ in zip(*params))
elif args.grid == 'disp_gc':
    folder = f'r{args.resolution}' + '_r_gc{:.2f}_disp_gc{:.2f}/'
    r_gc = (200, 600)
    disp_gc = (200, 300)
    param_keys = ('disp_core', 'r_core', )
    params = list(_.flatten() for _ in np.meshgrid(r_gc, disp_gc))
    param_folders = list(folder.format(*_) for _ in zip(*params))
elif args.grid == 'fuv_gc-disp_gc':
    folder = f'r{args.resolution}' + '_f_fuv_gc{:.2f}_r_fuv_gc{}_disp_gc{}_r_disp_gc{}/'
    r_gc = (0, 400, 800, 1200, 1600, 2000)
    fuv_gc = (10**np.array([0, 0.5, 1, 1.5, 2]), )
    r_core = (200, 600)
    disp_gc = (200, 300)
    grid_flag = (True, )
    param_keys = ('scale_gc', 'r_gc', 'disp_core', 'r_core', 'new_grid')
    params = list(_.flatten() for _ in np.meshgrid(fuv_gc, r_gc, disp_gc, r_core, grid_flag))
    param_folders = list(folder.format(*_) for _ in zip(*params))
else:
    params = None
    param_keys = None
    param_folders = None
    exit('Invalid grid queried; exiting...')

if args.exec == 'model':
    run_model = True
    run_rt = False
elif args.exec == 'RT':
    run_model = False
    run_rt = True
elif args.exec == 'all':
    run_model = True
    run_rt = True
else:
    run_model = False
    run_rt = False

debug = args.debug


# KOSMATAU3D PARAMETERS
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
              'scale_gc': 1.0,
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

i_max = len(params[0])

for i, param in enumerate(zip(*params)):
    
    # Set the required parameters
    for _, p in enumerate(param_keys):
        parameters[p] = param[_]
    parameters['folder'] = param_folders[i] 

    # Initialise model instance
    kosma = models.Model(**parameters)

    # Print details
    print('\n\n   ==> Model {} of {}'.format(i, i_max))
    print('       ' + models.constants.history)
    print('       ' + '-'*len(models.constants.history))
    if args.verbose:
        print()
        pprint(parameters)
    print()
    print("resolution: ".rjust(25) + f"{parameters['resolution']}")
    print("clump mass range: ".rjust(25) + f"{parameters['clump_mass_range']}")
    print("clump mass number: ".rjust(25) + f"{parameters['clump_mass_number']}")
    print()
    print("R_CMZ: ".rjust(25) + f"{parameters['r_cmz']}")
    print("ensemble mass factor: ".rjust(25) + f"{parameters['ensemble_mass_factor']}")
    print("H2 mass factor: ".rjust(25) + f"{parameters['h2_mass_factor']}")
    print("HI mass factor: ".rjust(25) + f"{parameters['hi_mass_factor']}")
    print("HI interclump ratio: ".rjust(25) + f"{parameters['interclump_hi_ratio']}")
    print("clump density factor: ".rjust(25) + f"{parameters['density_factor']}")
    print("FUV factor: ".rjust(25) + f"{parameters['fuv_factor']}")
    print()
    print('    -- parameters changed')
    for _, p in enumerate(param_keys):
        print(f'  {p} ->'.rjust(25) + f' {param[_]}')
    print()

    # Run model if selected
    if run_model:
        kosma.calculateModel(kind='linear', multiprocessing=0)

    # Set model folder directory
    directory = parameters['history_path'] + parameters['directory'] \
                + '/' + parameters['folder']

    # Run radiative transfer if selected
    if run_rt:
        models.radiativeTransfer.calculateObservation(directory=directory, dim='spherical',
                                                      multiprocessing=args.mp, 
                                                      vel_pool=False, pencil_beam=True, 
                                                      slRange=[(-np.pi, np.pi), 
                                                               (-2*np.pi/180, 2*np.pi/180)],
                                                      nsl=[721, 9], terminal=True, debug=False)

