#!/home/yanitski/anaconda3/lib/python3.8

'''
Available directories:
 - '/r{}_hi{}_f_fuv{:.2f}_const_fuv{}/'
   [[1.0, 0.8, 0.6, 0.4], ['0.10', '1.00', '10.00', '100.00', '1000.00'], [1.0, 1.4, 1.8, 2.2]]
 - '/r400_cm{}-{}_d{}_uv{:.0f}/'
   [[0.25, 0.5, 1.0, 2.0, 4.0], [1.0, 2.0], [0.25, 0.5, 1.0, 2.0, 4.0], [10, 100]]
 - '/r400_fcm{}_ficm{}/'
   [[0.25, 0.5, 1.0, 2.0, 4.0], [0.25, 0.5, 1.0, 2.0, 4.0]]
 - '/r400_fcm{}_cm{}/'
   [[0.25, 0.5, 1.0, 2.0, 4.0], ['0_2', '0_3', '-1_2', '-1_3']]
 - '/r400_ficm{}_icm{}/'
   [[0.25, 0.5, 1.0, 2.0, 4.0], ['-2', '-2_-1', '-3_-1', '-3_-2']]
 - '/r400_rcmz{:.0f}_uv{}/'
   [np.arange(0, 3001, 200), [10, 50, 100]]
'''


import numpy as np

from copy import deepcopy
from kosmatau3d import comparison as comp


# Folders and parameters of models
folders = ['/r400_hi{}_f_fuv{}_const_fuv1.8/']
          # ['/r400_fhi{}_f_uv{:.0f}_const_fuv{}/']
          # ['/r400_fhi{}_fuv{:.0f}/']
          # ['/r500_rcmz{}_fuv{}/',
          #  '/r500_fcm{}_fim{}/']
          # ['/r400_cm{}-{}_d{}_uv{:.0f}/',
          #  '/r400_fcm{}_ficm{}/',
          #  '/r400_fcm{}_cm{}/',
          #  '/r400_ficm{}_icm{}/',
          #  '/r400_rcmz{:.0f}_uv{}/']
parameters = [[[1.0, 0.8, 0.6, 0.4], ['0.10', '1.00', '10.00', '100.00', '1000.00']]]
             # [[[1.0, 0.8, 0.6, 0.4], (10**np.linspace(1, 3, num=5)).astype(int)]]
             # [[[np.arange(0, 3001, 500), (10**np.linspace(1.5, 2.5, num=5)).astype(int)],
             #  [[0.25, 0.5, 1.0, 2.0, 4.0], [0.25, 0.5, 1.0, 2.0, 4.0]]]
             # [[[0.5, 1.0, 2.0, 4.0], [1.0, 2.0], [0.5, 1.0, 2.0, 4.0], [10, 100]],
             #  [[0.25, 0.5, 1.0, 2.0, 4.0], [0.25, 0.5, 1.0, 2.0, 4.0]],
             #  [[0.25, 0.5, 1.0, 2.0, 4.0], ['-1_2', '-1_3', '0_2', '0_3']],
             #  [[0.25, 0.5, 1.0, 2.0, 4.0], ['-2', '-2_-1', '-3_-1', '-3_-2']],
             #  [np.arange(0, 3001, 200), [10, 50, 100]]]
plot_kwargs = [{'xlabel' : r'$f_{HI}$',
                'ylabel' : r'$f_{FUV}$',
                'output_file' : 'r400-fHI-fFUV-constFUV1.8'}]
               # {'xlabel' : r'$f_{HI}$',
               #  'ylabel' : r'$f_{FUV}$',
               #  'output_file' : 'r400-fHI-fFUV'}]
               # {'xlabel' : r'$R_{CMZ}$',
               #  'ylabel' : r'$f_{FUV}$',
               #  'output_file' : 'r500-Rcmz-fFUV'},
               # {'xlabel' : r'$f_{M_{cl}}$',
               #  'ylabel' : r'$f_{M_{icl}}$',
               #  'output_file' : 'r500-fcl-ficl'}]
               # {'i_x' : 0,
               #  'i_y' : 2,
               #  'log' : False,
               #  'normalise' : False,
               #  'figsize' : (10, 10),
               #  'xlabel' : r'$f_{clump}$',
               #  'ylabel' : r'$f_{\rho}$',
               #  'supxlabel' : r'$f_{FUV}$',
               #  'supylabel' : r'$f_{interclump}$',
               #  'clabel_xa' : 0.96,
               #  'pad_left' : 0.03,
               #  'pad_right' : 0.05,
               #  'pad_bottom' : 0.05,
               #  'pad_top' : 0.1,
               #  'wspace' : 0.5,
               #  'hspace' : 0.5,
               #  'fraction' : 0.031,
               #  'output_file' : 'Mcm-Micm-rho-FUV'},
               # {'xlabel' : r'$f_{clump}$',
               #  'ylabel' : r'$f_{interclump}$',
               #  'output_file' : 'Mcm-Micm'},
               # {'xlabel' : r'$f_{clump}$',
               #  'ylabel' : r'$m_{clump}$',
               #  'output_file' : 'Mcm-mcl'},
               # {'xlabel' : r'$f_{interclump}$',
               #  'ylabel' : r'$m_{interclump}$',
               #  'output_file' : 'Micm-micl_lat0'},
               # {'xlabel' : r'$R_{CMZ}$',
               #  'ylabel' : r'$f_{FUV}$',
               #  'output_file' : 'Rcmz-FUV'}]

# Path to observational data survey folders
path = '/mnt/hpc_backup/yanitski/projects/pdr/observational_data/MilkyWay/'

# target header and kernel for regridding
target_header = {
    'NAXIS': 3,
    'NAXIS1': 721,
    'NAXIS2': 361,
    'NAXIS3': 701,
    'CTYPE1': 'GLON-CAR',
    'CRVAL1': 0,
    'CDELT1': 0.5,
    'CRPIX1': 361,
    'CTYPE2': 'GLAT-CAR',
    'CRVAL2': 0,
    'CDELT2': 0.5,
    'CRPIX2': 181,
    'CTYPE3': 'VELO',
    'CRVAL3': 0,
    'CDELT3': 1,
    'CRPIX3': 351,
}
kernel_sigma = 0.25
target_kernel = (
    'gauss1d',
    (kernel_sigma,),
    kernel_sigma*3,
    kernel_sigma/2,
)

# List of survey missions to work on (None select all available surveys)
missions = None#'SEDIGISM'#['COBE-FIRAS', 'COGAL', 'Mopra', 'ThrUMMS', 'SEDIGISM']

# Colour map for the plots
cmap = 'magma'

# The type of comparison for the spectroscopic surveys
comp_type = 'pv'

# Specify whether the logged intensity is compared
log_comp = False

# Dictionary with the default kwargs for the comparison functions
default_comp_kwargs = {'lat' : 0,
                       'PRINT' : True,
                       'PLOT' : False,
                       'debug' : False}
default_plot_kwargs = {'normalise' : False,
                       'log' : False,
                       'likelihood' : False,
                       'levels' : 100,
                       'clabel_xa' : 0.95,
                       'fraction' : 0.045,
                       'pad' : 1,
                       'pad_left' : 0.1,
                       'pad_right' : 0.1,
                       'pad_bottom' : 0.05,
                       'pad_top' : 0.05,
                       'fontsize' : 24,
                       'save_plot' : True,
                       'output_format' : 'png',
                       'transparent' : True,
                       'debug' : False,
                       'COBE-FIRAS_files': ['FIRAS_LINE_EMISSION_MAP_HIGH_regridded',
                                            'FIRAS_LINE_EMISSION_MAP_HRES_regridded'],
                       'Planck_files' : ['planck_dust-commander_regridded']}
                                         # 'planck_Dust-GNILC-F353_regridded',
                                         # 'planck_Dust-GNILC-F545_regridded',
                                         # 'planck_Dust-GNILC-F857_regridded']}


for folder,params,plot_kwargs_adj in zip(folders, parameters, plot_kwargs):

    # Uncomment this to avoid the grids with more than two parameters
    # print(len(params))
    # if len(params) > 2:
    #     continue

    print('\nComparing {}'.format(folder))

    # kwargs that are common between the comparison functions
    kwargs_common = {'missions' : missions,
                     'model_param' : params,
                     'cmap' : cmap,
                     'comp_type' : comp_type,
                     'log_comp' : log_comp}

    # Compare models to observations
    kwargs = deepcopy(kwargs_common)
    kwargs.update({'model_dir':folder})
    kwargs.update(default_comp_kwargs)
    comp.model_selection(**kwargs)

    print('\n\nPlotting {}'.format(folder))

    # Plot the results
    kwargs = deepcopy(kwargs_common)
    kwargs.update({'file_format':folder.replace('/', '')})
    kwargs.update(default_plot_kwargs)
    kwargs.update(plot_kwargs_adj)
    comp.plot_comparison(**kwargs)

print('\nComparison finished...\n')
