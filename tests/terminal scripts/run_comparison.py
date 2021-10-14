#!/home/yanitski/anaconda3/lib/python3.8

'''
Available directories:
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
folders = ['/r400_cm{}-{}_d{}_uv{:.0f}/',
           '/r400_fcm{}_ficm{}/',
           '/r400_fcm{}_cm{}/',
           '/r400_ficm{}_icm{}/',
           '/r400_rcmz{:.0f}_uv{}/']
parameters = [[[0.5, 1.0, 2.0, 4.0], [1.0, 2.0], [0.5, 1.0, 2.0, 4.0], [10, 100]],
              [[0.25, 0.5, 1.0, 2.0, 4.0], [0.25, 0.5, 1.0, 2.0, 4.0]],
              [[0.25, 0.5, 1.0, 2.0, 4.0], ['-1_2', '-1_3', '0_2', '0_3']],
              [[0.25, 0.5, 1.0, 2.0, 4.0], ['-2', '-2_-1', '-3_-1', '-3_-2']],
              [np.arange(0, 3001, 200), [10, 50, 100]]]
plot_kwargs = [{'log' : True,
                'normalise' : False,
                'figsize' : (30, 10),
                'xlabel' : r'$f_{clump}$',
                'ylabel' : r'$f_{interclump}$',
                'supxlabel' : r'$f_{FUV}$',
                'supylabel' : r'$f_{\rho}$',
                'clabel_xa' : 0.98,
                'pad_left' : 0.03,
                'pad_right' : 0.05,
                'pad_bottom' : 0.03,
                'pad_top' : 0.05,
                'wspace' : 0.3,
                'fraction' : 0.031,
                'output_file' : 'Mcm-Micm-rho-FUV'},
               {'xlabel' : r'$f_{clump}$',
                'ylabel' : r'$f_{interclump}$',
                'output_file' : 'Mcm-Micm'},
               {'xlabel' : r'$f_{clump}$',
                'ylabel' : r'$m_{clump}$',
                'output_file' : 'Mcm-mcl'},
               {'xlabel' : r'$f_{interclump}$',
                'ylabel' : r'$m_{interclump}$',
                'output_file' : 'Micm-micl'},
               {'xlabel' : r'$R_{CMZ}$',
                'ylabel' : r'$f_{FUV}$',
                'output_file' : 'Rcmz-FUV'}]

# Path to observational data survey folders
path = '/mnt/hpc_backup/yanitski/projects/pdr/observational_data/MilkyWay/'

# target header and kernel for regridding
target_header = {
    'NAXIS': 3,
    'NAXIS1': 361,
    'NAXIS2': 181,
    'NAXIS3': 701,
    'CTYPE1': 'GLON-CAR',
    'CRVAL1': 0,
    'CDELT1': 1,
    'CRPIX1': 181,
    'CTYPE2': 'GLAT-CAR',
    'CRVAL2': 0,
    'CDELT2': 1,
    'CRPIX2': 91,
    'CTYPE3': 'VELO',
    'CRVAL3': 0,
    'CDELT3': 1,
    'CRPIX3': 351,
}
kernel_sigma = 0.5
target_kernel = (
    'gauss1d',
    (kernel_sigma,),
    kernel_sigma*3,
    kernel_sigma/2,
)

# List of survey missions to work on
missions = None#['COGAL', 'SEDIGISM', 'Mopra', 'ThrUMMS']

# Colour map for the plots
cmap = 'gist_ncar'

# Dictionary with the default kwargs for the comparison functions
default_comp_kwargs = {'lat' : None,
                       'PRINT' : True,
                       'PLOT' : False,
                       'debug' : False}
default_plot_kwargs = {'normalise' : True,
                       'levels' : 100,
                       'clabel' : r'$log_{10}(\mathcal{L})$',
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
                       'transparent' : True}


for folder,params,plot_kwargs_adj in zip(folders, parameters, plot_kwargs):
    # print(len(params))
    # if len(params) > 2:
    #     continue
    print('Comparing {}'.format(folder))

    # kwargs that are common between the comparison functions
    kwargs_common = {'missions':missions, 'model_param':params, 'cmap':cmap}

    # Compare models to observations
    kwargs = deepcopy(kwargs_common)
    kwargs.update({'model_dir':folder})
    kwargs.update(default_comp_kwargs)
    # comp.model_selection(**kwargs)

    # Plot the results
    kwargs = deepcopy(kwargs_common)
    kwargs.update({'file_format':folder.replace('/', '')})
    kwargs.update(default_plot_kwargs)
    kwargs.update(plot_kwargs_adj)
    comp.plot_comparison(**kwargs)

print('\nComparison finished...\n')
