#!/home/yanitski/anaconda3/lib/python3.8

import numpy as np

from kosmatau3d import comparison as comp


# Folders and parameters of models
folders = ['/r400_cm{}-{}_d{}_uv{:.0f}/',
           '/r400_fcm{}_ficm{}/',
           '/r400_fcm{}_cm{}/',
           '/r400_ficm{}_icm{}/',
           '/r400_rcmz{:.0f}_uv{}/']
parameters = [[[0.25, 0.5, 1.0, 2.0, 4.0], [1.0, 2.0], [0.25, 0.5, 1.0, 2.0, 4.0], [10, 100]],
              [[0.25, 0.5, 1.0, 2.0, 4.0], [0.25, 0.5, 1.0, 2.0, 4.0]],
              [[0.25, 0.5, 1.0, 2.0, 4.0], ['0_2', '0_3', '1_2', '1_3']],
              [[0.25, 0.5, 1.0, 2.0, 4.0], ['-2', '-2_-1', '-3_-1', '-3_-2']],
              [np.arange(0, 3001, 200), [10, 50, 100]]]

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
mission = None#['COGAL', 'SEDIGISM', 'Mopra', 'ThrUMMS']


for folder,params in zip(folders, parameters):
    comp.model_selection(
        mission=mission,
        lat=None,
        PLOT=False,
        debug=False,
        cmap='viridis',
        model_dir=folder,
        model_param=params,
        PRINT=True,
        )

    comp.plot_comparison(
        file_format=folder.replace('/', ''),
        missions=mission,
        model_param=params,
        xlabel=r'$R_{CMZ}$',
        ylabel=r'$f_{FUV}$',
        clabel=r'$log_{10}(\mathcal{L})$'
    )
