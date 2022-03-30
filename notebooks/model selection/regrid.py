import numpy as np
import matplotlib.pyplot as plt

from kosmatau3d import comparison as comp

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
kernel_sigma = 1
target_kernel = (
    'gauss1d',
    (kernel_sigma,),
    kernel_sigma*3,
    kernel_sigma/2,
)

# List of survey missions to work on
mission = ['COBE-FIRAS']#['COGAL', 'SEDIGISM', 'Mopra', 'ThrUMMS']

comp.regrid_observations(
    path='/mnt/hpc_backup/yanitski/projects/pdr/observational_data/MilkyWay/',
    mission=['Planck'],
    target_header=target_header,
    target_kernel=target_kernel,
)

