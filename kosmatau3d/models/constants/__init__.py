'''
A subpackage to contain all of the model constants used in the calculations
involved in :code:`kosmatau3d`.

This subpackage contains all of the constants and parameters used throughout
the program.
There are definitions to change the model parameters for when this is needed.

.. note::

   The dust wavelengths are hard-coded since these appear in Weingartner &
   Draine (2001).
   This will need to be updated if the dust calculations in KOSMA-:math:`\\tau`
   change.
'''


import inspect
import os

from numba import jit_module
import numpy as np
from numba.core.errors import NumbaDeprecationWarning, NumbaPendingDeprecationWarning
import warnings

from .change import *


warnings.simplefilter('ignore', category=NumbaDeprecationWarning)
warnings.simplefilter('ignore', category=NumbaPendingDeprecationWarning)


# Directory information
directory = ''
history = ''
filename = inspect.getframeinfo(inspect.currentframe()).filename
KOSMAPATH = os.path.abspath(os.path.dirname(filename)+'/../../')  # directory to kosmatau package
INPUTPATH = KOSMAPATH + '/input/'  # for model input properties
GRIDPATH = KOSMAPATH + '/grid/'  # for KOSMA-tau grid intensities, optical depths, and absorption
MOLECULARPATH = KOSMAPATH + '/molecular_data/'  # for KOSMA-tau ONION molecular data
HISTORYPATH = '~/KOSMA-tau^3/history/'  # history path

# Grid files
tb_grid_file = 'Tmb_LineCenter.dat'
tau_grid_file = 'tau_LineCenter.dat'
tau_fuv_grid_file = 'RhoMassAfuv.dat'

# Model input file information
h2_mass_file = ''
hi_mass_file = ''
h2_surface_density_file = ''
hi_surface_density_file = ''
h2_scale_height_file = ''
hi_scale_height_file = ''
h_number_density_file = ''
fuv_file = ''
r_gc = 4400
fuv_scale_gc = 1.0
mhi_scale_gc = 1.0
mh2_scale_gc = 1.0
disp_gmc = 1.1*1000**0.38
disp_gc = None
disp_r_gc = 4400
average_fuv = False
l_range = (912, 2066)  # for integration
velocity_file = ''
like_clumps = False
all_full = False

# Interpolation style (it accepts 'linear' or 'cubic'/'radial')
interpolation = 'linear'

# Data type for the `ensemble` module calculations
dtype = np.float64

# Factors (and constant) for fine-tuning the input data
hi_mass_factor = 1  # adjust the neutral hydrogen mass in the model
h2_mass_factor = 1  # adjust the molecular hydrogen mass in the model
interclump_hi_ratio = 1
interclump_wnm_ratio = 0.2  # adjust to change how much of the interclump medium is in the WNM
ensemble_mass_factor = [1, 1]  # adjust the mass factor for each ensemble
density_factor = 1
fuv_factor = 1

# CMZ cosmic ray ionisation rate (zeta_H2)
r_cmz = 0
zeta_cmz = 1e-14
zeta_sol = 2e-16
fuv_ism = 1

# Grid boundaries
density_limits = [3, 7]
mass_limits = [-3, 3]
uv_limits = [0, 6]

# UV adjustment
u_draine0 = 8.93717e-14  # erg/cm^-3 #2.89433*10**39 # erg/pc^3 #
# normUV = 1  # used to normalise the FUV distribution used
# globalUV = 10  # no longer used..

# Standard constants
mass_h = 1.007276466812*1.6605*10**-27  # in [kg]
mass_solar = 1.98852*10**30  # in [kg]
pc = 3.08567758149*10**16  # in [m]
pc_to_cm = 3.08567758149*10**18 #in [cm]

h = 6.62606957*10**-34  # in [J s]
c = 2.99792458*10**8  # in [m/s]
kB = 1.3806488*10**-23  # in [J/K]

# For observing the Milky Way from Earth (for use in a 3D model)
from_earth = True
rgal_earth = 8178  # from Abuter and the GRAVITY collaboration (2019)
rgal = 18000
scale_height_h2 = 70
scale_height_hi = 86
hd = 1000

# Integrated map properties
map_shape = (150, 100)  # shape in pixels of the desired map
map_size = (360, 180)  # size in degrees of the desired map
map_center = (0, 0)  # center of the desired map in degrees

# Model species
transitions = []
transition_number = 0
interclump_transitions = []
interclump_transition_number = 0
dust = 'all'
abundances = ['ELECTR', 'H', 'H2', 'H+']

# Model characteristics
n_sigma = 3
velocity_resolution = 3
voxel_size = 1000
velocity_number = 1000
velocity_bin = [-300, 300]
velocity_range = np.linspace(velocity_bin[0], velocity_bin[-1], num=velocity_number)
velocity_step = velocity_range[1] - velocity_range[0]

# Clump characteristics
clump_dispersion = 1.67/2.3548  # conversion from Gaussian FWHM to velocity dispersion (2.3548 = 2*sqrt(2*ln2))

ensembles = 2
clump_mass_number = [3, 1]
clump_log_mass_range = [[0, 2], [-2]]
clump_log_mass = [np.resize(np.linspace(clump_log_mass_range[i][0],
                                        clump_log_mass_range[i][-1],
                                        num=clump_mass_number[i]),
                          (1, clump_mass_number[i])) for i in range(len(clump_mass_number))]
# clumpLogMass[0].resize(1,clump_mass_number[0])
clump_max_indeces = [0, 0]
clump_n_max = [1, 100]
clump_fillingfactor = None
clump_density = [None, 1911]  # use None to take the voxel density
clump_log_fuv = None

interclump_idx = (False, True)
interclump_wnm_idx = (False, False)
interclump_hifactor = 1
interclump_fillingfactor = None
interclump_density = 1911  # as defined in the old version of this code
interclump_log_fuv = None
interclump_wnm_log_fuv = None
interclump_f_fuv_wnm = 1e4
interclump_mass_number = 1
interclump_log_mass_range = [-2]
interclump_log_mass = np.linspace(interclump_log_mass_range[0], interclump_log_mass_range[-1], num=interclump_mass_number)
interclump_log_mass.resize(1, interclump_mass_number)
interclump_max_indeces = 0
interclump_n_max = 100

ensemble_dispersion = 10./2.3548

# Initial mass function parameters (Heithausen et al. (1998))
alpha = 1.84
gamma = 2.31

# Statistics
scipy_probability = True
gauss = False
probability = 'binomial'
n_gauss = 1000
pn_gauss = 5

# Interpolation
log_encoded = True

# Finally, these are the wavelengths at which the dust emission will be computed
dust_names = np.array(['3.1mm', '2.4mm', '1.8mm', '1.3mm', '1.0mm', '850um', '700um', '550um', '420um', '300um',
                       '240um', '188.4um', '177.8um', '167.9um', '158.5um', '149.6um', '141.3um', '133.4um', '125.9um',
                       '118.9um', '112.2um', '105.9um', '100.0um', '94.41um', '89.13um', '84.14um', '79.43um', '74.99um',
                       '70.79um', '66.83um', '63.1um', '59.57um', '56.23um', '53.09um', '50.12um', '47.32um', '44.67um',
                       '44.28um', '42.75um', '41.33um', '39.99um', '38.74um', '37.57um', '36.47um', '35.42um', '34.44um',
                       '33.51um', '32.63um', '31.79um', '31.um', '30.24um', '29.52um', '28.83um', '28.18um', '27.55um',
                       '26.95um', '26.38um', '25.83um', '25.3um', '24.8um', '24.31um', '23.84um', '23.39um', '22.96um',
                       '22.54um', '22.14um', '21.75um', '21.38um', '21.01um', '20.66um', '20.32um', '20.um', '19.68um',
                       '19.37um', '19.07um', '18.79um', '18.5um', '18.23um', '17.97um', '17.71um', '17.46um', '17.22um',
                       '16.98um', '16.75um', '16.53um', '16.31um', '16.1um', '15.89um', '15.69um', '15.5um', '15.31um',
                       '15.12um', '14.94um', '14.42um', '13.93um', '13.78um', '13.62um', '13.48um', '13.33um', '13.19um',
                       '13.05um', '12.91um', '12.78um', '12.65um', '12.52um', '12.4um', '12.28um', '12.15um', '12.04um',
                       '11.92um', '11.81um', '11.7um', '11.59um', '11.48um', '11.37um', '11.27um', '11.17um', '11.07um',
                       '10.97um', '10.88um', '10.78um', '10.69um', '10.6um', '10.51um', '10.42um', '10.33um', '10.16um',
                       '10.08um', '9.998um', '9.918um', '9.84um', '9.762um', '9.686um', '9.611um', '9.537um', '9.464um',
                       '9.392um', '9.322um', '9.252um', '9.184um', '9.116um', '9.05um', '8.984um', '8.919um', '8.856um',
                       '8.793um', '8.731um', '8.67um', '8.61um', '8.55um', '8.492um', '8.434um', '8.377um', '8.321um',
                       '8.265um', '8.211um', '8.157um', '8.103um', '8.051um', '7.999um', '7.947um', '7.897um', '7.847um',
                       '7.797um', '7.749um', '7.701um', '7.653um', '7.606um', '7.56um', '7.514um', '7.424um', '7.38um',
                       '7.336um', '7.293um', '7.25um', '7.208um', '7.166um', '7.125um', '7.085um', '7.044um', '7.005um',
                       '6.965um', '6.926um', '6.888um', '6.85um', '6.812um', '6.775um', '6.738um', '6.702um', '6.595um',
                       '6.491um', '6.391um', '6.293um', '6.199um', '6.107um', '5.989um', '5.904um', '5.793um', '5.61um',
                       '5.39um', '5.209um', '4.999um', '4.805um', '4.592um', '4.396um', '4.203um', '3.999um', '3.936um',
                       '3.874um', '3.815um', '3.757um', '3.701um', '3.647um', '3.594um', '3.542um', '3.492um', '3.444um',
                       '3.397um', '3.342um', '3.324um', '3.306um', '3.297um', '3.29um', '3.28um', '3.246um', '3.204um',
                       '3.099um', '3.002um', '2.8um', '2.661um', '2.512um', '2.371um', '2.239um', '2.113um', '1.995um',
                       '1.884um', '1.778um', '1.679um', '1.585um', '1.496um', '1.413um', '1.334um', '1.259um', '1.189um',
                       '1.122um', '1.059um', '1.um', '944.1nm', '891.3nm', '841.4nm', '794.3nm', '749.9nm', '707.9nm',
                       '668.3nm', '631nm', '595.7nm', '562.3nm', '530.9nm', '501.2nm', '473.2nm', '446.7nm', '421.7nm',
                       '398.1nm', '375.8nm', '354.8nm', '335nm', '316.2nm', '298.5nm', '281.8nm', '266.1nm', '251.2nm',
                       '237.1nm', '223.9nm', '211.3nm', '199.5nm', '188.4nm', '177.8nm', '167.9nm', '158.5nm', '149.6nm',
                       '141.3nm', '133.4nm', '125.9nm', '121.6nm', '118.9nm', '112.2nm', '105.9nm', '100nm', '94.41nm',
                       '89.13nm', '84.14nm', '79.43nm', '74.99nm', '70.79nm', '66.83nm', '63.1nm', '59.57nm', '58.4nm',
                       '56.23nm', '53.09nm', '50.12nm', '47.32nm', '44.67nm', '42.17nm', '39.81nm', '37.58nm', '35.48nm',
                       '33.5nm', '31.62nm', '30.4nm', '29.85nm', '28.18nm', '26.61nm', '25.12nm', '23.71nm', '22.39nm',
                       '19.95nm', '17.78nm', '15.85nm', '14.13nm', '12.59nm', '11.22nm', '10.6nm', '10nm', '9nm', '8nm',
                       '7nm', '6nm', '5nm', '4nm', '3nm', '2nm', '1nm'])
wavelengths = np.array([3100., 2400., 1800., 1300., 1000., 850., 700., 550., 420., 300., 240., 188.4, 177.8, 167.9,
                        158.5, 149.6, 141.3, 133.4, 125.9, 118.9, 112.2, 105.9, 100., 94.41, 89.13, 84.14, 79.43, 74.99,
                        70.79, 66.83, 63.1, 59.57, 56.23, 53.09, 50.12, 47.32, 44.67, 44.28, 42.75, 41.33, 39.99, 38.74,
                        37.57, 36.47, 35.42, 34.44, 33.51, 32.63, 31.79, 31., 30.24, 29.52, 28.83, 28.18, 27.55, 26.95,
                        26.38, 25.83, 25.3, 24.8, 24.31, 23.84, 23.39, 22.96, 22.54, 22.14, 21.75, 21.38, 21.01, 20.66,
                        20.32, 20., 19.68, 19.37, 19.07, 18.79, 18.5, 18.23, 17.97, 17.71, 17.46, 17.22, 16.98, 16.75,
                        16.53, 16.31, 16.1, 15.89, 15.69, 15.5, 15.31, 15.12, 14.94, 14.42, 13.93, 13.78, 13.62, 13.48,
                        13.33, 13.19, 13.05, 12.91, 12.78, 12.65, 12.52, 12.4, 12.28, 12.15, 12.04, 11.92, 11.81, 11.7,
                        11.59, 11.48, 11.37, 11.27, 11.17, 11.07, 10.97, 10.88, 10.78, 10.69, 10.6, 10.51, 10.42, 10.33,
                        10.16, 10.08, 9.998, 9.918, 9.84, 9.762, 9.686, 9.611, 9.537, 9.464, 9.392, 9.322, 9.252, 9.184,
                        9.116, 9.05, 8.984, 8.919, 8.856, 8.793, 8.731, 8.67, 8.61, 8.55, 8.492, 8.434, 8.377, 8.321,
                        8.265, 8.211, 8.157, 8.103, 8.051, 7.999, 7.947, 7.897, 7.847, 7.797, 7.749, 7.701, 7.653,
                        7.606, 7.56, 7.514, 7.424, 7.38, 7.336, 7.293, 7.25, 7.208, 7.166, 7.125, 7.085, 7.044, 7.005,
                        6.965, 6.926, 6.888, 6.85, 6.812, 6.775, 6.738, 6.702, 6.595, 6.491, 6.391, 6.293, 6.199, 6.107,
                        5.989, 5.904, 5.793, 5.61, 5.39, 5.209, 4.999, 4.805, 4.592, 4.396, 4.203, 3.999, 3.936, 3.874,
                        3.815, 3.757, 3.701, 3.647, 3.594, 3.542, 3.492, 3.444, 3.397, 3.342, 3.324, 3.306, 3.297, 3.29,
                        3.28, 3.246, 3.204, 3.099, 3.002, 2.8, 2.661, 2.512, 2.371, 2.239, 2.113, 1.995, 1.884, 1.778,
                        1.679, 1.585, 1.496, 1.413, 1.334, 1.259, 1.189, 1.122, 1.059, 1., 0.9441, 0.8913, 0.8414,
                        0.7943, 0.7499, 0.7079, 0.6683, 0.631, 0.5957, 0.5623, 0.5309, 0.5012, 0.4732, 0.4467, 0.4217,
                        0.3981, 0.3758, 0.3548, 0.335, 0.3162, 0.2985, 0.2818, 0.2661, 0.2512, 0.2371, 0.2239, 0.2113,
                        0.1995, 0.1884, 0.1778, 0.1679, 0.1585, 0.1496, 0.1413, 0.1334, 0.1259, 0.1216, 0.1189, 0.1122,
                        0.1059, 0.1, 0.09441, 0.08913, 0.08414, 0.07943, 0.07499, 0.07079, 0.06683, 0.0631, 0.05957,
                        0.0584, 0.05623, 0.05309, 0.05012, 0.04732, 0.04467, 0.04217, 0.03981, 0.03758, 0.03548, 0.0335,
                        0.03162, 0.0304, 0.02985, 0.02818, 0.02661, 0.02512, 0.02371, 0.02239, 0.01995, 0.01778,
                        0.01585, 0.01413, 0.01259, 0.01122, 0.0106, 0.01, 0.009, 0.008, 0.007, 0.006, 0.005, 0.004,
                        0.003, 0.002, 0.001]) * 10**-6
# These are the flag and indeces to limit how many dust emission lines are calculated
dust_wavelengths = ''
limit_molecular = 5*10**-5
limit_pah = 5*10**-6
n_dust = wavelengths > 0
# [DEACTIVATED] These are variables to sort the wavelengths
sorted_wavelengths = []
sorted_indeces = wavelengths.argsort()

hclambda = h*c/wavelengths  # to aid in calculations

# jit_module(nopython=False)
