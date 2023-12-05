'''
A subpackage to load all of the KOSMA-:math:`\\tau` and input data required
to run :code:`kosmatau3d`.
This subpackage will contain the input data needed to properly simulate the PDR. 
All of the information specific to the object being simulated should be in 
their own folder in INPUTPATH. 
The KOSMA-tau grid data is located in the folder 'grid'.
'''


import numpy as np

from kosmatau3d.models import constants
from .methods import *


# Flags
grid_initialised = False
model_initialised = False

# Grid
clump_tb_centerline = None
interclump_tb_centerline = None
interclump_dust_tb_centerline = None
clump_tau_centerline = None
interclump_tau_centerline = None
interclump_dust_tau_centerline = None
clump_taufuv = None
interclump_taufuv = None
clump_column_density = None
interclump_column_density = None
clump_temperature = None
interclump_temperature = None
species_data = None
e_tilde_real = None
e_tilde_imaginary = None

# Model
h2_surface_mass_profile = None
hi_surface_mass_profile = None
h2_scale_height_profile = None
hi_scale_height_profile = None
number_density_profile = None
fuv_profile = None
galaxy_rotation_profile = None
