import numpy as np

from kosmatau3d.models import constants
from .methods import *


'''
This module will contain the input data needed to properly simulate the PDR. All of the information specific
to the object being simulated should be in their own folder in INPUTPATH. The KOSMA-tau grid data is located
in the folder 'grid'.
'''


# Flags
grid_initialised = False
model_initialised = False

# Grid
tb_centerline = None
tau_centerline = None
rho_mass_taufuv = None
column_density = None
temperature = None
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
