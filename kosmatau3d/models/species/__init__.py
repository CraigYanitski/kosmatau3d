"""
A subpackage to read and contain all of the species spectral line transition
data.
This relies on transition data used in KOSMA-:math:`\\tau`, and is mainly
used to plotting and to interpolate the dust continuum when calculating the 
voxel-averaged emissivity and absorption.
"""


#import numpy as np

from .molecules import *
#from kosmatau3d.models import observations


clump_transitions = []
clump_transition_indeces = []
clump_transition_frequencies = []
clump_transition_wavelengths = []
interclump_transitions = []
interclump_transition_indeces = []
interclump_transition_frequencies = []
interclump_transition_wavelengths = []
