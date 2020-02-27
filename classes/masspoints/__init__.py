from numba import jit_module
import numpy as np

from .masspoint import *

# Properties

logFUV = 0

clumpLogDensity = []
clumpRadius = []

interclumpLogDensity = []
interclumpRadius = []

# Emission

clumpIntensity = []
clumpOpticalDepth = []

interclumpIntensity = []
interclumpOpticalDepth = []

jit_module(nopython=False)
