import numpy as np
from numba import jit_module

from .combination import *

clumpCombination = []
interclumpCombination = []

clumpIntensity = []
clumpOpticalDepth = []
interclumpIntensity = []
interclumpOpticalDepth = []

jit_module(nopython=False)
