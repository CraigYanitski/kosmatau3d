import numpy as np
from astropy.io import fits
from astropy.wcs import WCS

from .cygrid import *

import constants

v_min = 0
v_max = 0
header = {
  'NAXIS': 3,
  'NAXIS1': constants.mapShape[0],
  'NAXIS2': constants.mapShape[1],
  'NAXIS3': constants.sortedWavelengths.size,
  'CTYPE1': 'GLON',
  'CTYPE2': 'GLAT',
  'CTYPE3': 'INTENSITY',
  'CUNIT1': 'deg',
  'CUNIT2': 'deg',
  'CUNIT3': 'K',
  'CDELT1': constants.mapSize[0]/constants.mapShape[0],
  'CDELT2': constants.mapSize[1]/constants.mapShape[1],
  'CRPIX1': constants.mapShape[0]/2.,
  'CRPIX2': constants.mapShape[1]/2.,
  'CRVAL1': constants.mapCenter[0],
  'CRVAL2': constants.mapCenter[1],
  'VELOCITY': 0,
  }

r2d = 180./np.pi
d2r = np.pi/180.