import numpy as np
from astropy.io import fits
from astropy.wcs import WCS

from .cygrid import *

import constants

v_min = 0
v_max = 0
header = {
  'BUNIT': 'K km/s',
  'NAXIS': 3,
  'NAXIS1': constants.mapShape[0],
  'NAXIS2': constants.mapShape[1],
  'NAXIS3': constants.sortedWavelengths.size,
  'CTYPE1': 'GLON-MOL',
  'CUNIT1': 'deg',
  'CDELT1': constants.mapSize[0]/constants.mapShape[0],
  'CRPIX1': constants.mapShape[0]/2.,
  'CRVAL1': constants.mapCenter[0],
  'CTYPE2': 'GLAT-MOL',
  'CUNIT2': 'deg',
  'CDELT2': constants.mapSize[1]/constants.mapShape[1],
  'CRPIX2': constants.mapShape[1]/2.,
  'CRVAL2': constants.mapCenter[1],
  # 'CTYPE3': 'FREQ',
  # 'CUNIT3': 'N/A',
  # 'CRPIX3': 'N/A',
  # 'CRVAL3': 'N/A',
  'VELOCITY': 0,
  }

r2d = 180./np.pi
d2r = np.pi/180.