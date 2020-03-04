import numpy as np
from astropy.io import fits

from ./cygrid import *

import constants

Cdelt = 0.1
kernelsize_sigma = 0.2
kernel_type = 'gauss1d'
kernel_params = (kernelsize_sigma,)
kernel_support = 3*kernelsize_sigma
hpx_maxres = kernelsize_sigma/2

header = {
  'NAXIS': 3,
  'NAXIS1': nPositions,
  'NAXIS2': nPositions,
  'NAXIS3': nWavelengths,
  'CTYPE1': 'GLON',
  'CTYPE2': 'GLAT',
  'CTYPE3': 'INTENSITY',
  'CUNIT1': 'DEG',
  'CUNIT2': 'DEG',
  'CUNIT3': 'K',
  'CRVAL1': mapcenter[0],
  'CRVAL2': mapcenter[1],
  'VELOCITY': velocity,
  }

r2d = 180./np.pi
d2r = np.pi/180.

def convertMap():
  
  # Retreive HDUs
  hdul = fits.open(constants.HISTORYPATH+constants.directory+directory+'integrated_intensity.fits')
  
  # List of grids
  grids = []

  for i in range(0, len(hdul), 4):
    positions = []
    intensityMap = []
    
    # Convert the Cartesian map to galactic coordinates for the input to cygrid
    tempPositions = hdul[i+0].data
    intensityMap.append(hdul[i+1].data)
    galPosition = [r2d*np.arcsin(tempPositions[:,0]/constants.rGalEarth),
                   r2d*np.arctan2(tempPositions[:,1], np.sqrt(constants.rGalEarth**2-tempPositions[:,0]))]
    galPosition = np.array(galPosition).T
    positions.append(galPosition)
    
    # Convert the Cartesian map to galactic coordinates for the input to cygrid
    tempPositions = hdul[i+2].data
    intensityMap.append(hdul[i+3].data)
    galPosition = [r2d*np.arcsin(tempPositions[:,0]/constants.rGalEarth),
                   r2d*np.arctan2(tempPositions[:,1], np.sqrt(constants.rGalEarth**2-tempPositions[:,0]))]
    galPosition = np.array(galPosition).T
    positions.append(galPosition)

    cg = convert2cygrid(positions, intensityMap, velocity=hdul[i].header['VELOCITY'])

    gridHDUL.append(cg)

  grids = fits.ImageHDU()