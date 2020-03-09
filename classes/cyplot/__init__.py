import numpy as np
from astropy.io import fits
from astropy.wcs import WCS

from .cygrid import *

import constants


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

def convertMap(directory=''):
  
  # Retreive HDUs
  hdul = fits.open(constants.HISTORYPATH+constants.directory+directory+'integrated_intensity.fits')

  # List of image grids
  grids = []

  # initialise the observing velocity
  velocity = None

# Cycle through observing velocities
  for i in range(0, len(hdul), 2):

    # Initialise the intensity map that will be gridded
    positions = []
    intensityMap = []

    # Move to the next HDU if this one has already been processed
    if hdul[i].header['VELOCITY']==velocity: continue

    print(hdul[i].header['VELOCITY'])

    # Specify the observing velocity
    velocity = hdul[i+0].header['VELOCITY']
    
    # Check if the data to be added contains the galactic center
    if hdul[i].header['DIREC']=='Inner disk':
      positions,intensityMap = convertData(hdul[i:i+2], direction='inner')

      # Check if the following data has the same observing velocity
      try:
        if hdul[i+2].header['VELOCITY']==velocity and hdul[i+2].header['DIREC']=='Outer disk':
          tempPositions,tempIntensity = convertData(hdul[i+2:i+4], direction='outer')
          positions = np.append(positions, tempPositions, axis=0)
          intensityMap = np.append(intensityMap, tempIntensity, axis=0)
      except: pass

    # Otherwise add the data for an observation towards the edge
    else:
      positions,intensityMap = convertData(hdul[i+2:i+4], direction='outer')

    # Grid the model emissions using cygrid
    cg = convert2cygrid(positions, intensityMap, velocity=hdul[i].header['VELOCITY'])

    # Add the image to the array
    grids.append(cg)

  # Convert the array to a numpy array with the zeroith dimension being the velocity
  grids = np.array(grids)

  # Define the target WCS object (used for all channels) to properly display the images
  wcs = WCS(header)

  return grids,wcs

def convertData(hdul, direction='inner'):
  
  # Extract emission data
  positions = hdul[0].data
  intensity = hdul[1].data
  
  if direction=='inner':
    # Convert Cartesian positions to spherical
    galPosition = [r2d*np.arcsin(positions[:,0]/constants.rGalEarth),
                   r2d*np.arctan2(positions[:,1], (np.sqrt(constants.rGal**2-positions[:,0]**2)+np.sqrt(constants.rGalEarth**2-positions[:,0]**2)))]
    galPosition = np.array(galPosition).T

  if direction=='outer':
    # Convert Cartesian positions to spherical
    galPosition = [r2d*np.arcsin(positions[:,0]/constants.rGalEarth),
                   r2d*np.arctan2(positions[:,1], (np.sqrt(constants.rGal**2-positions[:,0]**2)-np.sqrt(constants.rGalEarth**2-positions[:,0]**2)))]
    galPosition = np.array(galPosition).T

    # Correct position for observing the edge of the disk
    galPosition[galPosition[:,0]>=0,0] = 90 + galPosition[galPosition[:,0]>=0,0][::-1]
    galPosition[galPosition[:,0]<0,0]  = -90 + galPosition[galPosition[:,0]<0,0][::-1]
  
  return galPosition,intensity