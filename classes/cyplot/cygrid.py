import matplotlib.pyplot as plt
import numpy as np
import cygrid as cg
from astropy.io import fits
from astropy.wcs import WCS

import constants

import cyplot

def convertMap(directory='', verbose=False):
  
  # Retreive HDUs
  hdul = fits.open(constants.HISTORYPATH+constants.directory+directory+'integrated_intensity.fits')

  # List of image grids
  grids = []

  # initialise the observing velocity
  velocity = None

# Cycle through observing velocities
  for i in range(0, len(hdul), 2):

    # Skip if there is no data in the HDU
    if not np.size(hdul[i].data): continue

    # Initialise the intensity map that will be gridded
    positions = []
    intensityMap = []

    # Move to the next HDU if this one has already been processed
    if hdul[i].header['VELOCITY']==velocity: continue

    if verbose: print(hdul[i].header['VELOCITY'])

    # Specify the observing velocity
    velocity = hdul[i+0].header['VELOCITY']

    # Save this velocity if it is at the end of the range
    if velocity<cyplot.v_min or i==0:
      cyplot.v_min = velocity
    elif velocity>cyplot.v_max or i==len(hdul)-1:
      cyplot.v_max = velocity
    
    # Check if the data to be added contains the galactic center
    if hdul[i].header['DIREC']=='Inner disk':
      positions,intensityMap = convertData(hdul[i:i+2], direction='inner')

      # Check if the following data has the same observing velocity
      #a try-except routine is used in event i is pointing to the last HDU
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
  wcs = WCS(cyplot.header)

  return grids,wcs

def convertData(hdul, direction='inner'):
  
  # Extract emission data
  positions = hdul[0].data
  intensity = hdul[1].data
  
  if direction=='inner':
    # Convert Cartesian positions to spherical
    galPosition = [cyplot.r2d*np.arcsin(positions[:,0]/constants.rGalEarth),
                   cyplot.r2d*np.arctan2(positions[:,1], (np.sqrt(constants.rGal**2-positions[:,0]**2)+np.sqrt(constants.rGalEarth**2-positions[:,0]**2)))]
    galPosition = np.array(galPosition).T

  if direction=='outer':
    # Convert Cartesian positions to spherical
    galPosition = [cyplot.r2d*np.arcsin(positions[:,0]/constants.rGalEarth),
                   cyplot.r2d*np.arctan2(positions[:,1], (np.sqrt(constants.rGal**2-positions[:,0]**2)-np.sqrt(constants.rGalEarth**2-positions[:,0]**2)))]
    galPosition = np.array(galPosition).T

    # Correct position for observing the edge of the disk
    galPosition[galPosition[:,0]>=0,0] = 90 + galPosition[galPosition[:,0]>=0,0][::-1]
    galPosition[galPosition[:,0]<0,0]  = -90 + galPosition[galPosition[:,0]<0,0][::-1]
  
  return galPosition,intensity[:,0,:]

def convert2cygrid(positions, intensityMap, velocity=0):

  # Define the kernel properties
  kernelsize_sigma = 1
  kernel_type = 'gauss1d'
  kernel_params = (kernelsize_sigma,)
  kernel_support = 3*kernelsize_sigma
  hpx_maxres = kernelsize_sigma/2
  
  # Change the header for the cygrid WCS grid for the observing velocity
  cyplot.header['VELOCITY'] = velocity

  # Set the kernel of the gridder
  grid = cg.WcsGrid(cyplot.header)
  grid.set_kernel(
    kernel_type,
    kernel_params,
    kernel_support,
    hpx_maxres
    )

  # grid the positions using cygrid
  grid.grid(positions[:,0], positions[:,1], intensityMap)

  # Extract datacube of the required image
  image = grid.get_datacube()

  return image

def integrateVelocity(grid, splitWavelengths=False):

  # Get the resolution of the velocity dimension
  v_resolution = (cyplot.v_max-cyplot.v_min) / (grid.shape[0]-1)

  # Recreate velocity dimension
  v_range = np.arange(cyplot.v_min, cyplot.v_max+v_resolution, v_resolution)

  # Initialise interpolated grid
  interpGrid = np.zeros_like(grid[0,:,:,:])

  if splitWavelengths:

    # Cycle through wavelengths in model
    for w in range(grid.shape[1]):

      # Cycle through pixel positions
      for i in range(grid.shape[2]):
        for j in range(grid.shape[3]):

          # Integrate over velocities using the trapezoidal method
          interpGrid[w,i,j] = np.trapz(grid[:,w,i,j], v_range, axis=0)

  else:

    # Cycle through pixel positions
    for i in range(grid.shape[2]):
      for j in range(grid.shape[3]):

        # Integrate over velocities using the trapezoidal method
        interpGrid[:,i,j] = np.trapz(grid[:,:,i,j], v_range, axis=0)

  return interpGrid
