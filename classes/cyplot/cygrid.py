import matplotlib.pyplot as plt
import numpy as np
import cygrid as cg
from astropy.io import fits
from astropy.wcs import WCS

import constants

import cyplot

def convertMap(directory='', input_coord='spherical', integrate=True, verbose=False):
  
  # Retreive HDUs
  intensity = fits.open(constants.HISTORYPATH+constants.directory+directory+'integrated_intensity.fits')[1]

  # Re-adjust cygrid input header
  cyplot.header['NAXIS3'] = intensity.shape[-1]

  if integrate:

    lon = np.linspace(intensity.header['CRVAL2']-intensity.header['CDELT2']*intensity.header['CRPIX2'], \
                      intensity.header['CRVAL2']+intensity.header['CDELT2']*intensity.header['CRPIX2'], \
                      intensity.header['NAXIS2'])
    lat = np.linspace(intensity.header['CRVAL3']-intensity.header['CDELT3']*intensity.header['CRPIX3'], \
                      intensity.header['CRVAL3']+intensity.header['CDELT3']*intensity.header['CRPIX3'], \
                      intensity.header['NAXIS3'])
    positions = np.meshgrid(lon, lat)
    lon = positions[0].flatten()
    lat = positions[1].flatten()
    positions = np.array([lon, lat]).T
    
    intensity = integrateVelocity(intensity)

    print(intensity.data.shape)
    
    grid = convert2cygrid(positions, intensity)

    # grids.append(cg)

    # Convert the array to a numpy array with the zeroith dimension being the velocity
    # grids = np.array(grids)

  else:

    # List of image grids
    grid = []

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
      
      if input_coord=='spherical':
        positions = hdul[i].data
        positions *= cyplot.r2d #this could not be merged with the previous line due to a contiguous error
        intensityMap = hdul[i+1].data[:,0,:]

      elif input_coord=='Cartesian':
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
      cg = convert2cygrid(positions, intensityMap, velocity=hdul[i].header['VELOCITY'], grid_type='wcs')

      # Add the image to the array
      grid.append(cg)

    # Convert the array to a numpy array with the zeroith dimension being the velocity
    grid = np.array(grid)

  return grid

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

def convert2cygrid(positions, intensityMap, velocity=None, grid_type='wcs'):

  if grid_type=='wcs':
    # Define the kernel properties
    kernelsize_sigma = 1
    kernel_type = 'gauss1d'
    kernel_params = (kernelsize_sigma,)
    kernel_support = 3*kernelsize_sigma
    hpx_maxres = kernelsize_sigma/2
    
    # Change the header for the cygrid WCS grid for the observing velocity
    if velocity!=None: cyplot.header['VELOCITY'] = velocity

    # Set the kernel of the gridder
    grid = cg.WcsGrid(cyplot.header)

  elif grid_type=='sightline':
    # Define the kernel properties
    kernelsize_sigma = 1
    kernel_type = 'gauss1d'
    kernel_params = (kernelsize_sigma,)
    kernel_support = 3*kernelsize_sigma
    hpx_maxres = kernelsize_sigma/2
    
    # Change the header for the cygrid WCS grid for the observing velocity
    if velocity!=None: cyplot.header['VELOCITY'] = velocity

    # Set the kernel of the gridder
    grid = cg.SlGrid(positions[:,0], positions[:,1])


  grid.set_kernel(
    kernel_type,
    kernel_params,
    kernel_support,
    hpx_maxres
    )

  # Flatten the positional axes in intensity map to grid with cygrid. Note that there needs to be three dimensions to the data.
  lon = np.unique(positions[:,0])
  lat = np.unique(positions[:,1])
  # print(lon, lat)
  # print(positions[:,0], positions[:,1])
  shape = intensityMap.shape
  flatMap = np.zeros((shape[0]*shape[1], shape[2]))
  for j in range(shape[0]):
    for i in range(shape[1]):
      # print(j*shape[1]+i)
      flatMap[j*shape[1]+i,:] = intensityMap[j,i,:]

  # Grid the positions using cygrid
  # positions[positions[:,0]<0,0] = positions[positions[:,0]<0,0] + np.pi
  # positions[positions[:,0]>=0,0] = positions[positions[:,0]>=0,0] - np.pi
  # print(positions, flatMap.shape)
  # print(positions.shape, intensityMap.shape, flatMap.shape)
  grid.grid(positions[:,0]*cyplot.r2d, positions[:,1]*cyplot.r2d, flatMap)

  # Extract datacube of the required image
  image = grid.get_datacube()

  print(np.where(np.isnan(image[0,:,:]))[0].size)

  return grid

def integrateVelocity(grid, splitWavelengths=False, old=False):

  # Get the resolution of the velocity dimension
  vmin = grid.header['CRVAL4'] - grid.header['CDELT4']*grid.header['CRPIX4']
  vmax = grid.header['CRVAL4'] + grid.header['CDELT4']*grid.header['CRPIX4']
  vrange = np.linspace(vmin, vmax, num=grid.header['NAXIS4'])

  # # Recreate velocity dimension
  # v_range = np.arange(cyplot.v_min, cyplot.v_max+v_resolution, v_resolution)

  # # Initialise interpolated grid
  # interpGrid = np.zeros_like(grid[0,:,:,:])

  if splitWavelengths:

    # Cycle through wavelengths in model
    for w in range(grid.shape[1]):

      # Cycle through pixel positions
      for i in range(grid.shape[2]):
        for j in range(grid.shape[3]):

          # Integrate over velocities using the trapezoidal method
          interpGrid[w,i,j] = np.trapz(grid[:,w,i,j], v_range, axis=0)

  elif old:

    # Cycle through pixel positions
    for i in range(grid.shape[2]):
      for j in range(grid.shape[3]):

        # Integrate over velocities using the trapezoidal method
        interpGrid[:,i,j] = np.trapz(grid[:,:,i,j], v_range, axis=0)

  else:

    interpGrid = np.trapz(grid.data, vrange, axis=0)

  return interpGrid
