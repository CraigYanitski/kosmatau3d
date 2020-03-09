import matplotlib.pyplot as plt
import numpy as np
import cygrid as cg
from astropy.io import fits

import cyplot

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
