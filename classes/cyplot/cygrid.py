import matplotlib.pyplot as plt
import numpy as np
import cygrid as cg
from astropy.io import fits

import cyplot

def convert2cygrid(positions, intensityMap, velocity=0):

  nPositions = positions.shape[0]
  nWavelengths = intensityMap.shape[1]
  
  # create a header for the cygrid WCS grid
  #it has dimension (lon, lat, emission)
  cyplot.header['NAXIS1'] = nPositions
  cyplot.header['NAXIS2'] = nPositions
  cyplot.header['NAXIS3'] = nWavelengths
  cyplot.header['VELOCITY'] = velocity

  grid = cg.WcsGrid(header)
  grid.set_kernel(
    kernel_type,
    kernel_params,
    kernel_support,
    hpx_maxres
    )

  grid.grid(positions[:,0], positions[:,1], np.reshape(sig, (-1, 1)))
