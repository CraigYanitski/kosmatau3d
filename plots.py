# should at some point relpace plot_geo....

# different matplotlib plotting functions

import matplotlib.pyplot as plt
import numpy as np
# edit: Craig, 14.10.2019
import globals as gbl
#end edit

class plots:
  def __init__(self):
    #from PDR import gbl._globals
    self.ISum_map = gbl._globals['ISum_map']
    self.xy = gbl._globals['xy']
    self.specrange = gbl._globals['plots']['velocity_range']

  def plot_spec_simple(self, offsets = [0,0], nstep = 1000 ):
    """
    plot one, two or three specs, depending on offsets given. For GUI.
    """
    xr = self.specrange
    offsets = np.array(offsets)

    # array with x-values:
    xval = np.linspace(xr[0], xr[1], nstep)
    yval = np.zeros([nstep])# initialise array to contain y-values for different lines of sight
    print('Variable xy: ', self.xy)
    print('Variable xy size: ', self.xy.size)

    for k in range(int(self.xy.size/2)):
      print('k', k)
      print('Variable xy[k]: ', self.xy[k])
      if np.alltrue(np.equal(offsets, self.xy[k])):
        print('Variable k_found: ', k)
        for n in range(nstep):
          yval[n] = self.ISum_map[k](xval[n])

    #from PDR import gbl._globals
    gbl._globals['xval'] = xval
    gbl._globals['yval'] = yval
