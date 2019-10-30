
# create channel map from spectra and convolve with beam
from scipy.integrate import quad
from scipy import signal
import numpy as np
import copy
#from PDR import gbl._globals
# edit: Craig, 14.10.2019
import globals as gbl
#end edit

class maps:
    def __init__(self):
      #from PDR import gbl._globals
      self.xy = np.array(gbl._globals['xy'], float) # x-y positions
      self.ISum_map = gbl._globals['ISum_map'] 
      # Spectra at x-y positions (only y-values)
      self.distance = gbl._globals['compound']['distance']
      self.dint = gbl._globals['compound']['dint']
      self.pc = gbl._globals['constants']['pc']

    def gauss_kern(self, beamsize):
      """ Returns a normalized 2D gauss kernel array for convolutions. 
      FWHM is defined by beamsize, total size of matrix covers at 
      least 3 sigma of gaussian beam """
      #from PDR import gbl._globals
      # beam approximated by gauss. calculate FWHM from beamsize
      # (=Position where Power is halved)
      # beamsize in arcsec
      FWHM = beamsize/3600. * np.pi/180. * self.distance/\
             (float(self.dint) / float(self.pc)) 
      # FWHM measured in pixels
      print('Variable FWHM: ', FWHM)
      # pause = input('FWHM...')     
      sigma = FWHM/(2*np.sqrt(2*np.log(2))) 
      # standard deviation, measured in pixels

      threeSigma = np.ceil(gbl._globals['compound']['g_sigma'] * sigma) 
      # three sigma, rounded up, measured in pixels
      print('Variable 3sigma: ', threeSigma)
      # pause = input('3 sigma ...')   

      with open(gbl.KOSMAPATH+'maps/n-sigma.dat',"a") as sig: # open and append
        sig.write('3sigma' + str(threeSigma) + '\n')

      gbl._globals['threeSigma'] = threeSigma

      x, y = np.mgrid[-threeSigma:threeSigma + 1, -threeSigma:threeSigma + 1]
      # print 'x',x
      # print 'y',y

      g = np.exp(-(x**2/float(2*sigma**2) + y**2/float(2*sigma**2)))
      # print 'g', g
      # pause = input('gauss kernel ...')  
      return g / g.sum()

    def blur_image(self, im, dbeam) :
      """
        blurs the image by convolving with a gaussian kernel of typical
        size n. The optional keyword argument ny allows for a different
        size in the y direction.
      """
      g = self.gauss_kern(dbeam)
      improc = signal.convolve2d(im, g, mode = 'same') 

      # mode. 'full', 'same' or valid
      # same: The output is the same size as im, centered with respect 
      # to the 'full' output.
      # zero-padding is done by default

      # 'valid' does not work for 1D cut. use 'same'

      return(improc)   

    def create_map(self): 
      # create map without convolution
      #from PDR import gbl._globals
      
      no_los = 0 # count finished pixels
      maxx = max(self.xy[:,0])
      minx = min(self.xy[:,0])
      maxy = max(self.xy[:,1])
      miny = min(self.xy[:,1])
      gbl._globals['map_limits'] = {}
      gbl._globals['map_limits']['maxx'] = maxx
      gbl._globals['map_limits']['minx'] = minx
      gbl._globals['map_limits']['maxy'] = maxy
      gbl._globals['map_limits']['miny'] = miny

      # xvalues = np.linspace(vel - d_vel, vel + d_vel, nstep)
      xvalues = np.linspace(gbl._globals['compound']['vel'] - gbl._globals['compound']['d_vel'], \
                            gbl._globals['compound']['vel'] + gbl._globals['compound']['d_vel'], \
                            gbl._globals['compound']['nstep'])

      # print "maxx", maxx
      # print "minx", minx
      # print "maxy", maxy
      # print "miny", miny        
      # pause = input("max/min in maps") 

      mapPixels = gbl._globals['compound']['mapPixels']
      # number of pixels in map
      map_int = np.zeros([int(abs( maxx - minx ) + 1 ), \
                int(abs( maxy - miny ) + 1 )]) 
      # initialize map for integrated intensities
      map_int[:] = -1
      for x in np.array(np.linspace(minx, maxx, maxx - minx + 1), int):
        # loop over pixel, x-direction
        for y in np.array(np.linspace(miny, maxy, maxy - miny + 1), int):
          # loop over pixel, y-direction
          no_los = no_los + 1
          for k in range(int(self.xy.size/2)):
            # loop over number of pixels
            if self.xy[k][0] == x and self.xy[k][1] == y:
              # print 'x', x
              # print 'y', y
              """
              # integrate function (can take a lot of time)
              map_int[x - minx, y - miny] = (quad(self.ISum_map[k], \
    #              vel - d_vel, vel + d_vel, epsabs = 0.01, \
      #                                      epsrel = 10**(-4))[0]) 
              """
              # sample function and integrate
              print('\nIntegrating . . . ')
              yvalues = self.ISum_map[k]
              # for element in xvalues:
              #     yvalues.append(self.ISum_map[k](element))
              # yvalues = np.array(yvalues)
              print('Variable yvalues: ', yvalues)
              # pause = input('...')    
              map_int[x - minx, y - miny] = (np.trapz(yvalues, xvalues))
              print('Integrating . . . done :-)\n')

                # integrate spectra within the interval 
                # [vel - d_vel, vel + d_vel]
                # originial coordinates are shifted so that [xmin,ymin]
                # is the [0,0] in the output array

            print('Total pixels in map: ', mapPixels)
            print(no_los, ' pixels finished')
                      
      # carful: smallest x and y values have to go to the lower left corner 
      # of the array. otherwise map will me mirrored
      map_int = map_int.transpose()[::-1] 
      gbl._globals['map_int'] = map_int 
      # export map for plotting

    def create_convolved_map(self, dbeam): # convolve map
      #from PDR import gbl._globals
      print('dbeam: ', dbeam)
      map_int1 = copy.copy(gbl._globals['map_int'])
      # non-convolved map
      map_int2 = copy.copy(gbl._globals['map_int'])
      # copy of map needed to provide correct normalisation at boundaries
      for x in np.arange(map_int1.shape[0]):
        for y in np.arange(map_int1.shape[1]):
          if map_int1[x, y] == -1: map_int1[x, y] = 0 
          if map_int2[x, y] == -1: map_int2[x, y] = 0 
          else: map_int2[x, y] = 1 
      map_conv = self.blur_image(map_int1, dbeam)
      map_norm = self.blur_image(map_int2, dbeam)
      gbl._globals['map_conv'] = map_conv / map_norm
      # export COMPLETE map for plotting via GUI    
