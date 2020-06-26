import os
import numpy as np
from astropy.io import fits
import matplotlib.pyplot as plt
from tqdm import tqdm
#import progressbar as pb
from numba import jit
import importlib as il
import gc

import constants
import species
from Voxel import *
import interpolations
# import species

class VoxelGrid(object):
  '''
  This is a class to handle all of the voxels in KOSMA-tau^3. It contains a
  specified arrangement of voxels, and must coordinate with the Dimensions class
  to make the Shape class functional.
  '''
  # PRIVATE

  def __init__(self, shape):

    self.__shape = shape

    self.__voxelNumber = self.__shape.voxelNumber()

    self.__voxels = []

    self.__map = {}       #dictionary object to map the voxel indeces to the correct location

    #self.__species = None

    self.__voxelIntensity = []
    self.__voxelOpticalDepth = []

    self.__voxelFUV = []
    self.__voxelAfuv = []

    self.__x = []
    self.__y = []
    self.__z = []

    constants.history = 'r{}_im{}_cm{}_d{}_uv{}/'.format(int(constants.resolution), constants.interclumpMassFactor, constants.clumpMassFactor, constants.densityFactor, constants.globalUV)

    return

  def __initialiseGrid(self):
    #self.__species = species
    for i in range(self.__voxelNumber): self.__voxels.append(Voxel(i))
    return

  def __str__(self):
    return 'VoxelGrid\n  ->{} voxels\n  ->intensity {}\n  ->optical depth {}'.format(self.__voxelNumber, sum(self.__voxelIntensity), -np.log(np.exp(-self.__voxelOpticalDepth)))

  # PUBLIC
  #def createGrid(self, indeces):
  #  for i in indeces: self.__voxels.append(Voxel(i))
  #  return
  # def reloadModules(self):
  #   il.reload(Voxel)
  #   il.reload(Interpolate)
  #   for voxel in self.__grid:
  #     voxel.reloadModules()
  #   self.__interpolations.reloadModules()
  #   return

  def getDimensions(self):
    return self.__shape.getDimensions()

  def calculateEmission(self, index=0, debug=False, verbose=False):
    #print(observations.tauCenterline)
    interpolations.initialise()
    self.__initialiseGrid()

    print('\nCalculating Grid Emission...')

    x,y,z = self.__shape.voxelCartesianPositions()
    r,phi = self.__shape.voxelPolarPositions()
    #self.__unusedVoxels = []
    with tqdm(total=len(self.__voxels), desc='Voxels initialised', miniters=1, dynamic_ncols=True) as progress:
      for i,voxel in enumerate(self.__voxels):

        if verbose:
            print('\nMax X, Radius:', max(x), r[i], '\n')

        self.__x.append(x[i])
        self.__y.append(y[i])
        self.__z.append(z[i])

        voxel.setIndex(i)#-len(self.__unusedVoxels))
        voxel.setPosition(x[i], y[i], z[i], r[i], phi[i])
        voxel.setProperties(debug=debug)

        self.__voxelAfuv.append(voxel.getAfuv())
        voxel.calculateEmission()

        progress.update()

      progress.close()

      print('\nEmission calculated successfully.')
    
    return

  def writeEmission(self, verbose=False, debug=False):
    
    print('\nStreaming to fits files...')

    if debug:

      dim = [1, self.__voxelNumber]
      shdu_mass = self.shdu_header(name='Mass', units='Msol', filename='voxel_mass', dim=dim)

      dim = [1, self.__voxelNumber]
      shdu_clump_mass = self.shdu_header(name='Clump Mass', units='Msol', filename='voxel_clump_mass', dim=dim)

      dim = [1, self.__voxelNumber]
      shdu_interclump_mass = self.shdu_header(name='Interclump Mass', units='Msol', filename='voxel_interclump_mass', dim=dim)

      dim = [1, self.__voxelNumber]
      shdu_density = self.shdu_header(name='Density', units='cm^-3', filename='voxel_density', dim=dim)

    dim = [3, self.__voxelNumber]
    shdu_position = self.shdu_header(name='Position', units='pc', filename='voxel_position', dim=dim)

    dim = [1, self.__voxelNumber]
    shdu_velocity = self.shdu_header(name='Velocity', units='km/s', filename='voxel_velocity', dim=dim)
    
    dim = [constants.clumpMaxIndeces, self.__voxelNumber]
    shdu_clump_velocity = self.shdu_header(name='Clump velocity', units='km/s', filename='voxel_clump_velocity', dim=dim)
    
    dim = [constants.interclumpMaxIndeces, self.__voxelNumber]
    shdu_interclump_velocity = self.shdu_header(name='Interclump Velocity', units='km/s', filename='voxel_interclump_velocity', dim=dim)
    
    dim = [1,self.__voxelNumber]
    shdu_FUV = self.shdu_header(name='FUV', units='Draine', filename='voxel_fuv', dim=dim)
    shdu_Afuv = self.shdu_header(name='A_FUV', units='mag', filename='voxel_Afuv', dim=dim)

    # This is for a test of the voxel Emissions before streaming
    wav = np.append(constants.wavelengths[constants.nDust], species.moleculeWavelengths)
    nDust = constants.wavelengths[constants.nDust].size

     # The dimensions of the emissions are the expected velocity range (5 for clumps, 9 for interclumps),
    #the number of wavelengths at which the emission is calculated (#molecules + 333 dust wavelengths),
    #and finally the number of voxels.
    dim = [len(species.moleculeWavelengths)+nDust, constants.clumpMaxIndeces, self.__voxelNumber]
    shdu_clump_intensity = self.shdu_header(name='Clump intensity', units='K', filename='intensity_clump', dim=dim)
    shdu_clump_tau = self.shdu_header(name='Clump optical depth', units='1/cm', filename='opticalDepth_clump', dim=dim)
    
    dim = [len(species.moleculeWavelengths)+nDust, constants.interclumpMaxIndeces, self.__voxelNumber]
    shdu_interclump_intensity = self.shdu_header(name='Clump intensity', units='K', filename='intensity_interclump', dim=dim)
    shdu_interclump_tau = self.shdu_header(name='Clump optical depth', units='1/cm', filename='opticalDepth_interclump', dim=dim)

    with tqdm(total=len(self.__voxels), desc='Voxel emissions', miniters=1, dynamic_ncols=True) as progress:
      
      for i,voxel in enumerate(self.__voxels):
        
        gc.collect()
        
        emission = voxel.getEmission()

        # Optain the voxel emission data
        clumpIntensity = (emission[0])[0]

        if False:
          print()
          print(voxel.getPosition())
          print()
          plt.loglog(wav[nDust:], clumpIntensity[0,nDust:], marker='x', ms=2, ls='')
          plt.loglog(wav[:nDust], clumpIntensity[0,:nDust], marker='', ls='-', lw=1)
          plt.show()

        clumpTau = (emission[1])[0]
        clumpVelocity = voxel.getClumpVelocity()
        interclumpIntensity = (emission[0])[1]
        interclumpTau = (emission[1])[1]
        interclumpVelocity = voxel.getInterclumpVelocity()

        # Transform the voxel emission data to the maximum size in the model (to deal with numpy nd arrays)
        while clumpIntensity[:,0].size<constants.clumpMaxIndeces:
          clumpIntensity = np.append(clumpIntensity, np.zeros((1,clumpIntensity[0,:].size)), axis=0)
          clumpTau = np.append(clumpTau, np.zeros((1,clumpTau[0,:].size)), axis=0)
          clumpVelocity = np.append(clumpVelocity, [np.nan], axis=0)

        while interclumpIntensity[:,0].size<constants.interclumpMaxIndeces:
          interclumpIntensity = np.append(interclumpIntensity, np.zeros((1,interclumpIntensity[0,:].size)), axis=0)
          interclumpTau = np.append(interclumpTau, np.zeros((1,interclumpTau[0,:].size)), axis=0)
          interclumpVelocity = np.append(interclumpVelocity, [np.nan], axis=0)

        shdu_position.write(voxel.getPosition())
        velocity = voxel.getVelocity()
        if isinstance(velocity, float):
          shdu_velocity.write(np.array([velocity]))
        else:
          shdu_velocity.write(np.linalg.norm(velocity))
        shdu_clump_velocity.write(clumpVelocity)
        shdu_interclump_velocity.write(interclumpVelocity)
        shdu_FUV.write(np.array([voxel.getFUV()]))
        shdu_Afuv.write(np.array([voxel.getAfuv()]))

        if debug:

          shdu_mass.write(np.array([voxel.getMass()]))
          shdu_clump_mass.write(np.array([voxel.getClumpMass()]))
          shdu_interclump_mass.write(np.array([voxel.getInterclumpMass()]))
          shdu_density.write(np.array([voxel.getDensity()]))

        shdu_clump_intensity.write(np.array(clumpIntensity))
        shdu_clump_tau.write(np.array(clumpTau))

        shdu_interclump_intensity.write(interclumpIntensity)
        shdu_interclump_tau.write(interclumpTau)
        
        progress.update()
      
      progress.close()
    
    print('\nData files have been written successfully.\n')
    
    return

  def getVoxelNumber(self):
    return self.__voxelNumber

  def getVoxelPositions(self):
    return np.array([self.__x, self.__y, self.__z])

  def allVoxels(self):
    # Just in case all of the Voxel() instances need to be retrieved
    return self.__voxels

  def totalEmission(self):
    # Return the emission from all of the voxels, separated by observed velocity
    return np.array([self.__voxelIntensity.sum(1),self.__voxelOpticalDepth.sum(1)])

  def clumpEmission(self):
    # Return the emission from all of the voxels, separated by observed velocity
    return np.array([self.__voxelIntensity[:,0],self.__voxelOpticalDepth[:,0]])

  def interclumpEmission(self):
    # Return the emission from all of the voxels, separated by observed velocity
    return np.array([self.__voxelIntensity[:,1],self.__voxelOpticalDepth[:,1]])

  def getFUV(self):
    return self.__voxelFUV

  def getAfuv(self):
    return self.__voxelAfuv

  def printVoxels(self):
    for voxel in self.__voxels: voxel.printVoxel()
    return

  def shdu_header(self, name='', units='', filename=None, dim=None):

    if filename==None or dim==None: return

    header = fits.Header()

    header['SIMPLE'] = (True, 'conforms to FITS standard')
    header['BITPIX'] = (-64, 'element size')
    header['NAXIS'] = (len(dim), 'number of axes')
    header['EXTEND'] = True

    for i in range(len(dim)):
      header['NAXIS{}'.format(i+1)] = dim[i]

    # header['NAXIS1'] = self.__constants.velocityBins.size#emission[0].data[0,0,0,:].size
    # header['NAXIS2'] = len(self.__species[0].getMolecules())+len(self.__species[1].getDust())#emission[0].data[0,0,:,0].size
    # header['NAXIS3'] = 2#emission[0].data[:,0,0,0].size
    # header['NAXIS4'] = self.__voxelNumber#emission[0].data[0,:,0,0].size

    header['NAME'] = name
    header['UNITS'] = units

    directory = constants.HISTORYPATH + constants.directory + constants.history
    
    if not os.path.exists(directory): os.mkdir(directory)
    
    if not '.fits' in filename: filename = directory + filename + '.fits'
    else: filename = directory + filename
    
    if os.path.exists(filename): os.remove(filename)
    
    shdu = fits.StreamingHDU(filename, header)
    
    return shdu
