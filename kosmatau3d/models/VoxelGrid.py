import os
import numpy as np
from astropy.io import fits
import matplotlib.pyplot as plt
from tqdm import tqdm
#import progressbar as pb
from numba import jit
import importlib as il
import gc

from . import constants
from . import species
from .Voxel import *
from . import interpolations
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
    self.__voxelFUVabsorption = []

    self.__x = []
    self.__y = []
    self.__z = []

    constants.history = 'r{}_cm{}_d{}_uv{}/'.format(int(constants.resolution), '-'.join(str(f) for f in constants.clumpMassFactor), constants.densityFactor, constants.globalUV)

    return

  def __initialiseGrid(self):
    self.__voxels = []
    for i in range(self.__voxelNumber): self.__voxels.append(Voxel(i))
    return

  def __str__(self):
    return 'VoxelGrid\n  ->{} voxels\n  ->intensity {}\n  ->optical depth {}'.format(self.__voxelNumber, sum(self.__voxelIntensity), -np.log(np.exp(-self.__voxelOpticalDepth)))

  def __calculateProperties(self, X, Y, Z):
    # This is a method to calculate the dict to unpack into the argument for Voxel.setProperties().

    x,y = np.meshgrid(np.linspace(X-.5*constants.voxel_size, X+.5*constants.voxel_size,3), \
                      np.linspace(Y-.5*constants.voxel_size, Y+.5*constants.voxel_size,3))
    rPol = np.array([x.flatten(), y.flatten()]).T
    rPol = np.linalg.norm(rPol, axis=1)

    # Mass
    clumpMass = [interpolations.interpolateClumpMass(rPol), interpolations.interpolateInterclumpMass(rPol)]
    clumpMass = [constants.clumpMassFactor[ens]*np.asarray(clumpMass).mean(1)[ens] for ens in range(len(constants.clumpMassNumber))]

    # Velocity
    velocity = interpolations.interpolateRotationalVelocity(rPol)
    
    if constants.fromEarth:

      # Calculate the correction to the voxel velocity vectors
      relativeRpol = np.sqrt((x.flatten()-constants.rGalEarth)**2+y.flatten()**2)
      relativePhi = np.arctan2(y.flatten(), x.flatten()-constants.rGalEarth)
      relativeSigma = np.arccos((rPol**2+relativeRpol**2-constants.rGalEarth**2)/(2*rPol*relativeRpol))
      sigma = np.arctan2(Z, abs(x.flatten()-constants.rGalEarth))

      # Correct the relative velocity of the voxel
      velocityEarth = interpolations.interpolateRotationalVelocity(constants.rGalEarth)
      velocityCirc = velocity - velocityEarth*rPol/constants.rGalEarth

      velocity = (np.sign(relativePhi) * velocityCirc * np.sin(relativeSigma) * np.cos(sigma)).mean()

      if (rPol==0).any(): velocity = 0
      #self.__velocity = (velocity.mean()) * np.sin(self.__phi)
    
    else:
      velocity = np.array(velocity)

    # Use this to check the evaluation of the velocity field. It is still not working correctly...
    #print(self.__velocity)

    ensembleDispersion = interpolations.interpolateVelocityDispersion(rPol)
    
    ensembleDispersion = ensembleDispersion.mean()

    # Ensemble density
    ensembleDensity = interpolations.interpolateDensity(rPol)
    ensembleDensity = constants.densityFactor*ensembleDensity.mean()

    # FUV
    FUV = interpolations.interpolateFUVfield(rPol, Z)/constants.normUV*constants.globalUV
    FUV = np.clip(FUV, 1, None).mean()

    self.__properties = {  \
                        # Model parameters
                          'fromGrid' : True, \

                        # Voxel properties
                          'velocity' : velocity, \
                'ensembleDispersion' : ensembleDispersion, \
                      'ensembleMass' : clumpMass, \
                   'ensembleDensity' : [ensembleDensity,1911], \
                               'FUV' : [FUV,1] \
                        }

    return

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

  def calculateEmission(self, index=0, debug=False, timed=False, verbose=False):
    # This will initialise the grid of voxels and calculate their emission. This has to be
    #done in succession for each voxel since the calculations are modular (the temporary
    #voxel terms are changed when calculating different voxels). This can be rerun and it
    #will reinitialise the grid.
    
    if timed:
      t0 = time()

    #print(observations.tauCenterline)
    interpolations.initialise()
    self.__initialiseGrid()
    
    if timed:
      t1 = time()
      print('Grid initialised: {:.4f} s'.format(t1-t0))

    print('\nCalculating Grid Emission...')

    x,y,z = self.__shape.voxelCartesianPositions()
    r,phi = self.__shape.voxelPolarPositions()
    print()
    #self.__unusedVoxels = []
    with tqdm(total=len(self.__voxels), desc='Voxels initialised', miniters=1, dynamic_ncols=True) as progress:
      for i,voxel in enumerate(self.__voxels):
        
        if timed:
          t2 = time()

        if verbose:
            print('\nMax X, Radius:', max(x), r[i], '\n')

        self.__x.append(x[i])
        self.__y.append(y[i])
        self.__z.append(z[i])

        voxel.setIndex(i)#-len(self.__unusedVoxels))
        voxel.setPosition(x[i], y[i], z[i], r[i], phi[i])
        self.__calculateProperties(x[i], y[i], z[i])
        voxel.setProperties(**self.__properties)
        
        if timed:
          print('\nVoxel initialised: {:.4f} s'.format(time()-t2))

        self.__voxelFUVabsorption.append(voxel.getFUVabsorption())
        voxel.calculateEmission()
        
        if timed:
          print('Voxel calculated: {:.4f} s / {:.4f} s'.format(time()-t2, time()-t1))
          print()

        progress.update()

      progress.close()

      print('\nEmission calculated successfully.')
    
    return

  def writeEmission(self, verbose=False, debug=False):
    # This will stream the model parameters and emission to FITS files in the corresponding
    #directory for the model.
    
    print('\nStreaming to fits files...')

    if debug:

      # dim = [1, self.__voxelNumber]
      # shdu_mass = self.shdu_header(name='Mass', units='Msol', filename='voxel_mass', dim=dim)

      dim = [len(constants.clumpMassNumber), self.__voxelNumber]
      shdu_clump_mass = self.shdu_header(name='Ensemble Mass', units='Msol', filename='voxel_ensemble_mass', dim=dim)

      # dim = [1, self.__voxelNumber]
      # shdu_interclump_mass = self.shdu_header(name='Interclump Mass', units='Msol', filename='voxel_interclump_mass', dim=dim)

      dim = [len(constants.clumpMassNumber), self.__voxelNumber]
      shdu_density = self.shdu_header(name='Density', units='cm^-3', filename='voxel_density', dim=dim)

    dim = [3, self.__voxelNumber]
    shdu_position = self.shdu_header(name='Position', units='pc', filename='voxel_position', dim=dim)

    dim = [1, self.__voxelNumber]
    shdu_velocity = self.shdu_header(name='Velocity', units='km/s', filename='voxel_velocity', dim=dim)
    
    # dim = [constants.clumpMaxIndeces[0], self.__voxelNumber]
    # shdu_clump_velocity = self.shdu_header(name='Clump velocity', units='km/s', filename='voxel_clump_velocity', dim=dim)
    
    # dim = [constants.clumpMaxIndeces[1], self.__voxelNumber]
    # shdu_interclump_velocity = self.shdu_header(name='Interclump Velocity', units='km/s', filename='voxel_interclump_velocity', dim=dim)
    
    dim = [len(constants.clumpMassNumber),self.__voxelNumber]
    shdu_FUV = self.shdu_header(name='FUV', units='Draine', filename='voxel_fuv', dim=dim)
    shdu_FUVabsorption = self.shdu_header(name='tau_FUV', units='mag', filename='voxel_FUVabsorption', dim=dim)

    # This is for a test of the voxel Emissions before streaming
    wav = np.append(constants.wavelengths[constants.nDust], species.moleculeWavelengths)
    nDust = constants.wavelengths[constants.nDust].size

     # The dimensions of the emissions are the expected velocity range (5 for clumps, 9 for interclumps),
    #the number of wavelengths at which the emission is calculated (#molecules + 333 dust wavelengths),
    #and finally the number of voxels.
    # dim = [len(species.moleculeWavelengths)+nDust, constants.clumpMaxIndeces[0], self.__voxelNumber]
    dim = [len(species.moleculeWavelengths)+nDust, constants.velocityRange.size, self.__voxelNumber]
    shdu_clump_emissivity = self.shdu_header(name='Clump emissivity', units='K/pc', filename='emissivity_clump', dim=dim)
    shdu_clump_absorption = self.shdu_header(name='Clump absorption', units='1/pc', filename='absorption_clump', dim=dim)
    
    # dim = [len(species.moleculeWavelengths)+nDust, constants.clumpMaxIndeces[1], self.__voxelNumber]
    # shdu_interclump_intensity = self.shdu_header(name='Clump intensity', units='K', filename='intensity_interclump', dim=dim)
    # shdu_interclump_tau = self.shdu_header(name='Clump optical depth', units='1/cm', filename='opticalDepth_interclump', dim=dim)

    with tqdm(total=len(self.__voxels), desc='Voxel emissions', miniters=1, dynamic_ncols=True) as progress:
      
      for i,voxel in enumerate(self.__voxels):
        
        gc.collect()
        
        epsilon = np.sum(voxel.getEmissivity(), axis=0)
        kappa = np.sum(voxel.getAbsorption(), axis=0)

        # Optain the voxel emission data
        # clumpIntensity = intensity[0]
        # clumpTau = opticalDepth[0]
        # clumpVelocity = voxel.getClumpVelocity()[0]
        # # print(clumpVelocity.size)
        # interclumpIntensity = intensity[1]
        # interclumpTau = opticalDepth[1]
        # interclumpVelocity = voxel.getClumpVelocity()[1]

        if False:
          print()
          print(voxel.getPosition())
          print()
          plt.loglog(wav[nDust:], clumpIntensity[0,nDust:], marker='x', ms=2, ls='')
          plt.loglog(wav[:nDust], clumpIntensity[0,:nDust], marker='', ls='-', lw=1)
          plt.show()

        # Transform the voxel emission data to the maximum size in the model (to deal with numpy nd arrays)
        # while clumpIntensity[:,0].size<constants.clumpMaxIndeces[0]:
        #   clumpIntensity = np.append(clumpIntensity, np.zeros((1,clumpIntensity[0,:].size)), axis=0)
        #   clumpTau = np.append(clumpTau, np.zeros((1,clumpTau[0,:].size)), axis=0)
        #   clumpVelocity = np.append(clumpVelocity, [np.nan], axis=0)

        # while interclumpIntensity[:,0].size<constants.clumpMaxIndeces[1]:
        #   interclumpIntensity = np.append(interclumpIntensity, np.zeros((1,interclumpIntensity[0,:].size)), axis=0)
        #   interclumpTau = np.append(interclumpTau, np.zeros((1,interclumpTau[0,:].size)), axis=0)
        #   interclumpVelocity = np.append(interclumpVelocity, [np.nan], axis=0)

        shdu_position.write(voxel.getPosition())
        velocity = voxel.getVelocity()[0]
        if isinstance(velocity, float):
          shdu_velocity.write(np.array([velocity]))
        else:
          shdu_velocity.write(np.linalg.norm(velocity))
        # shdu_clump_velocity.write(clumpVelocity)
        # shdu_interclump_velocity.write(interclumpVelocity)
        shdu_FUV.write(np.array([voxel.getFUV()]))
        shdu_FUVabsorption.write(np.asarray([voxel.getFUVabsorption()]))

        if debug:

          # shdu_mass.write(np.array([voxel.getMass()]))
          shdu_clump_mass.write(np.array([voxel.getEnsembleMass()]))
          # shdu_interclump_mass.write(np.array([voxel.getInterclumpMass()]))
          shdu_density.write(np.array([voxel.getDensity()]))

        shdu_clump_emissivity.write(np.array(epsilon))
        shdu_clump_absorption.write(np.array(kappa))

        # shdu_interclump_intensity.write(interclumpIntensity)
        # shdu_interclump_tau.write(interclumpTau)
        
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
    
    if not os.path.exists(directory): os.makedirs(directory)
    
    if not '.fits' in filename: filename = directory + filename + '.fits'
    else: filename = directory + filename
    
    if os.path.exists(filename): os.remove(filename)
    
    shdu = fits.StreamingHDU(filename, header)
    
    return shdu
