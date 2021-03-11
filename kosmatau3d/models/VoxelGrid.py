import os
import numpy as np
from astropy.io import fits
import matplotlib.pyplot as plt
from tqdm import tqdm
#import progressbar as pb
from numba import jit
import importlib as il
import gc
from multiprocessing import Pool

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

    self.__voxelSpeciesEmissivity = []
    self.__voxelSpeciesAbsorption = []
    self.__voxelDustEmissivity = []
    self.__voxelDustAbsorption = []

    self.__voxelMass = []
    self.__voxelDensity = []
    self.__voxelFUV = []
    self.__voxelFUVabsorption = []
    
    self.__voxelVelocities = []
    self.__voxelDispersion = []

    self.__x = []
    self.__y = []
    self.__z = []

    constants.history = 'r{}_cm{}_d{}_uv{}/'.format(int(constants.voxel_size), '-'.join(str(f) for f in constants.clumpMassFactor), constants.densityFactor, constants.globalUV)

    return

  def __initialiseGrid(self):
    self.__voxels = []
    for i in range(self.__voxelNumber): self.__voxels.append(Voxel(i))
    return

  def __str__(self):
    return 'VoxelGrid\n  ->{} voxels\n'.format(self.__voxelNumber)

  def __calculateProperties(self, X, Y, Z, average=False):
    # This is a method to calculate the dict to unpack into the argument for Voxel.setProperties().

    if average:
      x,y = np.meshgrid(np.linspace(X-.5*constants.voxel_size, X+.5*constants.voxel_size,2), \
                        np.linspace(Y-.5*constants.voxel_size, Y+.5*constants.voxel_size,2))
    else:
      x = np.array([X])
      y = np.array([Y])
    rPol = np.array([x.flatten(), y.flatten()]).T
    rPol = np.linalg.norm(rPol, axis=1)

    # Mass
    ensembleMass = [interpolations.interpolateClumpMass(rPol), interpolations.interpolateInterclumpMass(rPol)]
    ensembleMass = [constants.clumpMassFactor[ens]*np.asarray(ensembleMass).mean(1)[ens] for ens in range(len(constants.clumpMassNumber))]

    # Velocity
    velocity = interpolations.interpolateRotationalVelocity(rPol)
    
    if constants.fromEarth:

      # Calculate the correction to the voxel velocity vectors
      relativeRpol = np.sqrt((x.flatten()-constants.rGalEarth)**2+y.flatten()**2)
      relativePhi = np.arctan2(y.flatten(), constants.rGalEarth-x.flatten())
      relativeSigma = np.arccos((rPol**2+relativeRpol**2-constants.rGalEarth**2)/(2*rPol*relativeRpol))
      relativeTheta = np.arctan(Z / relativeRpol)

      # Correct the relative velocity of the voxel
      velocityEarth = interpolations.interpolateRotationalVelocity(constants.rGalEarth)
      velocityCirc = velocity - velocityEarth*rPol/constants.rGalEarth

      velocity = (np.sign(relativePhi) * velocityCirc * np.sin(relativeSigma) * np.cos(relativeTheta))
      # velocity = (np.sign(relativePhi) * velocityCirc * np.sin(relativeSigma))
      # velocity = np.sign(np.arctan2(Y,X))*velocity*np.sin(relativeSigma) - velocityEarth*np.sin(relativePhi)
      ensembleDispersion = np.sqrt(velocity.std()**2+(constants.ensembleDispersion)**2) #this is leftover from Silke's version
      velocity = velocity.mean()

      if (rPol==0).any(): velocity = 0
      #self.__velocity = (velocity.mean()) * np.sin(self.__phi)
    
    else:
      velocity = velocity.mean()
      ensembleDispersion = constants.ensembleDispersion

    # Use this to check the evaluation of the velocity field. It is still not working correctly...
    #print(self.__velocity)

    # ensembleDispersion = interpolations.interpolateVelocityDispersion(rPol)
    #
    # ensembleDispersion = ensembleDispersion.mean()

    # Ensemble density
    ensembleDensity = interpolations.interpolateDensity(rPol)
    ensembleDensity = [constants.densityFactor*ensembleDensity.mean(),1911]

    # FUV
    FUV = interpolations.interpolateFUVfield(rPol, Z)/constants.normUV*constants.globalUV
    FUV = [np.clip(FUV, 1, None).mean(),1]
    
    # Save the properties in private lists
    self.__voxelVelocities.append(velocity)
    self.__voxelDispersion.append(ensembleDispersion)
    self.__voxelMass.append(ensembleMass)
    self.__voxelDensity.append(ensembleDensity)
    self.__voxelFUV.append(FUV)

    self.__properties = {
                        # Model parameters
                          'fromGrid' : True,

                        # Voxel properties
                          'velocity' : velocity,
                'ensembleDispersion' : ensembleDispersion,
                      'ensembleMass' : ensembleMass,
                   'ensembleDensity' : ensembleDensity,
                               'FUV' : FUV
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

  def calculateEmission(self, index=0, debug=False, timed=False, verbose=False, multiprocessing=0):
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
    
    # Setup fits files to stream the voxel emissivity and absorption.
    #species
    dim = [len(species.moleculeWavelengths), constants.velocityRange.size, self.__voxelNumber]
    shdu_voxel_emissivity_species = self.shdu_header(name='Clump species emissivity', units='K/pc', molecules=True, filename='species_emissivity', velocity=True, dim=dim)
    shdu_voxel_absorption_species = self.shdu_header(name='Clump species absorption', units='1/pc', molecules=True, filename='species_absorption', velocity=True, dim=dim)
    #dust
    nDust = constants.wavelengths[constants.nDust].size
    dim = [nDust, constants.velocityRange.size, self.__voxelNumber]
    shdu_voxel_emissivity_dust = self.shdu_header(name='Clump dust emissivity', units='K/pc', dust=True, filename='dust_emissivity', velocity=True, dim=dim)
    shdu_voxel_absorption_dust = self.shdu_header(name='Clump dust absorption', units='1/pc', dust=True, filename='dust_absorption', velocity=True, dim=dim)

    dim = [len(constants.clumpMassNumber), self.__voxelNumber]
    shdu_ensemble_mass = self.shdu_header(name='Ensemble Mass', units='Msol', filename='voxel_ensemble_mass', dim=dim)
    dim = [len(constants.clumpMassNumber), self.__voxelNumber]
    shdu_ensemble_density = self.shdu_header(name='Density', units='cm^-3', filename='voxel_density', dim=dim)
    dim = [3, self.__voxelNumber]
    shdu_voxel_position = self.shdu_header(name='Position', units='pc', filename='voxel_position', dim=dim)
    dim = [1, self.__voxelNumber]
    shdu_voxel_velocity = self.shdu_header(name='Velocity', units='km/s', filename='voxel_velocity', dim=dim)
    dim = [len(constants.clumpMassNumber),self.__voxelNumber]
    shdu_FUV = self.shdu_header(name='FUV', units='Draine', filename='voxel_fuv', dim=dim)
    shdu_FUVabsorption = self.shdu_header(name='tau_FUV', units='mag', filename='voxel_FUVabsorption', dim=dim)
    
    #test of multiprocessing
    def getProperties(voxel):
      i,voxel = voxel

      self.__x.append(x[i])
      self.__y.append(y[i])
      self.__z.append(z[i])

      voxel.setIndex(i)#-len(self.__unusedVoxels))
      voxel.setPosition(x[i], y[i], z[i], r[i], phi[i])
      self.__calculateProperties(x[i], y[i], z[i])
      voxel.setProperties(**self.__properties)
      
      voxel.calculateEmission()
      
      return voxel
      
    if multiprocessing:
      pool = Pool(multiprocessing)
      chunksize = int(len(self.__voxels)/1000/multiprocessing)
      voxels = pool.imap(getProperties, list(enumerate(self.__voxels)), chunksize)
    else:
      voxels = self.__voxels
    
    with tqdm(total=len(self.__voxels), desc='Voxels initialised', miniters=1, dynamic_ncols=True) as progress:
      
      # for i,voxel in enumerate(self.__voxels):
      for i,voxel in enumerate(voxels):
        
        if timed:
          t2 = time()

        if verbose:
          print('\nMax X, Radius:', max(x), r[i], '\n')
          
        if not multiprocessing:
          voxel = getProperties((i,voxel))

        ##modified for multiprocessing
        #------------------------------
        # self.__x.append(x[i])
        # self.__y.append(y[i])
        # self.__z.append(z[i])
        #------------------------------
        X,Y,Z = voxel.getPosition()
        self.__x.append(X)
        self.__y.append(Y)
        self.__z.append(Z)

        ## moved to getProperties
        #------------------------------------------------------
        # voxel.setIndex(i)#-len(self.__unusedVoxels))
        # voxel.setPosition(x[i], y[i], z[i], r[i], phi[i])
        # self.__calculateProperties(x[i], y[i], z[i])
        # voxel.setProperties(**self.__properties)
        #------------------------------------------------------
        
        if timed:
          print('\nVoxel initialised: {:.4f} s'.format(time()-t2))

        #this is likely unneeded now that I am directly writing to a fits file
        self.__voxelFUVabsorption.append(voxel.getFUVabsorption())
        
        # Save voxel properties
        # print(np.asarray([self.__properties['velocity']]))
        # shdu_voxel_position.write(np.asarray([x[i],y[i],z[i]]))
        # shdu_voxel_velocity.write(np.asarray([self.__properties['velocity']]))
        # shdu_ensemble_mass.write(np.asarray(self.__properties['ensembleMass']))
        # shdu_ensemble_density.write(np.asarray(self.__properties['ensembleDensity']))
        # shdu_FUV.write(np.asarray(self.__properties['FUV']))
        # shdu_FUVabsorption.write(np.asarray(self.__voxelFUVabsorption[-1]))
        
        shdu_voxel_position.write(voxel.getPosition())
        velocity = voxel.getVelocity()[0]
        if isinstance(velocity, float):
          shdu_voxel_velocity.write(np.array([velocity]))
        else:
          shdu_voxel_velocity.write(np.linalg.norm(velocity))
        shdu_ensemble_mass.write(np.asarray(voxel.getEnsembleMass()))
        shdu_ensemble_density.write(np.asarray(voxel.getDensity()))
        shdu_FUV.write(np.array([voxel.getFUV()]))
        shdu_FUVabsorption.write(np.asarray([voxel.getFUVabsorption()]))
        
        #moved to getProperties
        #------------------------------------
        # Calculate emissivity and absorption
        # voxel.calculateEmission()
        #------------------------------------
        
        # Save emissivity and absorption
        shdu_voxel_emissivity_species.write(np.sum(voxel.getSpeciesEmissivity(), axis=0))
        shdu_voxel_absorption_species.write(np.sum(voxel.getSpeciesAbsorption(), axis=0))
        shdu_voxel_emissivity_dust.write(np.sum(voxel.getDustEmissivity(), axis=0))
        shdu_voxel_absorption_dust.write(np.sum(voxel.getDustAbsorption(), axis=0))
        
        voxels[i] = None
        # self.__voxels[i] = None
        # print(len(self.__voxels))
        
        if timed:
          print('Voxel calculated: {:.4f} s / {:.4f} s'.format(time()-t2, time()-t1))
          print()

        progress.update()

      progress.close()

      print('\nEmission calculated successfully.')
    
    return
  
  def writeDataAlt(self):
    # Not yet working.
    
    voxel_species_emissivity = fits.ImageHDU(self.__voxelSpeciesEmissivity, header_species)
    voxel_species_absorption = fits.ImageHDU(self.__voxelSpeciesAbsorption, header_species)
    voxel_dust_emissivity = fits.ImageHDU(self.__voxelDustEmissivity, header_dust)
    voxel_dust_absorption = fits.ImageHDU(self.__voxelDustAbsorption, header_dust)
    
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
    dim = [len(species.moleculeWavelengths), constants.velocityRange.size, self.__voxelNumber]
    shdu_clump_emissivity_species = self.shdu_header(name='Clump species emissivity', units='K/pc', molecules=True, filename='species_emissivity_clump', dim=dim)
    shdu_clump_absorption_species = self.shdu_header(name='Clump species absorption', units='1/pc', molecules=True, filename='species_absorption_clump', dim=dim)
    
    dim = [nDust, constants.velocityRange.size, self.__voxelNumber]
    shdu_clump_emissivity_dust = self.shdu_header(name='Clump dust emissivity', units='K/pc', dust=True, filename='dust_emissivity_clump', dim=dim)
    shdu_clump_absorption_dust = self.shdu_header(name='Clump dust absorption', units='1/pc', dust=True, filename='dust_absorption_clump', dim=dim)
    
    # dim = [len(species.moleculeWavelengths)+nDust, constants.clumpMaxIndeces[1], self.__voxelNumber]
    # shdu_interclump_intensity = self.shdu_header(name='Clump intensity', units='K', filename='intensity_interclump', dim=dim)
    # shdu_interclump_tau = self.shdu_header(name='Clump optical depth', units='1/cm', filename='opticalDepth_interclump', dim=dim)

    with tqdm(total=len(self.__voxels), desc='Voxel emissions', miniters=1, dynamic_ncols=True) as progress:
      
      for i,voxel in enumerate(self.__voxels):
        
        gc.collect()
        
        epsilon_species = np.sum(voxel.getSpeciesEmissivity(), axis=0)
        kappa_species = np.sum(voxel.getSpeciesAbsorption(), axis=0)
        epsilon_dust = np.sum(voxel.getDustEmissivity(), axis=0)
        kappa_dust = np.sum(voxel.getDustAbsorption(), axis=0)

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

        shdu_clump_emissivity_species.write(np.array(epsilon_species))
        shdu_clump_absorption_species.write(np.array(kappa_species))
        shdu_clump_emissivity_dust.write(np.array(epsilon_dust))
        shdu_clump_absorption_dust.write(np.array(kappa_dust))

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

  def getFUVabsorption(self):
    return self.__voxelFUVabsorption

  def printVoxels(self):
    for voxel in self.__voxels: voxel.printVoxel()
    return

  def shdu_header(self, name='', units='', molecules=False, dust=False, filename=None, velocity=False, dim=None):

    if filename==None or dim==None: return

    header = fits.Header()

    header['SIMPLE'] = (True, 'conforms to FITS standard')
    header['BITPIX'] = (-64, 'element size')
    header['NAXIS'] = (len(dim), 'number of axes')
    header['EXTEND'] = True
    if molecules:
      header['SPECIES'] = ', '.join(species.molecules)
    if dust:
      header['DUST'] = ', '.join(constants.dustNames[constants.nDust])

    for i in range(len(dim)):
      header['NAXIS{}'.format(i+1)] = dim[i]
      if velocity&i==0:
        header['CTYPE1'] = 'Transition/wavelength'
        header['CUNIT1'] = 'm'
        header['CRPIX1'] = 'N/A'
        header['CRVAL1'] = 'N/A'
        header['CDELT1'] = 'N/A'
      elif velocity&i==1:
        header['CTYPE2'] = 'Velocity'
        header['CUNIT2'] = 'km\s'
        header['CRPIX2'] = (header['NAXIS2']-1)/2
        header['CRVAL2'] = 0
        header['CDELT2'] = (constants.velocityRange[-1]-constants.velocityRange[0])/(header['NAXIS2']-1)
      elif velocity&i==2:
        header['CTYPE3'] = 'Voxel'
        header['CUNIT1'] = 'N/A'
        header['CRPIX1'] = 'N/A'
        header['CRVAL1'] = 'N/A'
        header['CDELT1'] = 'N/A'

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
