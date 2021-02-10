import numpy as np
from numba import jit_module
import multiprocessing
from multiprocessing import Pool
from functools import partial
import matplotlib.pyplot as plt
# import matplotlib as mpl
# mpl.use('Qt4Agg')
import copy as c
import cmath
import scipy.interpolate as interpolate
import scipy.optimize as op
from astropy.io import fits
import importlib as il
from tqdm import tqdm
from time import time

import warnings

from .. import constants
from .. import interpolations
from .. import observations
from .. import radiativeTransfer

def eTildeReal(file='Ereal.dat'):
  eReal = np.genfromtxt(constants.GRIDPATH+file, names=['x', 'Ereal'])
  return (eReal['x'],eReal['Ereal'])

def eTildeImaginary(file='Eimag.dat'):
  eImaginary = np.genfromtxt(constants.GRIDPATH+file, names=['x', 'Eimaginary'])
  return (eImaginary['x'],eImaginary['Eimaginary'])

def calculateObservation(directory='', dim='xy', slRange=[(-np.pi,np.pi), (-np.pi/2,np.pi/2)], nsl=[50,25], terminal=True, plotV=False, multiprocessing=0, debug=False, verbose=False):

  if debug:
    sl = [5,5]
  
  constants.velocityRange = np.linspace(-300, 300, 500)
    
  # print('Load data')

  radiativeTransfer.voxelPositions = fits.open(directory+'/voxel_position.fits', mode='denywrite')
  radiativeTransfer.voxelVelocities = fits.open(directory+'/voxel_velocity.fits', mode='denywrite')
  radiativeTransfer.tempSpeciesEmissivity = fits.open(directory+'/species_emissivity.fits', mode='denywrite')# in K/pc
  radiativeTransfer.tempSpeciesAbsorption = fits.open(directory+'/species_absorption.fits', mode='denywrite')# in 1/pc
  radiativeTransfer.tempDustEmissivity = fits.open(directory+'/dust_emissivity.fits', mode='denywrite')# in K/pc
  radiativeTransfer.tempDustAbsorption = fits.open(directory+'/dust_absorption.fits', mode='denywrite')# in 1/pc
  
  # print('Data loaded :-)')

  xPositions,yPositions,zPositions = radiativeTransfer.voxelPositions[0].data[:,0],radiativeTransfer.voxelPositions[0].data[:,1],radiativeTransfer.voxelPositions[0].data[:,2]
  r = np.sqrt((xPositions-constants.rGalEarth)**2 + yPositions**2)

  radGrid = np.sqrt((xPositions-constants.rGalEarth)**2 + yPositions**2 + zPositions**2)
  lonGrid = np.arctan2(yPositions, -(xPositions-constants.rGalEarth))
  rPolar  = np.sqrt((xPositions-constants.rGalEarth)**2+yPositions**2)
  latGrid = np.arctan2(zPositions, rPolar)

  np.set_printoptions(threshold=100000)
  # print('\nLongitude\n', lonGrid, '\nLattitude\n', latGrid)

  #print('\nx\n', np.unique(xArray), '\ny\n', np.unique(yArray), '\nz\n', np.unique(zArray))

  if constants.fromEarth:
    # For an observation from Earth, the data is modified to account for Earth's position at (8750, 0, 0) pc.
    #The integrated intensity is then calculated in the y-z plane to account for different viewing angles across
    #the galaxy. This can be post-processed to convert to galactic coordinates.

    hdul = fits.HDUList()

    # Define the boundaries separating the inner and outer disk
    xBoundary = (xPositions>0)&(r>constants.rGalEarth)
    yBoundary = (yPositions<constants.rGalEarth)&(yPositions>-constants.rGalEarth)

    # np.set_printoptions(threshold=1000000)

    # Define sightlines calculated
    longrid = np.linspace(-np.pi, np.pi, num=nsl[0])
    latgrid = np.linspace(-np.pi/2, np.pi/2, num=nsl[1])
    # grid = np.meshgrid(lon, lat)
    # grid = np.array([grid[0].flatten(), grid[1].flatten()])

    # radiativeTransfer.vTqdm = tqdm(total=constants.velocityRange.size, desc='Observing velocity', miniters=1, dynamic_ncols=True)
    # if terminal: radiativeTransfer.slTqdm = tqdm(total=longrid.size*latgrid.size, desc='Sightline', miniters=1, dynamic_ncols=True)

    VintensityMapSpecies = []
    VintensityMapDust = []
    Vpositions    = []

    radiativeTransfer.sightlines = np.zeros((longrid.size, latgrid.size))

    if debug: velocityRange = [0]
    else: velocityRange = constants.velocityRange

    result = multiprocessCalculation(slRange=slRange, nsl=nsl, multiprocessing=multiprocessing, dim=dim, debug=debug)
    
    Vpositions = result[0]
    VintensityMapSpecies = result[1]
    VintensityMapDust = result[2]
    radiativeTransfer.sightlines = np.asarray(result[3]).max(0)
    vmin,vmax = result[4]
    
    # Save sightline lengths
    np.savetxt(directory+'/sightlines.csv', radiativeTransfer.sightlines, delimiter=',')
      
    # Convert to numpy arrays
    Vpositions = np.array(Vpositions[0])
    VintensityMapSpecies = np.array(VintensityMapSpecies)
    VintensityMapDust = np.array(VintensityMapDust)

    if verbose:
      print('Map position shape:')
      if Vpositions.ndim>1:
        print(Vpositions.shape)
      else:
        for p in VmapPositions: print(p.shape)
        
      print('Map intensity shapes:')
      print('Species')
      if VintensityMapSpecies.ndim>1:
        print(VintensityMapSpecies.shape)
      else:
        for intensity in intensityMapSpecies: print(intensity.shape)
      print('Dust')
      if VintensityMapDust.ndim>1:
        print(VintensityMapDust.shape)
      else:
        for intensity in intensityMapDust: print(intensity.shape)

    # Setup the data to be saved in a FITS file. It will create an HDU list with position, species, and dust HDUs.
    if not debug:
      # Create HDUs for the map position and intensity and add the velocity in the headers
      PositionHDU = fits.ImageHDU(Vpositions)
      
      print(VintensityMapSpecies.shape, np.shape(VintensityMapSpecies[0,0,0]))
      radiativeTransfer.intensitySpecies = VintensityMapSpecies
      IntensityHDUSpecies = fits.ImageHDU(VintensityMapSpecies)
      
      print(VintensityMapDust.shape, np.shape(VintensityMapDust[0,0,0]))
      radiativeTransfer.intensityDust = VintensityMapDust
      IntensityHDUDust = fits.ImageHDU(VintensityMapDust)
      # input('Integrated intensity shape: {}'.format(VintensityMap.shape))

      PositionHDU.header['TYPE1'] = 'Angle'
      PositionHDU.header['TYPE2'] = 'Position'
      PositionHDU.header['DIREC'] = 'Radial'

      IntensityHDUSpecies.header['TYPE'] = 'Species transitions'
      IntensityHDUSpecies.header['BUNIT'] = 'K'
      IntensityHDUSpecies.header['CTYPE1'] = 'Wavelength'
      IntensityHDUSpecies.header['CUNIT1'] = 'm'
      IntensityHDUSpecies.header['CRVAL1'] = 'N/A'
      IntensityHDUSpecies.header['CDELT1'] = 'N/A'
      IntensityHDUSpecies.header['CRPIX1'] = 'N/A'
      IntensityHDUSpecies.header['CTYPE2'] = 'GLON'
      IntensityHDUSpecies.header['CUNIT2'] = 'rad'
      IntensityHDUSpecies.header['CRVAL2'] = 0.
      IntensityHDUSpecies.header['CDELT2'] = 2*np.pi/(IntensityHDUSpecies.header['NAXIS2']-1)
      IntensityHDUSpecies.header['CRPIX2'] = (IntensityHDUSpecies.header['NAXIS2']-1)/2.
      IntensityHDUSpecies.header['CTYPE3'] = 'GLAT'
      IntensityHDUSpecies.header['CUNIT3'] = 'rad'
      IntensityHDUSpecies.header['CRVAL3'] = 0.
      IntensityHDUSpecies.header['CDELT3'] = np.pi/(IntensityHDUSpecies.header['NAXIS3']-1)
      IntensityHDUSpecies.header['CRPIX3'] = (IntensityHDUSpecies.header['NAXIS3']-1)/2.
      IntensityHDUSpecies.header['CTYPE4'] = 'Velocity'
      IntensityHDUSpecies.header['CUNIT4'] = 'km/s'
      IntensityHDUSpecies.header['CRVAL4'] = (vmax+vmin)/2.
      IntensityHDUSpecies.header['CDELT4'] = (vmax-vmin)/(IntensityHDUSpecies.header['NAXIS4']-1)
      IntensityHDUSpecies.header['CRPIX4'] = (IntensityHDUSpecies.header['NAXIS4'])/2.
      IntensityHDUSpecies.header['DIREC'] = 'Radial'
      IntensityHDUSpecies.header['SPECIES'] = radiativeTransfer.tempSpeciesEmissivity[0].header['SPECIES']

      IntensityHDUDust.header['TYPE'] = 'Dust continuum'
      IntensityHDUDust.header['BUNIT'] = 'K'
      IntensityHDUDust.header['CTYPE1'] = 'Wavelength'
      IntensityHDUDust.header['CUNIT1'] = 'm'
      IntensityHDUDust.header['CRVAL1'] = 'N/A'
      IntensityHDUDust.header['CDELT1'] = 'N/A'
      IntensityHDUDust.header['CRPIX1'] = 'N/A'
      IntensityHDUDust.header['CTYPE2'] = 'GLON'
      IntensityHDUDust.header['CUNIT2'] = 'rad'
      IntensityHDUDust.header['CRVAL2'] = 0.
      IntensityHDUDust.header['CDELT2'] = 2*np.pi/(IntensityHDUDust.header['NAXIS2']-1)
      IntensityHDUDust.header['CRPIX2'] = (IntensityHDUDust.header['NAXIS2']-1)/2.
      IntensityHDUDust.header['CTYPE3'] = 'GLAT'
      IntensityHDUDust.header['CUNIT3'] = 'rad'
      IntensityHDUDust.header['CRVAL3'] = 0.
      IntensityHDUDust.header['CDELT3'] = np.pi/(IntensityHDUDust.header['NAXIS3']-1)
      IntensityHDUDust.header['CRPIX3'] = (IntensityHDUDust.header['NAXIS3']-1)/2.
      IntensityHDUDust.header['CTYPE4'] = 'Velocity'
      IntensityHDUDust.header['CUNIT4'] = 'km/s'
      IntensityHDUDust.header['CRVAL4'] = (vmax+vmin)/2.
      IntensityHDUDust.header['CDELT4'] = (vmax-vmin)/(IntensityHDUDust.header['NAXIS4']-1)
      IntensityHDUDust.header['CRPIX4'] = (IntensityHDUDust.header['NAXIS4'])/2.
      IntensityHDUDust.header['DIREC'] = 'Radial'
      IntensityHDUDust.header['DUST'] = radiativeTransfer.tempDustEmissivity[0].header['DUST']

      hdul.append(PositionHDU)
      hdul.append(IntensityHDUSpecies)
      hdul.append(IntensityHDUDust)

      hdul.writeto(directory+'/channel_intensity.fits', overwrite=True)

      print('Intensity map written successfully :-)')

    radiativeTransfer.voxelPositions.close()
    radiativeTransfer.voxelVelocities.close()
    radiativeTransfer.tempSpeciesEmissivity.close()
    radiativeTransfer.tempSpeciesAbsorption.close()
    radiativeTransfer.tempDustEmissivity.close()
    radiativeTransfer.tempDustAbsorption.close()
  
  return

def multiprocessCalculation(slRange=[(-np.pi,np.pi), (-np.pi/2,np.pi/2)], nsl=[50,25], multiprocessing=0, dim='spherical', debug=False):
  
  Vpositions = []
  VintensityMapSpecies = []
  VintensityMapDust = []
  sightlines = []
  
  velChannel = partial(calculateVelocityChannel, slRange=slRange, nsl=nsl, dim=dim, debug=debug)
  vNum = constants.velocityRange.size
  # print('Pause')
  # input()
  t0 = time()
  
  if multiprocessing:
    pool = Pool(processes=multiprocessing)
    intensity = pool.imap(velChannel, list(enumerate(constants.velocityRange)), vNum/multiprocessing)
  else:
    intensity = []
    vTqdm = tqdm(total=constants.velocityRange.size, desc='Observing velocity', miniters=1, dynamic_ncols=True)
    # t0 = time()
    # print('{} of {} velocities'.format(0, constants.velocityRange.size))
    for iv in enumerate(constants.velocityRange):
      intensity.append(velChannel(iv))
      # print('\r{} of {} velocities: {:.2f}'.format(iv[0], constants.velocityRange.size, time()-t0))
      vTqdm.update()
      
  print(time()-t0)
      
  vmin = constants.velocityRange.max()
  vmax = constants.velocityRange.min()
  
  for n,i in enumerate(intensity):
  #   i.wait()
  #   print(i)
    if i[3].any():
      if constants.velocityRange[n]<vmin:
        vmin = constants.velocityRange[n]
      if constants.velocityRange[n]>vmax:
        vmax = constants.velocityRange[n]
  
      Vpositions.append(i[0])
      VintensityMapSpecies.append(i[1])
      VintensityMapDust.append(i[2])
      sightlines.append(i[3])
  
  return (Vpositions,VintensityMapSpecies,VintensityMapDust,sightlines,(vmin,vmax))

# for i_vel,velocity in enumerate(constants.velocityRange):
def calculateVelocityChannel(ivelocity, slRange=[(-np.pi,np.pi), (-np.pi/2,np.pi/2)], nsl=[50,25], dim='spherical', debug=False):

  # Convert the tuple to the desired index and velocity
  i_vel = ivelocity[0]
  velocity = ivelocity[1]
  print('Calculating velocity channel at {:.2f} km\s'.format(velocity))
  
  # Setup the sightlines that are calculated
  longrid = np.linspace(slRange[0][0], slRange[0][1], num=nsl[0])
  latgrid = np.linspace(slRange[1][0], slRange[1][1], num=nsl[1])
  sightlines = np.zeros((longrid.size, latgrid.size))

  # Find the voxels that exist at the observing velocity
  i_vox =  radiativeTransfer.tempSpeciesEmissivity[0].data[:,i_vel,:].any(1)#(clumpVelocity==velocity)
  # print(iV)

  # Update velocity progress bar
  # radiativeTransfer.vTqdm.update()
  
  # Reset sightline progress bar
  # if terminal: radiativeTransfer.slTqdm.reset()

  # Initialise the intensities and map
  position = []
  intensityMapSpecies = []
  intensityMapDust = []

  # Get indeces of voxels at the correct observing velocity
  # iClumpV = np.where(iCV)
  # iInterclumpV = np.where(iIV)

  if i_vox.any()==False:
    return 0,0,0,0

  # The voxel positions can be any of the voxels
  # radiativeTransfer.tempSpeciesEmissivity = tempSpeciesEmissivity#[iV,i,:] / constants.pc/100
  # radiativeTransfer.tempSpeciesAbsorption = tempSpeciesAbsorption#[iV,i,:] / constants.pc/100
  # radiativeTransfer.tempDustEmissivity = tempDustEmissivity#[iV,i,:] / constants.pc/100
  # radiativeTransfer.tempDustAbsorption = tempDustAbsorption#[iV,i,:] / constants.pc/100
  # radiativeTransfer.tempPosition = radiativeTransfer.voxelPositions[0].data[iV,:]
  # radiativeTransfer.tempClumpVelocity = clumpVelocity[iClumpV[0],:]
  # radiativeTransfer.tempInterclumpEmission = interclumpEmission[:,iInterclumpV[0],iInterclumpV[1],:]
  # radiativeTransfer.tempInterclumpPosition = voxelPositions[iInterclumpV[0],:]
  # radiativeTransfer.tempInterclumpVelocity = interclumpVelocity[iInterclumpV[0],:]

  velTest = []

  def sightlength(x, l):
    return constants.rGalEarth**2 - constants.rGal**2 + x**2 - 2*constants.rGalEarth*x*np.cos(l)

  for j,lat in enumerate(latgrid):

    # Initialise a list containing the intensities at each lattitude
    positionIntensitySpecies = []
    positionIntensityDust = []

    for i,lon in enumerate(longrid):

      # Calculate sightline length
      Rslh = op.root_scalar(sightlength, args=lon, x0=constants.rGalEarth, x1=constants.rGal).root
      thetaC = np.arctan(constants.hd/Rslh)
      # if verbose:
      #   print('\n',thetaC,'\n')
      if abs(lat)<thetaC:
        Rsl = Rslh/np.cos(lat)
      else:
        Rsl = constants.hd/np.sin(abs(lat))

      # Try to isolate voxels in LoS and work with all transitions, else do them one at a time
      try:
        result,vox,emission = setLOS(lon=lon, lat=lat, i_vox=i_vox, i_vel=i_vel, dim=dim, debug=debug)
      except OSError:
        print('OSError!!')
        intensity = []
        for spe in range(radiativeTransfer.tempSpeciesEmissivity[0].shape[-1]):
          result,vel = setLOS(lon=lon, lat=lat, i_vel=i_vel, i_spe=spe, dim=dim, debug=debug)
          if result:
            normfactor = Rsl/constants.voxel_size/vel  #to normalise to disk shape...
            radiativeTransfer.sightlines[i,j] = normfactor*vel
            if vel==1:
              radiativeTransfer.intensitySpecies *= normfactor
            else:
              calculateRadiativeTransfer(scale=normfactor, dust=False)
            intensity.append(radiativeTransfer.intensitySpecies[0])
            radiativeTransfer.intensitySpecies = []
          else:
            intensity.append(np.zeros(radiativeTransfer.tempSpeciesEmissivity[0].shape[-1]))
        positionIntensitySpecies.append(np.asarray(intensity))
        intensity = []
        for dust in range(radiativeTransfer.tempDustEmissivity[0].shape[-1]):
          result,vel = setLOS(lon=lon, lat=lat, i_vel=i_vel, i_dust=dust, dim=dim, debug=debug)
          if result:
            normfactor = Rsl/constants.voxel_size/vel  #to normalise to disk shape...
            if radiativeTransfer.sightlines[i,j]<vel: radiativeTransfer.sightlines[i,j] = normfactor*vel
            if vel==1:
              radiativeTransfer.intensityDust *= normfactor
            else:
              calculateRadiativeTransfer(scale=normfactor, species=False)
            intensity.append(radiativeTransfer.intensityDust[0])
            radiativeTransfer.intensityDust = []
          else:
            intensity.append(np.zeros(radiativeTransfer.tempDustEmissivity[0].shape[-1]))
        positionIntensityDust.append(np.asarray(intensity))
        result = 2
      # if vel>1: print('\n\n', vel, '\n', (radiativeTransfer.epsilonSpecies>0).any(), '\n\n')

      position.append([lon,lat])
      
      if vox:
        normfactor = Rsl/constants.voxel_size/vox  #to normalise to disk shape...
        if sightlines[i,j]<vox: sightlines[i,j] = normfactor*vox
      
      # Integrate along LoS
      if vox==1:
        positionIntensitySpecies.append(emission[0] * normfactor)
        positionIntensityDust.append(emission[1] * normfactor)
      elif vox>1:
        intensitySpecies,intensityDust = calculateRadiativeTransfer(*emission, scale=normfactor)
        positionIntensitySpecies.append(intensitySpecies)
        positionIntensityDust.append(intensityDust)
        # intensitySpecies = []
        # intensityDust = []
      else:
        positionIntensitySpecies.append(np.zeros(radiativeTransfer.tempSpeciesEmissivity[0].shape[-1]))
        positionIntensityDust.append(np.zeros(radiativeTransfer.tempDustEmissivity[0].shape[-1]))
      
      # if terminal: radiativeTransfer.slTqdm.update()
      
      if len(np.shape(positionIntensitySpecies[-1]))>1:
        print('Error', np.shape(positionIntensitySpecies[-1]))
        input()
    
    # Save intensities for each latitude
    intensityMapSpecies.append(positionIntensitySpecies)
    intensityMapDust.append(positionIntensityDust)
  
  # Save map for velocity channel
  # VintensityMapSpecies.append(intensityMapSpecies)
  # VintensityMapDust.append(intensityMapDust)
  # Vpositions.append(position)

  # if verbose:
  #   print('Evaluating {} km/s HDU'.format(velocity))


  # if plotV:
  #   fig = plt.figure()
  #   ax = fig.add_subplot(111, projection='mollweide')
  #   cb = ax.scatter(np.array(position)[:,0], np.array(position)[:,1], c=np.array(intensityMapSpecies)[:,:,0,23].flatten(), s=64, marker='s')
  #   plt.ion()
  #   fig.colorbar(cb)
  #   ax.grid(True)
  #   ax.set_xticklabels([])
  #   ax.set_yticklabels([])
  #   plt.show(block=True)
  
  return (position,intensityMapSpecies,intensityMapDust,sightlines)

def setLOS(x=0, y=0, z=0, lon=0, lat=0, i_vox=[], i_vel=0, i_spe=None, i_dust=None, dim='xy', reverse=True, debug=False, verbose=False):
  '''
  The emission dimensions should be velocity x species x 2 x voxels. Axis 1 should be split for intensity and optical depth.
  The positions dimensions should be 3 x voxels.
  '''
  LOSvoxels = []
  zPosition = []

  scale = constants.voxel_size*constants.pc*100   #pc should be in cm

  # #print(radiativeTransfer.tempClumpEmission.shape)
  #
  # # This block matches the voxel positions to add the ensembles belonging to the same voxel
  # nrows, ncols = radiativeTransfer.tempClumpPosition.shape
  # dtype={'names':['f{}'.format(i) for i in range(ncols)], \
  #        'formats':ncols * [radiativeTransfer.tempClumpPosition.dtype]}
  # common,iCommonClump,iCommonInterclump = np.intersect1d(radiativeTransfer.tempClumpPosition.view(dtype), radiativeTransfer.tempInterclumpPosition.view(dtype), return_indices=True)
  #
  # # print(iCommonClump.max())
  #
  # # Intensity and optical depth have shape (voxel, wavelength)
  # gridIntensity = radiativeTransfer.tempInterclumpEmission[0,:,:]
  # gridIntensity[iCommonInterclump,:] = gridIntensity[iCommonInterclump,:] + radiativeTransfer.tempClumpEmission[0,iCommonClump,:]# emission[0,:,:]#
  #
  # gridOpticalDepth = radiativeTransfer.tempInterclumpEmission[1,:,:]# emission[1,:,:]#
  # gridOpticalDepth[iCommonInterclump,:] += radiativeTransfer.tempClumpEmission[1,iCommonClump,:]
  
  # gridIntensity -> rt.tempSpeciesEmissivity/tempDustEmissivity, gridOpticalDepth -> rt.tempSpeciesAbsorption/tempDustAbsorption

  # Specify the voxel positions relative to Earth
  xGrid,yGrid,zGrid = radiativeTransfer.voxelPositions[0].data[:,0]-constants.rGalEarth,radiativeTransfer.voxelPositions[0].data[:,1],radiativeTransfer.voxelPositions[0].data[:,2]#grid.getVoxelPositions()
  
  if verbose:
    print('Centered at x={}, y={}, z={}'.format(x,y,z))

  if dim=='spherical':
    # Set sightline position
    x1LoS = lon
    x2LoS = lat

    # Convert voxel positions to spherical
    radGrid = np.sqrt((xGrid-constants.rGalEarth)**2 + yGrid**2 + zGrid**2)
    lonGrid = np.arctan2(yGrid, -(xGrid-constants.rGalEarth))
    if lon<0: lonGrid[lonGrid>0] = lonGrid[lonGrid>0] - 2*np.pi
    if lon>0: lonGrid[lonGrid<0] = lonGrid[lonGrid<0] + 2*np.pi
    rPolar  = np.sqrt((xGrid-constants.rGalEarth)**2+yGrid**2)
    latGrid = np.arctan2(zGrid, rPolar)
    if lat<0: latGrid[latGrid>0] = latGrid[latGrid>0] - np.pi
    if lat>0: latGrid[latGrid<0] = latGrid[latGrid<0] + np.pi

    # Choose the voxels in the sightline
    #adjustments for voxel orientation
    scaling = np.sqrt(2)# 2# 
    width = scaling*constants.voxel_size*np.max([np.sin(np.abs(lonGrid-np.pi/4)), np.sin(np.abs(lonGrid+np.pi/4))], axis=0)
    height = scaling*constants.voxel_size*np.max([np.sin(np.abs(latGrid-np.pi/4)), np.sin(np.abs(latGrid+np.pi/4))], axis=0)
    #angular size of voxels
    dLon = np.arctan(width/radGrid)
    dLon[(lat>1.0)|(lat<-1.0)] = np.pi
    dLat = np.arctan(height/radGrid)
    # dLat[(lat>1.0)|(lat<-1.0)] = np.pi/2.
    iLoS = np.where((abs(lonGrid-x1LoS)<=dLon)&(abs(latGrid-x2LoS)<=dLat)&(i_vox))[0]
  
  elif 'disk' in dim:
    x1LoS = y
    x2LoS = z
    iLoS = np.where((zGrid==z)&(yGrid==y)&(i_vox))[0]
  
  elif ('x' in dim) and ('y' in dim):
    x1LoS = x
    x2LoS = y
    iLoS = np.where((xGrid==x)&(yGrid==y)&(i_vox))[0]
  
  elif ('x' in dim) and ('z' in dim):
    x1LoS = x
    x2LoS = z
    iLoS = np.where((xGrid==x)&(zGrid==z)&(i_vox))[0]
  
  elif ('z' in dim) and ('y' in dim):
    x1LoS = z
    x2LoS = y
    iLoS = np.where((zGrid==z)&(yGrid==y)&(i_vox))[0]

  else:
    print('\nPlease enter valid dimensions.\n')
    return False
  
  if verbose:
    print(iLoS)
  
  if iLoS.size==1:

    intensitySpecies = scale*radiativeTransfer.tempSpeciesEmissivity[0].data[iLoS,i_vel,:][0,:] / constants.pc/100
    intensityDust = scale*radiativeTransfer.tempDustEmissivity[0].data[iLoS,i_vel,:][0,:] / constants.pc/100
    
    # print('\n\n', (radiativeTransfer.tempSpeciesEmissivity[0].data[iLOS,i_vel,:]>0).any(), '\n\n')

    if debug:

      radiativeTransfer.allIndeces.append(iLoS)
      radiativeTransfer.allK.append(np.zeros(iLoS.size))
      radiativeTransfer.allE.append(np.zeros(iLoS.size))
      radiativeTransfer.allKstep.append(np.zeros(iLoS.size-1))
      radiativeTransfer.allEstep.append(np.zeros(iLoS.size-1))

      with warnings.catch_warnings():
        warnings.simplefilter('ignore')
        a = np.zeros(iLoS.size-1)
        b = np.zeros(iLoS.size-1)

      # ar = []
      # br = []
      # ai = []
      # bi = []

      # for idx in range(len(i)-1):
      #   ar = np.append(ar, np.array(list(Ereal(a[idx,:]))))
      #   br = np.append(br, np.array(list(Ereal(b[idx,:]))))
      #   ai = np.append(ai, np.array(list(Eimag(a[idx,:]))))s
      #   bi = np.append(bi, np.array(list(Eimag(b[idx,:]))))

      radiativeTransfer.allAreal.append(a)
      radiativeTransfer.allBreal.append(b)
      # radiativeTransfer.allAimag.append(ai)
      # radiativeTransfer.allBimag.append(bi)

      radiativeTransfer.allIntensity.append(radiativeTransfer.intensity)

    vel = radiativeTransfer.voxelVelocities[0].data[iLoS]
    
    return 1,vel.size,[intensitySpecies,intensityDust]
  
  elif iLoS.size>1:
    
    if 'spherical' == dim:
      x3LoS = radGrid[iLoS]
    elif not 'x' in dim:
      x3LoS = xGrid[iLoS]
    elif not 'y' in dim:
      x3LoS = yGrid[iLoS]
    elif not 'z' in dim:
      x3LoS = zGrid[iLoS]
      
    if reverse:
      #order the voxels in the line-of-sight according to increasing radial distance
      i = np.argsort(x3LoS)
      iLoS_ordered = iLoS[i]
    else:
      #order the voxels in the line-of-sight according to decreasing radial distance
      i = np.argsort(x3LoS)[::-1]
      iLoS_ordered = iLoS[i]
    
    # print('\n\n', iLoS, i, iLoS_ordered, i_vel, (radiativeTransfer.tempSpeciesEmissivity[0].data[iLoS_ordered,i_vel,:]>0).any(), '\n\n')
    # print(radiativeTransfer.i_vox)
    # print(radiativeTransfer.x3LoS)
    # input()
    
    if (i_spe==None)&(i_dust==None):
      epsilonSpecies = c.copy(radiativeTransfer.tempSpeciesEmissivity[0].data[iLoS_ordered,i_vel,:]) / constants.pc/100
      epsilonStepSpecies = (epsilonSpecies[1:]-epsilonSpecies[:-1])/(scale)
      kappaSpecies = c.copy(radiativeTransfer.tempSpeciesAbsorption[0].data[iLoS_ordered,i_vel,:]) / constants.pc/100
      kappaStepSpecies = (kappaSpecies[1:]-kappaSpecies[:-1])/(scale)
      
      epsilonDust = c.copy(radiativeTransfer.tempDustEmissivity[0].data[iLoS_ordered,i_vel,:]) / constants.pc/100
      epsilonStepDust = (epsilonDust[1:]-epsilonDust[:-1])/(scale)
      kappaDust = c.copy(radiativeTransfer.tempDustAbsorption[0].data[iLoS_ordered,i_vel,:]) / constants.pc/100
      kappaStepDust = (kappaDust[1:]-kappaDust[:-1])/(scale)
    
    # elif i_spe!=None:
    #   radiativeTransfer.epsilonSpecies = (radiativeTransfer.tempSpeciesEmissivity[0].data[iLoS_ordered,i_vel,i_spe]) / constants.pc/100
    #   radiativeTransfer.epsilonStepSpecies = (radiativeTransfer.epsilonSpecies[1:]-radiativeTransfer.epsilonSpecies[:-1])/(scale)
    #   radiativeTransfer.kappaSpecies = (radiativeTransfer.tempSpeciesAbsorption[0].data[iLoS_ordered,i_vel,i_spe]) / constants.pc/100
    #   radiativeTransfer.kappaStepSpecies = (radiativeTransfer.kappaSpecies[1:]-radiativeTransfer.kappaSpecies[:-1])/(scale)
    # elif i_dust!=None:
    #   print(iLoS_ordered)
    #   radiativeTransfer.epsilonDust = (radiativeTransfer.tempDustEmissivity[0].data[iLoS_ordered,i_vel,i_dust]) / constants.pc/100
    #   radiativeTransfer.epsilonStepDust = (radiativeTransfer.epsilonDust[1:]-radiativeTransfer.epsilonDust[:-1])/(scale)
    #   radiativeTransfer.kappaDust = (radiativeTransfer.tempDustAbsorption[0].data[iLoS_ordered,i_vel,i_dust]) / constants.pc/100
    #   radiativeTransfer.kappaStepDust = (radiativeTransfer.kappaDust[1:]-radiativeTransfer.kappaDust[:-1])/(scale)

    if debug:

      radiativeTransfer.allIndeces.append(iLoS)
      
      radiativeTransfer.allKSpecies.append(np.array(radiativeTransfer.kappaSpecies))
      radiativeTransfer.allESpecies.append(np.array(radiativeTransfer.epsilonSpecies))
      radiativeTransfer.allKstepSpecies.append(np.array(radiativeTransfer.kappaStepSpecies))
      radiativeTransfer.allEstepSpecies.append(np.array(radiativeTransfer.epsilonStepSpecies))
      
      radiativeTransfer.allKDust.append(np.array(radiativeTransfer.kappaDust))
      radiativeTransfer.allEDust.append(np.array(radiativeTransfer.epsilonDust))
      radiativeTransfer.allKstepDust.append(np.array(radiativeTransfer.kappaStepDust))
      radiativeTransfer.allEstepDust.append(np.array(radiativeTransfer.epsilonStepDust))

      with warnings.catch_warnings():
        warnings.simplefilter('ignore')
        aSpecies = (radiativeTransfer.kappaSpecies[:-1,:]/np.sqrt(2*radiativeTransfer.kappaStepSpecies.astype(np.complex)))
        bSpecies = ((radiativeTransfer.kappaSpecies[:-1,:]+radiativeTransfer.kappaStepSpecies*scale)/np.sqrt(2*radiativeTransfer.kappaStepSpecies.astype(np.complex)))
        aDust = (radiativeTransfer.kappaDust[:-1,:]/np.sqrt(2*radiativeTransfer.kappaStepDust.astype(np.complex)))
        bDust = ((radiativeTransfer.kappaDust[:-1,:]+radiativeTransfer.kappaStepDust*scale)/np.sqrt(2*radiativeTransfer.kappaStepDust.astype(np.complex)))

      radiativeTransfer.allArealSpecies.append(aSpecies)
      radiativeTransfer.allBrealSpecies.append(bSpecies)
      radiativeTransfer.allArealDust.append(aDust)
      radiativeTransfer.allBrealDust.append(bDust)
      # radiativeTransfer.allAimagSpecies.append(aiSpecies)
      # radiativeTransfer.allBimagSpecies.append(biSpecies)
      # radiativeTransfer.allAimagDust.append(aiDust)
      # radiativeTransfer.allBimagDust.append(biDust)
  
  else:
    return 0,0,[]

  vel = radiativeTransfer.voxelVelocities[0].data[iLoS]
  
  if verbose:
    print('voxels:', i)
  
  return 2,vel.size,(epsilonSpecies,epsilonStepSpecies,kappaSpecies,kappaStepSpecies,epsilonDust,epsilonStepDust,kappaDust,kappaStepDust)

def calculateRadiativeTransfer(epsilonSpecies, epsilonStepSpecies, kappaSpecies, kappaStepSpecies, epsilonDust, epsilonStepDust, kappaDust, kappaStepDust, scale=1, backgroundI=0., species=True, dust=True, verbose=False, test=True):

  # Modify data and initialise observed intensity array
  # if not dust:
  #   radiativeTransfer.epsilonSpecies = radiativeTransfer.epsilonSpecies.reshape((-1,1))
  #   radiativeTransfer.epsilonStepSpecies = radiativeTransfer.epsilonStepSpecies.reshape((-1,1))
  #   radiativeTransfer.kappaSpecies = radiativeTransfer.kappaSpecies.reshape((-1,1))
  #   radiativeTransfer.kappaStepSpecies = radiativeTransfer.kappaStepSpecies.reshape((-1,1))
  #   intensitySpecies = np.array([backgroundI])
  #   nSteps = radiativeTransfer.epsilonStepSpecies.shape[0]
  # elif not species:
  #   radiativeTransfer.epsilonDust = radiativeTransfer.epsilonDust.reshape((-1,1))
  #   radiativeTransfer.epsilonStepDust = radiativeTransfer.epsilonStepDust.reshape((-1,1))
  #   radiativeTransfer.kappaDust = radiativeTransfer.kappaDust.reshape((-1,1))
  #   radiativeTransfer.kappaStepDust = radiativeTransfer.kappaStepDust.reshape((-1,1))
  #   intensityDust = np.array([backgroundI])
  #   nSteps = radiativeTransfer.epsilonStepDust.shape[0]
  # else:
  intensitySpecies = np.full(kappaSpecies.shape[1], backgroundI)
  intensityDust = np.full(kappaDust.shape[1], backgroundI)
  nSteps = epsilonStepSpecies.shape[0]
  
  # Adjust according to the average voxel depth
  if species:
    epsilonSpecies /= scale
    epsilonStepSpecies /= scale**2
    kappaSpecies /= scale
    kappaStepSpecies /= scale**2
  if dust:
    epsilonDust /= scale
    epsilonStepDust /= scale**2
    kappaDust /= scale
    kappaStepDust /= scale**2

  scale = scale*constants.voxel_size*constants.pc*100  #pc should be in cm
  # print(scale)
  
  np.set_printoptions(threshold=100000)
  warnings.filterwarnings('error')
  
  # # Boolean indeces to separate how the intensity is calculated
  # k0 = (radiativeTransfer.kappaStep==0)&(abs(radiativeTransfer.kappa[:-1]*constants.resolution)<10**-10)
  # kg = radiativeTransfer.kappa[:-1]>10**3*abs(radiativeTransfer.kappaStep)*constants.resolution
  # kE = ~(k0|kg)
  # kEg = ~(k0|kg)&(radiativeTransfer.kappaStep>0)
  # kEl = ~(k0|kg)&(radiativeTransfer.kappaStep<0)
  
  # Calculate the variables needed to utilise the E tilde tables
  with warnings.catch_warnings():
    warnings.simplefilter('ignore')
    if species:
      aSpecies = (kappaSpecies[:-1,:]/np.sqrt(2*kappaStepSpecies.astype(np.complex)))
      bSpecies = ((kappaSpecies[:-1,:]+kappaStepSpecies*scale)/np.sqrt(2*kappaStepSpecies.astype(np.complex)))
    if dust:
      aDust = (kappaDust[:-1,:]/np.sqrt(2*kappaStepDust.astype(np.complex)))
      bDust = ((kappaDust[:-1,:]+kappaStepDust*scale)/np.sqrt(2*kappaStepDust.astype(np.complex)))
  
  if verbose: print(kappaSpecies.shape[0])
  
  for i in range(nSteps):
  
    # Determine which form of integration the species transitions require
    if species:
      k0Species = (kappaStepSpecies[i,:]==0)&(abs(kappaSpecies[:-1,:][i,:]*scale)<10**-10)
      try:
        kgSpecies = ~k0Species&(kappaSpecies[:-1,:][i,:]>10**3*abs(kappaStepSpecies[i,:])*scale)
      except RuntimeWarning:
        # print(radiativeTransfer.x1LoS, radiativeTransfer.x2LoS)
        # print(radiativeTransfer.x3LoS)
        # print(radiativeTransfer.kappaStepSpecies[i,:])
        # print(radiativeTransfer.kappaSpecies[i,:])
        input()
        kgSpecies = ~k0&(kappaSpecies[:-1,:][i,:]>10**3*abs(kappaStepSpecies[i,:])*scale)
      kESpecies = ~(k0Species|kgSpecies)
      kEgSpecies = kESpecies&(kappaStepSpecies[i,:]>0)
      kElSpecies = kESpecies&(kappaStepSpecies[i,:]<0)

    # Determine which form of integration the dust continuum requires
    if dust:
      k0Dust = (kappaStepDust[i,:]==0)&(abs(kappaDust[:-1,:][i,:]*scale)<10**-10)
      try:
        kgDust = ~k0Dust&(kappaDust[:-1,:][i,:]>10**3*abs(kappaStepDust[i,:])*scale)
      except RuntimeWarning:
        # print(radiativeTransfer.x1LoS, radiativeTransfer.x2LoS)
        # print(radiativeTransfer.x3LoS)
        # print(radiativeTransfer.kappaStepDust[i,:])
        # print(radiativeTransfer.kappaDust[i,:])
        input()
        kgDust = ~k0&(kappaDust[:-1,:][i,:]>10**3*abs(kappaStepDust[i,:])*scale)
      kEDust = ~(k0Dust|kgDust)
      kEgDust = kEDust&(kappaStepDust[i,:]>0)
      kElDust = kEDust&(kappaStepDust[i,:]<0)
    
    # << Compute radiative transfer for the species transitions >>
    if species:
      if k0Species.any():
        intensitySpecies[k0Species] = intensitySpecies[k0Species] + epsilonSpecies[:-1,:][i,:][k0Species]*scale+0.5*epsilonStepSpecies[i,:][k0Species]*scale**2
    
      elif kgSpecies.any():
        try:
          intensitySpecies[kgSpecies] = np.exp(-kappaSpecies[i,:][kgSpecies]*scale) * intensitySpecies[kgSpecies] + \
                             ((epsilonSpecies[i,:][kgSpecies]*kappaSpecies[i,:][kgSpecies]+epsilonStepSpecies[i,:][kgSpecies]*(kappaSpecies[i,:][kgSpecies]*scale-1))/(kappaSpecies[i,:][kgSpecies]**2.)) - \
                             np.exp(-kappaSpecies[i,:][kgSpecies]*scale)*((epsilonSpecies[i,:][kgSpecies]*kappaSpecies[i,:][kgSpecies]-epsilonStepSpecies[i,:][kgSpecies])/(kappaSpecies[i,:][kgSpecies]**2.))
        except RuntimeWarning:
          print('ERROR: Integration option 2')
          intensitySpecies[kgSpecies] = np.exp(-kappaSpecies[i,:][kgSpecies]*scale) * (intensitySpecies[kgSpecies] + \
                             ((epsilonSpecies[i,:][kgSpecies]*kappaSpecies[i,:][kgSpecies]+epsilonStepSpecies[i,:][kgSpecies]*(kappaSpecies[i,:][kgSpecies]*scale-1))/(kappaSpecies[i,:][kgSpecies]**2.)) - \
                             ((epsilonSpecies[i,:][kgSpecies]*kappaSpecies[i,:][kgSpecies]-epsilonStepSpecies[i,:][kgSpecies])/(kappaSpecies[i,:][kgSpecies]**2.)))
  
      elif kEgSpecies.any():
        if verbose: print('\na, b:\n', aSpecies[i,:,:], '\n', bSpecies[i,:,:])
        aE = np.array(list(Ereal(aSpecies[i,:][kEgSpecies])))
        bE = np.array(list(Ereal(bSpecies[i,:][kEgSpecies])))
        intensitySpecies[kEgSpecies] = (epsilonStepSpecies[i,:][kEgSpecies]/kappaStepSpecies[i,:][kEgSpecies]*(1-np.exp(-kappaSpecies[:-1,:][i,:][kEgSpecies]*scale-kappaStepSpecies[i,:][kEgSpecies]/2.*scale**2.)) - \
                          (epsilonSpecies[:-1,:][i,:][kEgSpecies]*kappaStepSpecies[i,:][kEgSpecies]-epsilonStepSpecies[i,:][kEgSpecies]*kappaSpecies[:-1,:][i,:][kEgSpecies])/kappaStepSpecies[i,:][kEgSpecies] * \
                          np.sqrt(np.pi/(2.*abs(kappaStepSpecies[i][kEgSpecies].astype(np.complex)))) * \
                          (np.exp(aSpecies[i,:][kEgSpecies]**2.-bSpecies[i,:][kEgSpecies]**2.)*aE-bE) + \
                          intensitySpecies[kEgSpecies]*np.exp(-kappaSpecies[:-1][i,:][kEgSpecies]*scale-kappaStepSpecies[i,:][kEgSpecies]/2.*scale**2.)).real
  
      elif kElSpecies.any():
        if verbose: print('\na, b:\n', aSpecies[i,:,:], '\n', bSpecies[i,:])
        aE = np.array(list(Eimag(aSpecies[i,:][kElSpecies])))
        bE = np.array(list(Eimag(bSpecies[i,:][kElSpecies])))
        try:
          intensitySpecies[kElSpecies] = (epsilonStepSpecies[i,:][kElSpecies]/kappaStepSpecies[i,:][kElSpecies]*(1-np.exp(-kappaSpecies[:-1,:][i,:][kElSpecies]*scale-kappaStepSpecies[i,:][kElSpecies]/2.*scale**2.))\
                            -(epsilonSpecies[:-1,:][i,:][kElSpecies]*kappaStepSpecies[i,:][kElSpecies]-epsilonStepSpecies[i,:][kElSpecies]*kappaSpecies[:-1,:][i,:][kElSpecies])/kappaStepSpecies[i,:][kElSpecies] * \
                            np.sqrt(np.pi/(2.*abs(kappaStepSpecies[i][kElSpecies].astype(np.complex))))* \
                            (np.exp(aSpecies[i,:][kElSpecies]**2.-bSpecies[i,:][kElSpecies]**2.)*aE-bE) + \
                            intensitySpecies[kElSpecies]*np.exp(-kappaSpecies[:-1,:][i,:][kElSpecies]*scale-kappaStepSpecies[i,:][kElSpecies]/2.*scale**2.)).real
        except RuntimeWarning:
          print(kappaSpecies[-1,:][kElSpecies])
          print(kappaStepSpecies[-1,kElSpecies])
          print(aSpecies[i,kElSpecies]**2)
          print(bSpecies[i,kElSpecies]**2)
          input()
    
    # << Compute radiative transfer for the dust continuum >>
    if dust:
      if k0Dust.any():
        intensityDust[k0Dust] = intensityDust[k0Dust] + epsilonDust[:-1,:][i,:][k0Dust]*scale+0.5*epsilonStepDust[i,:][k0Dust]*scale**2
    
      elif kgDust.any():
        try:
          intensityDust[kgDust] = np.exp(-kappaDust[i,:][kgDust]*scale) * intensityDust[kgDust] + \
                             ((epsilonDust[i,:][kgDust]*kappaDust[i,:][kgDust]+epsilonStepDust[i,:][kgDust]*(kappaDust[i,:][kgDust]*scale-1))/(kappaDust[i,:][kgDust]**2.)) - \
                             np.exp(-kappaDust[i,:][kgDust]*scale)*((epsilonDust[i,:][kgDust]*kappaDust[i,:][kgDust]-epsilonStepDust[i,:][kgDust])/(kappaDust[i,:][kgDust]**2.))
        except RuntimeWarning:
          print('ERROR: Integration option 2')
          intensityDust[kgDust] = np.exp(-kappaDust[i,:][kgDust]*scale) * (intensityDust[kgDust] + \
                             ((epsilonDust[i,:][kgDust]*kappaDust[i,:][kgDust]+epsilonStepDust[i,:][kgDust]*(kappaDust[i,:][kgDust]*scale-1))/(kappaDust[i,:][kgDust]**2.)) - \
                             ((epsilonDust[i,:][kgDust]*kappaDust[i,:][kgDust]-epsilonStepDust[i,:][kgDust])/(kappaDust[i,:][kgDust]**2.)))
  
      elif kEgDust.any():
        if verbose: print('\na, b:\n', aDust[i,:,:], '\n', bDust[i,:,:])
        aE = np.array(list(Ereal(aDust[i,:][kEgDust])))
        bE = np.array(list(Ereal(bDust[i,:][kEgDust])))
        intensityDust[kEgDust] = (epsilonStepDust[i,:][kEgDust]/kappaStepDust[i,:][kEgDust]*(1-np.exp(-kappaDust[:-1,:][i,:][kEgDust]*scale-kappaStepDust[i,:][kEgDust]/2.*scale**2.)) - \
                          (epsilonDust[:-1,:][i,:][kEgDust]*kappaStepDust[i,:][kEgDust]-epsilonStepDust[i,:][kEgDust]*kappaDust[:-1,:][i,:][kEgDust])/kappaStepDust[i,:][kEgDust] * \
                          np.sqrt(np.pi/(2.*abs(kappaStepDust[i][kEgDust].astype(np.complex)))) * \
                          (np.exp(aDust[i,:][kEgDust]**2.-bDust[i,:][kEgDust]**2.)*aE-bE) + \
                          intensityDust[kEgDust]*np.exp(-kappaDust[:-1][i,:][kEgDust]*scale-kappaStepDust[i,:][kEgDust]/2.*scale**2.)).real
  
      elif kElDust.any():
        if verbose: print('\na, b:\n', aDust[i,:,:], '\n', bDust[i,:])
        aE = np.array(list(Eimag(aDust[i,:][kElDust])))
        bE = np.array(list(Eimag(bDust[i,:][kElDust])))
        try:
          intensityDust[kElDust] = (epsilonStepDust[i,:][kElDust]/kappaStepDust[i,:][kElDust]*(1-np.exp(-kappaDust[:-1,:][i,:][kElDust]*scale-kappaStepDust[i,:][kElDust]/2.*scale**2.))\
                            -(epsilonDust[:-1,:][i,:][kElDust]*kappaStepDust[i,:][kElDust]-epsilonStepDust[i,:][kElDust]*kappaDust[:-1,:][i,:][kElDust])/kappaStepDust[i,:][kElDust] * \
                            np.sqrt(np.pi/(2.*abs(kappaStepDust[i][kElDust].astype(np.complex))))* \
                            (np.exp(aDust[i,:][kElDust]**2.-bDust[i,:][kElDust]**2.)*aE-bE) + \
                            intensityDust[kElDust]*np.exp(-kappaDust[:-1,:][i,:][kElDust]*scale-kappaStepDust[i,:][kElDust]/2.*scale**2.)).real
        except RuntimeWarning:
          print(kappaDust[-1,:][kElDust])
          print(kappaStepDust[-1,kElDust])
          print(aDust[i,kElDust]**2)
          print(bDust[i,kElDust]**2)
          input()
  
  if verbose:
    print('Species intensity shape:', np.shape(intensitySpecies))
    print('Dust intensity shape:', np.shape(intensityDust))

  # if species:
  #   radiativeTransfer.intensitySpecies = intensitySpecies
  #   radiativeTransfer.allIntensitySpecies.append(intensitySpecies)
  # if dust:
  #   radiativeTransfer.intensityDust = intensityDust
  #   radiativeTransfer.allIntensityDust.append(intensityDust)

  return (intensitySpecies, intensityDust)

def Ereal(x, verbose=False):
  
  if verbose: print('E real input:', x)

  '''
  ## This is no longer needed if I have a complex matrix!
  ir = (x.imag)
  
  if (x.imag==0).any(): x = x.real
  # x should be a real number. remove imaginary party '0j' which
  # prevents ordering
  '''

  Er = np.zeros_like(x)

  il = x<0.01
  ig = x>8.0
  ib = ~(ig|il)

  if il.any():

    if verbose:
      print('x less than grid')

    Er[il] = (2*x[il]/np.sqrt(np.pi)).astype(np.complex)
  
  if ig.any():

    if verbose:
      print('x greater than grid')

    Er[ig] = (1/(np.sqrt(np.pi) * x[ig])).astype(np.complex)
  
  if ib.any():

    if verbose:
      print('x interpolated')

    Er[ib] = radiativeTransfer.eTildeReal(x[ib].real).astype(np.complex)

  return Er

def Eimag(x, verbose=False):
  
  if verbose: print('E imaginary input:', x)

  Ei = np.zeros_like(x)

  im = x==abs(x)*1j
  il = (abs(x)<0.01)&~im
  ig = (abs(x)>8.0)&~im
  ib = ~(ig|il|im)

  # Force x to be a positive real value
  x = abs(x)
  
  if im.any():

    if verbose:
      print('Maser case')

    # MASER case: treat with linear approximation
    Ei[im] = 1 + 2*x[im]/np.sqrt(np.pi)

  if il.any():

    if verbose:
      print('x less than grid')
    
    Ei[il] = (1 - 2*x[il]/np.sqrt(np.pi)).astype(np.complex)

  if ig.any():
    
    if verbose:
      print('x greater than grid')
    
    Ei[ig] = (1/(np.sqrt(np.pi) * x[ig])).astype(np.complex)

  if ib.any():
    
    if verbose:
      print('x interpolated')
    
    Ei[ib] = radiativeTransfer.eTildeImaginary(x[ib]).astype(np.complex)

  return Ei

def debugEmission(sl,i):

  k = radiativeTransfer.allK[sl][:,i]
  e = radiativeTransfer.allE[sl][:,i]
  kstep = radiativeTransfer.allKstep[sl][:,i]
  estep = radiativeTransfer.allEstep[sl][:,i]
  a = radiativeTransfer.allAreal[sl][:,i]
  b = radiativeTransfer.allBreal[sl][:,i]
  intensity = radiativeTransfer.allIntensity[sl][i]

  import pandas as pd

  d = {'kappa':k[:-1],'dk':kstep,'epsilon':e[:-1],'de':estep,'A':a,'B':b,'I':intensity,'dI':0}
  df = pd.DataFrame(data=d)

  scale = constants.resolution*constants.pc*100  #pc should be in cm

  intensity = 0

  for i in range(a.size):
    
    if kstep[i]==0 and (k[i]*scale)<10**-10:
      print('option 0')
      df.loc[[i],['dI']] = intensity + e[i]*scale + 0.5*estep[i]*scale**2
    
    elif k[i]>10**3*abs(kstep[i])*scale:
      print('option 1')
      df.loc[[i],['dI']] = intensity*np.exp(-k[i]*scale) + \
                           (e[i]*k[i]+estep[i]*(k[i]*scale-1))/k[i]**2 - \
                           (e[i]*k[i]-estep[i])/k[i]**2*np.exp(-k[i]*scale)
    
    elif kstep[i]>0:
      print('option 2')
      df.loc[[i],['dI']] = intensity*np.exp(-k[i]*scale-0.5*kstep[i]*scale**2) + \
                           estep[i]/kstep[i]*(1-np.exp(-k[i]*scale-0.5*kstep[i]*scale**2)) - \
                           (e[i]*kstep[i]-estep[i]*k[i])/kstep[i]*np.sqrt(np.pi/2/abs(kstep[i]))*\
                           (np.exp(a[i]**2-b[i]**2)*radiativeTransfer.Ereal(np.array([a[i]]), verbose=True)-radiativeTransfer.Ereal(np.array([b[i]]), verbose=True)).real
    
    elif kstep[i]<0:
      print('option 3')
      df.loc[[i],['dI']] = intensity*np.exp(-k[i]*scale-0.5*kstep[i]*scale**2) + \
                           estep[i]/kstep[i]*(1-np.exp(-k[i]*scale-0.5*kstep[i]*scale**2)) - \
                           (e[i]*kstep[i]-estep[i]*k[i])/kstep[i]*np.sqrt(np.pi/2/abs(kstep[i]))*\
                           (np.exp(a[i]**2-b[i]**2)*radiativeTransfer.Eimag(np.array([a[i]]), verbose=True)-radiativeTransfer.Eimag(np.array([b[i]]), verbose=True)).real

    intensity = c.copy(df['dI'][i])

  return df

# if __name__=='__main__':
#
#   print('spawned process')
#
#   directory = r'C:\Users\cyani\projects\pdr\KT3_history\MilkyWay\r200_cm1-1_d1_uv10'
#
#   constants.velocityRange = np.linspace(-300, 300, 500)
#
#   radiativeTransfer.voxelPositions = fits.open(directory+'/voxel_position.fits', mode='denywrite')
#   radiativeTransfer.voxelVelocities = fits.open(directory+'/voxel_velocity.fits', mode='denywrite')
#   radiativeTransfer.tempSpeciesEmissivity = fits.open(directory+'/species_emissivity.fits', mode='denywrite')# in K/pc
#   radiativeTransfer.tempSpeciesAbsorption = fits.open(directory+'/species_absorption.fits', mode='denywrite')# in 1/pc
#   radiativeTransfer.tempDustEmissivity = fits.open(directory+'/dust_emissivity.fits', mode='denywrite')# in K/pc
#   radiativeTransfer.tempDustAbsorption = fits.open(directory+'/dust_absorption.fits', mode='denywrite')# in 1/pc
#
#   multiprocessCalculation(slRange=[(-np.pi,np.pi), (-np.pi/2,np.pi/2)], nsl=[50,25], dim='spherical', multiprocessing=2)
  
  # multiprocessing.freeze_support()

# jit_module(nopython=False)
