import numpy as np
from numba import jit_module
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

def calculateObservation(directory='', dim='xy', sl=[50,50], terminal=True, plotV=False, debug=False, verbose=False):

  if debug:
    sl = [5,5]
    
  # print('Load data')

  voxelPositions = fits.open(directory+'/voxel_position.fits')
  radiativeTransfer.voxelVelocities = fits.open(directory+'/voxel_velocity.fits')
  radiativeTransfer.tempSpeciesEmissivity = fits.open(directory+'/species_emissivity_clump.fits')# in K/pc
  radiativeTransfer.tempSpeciesAbsorption = fits.open(directory+'/species_absorption_clump.fits')# in 1/pc
  radiativeTransfer.tempDustEmissivity = fits.open(directory+'/dust_emissivity_clump.fits')# in K/pc
  radiativeTransfer.tempDustAbsorption = fits.open(directory+'/dust_absorption_clump.fits')# in 1/pc
  
  # print('Data loaded :-)')

  if verbose:
    print(voxelPositions.shape, clumpEmission.shape, interclumpEmission.shape)

  xPositions,yPositions,zPositions = voxelPositions[0].data[:,0],voxelPositions[0].data[:,1],voxelPositions[0].data[:,2]
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
    longrid = np.linspace(-np.pi, np.pi, num=sl[0])
    latgrid = np.linspace(-np.pi/2, np.pi/2, num=sl[1])
    # grid = np.meshgrid(lon, lat)
    # grid = np.array([grid[0].flatten(), grid[1].flatten()])

    vTqdm = tqdm(total=constants.velocityRange.size, desc='Observing velocity', miniters=1, dynamic_ncols=True)
    if terminal: slTqdm = tqdm(total=longrid.size*latgrid.size, desc='Sightline', miniters=1, dynamic_ncols=True)

    VintensityMapSpecies = []
    VintensityMapDust = []
    Vpositions    = []
    vmin = constants.velocityRange.max()
    vmax = constants.velocityRange.min()

    radiativeTransfer.sightlines = np.zeros((longrid.size, latgrid.size))

    if debug: velocityRange = [0]
    else: velocityRange = constants.velocityRange

    for i_vel,velocity in enumerate(constants.velocityRange):

      # Find the voxels that exist at the observing velocity
      iV =  radiativeTransfer.tempSpeciesEmissivity[0].data[:,i_vel,:].any(1)#(clumpVelocity==velocity)
      #if np.any(iCV): input('{}'.format(iCV))

      # Update velocity progress bar
      vTqdm.update()
      # This is a test to arrange the sightlines in a sphere around the observer
      #if not (np.any(iCV) or np.any(iIV)): continue

      # Reset sightline progress bar
      if terminal: slTqdm.reset()

      # Initialise the intensities and map
      position = []
      intensityMapSpecies = []
      intensityMapDust = []

      # Get indeces of voxels at the correct observing velocity
      # iClumpV = np.where(iCV)
      # iInterclumpV = np.where(iIV)

      if iV.any()==False: continue

      if velocity<vmin: vmin = velocity
      if velocity>vmax: vmax = velocity

      # The voxel positions can be any of the voxels
      # radiativeTransfer.tempSpeciesEmissivity = tempSpeciesEmissivity#[iV,i,:] / constants.pc/100
      # radiativeTransfer.tempSpeciesAbsorption = tempSpeciesAbsorption#[iV,i,:] / constants.pc/100
      # radiativeTransfer.tempDustEmissivity = tempDustEmissivity#[iV,i,:] / constants.pc/100
      # radiativeTransfer.tempDustAbsorption = tempDustAbsorption#[iV,i,:] / constants.pc/100
      radiativeTransfer.tempPosition = voxelPositions[0].data[iV,:]
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
          if verbose:
            print('\n',thetaC,'\n')
          if abs(lat)<thetaC:
            Rsl = Rslh/np.cos(lat)
          else:
            Rsl = constants.hd/np.sin(abs(lat))

          # Isolate voxels in LoS
          result,vel = setLOS(lon=lon, lat=lat, i_vel=i, directory=directory, dim=dim, debug=debug)
          # print(vel)

          position.append([lon,lat])
          
          # Integrate along LoS
          if result:
            normfactor = Rsl/constants.voxel_size/vel  #to normalise to disk shape...
            if radiativeTransfer.sightlines[i,j]<vel: radiativeTransfer.sightlines[i,j] = normfactor*vel
            if vel==1:
              radiativeTransfer.intensitySpecies *= normfactor
              radiativeTransfer.intensityDust *= normfactor
            else:
              calculateRadiativeTransfer(scale=normfactor)
            positionIntensitySpecies.append(radiativeTransfer.intensitySpecies)
            positionIntensityDust.append(radiativeTransfer.intensityDust)
            radiativeTransfer.intensitySpecies = []
            radiativeTransfer.intensityDust = []
          else:
            positionIntensitySpecies.append(np.zeros(radiativeTransfer.tempSpeciesEmissivity[0].shape[-1]))
            positionIntensityDust.append(np.zeros(radiativeTransfer.tempDustEmissivity[0].shape[-1]))
          
          if terminal: slTqdm.update()
          
          if len(np.shape(positionIntensitySpecies[-1]))>1:
            print('Error', np.shape(positionIntensitySpecies[-1]))
            input()
        
        # Save intensities for each latitude
        intensityMapSpecies.append(positionIntensitySpecies)
        intensityMapDust.append(positionIntensityDust)
    
      VintensityMapSpecies.append(intensityMapSpecies)
      VintensityMapDust.append(intensityMapDust)
      Vpositions.append(position)

      if verbose:
        print('Evaluating {} km/s HDU'.format(velocity))


      if plotV:
        fig = plt.figure()
        ax = fig.add_subplot(111, projection='mollweide')
        cb = ax.scatter(np.array(position)[:,0], np.array(position)[:,1], c=np.array(intensityMapSpecies)[:,:,0,23].flatten(), s=64, marker='s')
        plt.ion()
        fig.colorbar(cb)
        ax.grid(True)
        ax.set_xticklabels([])
        ax.set_yticklabels([])
        plt.show(block=True)
    
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
      IntensityHDUSpecies.header['SPECIES'] = fits.open(directory+'/species_emissivity_clump.fits')[0].header['SPECIES']

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
      IntensityHDUDust.header['DUST'] = fits.open(directory+'/dust_emissivity_clump.fits')[0].header['DUST']

      hdul.append(PositionHDU)
      hdul.append(IntensityHDUSpecies)
      hdul.append(IntensityHDUDust)

      hdul.writeto(directory+'/channel_intensity.fits', overwrite=True)

      print('Intensity map written successfully :-)')

    voxelPositions.close()
    radiativeTransfer.voxelVelocities.close()
    radiativeTransfer.tempSpeciesEmissivity.close()
    radiativeTransfer.tempSpeciesAbsorption.close()
    radiativeTransfer.tempDustEmissivity.close()
    radiativeTransfer.tempDustAbsorption.close()
  
  return

def setLOS(emission=0, positions=0, x=0, y=0, z=0, lon=0, lat=0, i_vel=0, tempSpeciesEmissivity=[], tempSpeciesAbsorption=[], tempDustEmissivity=[], tempDustAbsorption=[], directory='', dim='xy', reverse=True, debug=False, verbose=False):
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
  xGrid,yGrid,zGrid = radiativeTransfer.tempPosition[:,0]-constants.rGalEarth,radiativeTransfer.tempPosition[:,1],radiativeTransfer.tempPosition[:,2]#grid.getVoxelPositions()
  
  if verbose:
    print('Centered at x={}, y={}, z={}'.format(x,y,z))

  if dim=='spherical':
    # Set sightline position
    radiativeTransfer.x1LoS = lon
    radiativeTransfer.x2LoS = lat

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
    iLOS = np.where((abs(lonGrid-radiativeTransfer.x1LoS)<=dLon)&(abs(latGrid-radiativeTransfer.x2LoS)<=dLat))[0]
  
  elif 'disk' in dim:
    radiativeTransfer.x1LoS = y
    radiativeTransfer.x2LoS = z
    iLOS = np.where((zGrid==z)&(yGrid==y))[0]
  
  elif ('x' in dim) and ('y' in dim):
    radiativeTransfer.x1LoS = x
    radiativeTransfer.x2LoS = y
    iLOS = np.where((xGrid==x)&(yGrid==y))[0]
  
  elif ('x' in dim) and ('z' in dim):
    radiativeTransfer.x1LoS = x
    radiativeTransfer.x2LoS = z
    iLOS = np.where((xGrid==x)&(zGrid==z))[0]
  
  elif ('z' in dim) and ('y' in dim):
    radiativeTransfer.x1LoS = z
    radiativeTransfer.x2LoS = y
    iLOS = np.where((zGrid==z)&(yGrid==y))[0]

  else:
    print('\nPlease enter valid dimensions.\n')
    return False
  
  if verbose:
    print(iLOS)
  
  if iLOS.size==1:
    
    # tempSpeciesEmissivity = fits.open(directory+'/species_emissivity_clump.fits')# / constants.pc/100
    # tempDustEmissivity = fits.open(directory+'/dust_emissivity_clump.fits')# / constants.pc/100

    radiativeTransfer.intensitySpecies = scale*radiativeTransfer.tempSpeciesEmissivity[0].data[iLOS,i_vel,:][0,:] / constants.pc/100
    radiativeTransfer.intensityDust = scale*radiativeTransfer.tempDustEmissivity[0].data[iLOS,i_vel,:][0,:] / constants.pc/100
    # print('Shape along LoS', radiativeTransfer.intensity.shape)
    
    # tempSpeciesEmissivity.close()
    # tempDustEmissivity.close()

    if debug:

      radiativeTransfer.allIndeces.append(iLOS)
      radiativeTransfer.allK.append(np.zeros(iLOS.size))
      radiativeTransfer.allE.append(np.zeros(iLOS.size))
      radiativeTransfer.allKstep.append(np.zeros(iLOS.size-1))
      radiativeTransfer.allEstep.append(np.zeros(iLOS.size-1))

      with warnings.catch_warnings():
        warnings.simplefilter('ignore')
        a = np.zeros(iLOS.size-1)
        b = np.zeros(iLOS.size-1)

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
  
  elif iLOS.size>1:
    
    if verbose:
      print('Intensity:', intensity, scale)
      print('Optical depth:', tau, scale)
    
    if 'spherical' == dim:
      radiativeTransfer.x3LoS = radGrid[iLOS]
    elif not 'x' in dim:
      radiativeTransfer.x3LoS = xGrid[iLOS]
    elif not 'y' in dim:
      radiativeTransfer.x3LoS = yGrid[iLOS]
    elif not 'z' in dim:
      radiativeTransfer.x3LoS = zGrid[iLOS]
      
    if reverse:
      i = np.argsort(radiativeTransfer.x3LoS)
    else:
      i = np.argsort(radiativeTransfer.x3LoS)[::-1]
    # if radiativeTransfer.x2LoS<-0.8: print(len(i))
    
    # tempSpeciesEmissivity = fits.open(directory+'/species_emissivity_clump.fits')# / constants.pc/100
    # tempSpeciesAbsorption = fits.open(directory+'/species_absorption_clump.fits')# / constants.pc/100
    # tempDustEmissivity = fits.open(directory+'/dust_emissivity_clump.fits')# / constants.pc/100
    # tempDustAbsorption = fits.open(directory+'/dust_absorption_clump.fits')# / constants.pc/100
    
    radiativeTransfer.epsilonSpecies = radiativeTransfer.tempSpeciesEmissivity[0].data[i,i_vel,:] / constants.pc/100
    radiativeTransfer.epsilonStepSpecies = (radiativeTransfer.epsilonSpecies[1:]-radiativeTransfer.epsilonSpecies[:-1])/(scale)
    radiativeTransfer.kappaSpecies = radiativeTransfer.tempSpeciesAbsorption[0].data[i,i_vel,:] / constants.pc/100
    radiativeTransfer.kappaStepSpecies = (radiativeTransfer.kappaSpecies[1:]-radiativeTransfer.kappaSpecies[:-1])/(scale)
    
    radiativeTransfer.epsilonDust = radiativeTransfer.tempDustEmissivity[0].data[i,i_vel,:] / constants.pc/100
    radiativeTransfer.epsilonStepDust = (radiativeTransfer.epsilonDust[1:]-radiativeTransfer.epsilonDust[:-1])/(scale)
    radiativeTransfer.kappaDust = radiativeTransfer.tempDustAbsorption[0].data[i,i_vel,:] / constants.pc/100
    radiativeTransfer.kappaStepDust = (radiativeTransfer.kappaDust[1:]-radiativeTransfer.kappaDust[:-1])/(scale)
    
    # tempSpeciesEmissivity.close()
    # tempSpeciesAbsorption.close()
    # tempDustEmissivity.close()
    # tempDustAbsorption.close()

    if debug:

      radiativeTransfer.allIndeces.append(iLOS)
      
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
    return False,0

  vel = radiativeTransfer.voxelVelocities[0].data[iLOS]
  
  if verbose:
    print('voxels:', i)
  
  return True,vel.size

def calculateRadiativeTransfer(scale=1, backgroundI=0., verbose=False, test=True):

  # if np.size(radiativeTransfer.intensity)>0: return
  
  intensitySpecies = np.full(radiativeTransfer.kappaSpecies.shape[1], backgroundI)
  intensityDust = np.full(radiativeTransfer.kappaDust.shape[1], backgroundI)
  
  # print(scale)
  # Rescale the LoS emission
  # scale = 1
  radiativeTransfer.kappaSpecies /= scale
  radiativeTransfer.kappaStepSpecies /= scale**2
  radiativeTransfer.epsilonSpecies /= scale
  radiativeTransfer.epsilonStepSpecies /= scale**2
  radiativeTransfer.kappaDust /= scale
  radiativeTransfer.kappaStepDust /= scale**2
  radiativeTransfer.epsilonDust /= scale
  radiativeTransfer.epsilonStepDust /= scale**2

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
  
  # Variables for utilising the E tilde tables
  with warnings.catch_warnings():
    warnings.simplefilter('ignore')
    aSpecies = (radiativeTransfer.kappaSpecies[:-1,:]/np.sqrt(2*radiativeTransfer.kappaStepSpecies.astype(np.complex)))
    bSpecies = ((radiativeTransfer.kappaSpecies[:-1,:]+radiativeTransfer.kappaStepSpecies*scale)/np.sqrt(2*radiativeTransfer.kappaStepSpecies.astype(np.complex)))
    aDust = (radiativeTransfer.kappaDust[:-1,:]/np.sqrt(2*radiativeTransfer.kappaStepDust.astype(np.complex)))
    bDust = ((radiativeTransfer.kappaDust[:-1,:]+radiativeTransfer.kappaStepDust*scale)/np.sqrt(2*radiativeTransfer.kappaStepDust.astype(np.complex)))
  
  if verbose: print(radiativeTransfer.kappaSpecies.shape[0])
  
  for i in range(len(radiativeTransfer.kappaStepSpecies)):
  
    # Determine which form of integration the species transitions require
    k0Species = (radiativeTransfer.kappaStepSpecies[i,:]==0)&(abs(radiativeTransfer.kappaSpecies[:-1,:][i,:]*scale)<10**-10)
    try:
      kgSpecies = ~k0Species&(radiativeTransfer.kappaSpecies[:-1,:][i,:]>10**3*abs(radiativeTransfer.kappaStepSpecies[i,:])*scale)
    except RuntimeWarning:
      print(radiativeTransfer.x1LoS, radiativeTransfer.x2LoS)
      print(radiativeTransfer.x3LoS)
      print(radiativeTransfer.kappaStepSpecies[i,:])
      print(radiativeTransfer.kappaSpecies[i,:])
      input()
      kgSpecies = ~k0&(radiativeTransfer.kappa[:-1,:][i,:]>10**3*abs(radiativeTransfer.kappaStep[i,:])*scale)
    kESpecies = ~(k0Species|kgSpecies)
    kEgSpecies = kESpecies&(radiativeTransfer.kappaStepSpecies[i,:]>0)
    kElSpecies = kESpecies&(radiativeTransfer.kappaStepSpecies[i,:]<0)
  
    # Determine which form of integration the dust continuum requires
    k0Dust = (radiativeTransfer.kappaStepDust[i,:]==0)&(abs(radiativeTransfer.kappaDust[:-1,:][i,:]*scale)<10**-10)
    try:
      kgDust = ~k0Dust&(radiativeTransfer.kappaDust[:-1,:][i,:]>10**3*abs(radiativeTransfer.kappaStepDust[i,:])*scale)
    except RuntimeWarning:
      print(radiativeTransfer.x1LoS, radiativeTransfer.x2LoS)
      print(radiativeTransfer.x3LoS)
      print(radiativeTransfer.kappaStepDust[i,:])
      print(radiativeTransfer.kappaDust[i,:])
      input()
      kgDust = ~k0&(radiativeTransfer.kappa[:-1,:][i,:]>10**3*abs(radiativeTransfer.kappaStep[i,:])*scale)
    kEDust = ~(k0Dust|kgDust)
    kEgDust = kEDust&(radiativeTransfer.kappaStepDust[i,:]>0)
    kElDust = kEDust&(radiativeTransfer.kappaStepDust[i,:]<0)
  
    if verbose:
      print('\ni_k0 Species\n', k0Species)
      print('i_kg Species\n', kgSpecies)
      print('i_kEg Species\n', kEgSpecies)
      print('i_kEl Species\n', kElSpecies)
      print('\ni_k0 Dust\n', k0Dust)
      print('i_kg Dust\n', kgDust)
      print('i_kEg Dust\n', kEgDust)
      print('i_kEl Dust\n', kElDust)
      
      print('\nSpecies kappa, kappa step:\n', radiativeTransfer.kappaSpecies[i], '\n', radiativeTransfer.kappaStepSpecies[i])
      print('\nSpecies epsilon, epsilon step, scale\n', radiativeTransfer.epsilonSpecies[i], radiativeTransfer.epsilonStepSpecies[i], scale)
      print('\nDust kappa, kappa step:\n', radiativeTransfer.kappaDust[i], '\n', radiativeTransfer.kappaStepDust[i])
      print('\nDust epsilon, epsilon step, scale\n', radiativeTransfer.epsilonDust[i], radiativeTransfer.epsilonStepDust[i], scale)
    
    # << Compute radiative transfer for the species transitions >>
  
    if k0Species.any():
      intensitySpecies[k0Species] = intensitySpecies[k0Species] + radiativeTransfer.epsilonSpecies[:-1,:][i,:][k0Species]*scale+0.5*radiativeTransfer.epsilonStepSpecies[i,:][k0Species]*scale**2
  
    elif kgSpecies.any():
      try:
        intensitySpecies[kgSpecies] = np.exp(-radiativeTransfer.kappaSpecies[i,:][kgSpecies]*scale) * intensitySpecies[kgSpecies] + \
                           ((radiativeTransfer.epsilonSpecies[i,:][kgSpecies]*radiativeTransfer.kappaSpecies[i,:][kgSpecies]+radiativeTransfer.epsilonStepSpecies[i,:][kgSpecies]*(radiativeTransfer.kappaSpecies[i,:][kgSpecies]*scale-1))/(radiativeTransfer.kappaSpecies[i,:][kgSpecies]**2.)) - \
                           np.exp(-radiativeTransfer.kappaSpecies[i,:][kgSpecies]*scale)*((radiativeTransfer.epsilonSpecies[i,:][kgSpecies]*radiativeTransfer.kappaSpecies[i,:][kgSpecies]-radiativeTransfer.epsilonStepSpecies[i,:][kgSpecies])/(radiativeTransfer.kappaSpecies[i,:][kgSpecies]**2.))
      except RuntimeWarning:
        print('ERROR: Integration option 2')
        intensitySpecies[kgSpecies] = np.exp(-radiativeTransfer.kappaSpecies[i,:][kgSpecies]*scale) * (intensitySpecies[kgSpecies] + \
                           ((radiativeTransfer.epsilonSpecies[i,:][kgSpecies]*radiativeTransfer.kappaSpecies[i,:][kgSpecies]+radiativeTransfer.epsilonStepSpecies[i,:][kgSpecies]*(radiativeTransfer.kappaSpecies[i,:][kgSpecies]*scale-1))/(radiativeTransfer.kappaSpecies[i,:][kgSpecies]**2.)) - \
                           ((radiativeTransfer.epsilonSpecies[i,:][kgSpecies]*radiativeTransfer.kappaSpecies[i,:][kgSpecies]-radiativeTransfer.epsilonStepSpecies[i,:][kgSpecies])/(radiativeTransfer.kappaSpecies[i,:][kgSpecies]**2.)))

    elif kEgSpecies.any():
      if verbose: print('\na, b:\n', aSpecies[i,:,:], '\n', bSpecies[i,:,:])
      aE = np.array(list(Ereal(aSpecies[i,:][kEgSpecies])))
      bE = np.array(list(Ereal(bSpecies[i,:][kEgSpecies])))
      intensitySpecies[kEgSpecies] = (radiativeTransfer.epsilonStepSpecies[i,:][kEgSpecies]/radiativeTransfer.kappaStepSpecies[i,:][kEgSpecies]*(1-np.exp(-radiativeTransfer.kappaSpecies[:-1,:][i,:][kEgSpecies]*scale-radiativeTransfer.kappaStepSpecies[i,:][kEgSpecies]/2.*scale**2.)) - \
                        (radiativeTransfer.epsilonSpecies[:-1,:][i,:][kEgSpecies]*radiativeTransfer.kappaStepSpecies[i,:][kEgSpecies]-radiativeTransfer.epsilonStepSpecies[i,:][kEgSpecies]*radiativeTransfer.kappaSpecies[:-1,:][i,:][kEgSpecies])/radiativeTransfer.kappaStepSpecies[i,:][kEgSpecies] * \
                        np.sqrt(np.pi/(2.*abs(radiativeTransfer.kappaStepSpecies[i][kEgSpecies].astype(np.complex)))) * \
                        (np.exp(aSpecies[i,:][kEgSpecies]**2.-bSpecies[i,:][kEgSpecies]**2.)*aE-bE) + \
                        intensitySpecies[kEgSpecies]*np.exp(-radiativeTransfer.kappaSpecies[:-1][i,:][kEgSpecies]*scale-radiativeTransfer.kappaStepSpecies[i,:][kEgSpecies]/2.*scale**2.)).real

    elif kElSpecies.any():
      if verbose: print('\na, b:\n', aSpecies[i,:,:], '\n', bSpecies[i,:])
      aE = np.array(list(Eimag(aSpecies[i,:][kElSpecies])))
      bE = np.array(list(Eimag(bSpecies[i,:][kElSpecies])))
      try:
        intensitySpecies[kElSpecies] = (radiativeTransfer.epsilonStepSpecies[i,:][kElSpecies]/radiativeTransfer.kappaStepSpecies[i,:][kElSpecies]*(1-np.exp(-radiativeTransfer.kappaSpecies[:-1,:][i,:][kElSpecies]*scale-radiativeTransfer.kappaStepSpecies[i,:][kElSpecies]/2.*scale**2.))\
                          -(radiativeTransfer.epsilonSpecies[:-1,:][i,:][kElSpecies]*radiativeTransfer.kappaStepSpecies[i,:][kElSpecies]-radiativeTransfer.epsilonStepSpecies[i,:][kElSpecies]*radiativeTransfer.kappaSpecies[:-1,:][i,:][kElSpecies])/radiativeTransfer.kappaStepSpecies[i,:][kElSpecies] * \
                          np.sqrt(np.pi/(2.*abs(radiativeTransfer.kappaStepSpecies[i][kElSpecies].astype(np.complex))))* \
                          (np.exp(aSpecies[i,:][kElSpecies]**2.-bSpecies[i,:][kElSpecies]**2.)*aE-bE) + \
                          intensitySpecies[kElSpecies]*np.exp(-radiativeTransfer.kappaSpecies[:-1,:][i,:][kElSpecies]*scale-radiativeTransfer.kappaStepSpecies[i,:][kElSpecies]/2.*scale**2.)).real
      except RuntimeWarning:
        print(radiativeTransfer.kappaSpecies[-1,:][kElSpecies])
        print(radiativeTransfer.kappaStepSpecies[-1,kElSpecies])
        print(aSpecies[i,kElSpecies]**2)
        print(bSpecies[i,kElSpecies]**2)
        input()
    
    # << Compute radiative transfer for the dust continuum >>
  
    if k0Dust.any():
      intensityDust[k0Dust] = intensityDust[k0Dust] + radiativeTransfer.epsilonDust[:-1,:][i,:][k0Dust]*scale+0.5*radiativeTransfer.epsilonStepDust[i,:][k0Dust]*scale**2
  
    elif kgDust.any():
      try:
        intensityDust[kgDust] = np.exp(-radiativeTransfer.kappaDust[i,:][kgDust]*scale) * intensityDust[kgDust] + \
                           ((radiativeTransfer.epsilonDust[i,:][kgDust]*radiativeTransfer.kappaDust[i,:][kgDust]+radiativeTransfer.epsilonStepDust[i,:][kgDust]*(radiativeTransfer.kappaDust[i,:][kgDust]*scale-1))/(radiativeTransfer.kappaDust[i,:][kgDust]**2.)) - \
                           np.exp(-radiativeTransfer.kappaDust[i,:][kgDust]*scale)*((radiativeTransfer.epsilonDust[i,:][kgDust]*radiativeTransfer.kappaDust[i,:][kgDust]-radiativeTransfer.epsilonStepDust[i,:][kgDust])/(radiativeTransfer.kappaDust[i,:][kgDust]**2.))
      except RuntimeWarning:
        print('ERROR: Integration option 2')
        intensityDust[kgDust] = np.exp(-radiativeTransfer.kappaDust[i,:][kgDust]*scale) * (intensityDust[kgDust] + \
                           ((radiativeTransfer.epsilonDust[i,:][kgDust]*radiativeTransfer.kappaDust[i,:][kgDust]+radiativeTransfer.epsilonStepDust[i,:][kgDust]*(radiativeTransfer.kappaDust[i,:][kgDust]*scale-1))/(radiativeTransfer.kappaDust[i,:][kgDust]**2.)) - \
                           ((radiativeTransfer.epsilonDust[i,:][kgDust]*radiativeTransfer.kappaDust[i,:][kgDust]-radiativeTransfer.epsilonStepDust[i,:][kgDust])/(radiativeTransfer.kappaDust[i,:][kgDust]**2.)))

    elif kEgDust.any():
      if verbose: print('\na, b:\n', aDust[i,:,:], '\n', bDust[i,:,:])
      aE = np.array(list(Ereal(aDust[i,:][kEgDust])))
      bE = np.array(list(Ereal(bDust[i,:][kEgDust])))
      intensityDust[kEgDust] = (radiativeTransfer.epsilonStepDust[i,:][kEgDust]/radiativeTransfer.kappaStepDust[i,:][kEgDust]*(1-np.exp(-radiativeTransfer.kappaDust[:-1,:][i,:][kEgDust]*scale-radiativeTransfer.kappaStepDust[i,:][kEgDust]/2.*scale**2.)) - \
                        (radiativeTransfer.epsilonDust[:-1,:][i,:][kEgDust]*radiativeTransfer.kappaStepDust[i,:][kEgDust]-radiativeTransfer.epsilonStepDust[i,:][kEgDust]*radiativeTransfer.kappaDust[:-1,:][i,:][kEgDust])/radiativeTransfer.kappaStepDust[i,:][kEgDust] * \
                        np.sqrt(np.pi/(2.*abs(radiativeTransfer.kappaStepDust[i][kEgDust].astype(np.complex)))) * \
                        (np.exp(aDust[i,:][kEgDust]**2.-bDust[i,:][kEgDust]**2.)*aE-bE) + \
                        intensityDust[kEgDust]*np.exp(-radiativeTransfer.kappaDust[:-1][i,:][kEgDust]*scale-radiativeTransfer.kappaStepDust[i,:][kEgDust]/2.*scale**2.)).real

    elif kElDust.any():
      if verbose: print('\na, b:\n', aDust[i,:,:], '\n', bDust[i,:])
      aE = np.array(list(Eimag(aDust[i,:][kElDust])))
      bE = np.array(list(Eimag(bDust[i,:][kElDust])))
      try:
        intensityDust[kElDust] = (radiativeTransfer.epsilonStepDust[i,:][kElDust]/radiativeTransfer.kappaStepDust[i,:][kElDust]*(1-np.exp(-radiativeTransfer.kappaDust[:-1,:][i,:][kElDust]*scale-radiativeTransfer.kappaStepDust[i,:][kElDust]/2.*scale**2.))\
                          -(radiativeTransfer.epsilonDust[:-1,:][i,:][kElDust]*radiativeTransfer.kappaStepDust[i,:][kElDust]-radiativeTransfer.epsilonStepDust[i,:][kElDust]*radiativeTransfer.kappaDust[:-1,:][i,:][kElDust])/radiativeTransfer.kappaStepDust[i,:][kElDust] * \
                          np.sqrt(np.pi/(2.*abs(radiativeTransfer.kappaStepDust[i][kElDust].astype(np.complex))))* \
                          (np.exp(aDust[i,:][kElDust]**2.-bDust[i,:][kElDust]**2.)*aE-bE) + \
                          intensityDust[kElDust]*np.exp(-radiativeTransfer.kappaDust[:-1,:][i,:][kElDust]*scale-radiativeTransfer.kappaStepDust[i,:][kElDust]/2.*scale**2.)).real
      except RuntimeWarning:
        print(radiativeTransfer.kappaDust[-1,:][kElDust])
        print(radiativeTransfer.kappaStepDust[-1,kElDust])
        print(aDust[i,kElDust]**2)
        print(bDust[i,kElDust]**2)
        input()
  
  if verbose:
    print('Species intensity shape:', np.shape(intensitySpecies))
    print('Dust intensity shape:', np.shape(intensityDust))

  radiativeTransfer.intensitySpecies = intensitySpecies
  radiativeTransfer.intensityDust = intensityDust
  radiativeTransfer.allIntensitySpecies.append(intensitySpecies)
  radiativeTransfer.allIntensityDust.append(intensityDust)

  return

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

# jit_module(nopython=False)
