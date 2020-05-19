import numpy as np
from numba import jit_module
import matplotlib as mpl
mpl.use('Qt4Agg')
import matplotlib.pyplot as plt
import copy as c
import cmath
import scipy.interpolate as interpolate
from astropy.io import fits
import importlib as il
from tqdm import tqdm

import warnings

import constants
import interpolations
import observations
import radiativeTransfer

def eTildeReal(file='Ereal.dat'):
  eReal = np.genfromtxt(constants.GRIDPATH+file, names=['x', 'Ereal'])
  return (eReal['x'],eReal['Ereal'])

def eTildeImaginary(file='Eimag.dat'):
  eImaginary = np.genfromtxt(constants.GRIDPATH+file, names=['x', 'Eimaginary'])
  return (eImaginary['x'],eImaginary['Eimaginary'])

def calculateObservation(directory='', dim='xy', sl=[50,50], terminal=True, plotV=False, debug=False, verbose=False):

  if debug:
    sl = [5,5]

  voxelPositions = fits.open(constants.HISTORYPATH+constants.directory+directory+'/voxel_position.fits')[0].data
  radiativeTransfer.voxelVelocities = fits.open(constants.HISTORYPATH+constants.directory+directory+'/voxel_velocity.fits')[0].data #test
  clumpIntensity = fits.open(constants.HISTORYPATH+constants.directory+directory+'/intensity_clump.fits')[0].data
  clumpOpticalDepth = fits.open(constants.HISTORYPATH+constants.directory+directory+'/opticalDepth_clump.fits')[0].data
  clumpVelocity = fits.open(constants.HISTORYPATH+constants.directory+directory+'/voxel_clump_velocity.fits')[0].data
  interclumpIntensity = fits.open(constants.HISTORYPATH+constants.directory+directory+'/intensity_interclump.fits')[0].data
  interclumpOpticalDepth = fits.open(constants.HISTORYPATH+constants.directory+directory+'/opticalDepth_interclump.fits')[0].data
  interclumpVelocity = fits.open(constants.HISTORYPATH+constants.directory+directory+'/voxel_interclump_velocity.fits')[0].data

  # Create emission arrays; the current shape is (emission type, voxel, wavelength, velocity)
  clumpEmission = np.array([clumpIntensity, clumpOpticalDepth])
  interclumpEmission = np.array([interclumpIntensity, interclumpOpticalDepth])

  if verbose:
    print(voxelPositions.shape, clumpEmission.shape, interclumpEmission.shape)

  xPositions,yPositions,zPositions = voxelPositions[:,0],voxelPositions[:,1],voxelPositions[:,2]
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

    VintensityMap = []
    Vpositions    = []
    vmin = constants.velocityRange.max()
    vmax = constants.velocityRange.min()

    radiativeTransfer.sightlines = np.zeros((longrid.size, latgrid.size))

    if debug: velocityRange = [0]
    else: velocityRange = constants.velocityRange

    print(sl)

    for i,velocity in enumerate(constants.velocityRange):

      # Find the voxels that exist at the observing velocity
      iCV = (clumpVelocity==velocity)
      iIV = (interclumpVelocity==velocity)
      #if np.any(iCV): input('{}'.format(iCV))

      # Update velocity progress bar
      vTqdm.update()

      if dim=='spherical':
        # This is a test to arrange the sightlines in a sphere around the observer
        #if not (np.any(iCV) or np.any(iIV)): continue

        # Reset sightline progress bar
        if terminal: slTqdm.reset()

        # Initialise the intensities and map
        position = []
        intensityMap = []

        # Get indeces of voxels at the correct observing velocity
        iClumpV = np.where(iCV)
        iInterclumpV = np.where(iIV)

        if (iCV.any()==False) | (iIV.any()==False): continue

        if velocity<vmin: vmin = velocity
        if velocity>vmax: vmax = velocity

        # The voxel positions can be any of the voxels
        radiativeTransfer.tempClumpEmission = clumpEmission[:,iClumpV[0],iClumpV[1],:]
        radiativeTransfer.tempClumpPosition = voxelPositions[iClumpV[0],:]
        radiativeTransfer.tempClumpVelocity = clumpVelocity[iClumpV[0],:]
        radiativeTransfer.tempInterclumpEmission = interclumpEmission[:,iInterclumpV[0],iInterclumpV[1],:]
        radiativeTransfer.tempInterclumpPosition = voxelPositions[iInterclumpV[0],:]
        radiativeTransfer.tempInterclumpVelocity = interclumpVelocity[iInterclumpV[0],:]

        velTest = []

        for j,lat in enumerate(latgrid):
          positionIntensity = []
          for i,lon in enumerate(longrid):
            result,vel = setLOS(lon=lon, lat=lat, dim=dim, debug=debug)
            if radiativeTransfer.sightlines[i,j]<vel: radiativeTransfer.sightlines[i,j] = vel
            position.append([lon,lat])
            if result:
              calculateRadiativeTransfer()
              positionIntensity.append(radiativeTransfer.intensity)
              radiativeTransfer.intensity = []
            else:
              positionIntensity.append(np.zeros(clumpEmission.shape[-1]))
            if terminal: slTqdm.update()
            if len(np.shape(positionIntensity[-1]))>1:
              print('Error', np.shape(positionIntensity[-1]))
              input()
          intensityMap.append(positionIntensity)
      
        VintensityMap.append(intensityMap)
        Vpositions.append(position)

      if verbose:
        print('Evaluating {} km/s HDU'.format(velocity))


      if plotV:
        fig = plt.figure()
        ax = fig.add_subplot(111, projection='mollweide')
        cb = ax.scatter(np.array(position)[:,0], np.array(position)[:,1], c=np.array(intensityMap)[:,:,0,23].flatten(), s=64, marker='s')
        plt.ion()
        fig.colorbar(cb)
        ax.grid(True)
        ax.set_xticklabels([])
        ax.set_yticklabels([])
        plt.show(block=True)
      
    # Convert to numpy arrays
    Vpositions = np.array(Vpositions[0])
    VintensityMap = np.array(VintensityMap)

    if verbose:
      print('Map position shape:')
      if Vpositions.ndim>1:
        print(Vpositions.shape)
      else:
        for p in VmapPositions: print(p.shape)

    if verbose:
      print('Map intensity shape:')
      if VintensityMap.ndim>1:
        print(VintensityMap.shape)
      else:
        for Vintensity in intensityMap: print(intensity.shape)

    if not debug:
      # Create HDUs for the map position and intensity and add the velocity in the headers
      PositionHDU = fits.ImageHDU(Vpositions)
      print(VintensityMap.shape, np.shape(VintensityMap[0,0,0]))
      radiativeTransfer.intensity = VintensityMap
      IntensityHDU = fits.ImageHDU(VintensityMap)
      # input('Integrated intensity shape: {}'.format(VintensityMap.shape))

      PositionHDU.header['TYPE1'] = 'Angle'
      PositionHDU.header['TYPE2'] = 'Position'
      PositionHDU.header['DIREC'] = 'Radial'

      IntensityHDU.header['BUNIT'] = 'K'
      IntensityHDU.header['CTYPE1'] = 'Wavelength'
      IntensityHDU.header['CUNIT1'] = 'm'
      IntensityHDU.header['CRVAL1'] = 'N/A'
      IntensityHDU.header['CDELT1'] = 'N/A'
      IntensityHDU.header['CRPIX1'] = 'N/A'
      IntensityHDU.header['CTYPE2'] = 'GLON'
      IntensityHDU.header['CUNIT2'] = 'rad'
      IntensityHDU.header['CRVAL2'] = 0.
      IntensityHDU.header['CDELT2'] = 2*np.pi/(IntensityHDU.header['NAXIS2']-1)
      IntensityHDU.header['CRPIX2'] = (IntensityHDU.header['NAXIS2']-1)/2.
      IntensityHDU.header['CTYPE3'] = 'GLAT'
      IntensityHDU.header['CUNIT3'] = 'rad'
      IntensityHDU.header['CRVAL3'] = 0.
      IntensityHDU.header['CDELT3'] = np.pi/(IntensityHDU.header['NAXIS3']-1)
      IntensityHDU.header['CRPIX3'] = (IntensityHDU.header['NAXIS3']-1)/2.
      IntensityHDU.header['CTYPE4'] = 'Velocity'
      IntensityHDU.header['CUNIT4'] = 'km/s'
      IntensityHDU.header['CRVAL4'] = (vmax+vmin)/2.
      IntensityHDU.header['CDELT4'] = (vmax-vmin)/(IntensityHDU.header['NAXIS4']-1)
      IntensityHDU.header['CRPIX4'] = (IntensityHDU.header['NAXIS4'])/2.
      IntensityHDU.header['DIREC'] = 'Radial'

      hdul.append(PositionHDU)
      hdul.append(IntensityHDU)

      hdul.writeto(constants.HISTORYPATH+constants.directory+directory+'/integrated_intensity.fits', overwrite=True)
  
  elif dim=='xy':
    # This is an observation from outside the galaxy, face-on

    # Initialise the intensities and map
    position = []
    intensityMap = []
    
    Array = np.unique([xPositions,yPositions], axis=1).T

    if verbose:
      print('\nx\n{}\n\ny\n{}\n'.format(Array[:,0],Array[:,1]))

    for x,y in Array:
      #for y in np.unique(yArray):
        result = orientation.setLOS(emission=clumpEmission+interclumpEmission, positions=voxelPositions, x=x, y=y, dim=dim)
        if result:
          position.append([x,y])
          intensity = orientation.calculateRadiativeTransfer(velocity)
          intensityMap.append(intensity)
  
    mapPositions = np.array(position)
    intensityMap = np.array(intensityMap)
  
  elif dim=='xz':
    # This is an observation from outside the galaxy, edge-on

    # Initialise the intensities and map
    position = []
    intensityMap = []

    xArray,zArray = np.unique([xPositions,zPositions], axis=1)
    for x in np.unique(xArray):
      for z in np.unique(zArray):
        orientation.setLOS(emission=clumpEmission+interclumpEmission, positions=voxelPositions, x=x, z=z, dim=dim)
        position.append([x,z])
        intensity = orientation.calculateRadiativeTransfer(velocity)
        intensityMap.append(intensity)
  
    mapPositions = np.array(position)
    intensityMap = np.array(intensityMap)
  
  elif dim=='yz':
    # This is an observation from outside the galaxy, edge-on

    # Initialise the intensities and map
    position = []
    intensityMap = []

    yArray,zArray = np.unique([yPositions,zPositions], axis=1)
    for y in np.unique(yArray):
      for z in np.unique(zArray):
        orientation.setLOS(emission=clumpEmission+interclumpEmission, positions=voxelPositions, y=y, z=z, dim=dim)
        position.append([y,z])
        intensity = orientation.calculateRadiativeTransfer(velocity)
        intensityMap.append(intensity)
  
    mapPositions = np.array(position)
    intensityMap = np.array(intensityMap)
  
  return

def setLOS(emission=0, positions=0, x=0, y=0, z=0, lon=0, lat=0, dim='xy', reverse=False, debug=False, verbose=False):
  '''
  The emission dimensions should be velocity x species x 2 x voxels. Axis 1 should be split for intensity and optical depth.
  The positions dimensions should be 3 x voxels.
  '''
  LOSvoxels = []
  zPosition = []
  epsilon = []
  kappa = []

  scale = constants.resolution*constants.pc

  #print(radiativeTransfer.tempClumpEmission.shape)

  # This block matches the voxel positions to add the ensembles belonging to the same voxel
  nrows, ncols = radiativeTransfer.tempClumpPosition.shape
  dtype={'names':['f{}'.format(i) for i in range(ncols)], \
         'formats':ncols * [radiativeTransfer.tempClumpPosition.dtype]}
  common,iCommonClump,iCommonInterclump = np.intersect1d(radiativeTransfer.tempClumpPosition.view(dtype), radiativeTransfer.tempInterclumpPosition.view(dtype), return_indices=True)
  
  # print(iCommonClump.max())

  # Intensity and optical depth have shape (voxel, wavelength)
  gridIntensity = radiativeTransfer.tempInterclumpEmission[0,:,:]
  gridIntensity[iCommonInterclump,:] = gridIntensity[iCommonInterclump,:] + radiativeTransfer.tempClumpEmission[0,iCommonClump,:]# emission[0,:,:]#

  gridOpticalDepth = radiativeTransfer.tempInterclumpEmission[1,:,:]# emission[1,:,:]#
  gridOpticalDepth[iCommonInterclump,:] += radiativeTransfer.tempClumpEmission[1,iCommonClump,:]

  # Specify the voxel positions relative to Earth
  xGrid,yGrid,zGrid = radiativeTransfer.tempInterclumpPosition[:,0]-constants.rGalEarth,radiativeTransfer.tempInterclumpPosition[:,1],radiativeTransfer.tempInterclumpPosition[:,2]#grid.getVoxelPositions()
  
  if verbose:
    print('Centered at x={}, y={}, z={}'.format(x,y,z))

  if dim=='spherical':
    # Set sightline position
    radiativeTransfer.x1LoS = lon
    radiativeTransfer.x2LoS = lat

    # Convert voxel positions to spherical
    radGrid = np.sqrt((xGrid-constants.rGalEarth)**2 + yGrid**2 + zGrid**2)
    lonGrid = np.arctan2(yGrid, -(xGrid-constants.rGalEarth))
    if lon<0: lonGrid[lonGrid>0] = lonGrid[lonGrid>0] + 2*np.pi
    if lon>0: lonGrid[lonGrid<0] = lonGrid[lonGrid<0] - 2*np.pi
    rPolar  = np.sqrt((xGrid-constants.rGalEarth)**2+yGrid**2)
    latGrid = np.arctan2(zGrid, rPolar)

    # Choose the voxels in the sightline
    #adjustments for voxel orientation
    width = np.sqrt(2)*constants.resolution*np.max([np.sin(np.abs(lonGrid-np.pi/4)), np.sin(np.abs(lonGrid+np.pi/4))], axis=0)
    height = np.sqrt(2)*constants.resolution*np.max([np.sin(np.abs(latGrid-np.pi/4)), np.sin(np.abs(latGrid+np.pi/4))], axis=0)
    #angular size of voxels
    dLon = np.arctan(width/2/radGrid)
    dLat = np.arctan(height/2/radGrid)
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
  #self.__losVoxels = grid.allVoxels()[iLOS]
  
  if iLOS.size==1:

    radiativeTransfer.intensity = gridIntensity[iLOS,:][0]
    # print('Shape along LoS', radiativeTransfer.intensity.shape)
  
  elif iLOS.size>1:
    #for i in iLOS:
      #LOSvoxels.append(voxels[i])
    intensity = gridIntensity[iLOS,:]
    tau = gridOpticalDepth[iLOS,:]
      #if z:
      #print(factor, scale, intensity)
    epsilon = (intensity/scale)
    kappa = (tau/scale)
      #else:
      #  epsilon.append(intensity/(constants.resolution))
      #  kappa.append(tau/(constants.resolution))
    if verbose:
      print('Intensity:', intensity, scale)
      print('Optical depth:', tau, scale)
    #epsilonStep = (epsilon[1:]-epsilon[:-1])/(scale)
    #kappaStep = (kappa[1:]-kappa[:-1])/(scale)
    #orientation.losVoxels = LOSvoxels
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
    radiativeTransfer.epsilon = np.array(epsilon, dtype=np.float)[i]
    radiativeTransfer.epsilonStep = (radiativeTransfer.epsilon[1:]-radiativeTransfer.epsilon[:-1])/(scale)
    radiativeTransfer.kappa = np.array(kappa, dtype=np.float)[i]
    radiativeTransfer.kappaStep = (radiativeTransfer.kappa[1:]-radiativeTransfer.kappa[:-1])/(scale)

    if debug:

      radiativeTransfer.allIndeces.append(iLOS)
      radiativeTransfer.allK.append(np.array(radiativeTransfer.kappa))
      radiativeTransfer.allE.append(np.array(radiativeTransfer.epsilon))
      radiativeTransfer.allKstep.append(np.array(radiativeTransfer.kappaStep))
      radiativeTransfer.allEstep.append(np.array(radiativeTransfer.epsilonStep))

      with warnings.catch_warnings():
        warnings.simplefilter('ignore')
        a = (radiativeTransfer.kappa[:-1,:]/np.sqrt(2*radiativeTransfer.kappaStep.astype(np.complex)))
        b = ((radiativeTransfer.kappa[:-1,:]+radiativeTransfer.kappaStep*scale)/np.sqrt(2*radiativeTransfer.kappaStep.astype(np.complex)))

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
  
  else:
    #print('WARNING: No LOS at position x={}, y={}, z={}.'.format(x,y,z))
    return False,0

  vel = radiativeTransfer.voxelVelocities[iLOS]
  
  if verbose:
    print('voxels:', i)
  
  return True,vel.size

def calculateRadiativeTransfer(backgroundI=0., verbose=False, test=True):

  if np.size(radiativeTransfer.intensity)>0: return
  
  intensity = np.full(radiativeTransfer.kappa.shape[1], backgroundI)
  scale = constants.resolution*constants.pc
  
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
    a = (radiativeTransfer.kappa[:-1,:]/np.sqrt(2*radiativeTransfer.kappaStep.astype(np.complex)))
    b = ((radiativeTransfer.kappa[:-1,:]+radiativeTransfer.kappaStep*scale)/np.sqrt(2*radiativeTransfer.kappaStep.astype(np.complex)))
  
  if verbose: print(radiativeTransfer.kappa.shape[0])
  
  for i in range(len(radiativeTransfer.kappaStep)):
  
    k0 = (radiativeTransfer.kappaStep[i,:]==0)&(abs(radiativeTransfer.kappa[:-1,:][i,:]*scale)<10**-10)
    try:
      kg = radiativeTransfer.kappa[:-1,:][i,:]>10**3*abs(radiativeTransfer.kappaStep[i,:])*scale
    except RuntimeWarning:
      print(radiativeTransfer.x1LoS, radiativeTransfer.x2LoS)
      print(radiativeTransfer.x3LoS)
      print(radiativeTransfer.kappaStep[i,:])
      print(radiativeTransfer.kappa[i,:])
      print(radiativeTransfer.tempInterclumpEmission[i:i+2,1,:])
      input()
      kg = radiativeTransfer.kappa[:-1,:][i,:]>10**3*abs(radiativeTransfer.kappaStep[i,:])*scale
    kE = ~(k0|kg)
    kEg = ~(k0|kg)&(radiativeTransfer.kappaStep[i,:]>0)
    kEl = ~(k0|kg)&(radiativeTransfer.kappaStep[i,:]<0)
  
    if verbose:
      print('i_k0\n', k0)
      print('i_kg\n', kg)
      print('i_kEg\n', kEg)
      print('i_kEl\n', kEl)
    #print(self.__epsilon[0][kg],self.__epsilonStep[0][kg],self.__kappa[0][kg])
  
    if verbose:
      print('\nkappa, kappa step:\n', radiativeTransfer.kappa[i], '\n', radiativeTransfer.kappaStep[i])
      print('\nepsilon, epsilon step, scale\n', radiativeTransfer.epsilon[i], radiativeTransfer.epsilonStep[i], scale)
  
    if test:
  
      if k0.any():
        intensity[k0] = intensity[k0] + radiativeTransfer.epsilon[:-1,:][i,:][k0]*scale+0.5*radiativeTransfer.epsilonStep[i,:][k0]*scale**2
    
      elif kg.any():
        try:
          intensity[kg] = np.exp(-radiativeTransfer.kappa[i,:][kg]*scale) * (intensity[kg] + \
                             ((radiativeTransfer.epsilon[i,:][kg]*radiativeTransfer.kappa[i,:][kg]+radiativeTransfer.epsilonStep[i,:][kg]*(radiativeTransfer.kappa[i,:][kg]*scale-1))/(radiativeTransfer.kappa[i,:][kg]**2.))*np.exp(radiativeTransfer.kappa[i,:][kg]*scale) - \
                             ((radiativeTransfer.epsilon[i,:][kg]*radiativeTransfer.kappa[i,:][kg]-radiativeTransfer.epsilonStep[i,:][kg])/(radiativeTransfer.kappa[i,:][kg]**2.)))
  
        except RuntimeWarning:
          intensity[kg] = np.exp(-radiativeTransfer.kappa[i,:][kg]*scale) * (intensity[kg] + \
                             ((radiativeTransfer.epsilon[i,:][kg]*radiativeTransfer.kappa[i,:][kg]+radiativeTransfer.epsilonStep[i,:][kg]*(radiativeTransfer.kappa[i,:][kg]*scale-1))/(radiativeTransfer.kappa[i,:][kg]**2.)) - \
                             ((radiativeTransfer.epsilon[i,:][kg]*radiativeTransfer.kappa[i,:][kg]-radiativeTransfer.epsilonStep[i,:][kg])/(radiativeTransfer.kappa[i,:][kg]**2.)))

      elif kEg.any():
        if verbose: print('\na, b:\n', a[i,:,:], '\n', b[i,:,:])
  
        aE = np.array(list(Ereal(a[i,:][kEg])))
        bE = np.array(list(Ereal(b[i,:][kEg])))
        intensity[kEg] = (radiativeTransfer.epsilonStep[i,:][kEg]/radiativeTransfer.kappaStep[i,:][kEg]*(1-np.exp(-radiativeTransfer.kappa[:-1,:][i,:][kEg]*scale-radiativeTransfer.kappaStep[i,:][kEg]/2.*scale**2.)) - \
                          (radiativeTransfer.epsilon[:-1,:][i,:][kEg]*radiativeTransfer.kappaStep[i,:][kEg]-radiativeTransfer.epsilonStep[i,:][kEg]*radiativeTransfer.kappa[:-1,:][i,:][kEg])/radiativeTransfer.kappaStep[i,:][kEg] * \
                          np.sqrt(np.pi/(2.*abs(radiativeTransfer.kappaStep[i][kEg].astype(np.complex)))) * \
                          (np.exp(a[i,:]**2.-b[i,:]**2.)*aE-bE) + \
                          intensity[kEg]*np.exp(-radiativeTransfer.kappa[:-1][i,:][kEg]*scale-radiativeTransfer.kappaStep[i,:][kEg]/2.*scale**2.)).real 
  
      elif kEl.any():
        if verbose: print('\na, b:\n', a[i,:,:], '\n', b[i,:])
  
        aE = np.array(list(Eimag(a[i,:][kEl])))
        bE = np.array(list(Eimag(b[i,:][kEl])))
        
        intensity[kEl] = (radiativeTransfer.epsilonStep[i,:][kEl]/radiativeTransfer.kappaStep[i,:][kEl]*(1-np.exp(-radiativeTransfer.kappa[:-1,:][i,:][kEl]*scale-radiativeTransfer.kappaStep[i,:][kEl]/2.*scale**2.))\
                          -(radiativeTransfer.epsilon[:-1][i,:][kEl]*radiativeTransfer.kappaStep[i,:][kEl]-radiativeTransfer.epsilonStep[i,:][kEl]*radiativeTransfer.kappa[:-1,:][i,:][kEl])/radiativeTransfer.kappaStep[i,:][kEl] * \
                          np.sqrt(np.pi/(2.*abs(radiativeTransfer.kappaStep[i][kEl].astype(np.complex))))* \
                          (np.exp(a[i,:]**2.-b[i,:]**2.)*aE-bE) + \
                          intensity[kEl]*np.exp(-radiativeTransfer.kappa[:-1,:][i,:][kEl]*scale-radiativeTransfer.kappaStep[i,:][kEl]/2.*scale**2.)).real
    else:
  
      if (radiativeTransfer.kappaStep[i,:]==0).any() & (abs(radiativeTransfer.kappa[i,:]*scale)<10**-10).any():
        intensity += radiativeTransfer.epsilon[i,:]*scale+0.5*radiativeTransfer.epsilonStep[i,:]*scale**2
  
      elif (radiativeTransfer.kappa[i,:]>10**3*abs(radiativeTransfer.kappaStep[i,:])*scale).any():
        intensity = np.exp(-radiativeTransfer.kappa[i,:]*scale) * (intensity + \
                           ((radiativeTransfer.epsilon[i,:]*radiativeTransfer.kappa[i,:]+radiativeTransfer.epsilonStep[i,:]*(radiativeTransfer.kappa[i,:]*scale-1))/(radiativeTransfer.kappa[i,:]**2.))*np.exp(radiativeTransfer.kappa[i,:]*scale) - \
                           ((radiativeTransfer.epsilon[i,:]*radiativeTransfer.kappa[i,:]-radiativeTransfer.epsilonStep[i,:])/(radiativeTransfer.kappa[i,:]**2.)))
  
      elif (radiativeTransfer.kappaStep[i,:]>0).any():
        if verbose: print('\na, b:\n', a[i,:], '\n', b[i,:])
        
        aE = np.array(list(Ereal(a[i,:])))
        bE = np.array(list(Ereal(b[i,:])))
        
        intensity = (radiativeTransfer.epsilonStep[i,:]/radiativeTransfer.kappaStep[i,:]*(1-np.exp(-radiativeTransfer.kappa[i,:]*scale-radiativeTransfer.kappaStep[:-1][i,:]/2.*scale**2.)) - \
                     (radiativeTransfer.epsilon[i,:]*radiativeTransfer.kappaStep[i,:]-radiativeTransfer.epsilonStep[i,:]*radiativeTransfer.kappa[i,:])/radiativeTransfer.kappaStep[i,:] * \
                     np.sqrt(np.pi/(2.*abs(radiativeTransfer.kappaStep[i,:].astype(np.complex)))) * \
                     (np.exp(a[i,:]**2.-b[i,:]**2.)*aE-bE) + \
                     intensity*np.exp(-radiativeTransfer.kappa[:-1][i,:]*scale-radiativeTransfer.kappaStep[i,:]/2.*scale**2.)).real 
  
      elif (radiativeTransfer.kappaStep[i]<0).any():
        if verbose: print('\na, b:\n', a[i], '\n', b[i])
        
        aE = np.array(list(Eimag(a[i,:])))
        bE = np.array(list(Eimag(b[i,:])))
        
        intensity = (radiativeTransfer.epsilonStep[i,:]/radiativeTransfer.kappaStep[i,:]*(1-np.exp(-radiativeTransfer.kappa[i,:]*scale-radiativeTransfer.kappaStep[i,:]/2.*scale**2.))\
                     -(radiativeTransfer.epsilon[i,:]*radiativeTransfer.kappaStep[i,:]-radiativeTransfer.epsilonStep[i,:]*radiativeTransfer.kappa[i,:])/radiativeTransfer.kappaStep[i,:] * \
                     np.sqrt(np.pi/(2.*abs(radiativeTransfer.kappaStep[i,:].astype(np.complex))))* \
                     (np.exp(a[i,:]**2.-b[i,:]**2.)*aE-bE) + \
                     intensity*np.exp(-radiativeTransfer.kappa[i,:]*scale-radiativeTransfer.kappaStep[i,:]/2.*scale**2.)).real
  
  if verbose:
    print(np.shape(intensity))

  radiativeTransfer.intensity = intensity

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
    Er[il] = (2*x[il]/np.sqrt(np.pi)).astype(np.complex)
  
  if ig.any():
    Er[ig] = (1/(np.sqrt(np.pi) * x[ig])).astype(np.complex)
  
  if ib.any():
    Er[ib] = radiativeTransfer.eTildeReal(x[ib].real).astype(np.complex)

  return Er

def Eimag(x, verbose=False):
  
  if verbose: print('E imaginary input:', x)

  Ei = np.zeros_like(x)

  im = x==abs(x)*1j
  il = (x<0.01)&~im
  ig = (x>8.0)&~im
  ib = ~(ig|il|im)

  # Force x to be a positive real value
  x = abs(x)
  
  if im.any():
    # MASER case: treat with linear approximation
    Ei[im] = 1 + 2*x[im]/np.sqrt(np.pi)

  if il.any():
    Ei[il] = (1 - 2*x[il]/np.sqrt(np.pi)).astype(np.complex)

  if ig.any():
    Ei[ig] = (1/(np.sqrt(np.pi) * x[ig])).astype(np.complex)

  if ib.any():
    Ei[ib] = radiativeTransfer.eTildeImaginary(x[ib]).astype(np.complex)

  return Ei

# jit_module(nopython=False)
