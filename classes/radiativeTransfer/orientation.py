import numpy as np
from numba import jit_module
import copy as c
import cmath
import scipy.interpolate as interpolate
from astropy.io import fits
import importlib as il

import constants
import interpolations
import observations
import radiativeTransfer

def eTildeReal(file='Ereal.dat'):
  eReal = np.genfromtxt('/home/craig/projects/pdr/KOSMA-tau^3/grid/'+file, names=['x', 'Ereal'])
  return (eReal['x'],eReal['Ereal'])

def eTildeImaginary(file='Eimag.dat'):
  eImaginary = np.genfromtxt('/home/craig/projects/pdr/KOSMA-tau^3/grid/'+file, names=['x', 'Eimaginary'])
  return (eImaginary['x'],eImaginary['Eimaginary'])

def calculateObservation(directory='', dim='xy', verbose=False):
  voxelPositions = fits.open(constants.HISTORYPATH+constants.directory+directory+'/voxel_position.fits')[0].data
  clumpIntensity = fits.open(constants.HISTORYPATH+constants.directory+directory+'/intensity_clump.fits')[0].data
  clumpOpticalDepth = fits.open(constants.HISTORYPATH+constants.directory+directory+'/opticalDepth_clump.fits')[0].data
  clumpVelocity = fits.open(constants.HISTORYPATH+constants.directory+directory+'/voxel_clump_velocity.fits')[0].data
  interclumpIntensity = fits.open(constants.HISTORYPATH+constants.directory+directory+'/intensity_interclump.fits')[0].data
  interclumpOpticalDepth = fits.open(constants.HISTORYPATH+constants.directory+directory+'/opticalDepth_interclump.fits')[0].data
  interclumpVelocity = fits.open(constants.HISTORYPATH+constants.directory+directory+'/voxel_interclump_velocity.fits')[0].data

  # Create emission arrays; the current shape is (emission type, voxel, wavelength, velocity)
  clumpEmission = np.array([clumpIntensity, clumpOpticalDepth])
  interclumpEmission = np.array([interclumpIntensity, interclumpOpticalDepth])

  print(voxelPositions[0].shape, clumpEmission.shape, interclumpEmission.shape)

  xPositions,yPositions,zPositions = voxelPositions[:,0],voxelPositions[:,1],voxelPositions[:,2]
  r = np.sqrt(xPositions**2 + yPositions**2)

  #print('\nx\n', np.unique(xArray), '\ny\n', np.unique(yArray), '\nz\n', np.unique(zArray))

  if constants.fromEarth:
    # For an observation from Earth, the data is modified to account for Earth's position at (8750, 0, 0) pc.
    #The integrated intensity is then calculated in the y-z plane to account for different viewing angles across
    #the galaxy. This can be post-processed to convert to galactic coordinates.

    hdul = fits.HDUList()

    for i,velocity in enumerate(constants.velocityRange):

      # Initialise the intensities and map
      position = []
      intensityMap = []

      # Find the voxels that exist at the observing velocity
      boundary = ((xPositions>0)&(r>constants.rGalEarth))|(yPositions>constants.rGalEarth)|(yPositions<-constants.rGalEarth)
      iCV = (clumpVelocity==velocity)
      iIV = (interclumpVelocity==velocity)

      # Separate the voxels in the direction of the Milky Way center (Inner) and the edge of the galaxy (Outer)
      iInner = np.where((~boundary)&(iCV.any(1)&iIV.any(1)))[0]
      iOuter = np.where(boundary&(iCV.any(1)&iIV.any(1)))[0]

      # Define the indeces contributing to the integrated intensity at this velocity
      iInnerClump = np.where(iCV[iInner,:])
      iInnerInterclump = np.where(iIV[iInner,:])
      iOuterClump = np.where(iCV[iOuter,:])
      iOuterInterclump = np.where(iIV[iOuter,:])

      if not (len(iInner)+len(iOuter)): continue

      # print(iInnerClump)
      # print(iInnerInterclump)

      if len(iInner):
        ArrayInner = np.unique([yPositions[~boundary], zPositions[~boundary]], axis=1).T
        radiativeTransfer.tempClumpEmission = clumpEmission[:,iInnerClump[0],:,iInnerClump[1]]
        radiativeTransfer.tempClumpPosition = voxelPositions[iInnerClump[0],:]
        radiativeTransfer.tempClumpVelocity = clumpVelocity[iInnerClump[0],:]
        radiativeTransfer.tempInterclumpEmission = interclumpEmission[:,iInnerInterclump[0],:,iInnerInterclump[1]]
        radiativeTransfer.tempInterclumpPosition = voxelPositions[iInnerInterclump[0],:]
        radiativeTransfer.tempInterclumpVelocity = interclumpVelocity[iInnerInterclump[0],:]

        # print(ArrayInner)

        for y,z in ArrayInner:
          result = setLOS(y=y, z=z, dim='inner disk')
          if result:
            position.append([y,z])
            calculateRadiativeTransfer()
            intensityMap.append(radiativeTransfer.intensity)
    
      # Convert to numpy arrays
      mapPositions = np.array(position)
      intensityMap = np.array(intensityMap)

      # Create HDUs for the map position and intensity and add the velocity in the headers
      PositionHDU = fits.ImageHDU(mapPositions)
      IntensityHDU = fits.ImageHDU(intensityMap)
      PositionHDU.header['Velocity'] = velocity
      PositionHDU.header['Direction'] = 'Inner disk'
      IntensityHDU.header['Velocity'] = velocity
      IntensityHDU.header['Direction'] = 'Inner disk'

      # Add the HDUs to the HDU list
      hdul.append(PositionHDU)
      hdul.append(IntensityHDU)

      positions = []
      intensityMap = []

      if len(iOuter):
        ArrayOuter = np.unique([yPositions[iOuter], zPositions[iOuter]], axis=1).T
        radiativeTransfer.tempClumpEmission = clumpEmission[:,iOuterClump[0],:,iOuterClump[1]]
        radiativeTransfer.tempClumpPosition = voxelPositions[iOuterClump[0],:]
        radiativeTransfer.tempClumpVelocity = clumpVelocity[iOuterClump[0],:]
        radiativeTransfer.tempInterclumpEmission = interclumpEmission[:,iOuterInterclump[0],:,iOuterInterclump[1]]
        radiativeTransfer.tempInterclumpPosition = voxelPositions[iOuterInterclump[0],:]
        radiativeTransfer.tempInterclumpVelocity = interclumpVelocity[iOuterInterclump[0],:]

        #print(ArrayInner)

        for y,z in ArrayOuter:
          result = setLOS(y=y, z=z, dim='outer disk', reverse=True)
          if result:
            position.append([y,z])
            calculateRadiativeTransfer()
            intensityMap.append(radiativeTransfer.intensity)
    
      # Convert to numpy arrays
      mapPositions = np.array(position)
      intensityMap = np.array(intensityMap)

      # Create HDUs for the map position and intensity and add the velocity in the headers
      PositionHDU = fits.ImageHDU(mapPositions)
      IntensityHDU = fits.ImageHDU(intensityMap)
      PositionHDU.header['Velocity'] = velocity
      PositionHDU.header['Direction'] = 'Outer disk'
      IntensityHDU.header['Velocity'] = velocity
      IntensityHDU.header['Direction'] = 'Outer disk'

      # Add the HDUs to the HDU list
      hdul.append(PositionHDU)
      hdul.append(IntensityHDU)

      print('Writing {} km/s HDU'.format(velocity))

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

def setLOS(emission=0, positions=0, x=0, y=0, z=0, dim='xy', reverse=False, verbose=False):
  '''
  The emission dimensions should be velocity x species x 2 x voxels. Axis 1 should be split for intensity and optical depth.
  The positions dimensions should be 3 x voxels.
  '''
  LOSvoxels = []
  zPosition = []
  epsilon = []
  kappa = []

  scale = constants.resolution*constants.pc

  # This block matches the voxel positions to add the ensembles belonging to the same voxel
  nrows, ncols = radiativeTransfer.tempClumpPosition.shape
  dtype={'names':['f{}'.format(i) for i in range(ncols)], \
         'formats':ncols * [radiativeTransfer.tempClumpPosition.dtype]}
  common,iCommonClump,iCommonInterclump = np.intersect1d(radiativeTransfer.tempClumpPosition.view(dtype), radiativeTransfer.tempInterclumpPosition.view(dtype), return_indices=True)
  
  # Intensity and optical depth have shape (voxel, wavelength)
  gridIntensity = radiativeTransfer.tempInterclumpEmission[:,0,:]
  gridIntensity[iCommonInterclump,:] = gridIntensity[iCommonInterclump,:] + radiativeTransfer.tempClumpEmission[iCommonClump,0,:]# emission[0,:,:]#

  gridOpticalDepth = radiativeTransfer.tempInterclumpEmission[:,1,:]# emission[1,:,:]#
  gridOpticalDepth[iCommonInterclump,:] += radiativeTransfer.tempClumpEmission[iCommonClump,1,:]

  xGrid,yGrid,zGrid = radiativeTransfer.tempInterclumpPosition[:,0],radiativeTransfer.tempInterclumpPosition[:,1],radiativeTransfer.tempInterclumpPosition[:,2]#grid.getVoxelPositions()
  
  if verbose:
    print('Centered at x={}, y={}, z={}'.format(x,y,z))
  
  if 'disk' in dim:
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
    radiativeTransfer.intensity = gridIntensity[iLOS, :]
    return True
  
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
    if not 'x' in dim:
      radiativeTransfer.x3LoS = xGrid[iLOS]
    elif not 'y' in dim:
      radiativeTransfer.x3LoS = yGrid[iLOS]
    elif not 'z' in dim:
      radiativeTransfer.x3LoS = zGrid[iLOS]
    if reverse:
      i = np.argsort(radiativeTransfer.x3LoS)
    else:
      i = np.argsort(radiativeTransfer.x3LoS)[::-1]
    radiativeTransfer.epsilon = np.array(epsilon, dtype=np.float)[i]
    radiativeTransfer.epsilonStep = (radiativeTransfer.epsilon[1:]-radiativeTransfer.epsilon[:-1])/(scale)
    radiativeTransfer.kappa = np.array(kappa, dtype=np.float)[i]
    radiativeTransfer.kappaStep = (radiativeTransfer.kappa[1:]-radiativeTransfer.kappa[:-1])/(scale)
  
  else:
    #print('WARNING: No LOS at position x={}, y={}, z={}.'.format(x,y,z))
    return False
  
  if verbose:
    print('voxels:', i)
  
  return True

def calculateRadiativeTransfer(backgroundI=0., verbose=False, test=False):
  
  if np.size(radiativeTransfer.intensity)==1: return
  
  intensity = np.full(radiativeTransfer.kappa[0].shape, backgroundI)
  scale = constants.resolution*constants.pc
  
  # Boolean indeces to separate how the intensity is calculated
  k0 = (radiativeTransfer.kappaStep==0)&(abs(radiativeTransfer.kappa[:-1]*constants.resolution)<10**-10)
  kg = radiativeTransfer.kappa[:-1]>10**3*abs(radiativeTransfer.kappaStep)*constants.resolution
  kE = ~(k0|kg)
  kEg = ~(k0|kg)&(radiativeTransfer.kappaStep>0)
  kEl = ~(k0|kg)&(radiativeTransfer.kappaStep<0)
  
  # Variables for utilising the E tilde tables
  a = (radiativeTransfer.kappa[:-1,:]/np.sqrt(2*radiativeTransfer.kappaStep.astype(np.complex)))
  b = ((radiativeTransfer.kappa[:-1,:]+radiativeTransfer.kappaStep*scale)/np.sqrt(2*radiativeTransfer.kappaStep.astype(np.complex)))
  
  if verbose: print(len(radiativeTransfer.losVoxels[:-1]))
  
  for i in range(len(radiativeTransfer.kappaStep)):
  
    k0 = (radiativeTransfer.kappaStep[i,:]==0)&(abs(radiativeTransfer.kappa[:-1,:][i,:]*constants.resolution)<10**-10)
    kg = radiativeTransfer.kappa[:-1,:][i,:]>10**3*abs(radiativeTransfer.kappaStep[i,:])*constants.resolution
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
        intensity[k0] += radiativeTransfer.epsilon[:-1,:][i,:][k0]*scale+0.5*radiativeTransfer.epsilonStep[i,:][k0]*scale**2
  
      if kg.any():
        intensity[kg] = np.exp(-radiativeTransfer.kappa[i,:][kg]*scale) * (intensity[kg] + \
                           ((radiativeTransfer.epsilon[i,:][kg]*radiativeTransfer.kappa[i,:][kg]+radiativeTransfer.epsilonStep[i,:][kg]*(radiativeTransfer.kappa[i,:][kg]*scale-1))/(radiativeTransfer.kappa[i,:][kg]**2.))*np.exp(radiativeTransfer.kappa[i,:][kg]*scale) - \
                           ((radiativeTransfer.epsilon[i,:][kg]*radiativeTransfer.kappa[i,:][kg]-radiativeTransfer.epsilonStep[i,:][kg])/(radiativeTransfer.kappa[i,:][kg]**2.)))
  
      if kEg.any():
        if verbose: print('\na, b:\n', a[i,:,:], '\n', b[i,:,:])
  
        aE = np.array(list(Ereal(a[i,:][kEg])))
        bE = np.array(list(Ereal(b[i,:][kEg])))
        intensity[kEg] = (radiativeTransfer.epsilonStep[i,:][kEg]/radiativeTransfer.kappaStep[i,:][kEg]*(1-np.exp(-radiativeTransfer.kappa[:-1,:][i,:][kEg]*scale-radiativeTransfer.kappaStep[i,:][kEg]/2.*scale**2.)) - \
                          (radiativeTransfer.epsilon[:-1,:][i,:][kEg]*radiativeTransfer.kappaStep[i,:][kEg]-radiativeTransfer.epsilonStep[i,:][kEg]*radiativeTransfer.kappa[:-1,:][i,:][kEg])/radiativeTransfer.kappaStep[i,:][kEg] * \
                          np.sqrt(np.pi/(2.*abs(radiativeTransfer.kappaStep[i][kEg].astype(np.complex)))) * \
                          (np.exp(a[i,:]**2.-b[i,:]**2.)*aE-bE) + \
                          intensity[kEg]*np.exp(-radiativeTransfer.kappa[:-1][i,:][kEg]*scale-radiativeTransfer.kappaStep[i,:][kEg]/2.*scale**2.)).real 
  
      if kEl.any():
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
  
  radiativeTransfer.intensity = intensity
  return

def Ereal(x, verbose=False):
  
  if verbose: print('E real input:', x)
  
  if (x.imag==0).any(): x = x.real
  # x should be a real number. remove imaginary party '0j' which
  # prevents ordering
  
  if (x<0.01).any():
    return 2*x/np.sqrt(np.pi)
  
  elif (x>8.0).any():
    return 1/(np.sqrt(np.pi) * x)
  
  else:
    return np.array(list(map(radiativeTransfer.eTildeReal, x)))

def Eimag(x, verbose=False):
  
  if verbose: print('E imaginary input:', x)
  
  if (x==abs(x)*1j).any():
    # maser case. treated in linear approximation
    x = abs(x)
    return 1 + 2*x/np.sqrt(np.pi)
  
  else:
    x = abs(x)
    # x needs to be real and positive, i.e. abs(a) or abs(b)
  
    if (x<0.01).any():
      return 1 - 2*x/np.sqrt(np.pi)
  
    elif (x>8.0).any():
      return 1/(np.sqrt(np.pi) * x)
  
    else:
      return np.array(list(map(radiativeTransfer.eTildeImaginary, x)))

# jit_module(nopython=False)
