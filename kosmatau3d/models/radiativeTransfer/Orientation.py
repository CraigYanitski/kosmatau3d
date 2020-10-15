import numpy as np
import copy as c
import cmath
import scipy.interpolate as interpolate
import importlib as il

import constants
import observations
import orientation

#from Dimensions import *
#from Observations import *
# class Orientation(object):
#   '''
#   This is a class to alter how the PDR is viewed. This will make it easy to
#   adjust the viewing angle and check the observed linewidths. For now, the
#   LoS is parallel to the z-axis, and the position in the x-y plane can be
#   specified as args. In the future, I plan to create a map of the radiative
#   transfer according to the projected size of the PDR. For the Milky Way example,
#   this is obviously not taken from the perspective of Earth.
#   '''

#   # PRIVATE
#   def __init__(self, dimensions, backgroundI=0., interpolation='linear'):
#     #self.__dimensions = dimensions
#     self.__xLoS = self.__yLoS = 0
#     zRange = np.unique(dimensions.voxelCartesianPosition()[2])
#     #self.__zLoS = np.arange(min(zRange),max(zRange)+1)
#     self.__losVoxels = []
#     self.__losZ = []
#     self.__losKappa = []
#     self.__losKappaStep = []
#     self.__losEpsilon = []
#     self.__losEpsilonStep = []
#     self.__backgroundI = backgroundI
#     self.__intensity = np.array([0.])
#     #observations = Observations()
#     eTildeReal = observations.eTildeReal
#     eTildeImaginary = observations.eTildeImaginary
#     self.__eTildeReal = interpolate.interp1d(eTildeReal[0], eTildeReal[1], kind=interpolation)
#     self.__eTildeImaginary = interpolate.interp1d(eTildeImaginary[0], eTildeImaginary[1], kind=interpolation)
#     return

  #PUBLIC
  # def reloadModules(self):
  #   constants.resolution.reloadModules()
  #   for voxel in self.__losVoxels: voxel.reloadModules()
  #   return
def setLOS(emission=0, positions=0, x=0, y=0, z=0, dim='xy', verbose=False):
  '''
  The emission dimensions should be velocity x species x 2 x voxels. Axis 1 should be split for intensity and optical depth.
  The positions dimensions should be 3 x voxels.
  '''
  LOSvoxels = []
  zPosition = []
  epsilon = []
  kappa = []
  factor = 1#np.exp(-z**2/500.**2)      #this is to make the disk more centrally-located
  scale = constants.resolution*constants.pc
  gridIntensity = emission[:,:,0,:]#voxels = grid.allVoxels()
  gridOpticalDepth = emission[:,:,1,:]#voxels = grid.allVoxels()
  xGrid,yGrid,zGrid = positions#grid.getVoxelPositions()
  if verbose:
    print('Centered at x={}, y={}, z={}'.format(x,y,z))
  if ('x' in dim) and ('y' in dim):
    orientation.x1LoS = x
    orientation.x2LoS = y
    iLOS = np.where((xGrid==x)&(yGrid==y))[0]
  elif ('x' in dim) and ('z' in dim):
    orientation.x1LoS = x
    orientation.x2LoS = z
    iLOS = np.where((xGrid==x)&(zGrid==z))[0]
  elif ('z' in dim) and ('y' in dim):
    orientation.x1LoS = z
    orientation.x2LoS = y
    iLOS = np.where((zGrid==z)&(yGrid==y))[0]
  else:
    print('\nPlease enter valid dimensions.\n')
    return
  if verbose:
    print(iLOS)
  #self.__losVoxels = grid.allVoxels()[iLOS]
  if iLOS.size==1:
    #LOSvoxels.append(voxels[iLOS[0]])
    orientation.intensity = gridIntensity[:,:,iLOS]
    #tau = gridOpticalDepth[:,:,iLOS]
    # epsilon.append(factor*intensity/(scale))
    # kappa.append(factor*tau/(scale))
    # epsilonStep = [0]
    # kappaStep = [0]
    #orientation.intensity = intensity
    return
  elif iLOS.size>1:
    #for i in iLOS:
      #LOSvoxels.append(voxels[i])
    intensity = gridIntensity[:,:,iLOS]
    tau = gridOpticalDepth[:,:,iLOS]
      #if z:
      #print(factor, scale, intensity)
    epsilon.append(factor*intensity/scale)
    kappa.append(factor*tau/scale)
      #else:
      #  epsilon.append(intensity/(constants.resolution))
      #  kappa.append(tau/(constants.resolution))
    if verbose:
      print('Intensity:', intensity, scale)
      print('Optical depth:', tau, scale)
    #epsilonStep = (epsilon[1:]-epsilon[:-1])/(scale)
    #kappaStep = (kappa[1:]-kappa[:-1])/(scale)
    #orientation.losVoxels = LOSvoxels
    orientation.zPosition = zGrid[iLOS]
    i = np.argsort(orientation.zPosition)[::-1]
    orientation.epsilon = np.array(epsilon, dtype=np.float)[i]
    orientation.epsilonStep = (orientation.epsilon[1:]-orientation.epsilon[:-1])/(scale)
    orientation.kappa = np.array(kappa, dtype=np.float)[i]
    orientation.kappaStep = (orientation.kappa[1:]-orientation.kappa[:-1])/(scale)
  # for voxel in grid.allVoxels():
  #   x,y,z = voxel.getPosition()
  #   if ('x' in dim) and ('y' in dim):
  #     x1 = x
  #     x2 = y
  #     z0 = c.copy(z)
  #   elif ('x' in dim) and ('z' in dim):
  #     x1 = x
  #     x2 = z
  #     z0 = c.copy(y)
  #   elif ('z' in dim) and ('y' in dim):
  #     x1 = z
  #     x2 = y
  #     z0 = c.copy(x)
  #   if (x1<(self.__x1LoS+constants.resolution/2.)) and (x1>(self.__x1LoS-constants.resolution/2.)) and (x2<(self.__x2LoS+constants.resolution/2.)) and (x2>(self.__x2LoS-constants.resolution/2.)):
  #     intensity,tau,fuv = voxel.getEmission()
  #     voxels.append(voxel)
  #     zPosition.append(z0)
  #     #if z:
  #     factor = 1#np.exp(-z**2/500.**2)      #this is to make the disk more centrally-located
  #     epsilon.append(factor*intensity/(scale))
  #     kappa.append(factor*tau/(scale))
  #     #else:
  #     #  epsilon.append(intensity/(constants.resolution))
  #     #  kappa.append(tau/(constants.resolution))
  #     if verbose: print(tau, scale)
  else: print('WARNING: No LOS at position x={}, y={}, z={}.'.format(x,y,z))
  if verbose:
    print('voxels:', i)
  return
def calculateRadiativeTransfer(velocity, verbose=False, test=True):
  if len(self.__losVoxels)==1: return orientation.intensity[0]
  intensity = np.full(orientation.kappa[0].shape, orientation.backgroundI)
  scale = constants.resolution*constants.pc
  # Boolean indeces to separate how the intensity is calculated
  k0 = (orientation.kappaStep==0)&(abs(orientation.kappa[:-1]*constants.resolution)<10**-10)
  kg = orientation.kappa[:-1]>10**3*abs(orientation.kappaStep)*constants.resolution
  kE = ~(k0|kg)
  kEg = ~(k0|kg)&(orientation.kappaStep>0)
  kEl = ~(k0|kg)&(orientation.kappaStep<0)
  # Variables for utilising the E tilde tables
  a = (orientation.kappa[:-1,:,:]/np.sqrt(2*orientation.kappaStep.astype(np.complex)))
  b = ((orientation.kappa[:-1,:,:]+orientation.kappaStep*scale)/np.sqrt(2*orientation.kappaStep.astype(np.complex)))
  if verbose: print(len(orientation.losVoxels[:-1]))
  for i in range(len(orientation.kappaStep)):
    k0 = (orientation.kappaStep[i,:,:]==0)&(abs(orientation.kappa[:-1,:,:][i,:,:]*constants.resolution)<10**-10)
    kg = orientation.kappa[:-1,:,:][i,:,:]>10**3*abs(orientation.kappaStep[i,:,:])*constants.resolution
    kE = ~(k0|kg)
    kEg = ~(k0|kg)&(orientation.kappaStep[i,:,:]>0)
    kEl = ~(k0|kg)&(orientation.kappaStep[i,:,:]<0)
    if verbose:
      print('i_k0\n', k0)
      print('i_kg\n', kg)
      print('i_kEg\n', kEg)
      print('i_kEl\n', kEl)
    #print(self.__epsilon[0][kg],self.__epsilonStep[0][kg],self.__kappa[0][kg])
    if verbose:
      print('\nkappa, kappa step:\n', orientation.kappa[i], '\n', orientation.kappaStep[i])
      print('\nepsilon, epsilon step, scale\n', orientation.epsilon[i], orientation.epsilonStep[i], scale)
    if test:
      if k0.any():
        intensity[k0] += orientation.epsilon[:-1,:,:][i,:,:][k0]*scale+0.5*orientation.epsilonStep[i,:,:][k0]*scale**2
      if kg.any():
        intensity[kg] = np.exp(-orientation.kappa[i,:,:][kg]*scale) * (intensity[kg] + \
                           ((orientation.epsilon[i,:,:][kg]*orientation.kappa[i,:,:][kg]+orientation.epsilonStep[i,:,:][kg]*(orientation.kappa[i,:,:][kg]*scale-1))/(orientation.kappa[i,:,:][kg]**2.))*np.exp(orientation.kappa[i,:,:][kg]*scale) - \
                           ((orientation.epsilon[i,:,:][kg]*orientation.kappa[i,:,:][kg]-orientation.epsilonStep[i,:,:][kg])/(orientation.kappa[i,:,:][kg]**2.)))
      if kEg.any():
        if verbose: print('\na, b:\n', a[i,:,:], '\n', b[i,:,:])
        aE = np.array(list(Ereal(a[i,:,:][kEg])))
        bE = np.array(list(Ereal(b[i,:,:][kEg])))
        intensity[kEg] = (orientation.epsilonStep[i,:,:][kEg]/orientation.kappaStep[i,:,:][kEg]*(1-np.exp(-orientation.kappa[:-1,:,:][i,:,:][kEg]*scale-orientation.kappaStep[i,:,:][kEg]/2.*scale**2.)) - \
                          (orientation.epsilon[:-1,:,:][i,:,:][kEg]*orientation.kappaStep[i,:,:][kEg]-orientation.epsilonStep[i,:,:][kEg]*orientation.kappa[:-1,:,:][i,:,:][kEg])/orientation.kappaStep[i,:,:][kEg] * \
                          np.sqrt(np.pi/(2.*abs(orientation.kappaStep[i][kEg].astype(np.complex)))) * \
                          (np.exp(a[i,:,:]**2.-b[i,:,:]**2.)*aE-bE) + \
                          intensity[kEg]*np.exp(-orientation.kappa[:-1][i,:,:][kEg]*scale-orientation.kappaStep[i,:,:][kEg]/2.*scale**2.)).real 
      if kEl.any():
        if verbose: print('\na, b:\n', a[i,:,:], '\n', b[i,:,:])
        aE = np.array(list(Eimag(a[i,:,:][kEl])))
        bE = np.array(list(Eimag(b[i,:,:][kEl])))
        intensity[kEl] = (orientation.epsilonStep[i,:,:][kEl]/orientation.kappaStep[i,:,:][kEl]*(1-np.exp(-orientation.kappa[:-1,:,:][i,:,:][kEl]*scale-orientation.kappaStep[i,:,:][kEl]/2.*scale**2.))\
                          -(orientation.epsilon[:-1][i,:,:][kEl]*orientation.kappaStep[i,:,:][kEl]-orientation.epsilonStep[i,:,:][kEl]*orientation.kappa[:-1,:,:][i,:,:][kEl])/orientation.kappaStep[i,:,:][kEl] * \
                          np.sqrt(np.pi/(2.*abs(orientation.kappaStep[i][kEl].astype(np.complex))))* \
                          (np.exp(a[i,:,:]**2.-b[i,:,:]**2.)*aE-bE) + \
                          intensity[kEl]*np.exp(-orientation.kappa[:-1,:,:][i,:,:][kEl]*scale-orientation.kappaStep[i,:,:][kEl]/2.*scale**2.)).real
    else:
      if (orientation.kappaStep[i,:,:]==0).any() & (abs(orientation.kappa[i]*scale)<10**-10).any():
        intensity += orientation.epsilon[i]*scale+0.5*orientation.epsilonStep[i]*scale**2
      elif (orientation.kappa[i]>10**3*abs(orientation.kappaStep[i])*scale).any():
        intensity = np.exp(-orientation.kappa[i]*scale) * (intensity + \
                           ((orientation.epsilon[i]*orientation.kappa[i]+orientation.epsilonStep[i]*(orientation.kappa[i]*scale-1))/(orientation.kappa[i]**2.))*np.exp(orientation.kappa[i]*scale) - \
                           ((orientation.epsilon[i]*orientation.kappa[i]-orientation.epsilonStep[i])/(orientation.kappa[i]**2.)))
      elif (orientation.kappaStep[i]>0).any():
        if verbose: print('\na, b:\n', a[i], '\n', b[i])
        aE = np.array(list(Ereal(a[i])))
        bE = np.array(list(Ereal(b[i])))
        intensity = (orientation.epsilonStep[i]/orientation.kappaStep[i]*(1-np.exp(-orientation.kappa[i]*scale-orientation.kappaStep[:-1][i]/2.*scale**2.)) - \
                     (orientation.epsilon[i]*orientation.kappaStep[i]-orientation.epsilonStep[i]*orientation.kappa[i])/orientation.kappaStep[i] * \
                     np.sqrt(np.pi/(2.*abs(orientation.kappaStep[i].astype(np.complex)))) * \
                     (np.exp(a[i]**2.-b[i]**2.)*aE-bE) + \
                     intensity*np.exp(-orientation.kappa[:-1][i]*scale-orientation.kappaStep[i]/2.*scale**2.)).real 
      elif (orientation.kappaStep[i]<0).any():
        if verbose: print('\na, b:\n', a[i], '\n', b[i])
        aE = np.array(list(Eimag(a[i])))
        bE = np.array(list(Eimag(b[i])))
        intensity = (orientation.epsilonStep[i]/orientation.kappaStep[i]*(1-np.exp(-orientation.kappa[i]*scale-orientation.kappaStep[i]/2.*scale**2.))\
                     -(orientation.epsilon[i]*orientation.kappaStep[i]-orientation.epsilonStep[i]*orientation.kappa[i])/orientation.kappaStep[i] * \
                     np.sqrt(np.pi/(2.*abs(orientation.kappaStep[i].astype(np.complex))))* \
                     (np.exp(a[i]**2.-b[i]**2.)*aE-bE) + \
                     intensity*np.exp(-orientation.kappa[i]*scale-orientation.kappaStep[i]/2.*scale**2.)).real
  return intensity
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
    return np.array(list(map(self.__eTildeReal, x)))
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
      return np.array(list(map(self.__dTildeImaginary, x)))
