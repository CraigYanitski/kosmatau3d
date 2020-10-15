import numpy as np
import copy as c
import cmath
import scipy.interpolate as interpolate
import importlib as il

import constants
import observations

#from Dimensions import *
#from Observations import *
class Orientation(object):
  '''
  This is a class to alter how the PDR is viewed. This will make it easy to
  adjust the viewing angle and check the observed linewidths. For now, the
  LoS is parallel to the z-axis, and the position in the x-y plane can be
  specified as args. In the future, I plan to create a map of the radiative
  transfer according to the projected size of the PDR. For the Milky Way example,
  this is obviously not taken from the perspective of Earth.
  '''

  # PRIVATE
  def __init__(self, dimensions, backgroundI=0., interpolation='linear'):
    #self.__dimensions = dimensions
    self.__xLoS = self.__yLoS = 0
    zRange = np.unique(dimensions.voxelCartesianPosition()[2])
    #self.__zLoS = np.arange(min(zRange),max(zRange)+1)
    self.__losVoxels = []
    self.__losZ = []
    self.__losKappa = []
    self.__losKappaStep = []
    self.__losEpsilon = []
    self.__losEpsilonStep = []
    self.__backgroundI = backgroundI
    self.__intensity = np.array([0.])
    #observations = Observations()
    eTildeReal = observations.eTildeReal
    eTildeImaginary = observations.eTildeImaginary
    self.__eTildeReal = interpolate.interp1d(eTildeReal[0], eTildeReal[1], kind=interpolation)
    self.__eTildeImaginary = interpolate.interp1d(eTildeImaginary[0], eTildeImaginary[1], kind=interpolation)
    return

  #PUBLIC
  def reloadModules(self):
    constants.resolution.reloadModules()
    for voxel in self.__losVoxels: voxel.reloadModules()
    return
  def setLOS(self, grid, x=0, y=0, z=0, dim='xy', verbose=False):
    LOSvoxels = []
    zPosition = []
    epsilon = []
    kappa = []
    factor = 1#np.exp(-z**2/500.**2)      #this is to make the disk more centrally-located
    scale = constants.resolution*constants.pc
    voxels = grid.allVoxels()
    xGrid,yGrid,zGrid = grid.getVoxelPositions()
    if verbose:
      print('Centered at x={}, y={}, z={}'.format(x,y,z))
    if ('x' in dim) and ('y' in dim):
      self.__x1LoS = x
      self.__x2LoS = y
      iLOS = np.where((xGrid==x)&(yGrid==y))[0]
    elif ('x' in dim) and ('z' in dim):
      self.__x1LoS = x
      self.__x2LoS = z
      iLOS = np.where((xGrid==x)&(zGrid==z))[0]
    elif ('z' in dim) and ('y' in dim):
      self.__x1LoS = z
      self.__x2LoS = y
      iLOS = np.where((zGrid==z)&(yGrid==y))[0]
    else:
      print('\nPlease enter valid dimensions.\n')
      return
    if verbose:
      print(iLOS)
    #self.__losVoxels = grid.allVoxels()[iLOS]
    if iLOS.size==1:
      LOSvoxels.append(voxels[iLOS[0]])
      intensity,tau,fuv = LOSvoxels[0].getEmission()
      # epsilon.append(factor*intensity/(scale))
      # kappa.append(factor*tau/(scale))
      # epsilonStep = [0]
      # kappaStep = [0]
      self.__intensity = intensity
    elif iLOS.size>1:
      for i in iLOS:
        LOSvoxels.append(voxels[i])
        intensity,tau,fuv = LOSvoxels[-1].getEmission()
        #if z:
        epsilon.append(factor*intensity/(scale))
        kappa.append(factor*tau/(scale))
        #else:
        #  epsilon.append(intensity/(constants.resolution))
        #  kappa.append(tau/(constants.resolution))
        if verbose:
          print('Intensity:', intensity, scale)
          print('Optical depth:', tau, scale)
      #epsilonStep = (epsilon[1:]-epsilon[:-1])/(scale)
      #kappaStep = (kappa[1:]-kappa[:-1])/(scale)
      self.__losVoxels = LOSvoxels
      zPosition = zGrid[iLOS]
      i = np.argsort(zPosition)[::-1]
      self.__epsilon = np.array(epsilon, dtype=np.float)[i]
      self.__epsilonStep = (self.__epsilon[1:]-self.__epsilon[:-1])/(scale)
      self.__kappa = np.array(kappa, dtype=np.float)[i]
      self.__kappaStep = (self.__kappa[1:]-self.__kappa[:-1])/(scale)
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
  def calculateRadiativeTransfer(self, velocity, verbose=False, test=True):
    if len(self.__losVoxels)>1: return self.__intensity
    intensity = np.full(self.__kappa[0].shape, self.__backgroundI)
    scale = constants.resolution*constants.pc
    # Boolean indeces to separate how the intensity is calculated
    k0 = (self.__kappaStep==0)&(abs(self.__kappa[:-1]*constants.resolution)<10**-10)
    kg = self.__kappa[:-1]>10**3*abs(self.__kappaStep)*constants.resolution
    kE = ~(k0|kg)
    kEg = ~(k0|kg)&(self.__kappaStep>0)
    kEl = ~(k0|kg)&(self.__kappaStep<0)
    # Variables for utilising the E tilde tables
    a = (self.__kappa[:-1,:,:]/np.sqrt(2*self.__kappaStep.astype(np.complex)))
    b = ((self.__kappa[:-1,:,:]+self.__kappaStep*scale)/np.sqrt(2*self.__kappaStep.astype(np.complex)))
    if verbose: print(len(self.__losVoxels[:-1]))
    for i in range(len(self.__losVoxels[:-1])):
      k0 = (self.__kappaStep[i,:,:]==0)&(abs(self.__kappa[:-1,:,:][i,:,:]*constants.resolution)<10**-10)
      kg = self.__kappa[:-1,:,:][i,:,:]>10**3*abs(self.__kappaStep[i,:,:])*constants.resolution
      kE = ~(k0|kg)
      kEg = ~(k0|kg)&(self.__kappaStep[i,:,:]>0)
      kEl = ~(k0|kg)&(self.__kappaStep[i,:,:]<0)
      if verbose:
        print('i_k0\n', k0)
        print('i_kg\n', kg)
        print('i_kEg\n', kEg)
        print('i_kEl\n', kEl)
      #print(self.__epsilon[0][kg],self.__epsilonStep[0][kg],self.__kappa[0][kg])
      if verbose:
        print('\nkappa, kappa step:\n', self.__kappa[i], '\n', self.__kappaStep[i])
        print('\nepsilon, epsilon step, scale\n', self.__epsilon[i], self.__epsilonStep[i], scale)
      if test:
        if k0.any():
          intensity[k0] += self.__epsilon[:-1,:,:][i,:,:][k0]*scale+0.5*self.__epsilonStep[i,:,:][k0]*scale**2
        if kg.any():
          intensity[kg] = np.exp(-self.__kappa[i,:,:][kg]*scale) * (intensity[kg] + \
                             ((self.__epsilon[i,:,:][kg]*self.__kappa[i,:,:][kg]+self.__epsilonStep[i,:,:][kg]*(self.__kappa[i,:,:][kg]*scale-1))/(self.__kappa[i,:,:][kg]**2.))*np.exp(self.__kappa[i,:,:][kg]*scale) - \
                             ((self.__epsilon[i,:,:][kg]*self.__kappa[i,:,:][kg]-self.__epsilonStep[i,:,:][kg])/(self.__kappa[i,:,:][kg]**2.)))
        if kEg.any():
          if verbose: print('\na, b:\n', a[i,:,:], '\n', b[i,:,:])
          aE = np.array(list(self.Ereal(a[i,:,:][kEg])))
          bE = np.array(list(self.Ereal(b[i,:,:][kEg])))
          intensity[kEg] = (self.__epsilonStep[i,:,:][kEg]/self.__kappaStep[i,:,:][kEg]*(1-np.exp(-self.__kappa[:-1,:,:][i,:,:][kEg]*scale-self.__kappaStep[i,:,:][kEg]/2.*scale**2.)) - \
                            (self.__epsilon[:-1,:,:][i,:,:][kEg]*self.__kappaStep[i,:,:][kEg]-self.__epsilonStep[i,:,:][kEg]*self.__kappa[:-1,:,:][i,:,:][kEg])/self.__kappaStep[i,:,:][kEg] * \
                            np.sqrt(np.pi/(2.*abs(self.__kappaStep[i][kEg].astype(np.complex)))) * \
                            (np.exp(a[i,:,:]**2.-b[i,:,:]**2.)*aE-bE) + \
                            intensity[kEg]*np.exp(-self.__kappa[:-1][i,:,:][kEg]*scale-self.__kappaStep[i,:,:][kEg]/2.*scale**2.)).real 
        if kEl.any():
          if verbose: print('\na, b:\n', a[i,:,:], '\n', b[i,:,:])
          aE = np.array(list(self.Eimag(a[i,:,:][kEl])))
          bE = np.array(list(self.Eimag(b[i,:,:][kEl])))
          intensity[kEl] = (self.__epsilonStep[i,:,:][kEl]/self.__kappaStep[i,:,:][kEl]*(1-np.exp(-self.__kappa[:-1,:,:][i,:,:][kEl]*scale-self.__kappaStep[i,:,:][kEl]/2.*scale**2.))\
                            -(self.__epsilon[:-1][i,:,:][kEl]*self.__kappaStep[i,:,:][kEl]-self.__epsilonStep[i,:,:][kEl]*self.__kappa[:-1,:,:][i,:,:][kEl])/self.__kappaStep[i,:,:][kEl] * \
                            np.sqrt(np.pi/(2.*abs(self.__kappaStep[i][kEl].astype(np.complex))))* \
                            (np.exp(a[i,:,:]**2.-b[i,:,:]**2.)*aE-bE) + \
                            intensity[kEl]*np.exp(-self.__kappa[:-1,:,:][i,:,:][kEl]*scale-self.__kappaStep[i,:,:][kEl]/2.*scale**2.)).real
      else:
        if (self.__kappaStep[i,:,:]==0).any() & (abs(self.__kappa[i]*scale)<10**-10).any():
          intensity += self.__epsilon[i]*scale+0.5*self.__epsilonStep[i]*scale**2
        elif (self.__kappa[i]>10**3*abs(self.__kappaStep[i])*scale).any():
          intensity = np.exp(-self.__kappa[i]*scale) * (intensity + \
                             ((self.__epsilon[i]*self.__kappa[i]+self.__epsilonStep[i]*(self.__kappa[i]*scale-1))/(self.__kappa[i]**2.))*np.exp(self.__kappa[i]*scale) - \
                             ((self.__epsilon[i]*self.__kappa[i]-self.__epsilonStep[i])/(self.__kappa[i]**2.)))
        elif (self.__kappaStep[i]>0).any():
          if verbose: print('\na, b:\n', a[i], '\n', b[i])
          aE = np.array(list(self.Ereal(a[i])))
          bE = np.array(list(self.Ereal(b[i])))
          intensity = (self.__epsilonStep[i]/self.__kappaStep[i]*(1-np.exp(-self.__kappa[i]*scale-self.__kappaStep[:-1][i]/2.*scale**2.)) - \
                       (self.__epsilon[i]*self.__kappaStep[i]-self.__epsilonStep[i]*self.__kappa[i])/self.__kappaStep[i] * \
                       np.sqrt(np.pi/(2.*abs(self.__kappaStep[i].astype(np.complex)))) * \
                       (np.exp(a[i]**2.-b[i]**2.)*aE-bE) + \
                       intensity*np.exp(-self.__kappa[:-1][i]*scale-self.__kappaStep[i]/2.*scale**2.)).real 
        elif (self.__kappaStep[i]<0).any():
          if verbose: print('\na, b:\n', a[i], '\n', b[i])
          aE = np.array(list(self.Eimag(a[i])))
          bE = np.array(list(self.Eimag(b[i])))
          intensity = (self.__epsilonStep[i]/self.__kappaStep[i]*(1-np.exp(-self.__kappa[i]*scale-self.__kappaStep[i]/2.*scale**2.))\
                       -(self.__epsilon[i]*self.__kappaStep[i]-self.__epsilonStep[i]*self.__kappa[i])/self.__kappaStep[i] * \
                       np.sqrt(np.pi/(2.*abs(self.__kappaStep[i].astype(np.complex))))* \
                       (np.exp(a[i]**2.-b[i]**2.)*aE-bE) + \
                       intensity*np.exp(-self.__kappa[i]*scale-self.__kappaStep[i]/2.*scale**2.)).real
    return intensity
  def Ereal(self, x, verbose=False):
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
  def Eimag(self, x, verbose=False):
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