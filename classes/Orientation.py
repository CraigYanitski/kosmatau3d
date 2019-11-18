import numpy as np
import cmath
import scipy.interpolate as interpolate
import importlib as il
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
  def __init__(self, dimensions, observations, backgroundI=0, interpolation='linear'):
    #self.__dimensions = dimensions
    self.__scale = dimensions.getResolution()
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
    self.__intensity = 0
    #observations = Observations()
    eTildeReal = observations.eTildeReal
    eTildeImaginary = observations.eTildeImaginary
    self.__eTildeReal = interpolate.interp1d(eTildeReal[0], eTildeReal[1], kind=interpolation)
    self.__eTildeImaginary = interpolate.interp1d(eTildeImaginary[0], eTildeImaginary[1], kind=interpolation)
    return

  #PUBLIC
  def reloadModules(self):
    self.__scale.reloadModules()
    for voxel in self.__losVoxels: voxel.reloadModules()
    return
  def setLOS(self, grid, x=0, y=0, verbose=True):
    voxels = []
    zPosition = []
    epsilon = []
    kappa = []
    scale = self.__scale*3.086*10**18
    self.__xLoS = x
    self.__yLoS = y
    for voxel in grid.allVoxels():
      x,y,z = voxel.getPosition()
      intensity,tau,fuv = voxel.getEmission()
      if (x<(self.__xLoS+self.__scale/2.)) and (x>(self.__xLoS-self.__scale/2.)) and (y<(self.__yLoS+self.__scale/2.)) and (y>(self.__yLoS-self.__scale/2.)):
        voxels.append(voxel)
        zPosition.append(z)
        #if z:
        factor = np.exp(-z**2/500.**2)
        epsilon.append(factor*intensity/(self.__scale))
        kappa.append(factor*tau/(self.__scale))
        #else:
        #  epsilon.append(intensity/(self.__scale))
        #  kappa.append(tau/(self.__scale))
        if verbose: print(tau, self.__scale)
    i = np.argsort(zPosition)
    print(i)
    for idx in i:
      self.__losVoxels.append(voxels[idx])
    self.__losZ = np.array(zPosition)[i[::-1]]
    self.__epsilon = np.array(epsilon, dtype=np.float)[i[::-1]]
    self.__epsilonStep = (self.__epsilon[1:]-self.__epsilon[:-1])/(scale)
    self.__kappa = np.array(kappa, dtype=np.float)[i[::-1]]
    self.__kappaStep = (self.__kappa[1:]-self.__kappa[:-1])/(scale)
    return
  def calculateRadiativeTransfer(self, velocity, verbose=True):
    intensity = np.full((1, 10), self.__backgroundI)
    scale = self.__scale*3.086*10**18
    # Boolean indeces to separate how the intensity is calculated
    k0 = (self.__kappaStep==0)&(abs(self.__kappa[:-1]*self.__scale)<10**-10)
    kg = self.__kappa[:-1]>10**3*abs(self.__kappaStep)*self.__scale
    kE = ~(k0|kg)
    kEg = ~(k0|kg)&(self.__kappaStep>0)
    kEl = ~(k0|kg)&(self.__kappaStep<0)
    # Variables for utilising the E tilde tables
    if verbose: print('\nkappa, kappa step:\n', self.__kappa, '\n', self.__kappaStep)
    a = (self.__kappa[:-1]/np.sqrt(2*self.__kappaStep.astype(np.complex))).astype(np.complex)
    b = ((self.__kappa[:-1]+self.__kappaStep*self.__scale)/np.sqrt(2*self.__scale)).astype(np.complex)
    for i in range(len(self.__losVoxels[:-1])):
      if (self.__kappaStep[i]==0).any() & (abs(self.__kappa[i]*scale)<10**-10).any():
        intensity += self.__epsilon[i]*scale+0.5*self.__epsilonStep[i]*scale**2
      elif (self.__kappa[i]>10**3*abs(self.__kappaStep[i])*scale).any():
        intensity = np.exp(-self.__kappa[i]*scale) * (intensity + \
                      ((self.__epsilon[i]*self.__kappa[i]+self.__epsilonStep[i]*(self.__kappa[i]*scale-1))/(self.__kappa[i]**2.))*np.exp(self.__kappa[i]*scale) - \
                      ((self.__epsilon[i]*self.__kappa[i]-self.__epsilonStep[i])/(self.__kappa[i]**2.)))
      elif (self.__kappaStep[i]>0).any():
        print('\na, b:\n', a[i], '\n', b[i])
        aE = np.array(list(self.Ereal(a[i])))
        bE = np.array(list(self.Ereal(b[i])))
        intensity = (self.__epsilonStep[i]/self.__kappaStep[i]*(1-np.exp(-self.__kappa[i]*scale-self.__kappaStep[:-1][i]/2.*scale**2.)) - \
                         (self.__epsilon[i]*self.__kappaStep[i]-self.__epsilonStep[i]*self.__kappa[i])/self.__kappaStep[i] * \
                         np.sqrt(np.pi/(2.*abs(self.__kappaStep[i]))) * \
                         (np.exp(a[i]**2.-b[i]**2.)*aE-bE) + \
                         intensity*np.exp(-self.__kappa[:-1][i]*scale-self.__kappaStep[i]/2.*scale**2.)).real 
      elif (self.__kappaStep[i]<0).any():
        aE = np.array(list(self.Eimag(a[i])))
        bE = np.array(list(self.Eimag(b[i])))
        intensity = (self.__epsilonStep[i]/self.__kappaStep[i]*(1-np.exp(-self.__kappa[i]*scale-self.__kappaStep[i]/2.*scale**2.))\
                         -(self.__epsilon[i]*self.__kappaStep[i]-self.__epsilonStep[i]*self.__kappa[i])/self.__kappaStep[i] * \
                         np.sqrt(np.pi/(2.*abs(self.__kappaStep[i])))* \
                         (np.exp(a[i]**2.-b[i]**2.)*aE-bE) + \
                         intensity*np.exp(-self.__kappa[i]*scale-self.__kappaStep[i]/2.*scale**2.)).real
    return intensity
  def Ereal(self, x, verbose=True):
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
  def Eimag(self, x, verbose=True):
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