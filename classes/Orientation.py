import numpy as np
import cmath
import scipy.interpolate as interpolate
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
  def setLOS(self, grid, x=0, y=0):
    voxels = []
    zPosition = []
    epsilon = []
    kappa = []
    scale = self.__scale*3.086*10**18
    for voxel in grid.allVoxels():
      self.__xLoS = x
      self.__yLoS = y
      x,y,z = voxel.getPosition()
      intensity,tau,fuv = voxel.getEmission()
      if (x<self.__xLoS+self.__scale/2.) and (x>self.__xLoS-self.__scale/2.) and (y<self.__yLoS+self.__scale/2.) and (y>self.__yLoS-self.__scale/2.):
        voxels.append(voxel)
        zPosition.append(z)
        epsilon.append(intensity/(self.__scale))
        kappa.append(tau/(self.__scale))
    i = np.argsort(zPosition)
    print(i)
    for idx in i:
      self.__losVoxels.append(voxels[idx])
    self.__losZ = np.array(zPosition)[i[::-1]]
    self.__epsilon = np.array(epsilon)[i[::-1]]
    self.__epsilonStep = (self.__epsilon[1:]-self.__epsilon[:-1])/(scale)
    self.__kappa = np.array(kappa)[i[::-1]]
    self.__kappaStep = (self.__kappa[1:]-self.__kappa[:-1])/(scale)
    return
  def calculateRadiativeTransfer(self, velocity):
    intensity = np.full((1, len(velocity)), self.__backgroundI)
    scale = self.__scale*3.086*10**18
    # Boolean indeces to separate how the intensity is calculated
    k0 = (self.__kappaStep==0)&(abs(self.__kappa[:-1]*self.__scale)<10**-10)
    kg = self.__kappa[:-1]>10**3*abs(self.__kappaStep)*self.__scale
    kE = ~(k0|kg)
    kEg = ~(k0|kg)&(self.__kappaStep>0)
    kEl = ~(k0|kg)&(self.__kappaStep<0)
    # Variables for utilising the E tilde tables
    a = self.__kappa[:-1]/cmath.sqrt(2*self.__kappaStep)
    b = (self.__kappa[:-1]+self.__kappaStep*self.__scale)/cmath.sqrt(2*self.__scale)
    for i in range(len(self.__losVoxels[:-1])):
      if (self.__kappaStep[i]==0) & (abs(self.__kappa[i]*scale)<10**-10):
        intensity += self.__epsilon[i]*scale+0.5*self.__epsilonStep[i]*scale**2
      elif self.__kappa[i]>10**3*abs(self.kappaStep[i])*scale:
        intensity = np.exp(-self.__kappa[i]*scale) * (intensity + \
                      ((self.__epsilon[i]*self.__kappa[i]+self.__epsilonStep[i]*(self.__kappa[i]*scale-1))/(self.__kappa[i]**2.))*np.exp(self.__kappa[i]*scale) \
                      ((self.__epsilon[i]*self.__kappa[i]-self.__epsilonStep[i])/(self.__kappa[i]**2.)))
      elif self.__kappaStep>0:
        intensity = (self.__epsilonStep[i]/self.__kappaStep[i]*(1-np.exp(-self.__kappa[i]*scale-self.__kappaStep[:-1][i]/2.*scale**2.)) - \
                         (self.__epsilon[i]*self.__kappaStep[i]-self.__epsilonStep[i]*self.__kappa[i])/self.__kappaStep[i] * \
                         cmath.sqrt(np.pi/(2.*abs(self.__kappaStep[i])))* \
                         (np.exp(a[self.__kappaStep[i]>0]**2.-b[self.__kappaStep[i]>0]**2.)*self.__eTildeReal(a[self.__kappaStep[i]>0])-self.__eTildeReal(b[self.__kappaStep[i]>0])) + \
                         intensity*np.exp(-self.__kappa[:-1][i]*scale-self.__kappaStep[i]/2.*scale**2.)).real 
      elif self.__kappaStep<0:
        intensity = (self.__epsilonStep[i]/self.__kappaStep[i]*(1-np.exp(-self.__kappa[i]*scale-self.__kappaStep[i]/2.*scale**2.))\
                         -(self.__epsilon[i]*self.__kappaStep[i]-self.__epsilonStep[i]*self.__kappa[i])/self.__kappaStep[i] * \
                         cmath.sqrt(np.pi/(2.*abs(self.__kappaStep[i])))* \
                         (np.exp(a[self.__kappaStep[i]<0]**2.-b[self.__kappaStep[i]<0]**2.)*self.__eTildeImaginary(a[self.__kappaStep[i]<0])-self.__eTildeImaginary(b[self.__kappaStep[i]<0])) + \
                         intensity*np.exp(-self.__kappa[i]*scale-self.__kappaStep[i]/2.*scale**2.)).real
    return intensity
  def Ereal(x):
    # print 'Ereal: x', x
    if x.imag == 0: x = x.real
    # x should be a real number. remove imaginary party '0j' which
    # prevents ordering
    if x < 0.01:
      return 2*x/np.sqrt(np.pi)
    elif x > 8.0:
      return 1/(np.sqrt(np.pi) * x)
    else:
      return self.__eTildeReal(x)
  def Eimag(x):
    if x == abs(x)*1j:
      # maser case. treated in linear approximation
      x = abs(x)
      return 1 + 2*x/np.sqrt(np.pi)
    else:
      x = abs(x)
      # x needs to be real and positive, i.e. abs(a) or abs(b)
      if x < 0.01:
        return 1 - 2*x/np.sqrt(np.pi)
      elif x > 8.0:
        return 1/(np.sqrt(np.pi) * x)
      else:
        return self.__dTildeImaginary(x)