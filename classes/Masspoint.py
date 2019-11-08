import numpy as np
from Molecules import *
from Dust import *
class Masspoint(object):
  '''
  This is a class to handle one fractal mass in a combination.
  It will have the associated emission and extinction information from the KOSMA-tau simulations,
  which will be used to contribute to the Combination class' intrinsic intensity and optical depth.
  At the moment, the emission interpolation is performed for each individual species. It will be
  an future update to perform the interpolation for all species at the same time.
  '''
  # PRIVATE
  def __init__(self, species, interpolations, density=0, mass=0, fuv=0, number=0):
    self.__species = species     #list of both moleculular and dust species
    self.__interpolations = interpolations
    self.__density = density
    self.__mass = mass 	#mass of this KOSMA-tau simulation
    self.__FUV = fuv
    self.__number = number    #number of this KOSMA-tau model in the combination
    #self.__FUV = FUVfield()     #the FUV field for this combination of mass points
    self.__intensity_xi = np.array([])       #velocity-averaged intensity of this combination of masspoints
    self.__opticalDepth_xi = np.array([])    #velocity-averaged optical depth of this combination of masspoints
    return
  def __str__(self):
    return 'Simulated KOSMA-tau clump of mass {}'.format(10**self.mass)

  # PUBLIC
  def calculateEmission(self, vRange, vDispersion):
    self.__intensity = np.zeros(vRange.size)
    self.__opticalDepth = np.zeros(vRange.size)
    if self.__number==0:
      self.__intensity_xi = np.zeros((len(vRange),len(vRange)))
      self.__opticalDepth_xi = np.zeros((len(vRange),len(vRange)))
    else:
      for element in self.__species:
        interpolationPoint = [self.__density, self.__mass, self.__FUV.getFUV()]
        print(element)
        if isinstance(element, Molecules):
          self.__intensity_xi = self.__interpolations.interpolateIntensity(interpolationPoint, np.array(element.getInterpolationIndeces())).sum()*self.__number*np.exp(-1/2.*((vRange-vRange.reshape(vRange.size,1))/(vDispersion))**2)
          self.__opticalDepth_xi = self.__interpolations.interpolateTau(interpolationPoint, element.getInterpolationIndeces()).sum()*self.__number*np.exp(-1/2.*((vRange-vRange.reshape(vRange.size,1))/(vDispersion))**2)
        elif isinstance(element, Dust):
          self.__intensity_xi = self.__interpolations.interpolateIntensity(interpolationPoint, element.getInterpolationIndeces()).sum()*self.__number
          self.__opticalDepth_xi = self.__interpolations.interpolateTau(interpolationPoint, element.getInterpolationIndeces()).sum()*self.__number
    return np.stack((self.__intensity_xi,self.__opticalDepth_xi))