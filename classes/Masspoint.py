import numpy
class Masspoint(object):
  '''
  This is a class to handle one fractal mass in a combination.
  It will have the associated emission and extinction information from the KOSMA-tau simulations,
  which will be used to contribute to the Combination class' intrinsic intensity and optical depth.
  At the moment, the emission interpolation is performed for each individual species. It will be
  an future update to perform the interpolation for all species at the same time.
  '''
  # PRIVATE
  def __init__(self, species, interpolations, density, mass, fuv, number):
    self.__species = species     #list of both moleculular and dust species
    self.__interpolations = interpolations
    self.__density = density
    self.__mass = mass 	#mass of this KOSMA-tau simulation
    self.__fuv = fuv
    self.__number = number    #number of this KOSMA-tau model
    #self.__FUV = FUVfield()     #the FUV field for this combination of mass points
    self.__intensity = np.array([])       #velocity-averaged intensity of this combination of masspoints
    self.__opticalDepth = np.array([])    #velocity-averaged optical depth of this combination of masspoints
    return

  # PUBLIC
  def calculateEmission(self, species, vRange, vDispersion):
    self.__intensity = np.zeros(vrange.size)
    self.__opticalDepth = np.zeros(vrange.size)
    for element in self.__species:
      interpolationPoint = [self.__density, element, self.fuv]
      if isinstance(element, Molecule):
        self.__intensity_xi += self.__interpolations.interpolateIntensity(interpolationPoint)*self.__number*np.exp(-1/2.*((vRange-vRange.reshape(vRange.size,1))/(vDispersion))**2)
  	    self.__opticalDepth_xi += self.__interpolations.interpolateTau(interpolationPoint)*self.__number*np.exp(-1/2.*((vRange-vRange.reshape(vRange.size,1))/(vDispersion))**2)
      elif isinstance(element, Dust):
        self.__intensity_xi += self.__interpolations.interpolateIntensity(interpolationPoint) * self.__number
        self.__opticalDepth_xi += self.__interpolations.interpolateTau(interpolationPoint) * self.__number
      return np.stack((self.__intensity,self.__opticalDepth))