import numpy
class Masspoint(object):
  '''
  This is a class to handle one fractal mass in a combination.
  It will have the associated emission and extinction information from the KOSMA-tau simulations,
  which will be used to contricute to the Combination class' intrinsic intensity and optical depth.
  '''
  # PRIVATE
  def __init__(self, mass, combination):
    self.__listSpecies = []     #list of both moleculular and dust species
    self.__mass = mass 	#mass of this KOSMA-tau simulation
    self.__number = number    #number of this KOSMA-tau model
    #self.__FUV = FUVfield()     #the FUV field for this combination of mass points
    self.__intensity = np.array([])       #velocity-averaged intensity of this combination of masspoints
    self.__opticalDepth = np.array([])    #velocity-averaged optical depth of this combination of masspoints
    return

  # PUBLIC
  def calculateEmission(self, species, vrange):
    self.__intensity = np.zeros(vrange.size)
    self.__opticalDepth = np.zeros(vrange.size)
    for element in self.__listSpecies:
      if isinstance(element, Molecule):
        self.__intensity_xi += constants.interpolateI(element.getName())*self.__number*np.exp(-1/2.*((vVox-vObs)/(sigmaV))**2)
  	    self.__opticalDepth_xi += constants.interpolateTau(element.getName())*self.__number*np.exp(-1/2.*((vVox-vObs)/(sigmaV))**2)
      elif isinstance(element, Dust):
        self.__intensity_xi += constants.interpolateI(element.getName()) * self.__number
        self.__opticalDepth_xi += constants.interpolateTau(element.getName()) * self.__number
      return np.stack((self.__intensity,self.__opticalDepth))

