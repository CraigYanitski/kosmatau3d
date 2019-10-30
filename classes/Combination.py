import Masspoint
class Combination(object):
  '''
  This is a class to handle a combination of fractal masses in an ensemble.
  It will have its associated probability, which will scale its intrinsic
  intensity and optical depth. It returns a tuple of the combination's
  probability, intensity, optical depth, and FUV field.
  '''
  # PRIVATE
  def __init__(self, combination=[], masses=[], probability=0):
    self.__combination = combination 	#list of the number of each masspoint
    self.__listMasspoints = masses    #list of masses in the combination
    self.__FUV = 0                    #the FUV field for this combination of mass points
    self.__intensity = []             #velocity-averaged intensity of this combination of masspoints
    self.__opticalDepth = []          #velocity-averaged optical depth of this combination of masspoints
    self.__probability = 0           #the probability of this combination of masspoints
    self.__listElements = elements    #list of the modelled elements; taken from the overall Model() classs
    return
  def __setFUV(self, fuvField):
    self.__FUV = fuvField
    return

  # PUBLIC
  #def addMolecule(self, element):
  #  self.__listMolecules.append(MolecularEmission(element))
  #  return
  #def addDust(self, element):
  #  self.__listDust.append(DustEmission(element))
  #  return
  def addFUV(self, fuvField):
    self.__setFUV(fuvField)
    return
  def addMasspoint(self, masspoint, number):
    self.__listMasspoints.append(masspoint)
    self.__combination.append(number)
    return
  def calculateEmission(self, vrange):
    for i,masspoint in enumerate(self.__listMasspoints):
      (intensity,opticalDepth) = masspoint.calculateEmission(species, self.__combination[i], vrange)
      self.__intensity.append(intensity)
      self.__opticalDepth.append(opticalDepth)
    return
    self.__intensity = np.array(self.__intensity)
    self.__opticalDepth = np.array(self.__opticalDepth)
  def getScaledCombinationEmission(self):
    return self.__probability*np.stack((self.__intensity, np.exp(-self.__opticalDepth), np.exp(-self.__FUV)))
