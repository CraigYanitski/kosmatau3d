class Molecules(object):
  '''
  This is a class to calculate and solve for the intensity emission by the
  molecules in a PDR.
  '''
  # PRIVATE
  def __init__(self):
    self.__numbers = []
    self.__molecules = []  #private variable for the name of the molecule
    self.__transitions = []
    self.__frequencies = []
   	#self.__intensity = 0      #private variable for the intensity emitted (I_x)
  	#self.__opticalDepth = 0   #private variable for the optical depth (tau_x)
    return
  #def __setI(self, I):
  #	self.__intensity = I
  #	return
  #def __setTau(self, tau):
  #	self.__opticalDepth = tau
  #	return

  #PUBLIC
  def addMolecule(self, molecule, transition, frequency, number):
    self.__numbers.append(number)
    self.__molecules.append(molecule + ' {}'.format(transition))
    self.__transitions.append(transition)
    self.__frequencies.append(frequency)
    return
  def addTransition(self, molecule, transition, frequency, number):
    self.__numbers.append(number)
    self.__molecules.append(molecule + ' {}'.format(transition))
    self.__transitions.append(transition)
    self.__frequencies.append(frequency)
    return
  def getNumbers(self):
    return self.__numbers
  def getMolecules(self):
    return self.__molecules
  def getFrequencies(self):
    return self.__frequencies
  def getTransitions(self):
    return self.__transitions
  #def calculateEmission(self, velocity):
  #	intensity_xi = inten_x_i[sp, vo]+intensityCl[sp][0][ma]*super()._Combination__combination*np.exp(-1/2.*((gbl._globals['compound']['ranges']['vbin'][v]-velocity)/(sigma_cl_j[ma]))**2)
  #	opticalDepth_xi = tau_x_i[sp, vo]+tauCl[sp][0][ma]*super()._Combination__combination*np.exp(-1/2.*((gbl._globals['compound']['ranges']['vbin'][v]-velocity)/(super()._Voxel__sigma))**2)
  #	self.__setI(intensity_xi)
  #	self.__setTau(opticalDepth_xi)
  #	return