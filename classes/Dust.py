import numpy as np
class Dust(object):
  '''
  This is a class to calculate and solve for the intensity emission by the
  dust in a PDR.
  '''
  # PRIVATE
  def __init__(self):
    self.__fileIndex = []             #a list for the dust index in the grid
    self.__interpolationIndex = []    #a list for the dust index in the respective interpolation list
    self.__dust = []  #private variable for the name of the dust element
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
  def addDust(self, dust, transition, frequency, index):
    self.__fileIndex.append(index)
    self.__interpolationIndex.append(len(self.__dust))
    self.__dust.append(dust + ' {}'.format(transition))
    self.__transitions.append(transition)
    self.__frequencies.append(frequency)
    return
  def addTransition(self, dust, transition, frequency, index):
    '''This function might not be needed. Initially, I wanted to be able to add
       separate dust with their respective transitions, but it might be more
       computationally-efficient to just add each molecule...'''
    self.__fileIndex.append(index)
    self.__interpolationIndex.append(len(self.__dust))
    self.__molecules.append(dust + ' {}'.format(transition))
    self.__transitions.append(transition)
    self.__frequencies.append(frequency)
    return
  def getFileIndeces(self):
    return np.array(self.__fileIndex, dtype=int)
  def getInterpolationIndeces(self):
    return np.array(self.__interpolationIndex, dtype=int)
  def getDust(self):
    return np.array(self.__dust)
  def getFrequencies(self):
    return np.array(self.__frequencies)
  def getTransitions(self):
    return np.array(self.__transitions)
  #def calculateEmission(self):
  #	for mass in super()._Combination__masspoints
  #    intensity_xi[vo] = inten_x_i[vo] + intensityCl[0][ma] * super()._Combination__combination[ma]
  #    tau_x_i[vo] = tau_x_i[vo] + tauCl[0][ma] * super()._Combination__combination[ma]
  #  return