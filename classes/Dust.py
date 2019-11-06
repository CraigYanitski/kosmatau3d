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
    self.__dust = []                  #list for the name of the dust elements being tracked
    self.__transitions = []           #list for the transitions being tracked
    self.__frequencies = []           #list for the frequencies being tracked
    return
  def __str__(self):
    '''Overloaded print routine to print useful information.'''
    printout = '  {} Dust'.format(len(self.__dust))
    for i in range(len(self.__dust)):
      printout += '\n    ->{} transition {}, at {} GHz'.format(self.__dust[i], self.__transitions[i], self.__frequencies[i])
    return printout

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
       computationally-efficient to just add each dust element...'''
    self.__fileIndex.append(index)
    self.__interpolationIndex.append(len(self.__dust))
    self.__dust.append(dust)
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
