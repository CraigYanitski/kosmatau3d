import numpy as np
import constants
class Molecules(object):
  '''
  This is a class to calculate and solve for the intensity emission by the
  molecules in a PDR.
  '''
  # PRIVATE
  def __init__(self):
    self.__fileIndex = []               #a list for the molecule index in the grid
    self.__interpolationIndex = []      #a list for the molecule index in the respective interpolation list
    self.__molecules = []               #list for the name of the molecule
    self.__transitions = []             #list for the transitions being tracked
    self.__frequencies = []             #list for the frequencies being tracked
    self.__wavelengths = []             #list for the wavelengths being tracked
    return
  def __str__(self):
    '''Overloaded print routine to print useful information.'''
    printout = '  {} Molecules:'.format(len(self.__molecules))
    for i in range(len(self.__molecules)):
      printout += '\n    ->{} transition {}, at {} GHz'.format(self.__molecules[i], self.__transitions[i], self.__frequencies[i])
    return printout

  #PUBLIC
  def reset(self):
    self.__fileIndex = []
    self.__interpolationIndex = []
    self.__molecules = []
    self.__transitions = []
    self.__frequencies = []
    self.__wavelengths = []

  def addMolecule(self, molecule, transition, frequency, index):
    self.__fileIndex.append(index)
    self.__interpolationIndex.append(len(self.__molecules))
    self.__molecules.append(molecule)
    self.__transitions.append(transition)
    self.__frequencies.append(frequency)
    self.__wavelengths.append(constants.c/frequency/10**9)
    return
  def addTransition(self, molecule, transition, frequency, index):
    '''This function might not be needed. Initially, I wanted to be able to add
       separate molecules with their respective transitions, but it might be more
       computationally-efficient to just add each molecule...'''
    self.__fileIndex.append(index)
    self.__interpolationIndex.append(len(self.__molecules))
    self.__molecules.append(molecule)
    self.__transitions.append(transition)
    self.__frequencies.append(frequency)
    self.__wavelengths.append(constants.c/frequency/10**9)
    return
  def getFileIndeces(self):
    return np.array(self.__fileIndex, dtype=int)
  def getInterpolationIndeces(self):
    return np.array(self.__interpolationIndex, dtype=int)
  def getMolecules(self):
    return self.__molecules
  def getFrequencies(self):
    return self.__frequencies
  def getWavelengths(self):
    return self.__wavelengths
  def getTransitions(self):
    return self.__transitions