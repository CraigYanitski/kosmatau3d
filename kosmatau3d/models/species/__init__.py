import numpy as np

from .. import observations

from .molecules import *

molecules = []
moleculeIndeces = []
moleculeFrequencies = []
moleculeWavelengths = []

def addMoleculeOld(molecule):
  # This is currently setup to accept an acsii input only, in the format:
  # '{molecule} {transition}'. This will be cross-referenced with the molecules in
  # the model to determine what needs to be interpolated, and an error will
  # be raised if the molecular transition is not available in the model.

  if molecule in molecules:
    pass

  elif molecule in constants.molecules:
    molecules.append(molecule)
    moleculeIndeces.append(np.where(constants.molecules==molecule)[0][0])

  else:
    print('Molecular transition {} not available. Please use a different grid or select a different molecule.')

  return