import numpy as np

import observations

from .Dust import Dust
from .Molecules import Molecules

molecules = Molecules()
dust = Dust()

moleculeWavelengths = []

speciesNames = None

def reset():
  molecules.reset()
  moleculeWavelengths = []
  speciesNames = None
  return

def addDust(dustElement, transition):
  (numbers,species,transitions,frequencies) = observations.speciesData
  i = (species==dustElement)&(transitions==transition)
  if transition in dust.getTransitions():
    dust.addTransition(dustElement, transition, frequencies[i], numbers[i])
  else:
    dust.addDust(dustElement, transition, frequencies[i], numbers[i])
  speciesNames = np.append(molecules.getMolecules(), dust.getDust())
  return

def addMolecule(molecule, transition):
  (numbers,species,transitions,frequencies) = observations.speciesData
  i = (species==molecule)&(transitions==transition)
  if molecule in molecules.getMolecules():
    molecules.addTransition(molecule, transition, frequencies[i], numbers[i])
  else:
    molecules.addMolecule(molecule, transition, frequencies[i], numbers[i])
  speciesNames = np.append(molecules.getMolecules(), dust.getDust())
  moleculeWavelengths = molecules.getWavelengths()
  return
