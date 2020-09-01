import numpy as np

import constants
import interpolations
import species

def addMolecules(molecules):
  # This is currently setup to accept an acsii input only, in the format:
  # '{molecule} {transition}'. This will be cross-referenced with the molecules in
  # the model to determine what needs to be interpolated, and an error will
  # be raised if the molecular transition is not available in the model.

  if not isinstance(molecules, list) and not isinstance(molecules, np.ndarray):
    molecules = [molecules]

  for molecule in molecules:

    if molecule in species.molecules:
      pass

    elif molecule in constants.molecules:
      species.molecules.append(molecule)
      species.moleculeIndeces.append(np.where(constants.molecules==molecule))

      if molecule.split()[0]=='13CO':
        file = open(constants.MOLECULARPATH+'13C16O.INP').readlines()

      elif molecule.split()[0]=='H3O+':
        file = open(constants.MOLECULARPATH+'p-H3O+.INP').readlines()

      else:
        file = open(constants.MOLECULARPATH+'{}.INP'.format(molecule.split(' ')[0])).readlines()

      info = file[1].split()
      energies = []
      levels = []
      for i in range(int(info[0])):
        levels.append(file[2+i].split()[0])
        energies.append(float(file[2+i].split()[2]))
      energies = np.array(energies)
      levels = np.array(levels)
      transition = file[1+int(info[0])+int(molecule.split(' ')[1])].split()[:2]
      nINI = levels==transition[0]
      nFIN = levels==transition[1]
      species.moleculeFrequencies.append(2.99792458*10**10*(energies[nINI]-energies[nFIN]))
      species.moleculeWavelengths.append(1./100/(energies[nINI]-energies[nFIN]))

      # interpolations.initialise()

    else:
      print('MODEL ERROR: Molecular transition {} not available. Please use a different grid or select a different molecule.')

  return