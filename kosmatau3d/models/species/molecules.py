import numpy as np

from .. import constants
from .. import interpolations
from .. import species

def addMolecules(molecules):
  # This is currently setup to accept an acsii input only, in the format:
  # '{molecule} {transition}'. This will be cross-referenced with the molecules in
  # the model to determine what needs to be interpolated, and an error will
  # be raised if the molecular transition is not available in the model.

  if not isinstance(molecules, list) and not isinstance(molecules, np.ndarray):
    molecules = [molecules]

  if molecules[0]=='all':
    molecules = constants.molecules

  species.molecules = []
  species.moleculeIndeces = []
  species.moleculeFrequencies = []
  species.moleculeWavelengths = []

  for molecule in molecules:

    if molecule in species.molecules:
      pass

    elif molecule in constants.molecules:
      species.molecules.append(molecule)
      species.moleculeIndeces.append(np.where(constants.molecules==molecule))

      if molecule.split()[0]=='13CO':
        file = open(constants.MOLECULARPATH+'C13O16.INP').readlines()

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

def speciesIndeces(molecules=[]):

  if isinstance(molecules, str):
    molecules = [molecules]

  if len(molecules)==0:
    return

  elif molecules[0]=='all':
    indeces = np.arange(len(species.molecules))

  else:
    indeces = []
    for mol in molecules:
      i = np.where(np.asarray(species.molecules)==mol)[0]
      if len(i)==0:
        mol = mol.split()
        if mol[0]=='dust':
          i = np.where(constants.dustNames==mol[1])[0]
          if len(i)==0:
            print('Dust wavelength {} not found.'.format(mol))
            continue
          indeces.append(i[0])
        else:
          print('Species transition {} not found.'.format(mol))
          continue
      else:
        indeces.append(constants.wavelengths[constants.nDust].size+i[0])

  return indeces