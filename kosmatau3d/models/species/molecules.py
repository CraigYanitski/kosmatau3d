import numpy as np

from kosmatau3d.models import constants
#from kosmatau3d.models import interpolations
from kosmatau3d.models import species


def add_transitions(transitions, interclump=False):
    # This is currently setup to accept an acsii input only, in the format:
    # '{molecule} {transition}'. This will be cross-referenced with the molecules in
    # the model to determine what needs to be interpolated, and an error will
    # be raised if the molecular transition is not available in the model.
  
    if not isinstance(transitions, list) and not isinstance(transitions, np.ndarray):
        transitions = [transitions]
  
    if transitions[0] == 'all':
        transitions = constants.transitions
  
    if interclump:
        species.interclump_transitions = []
        species.interclump_transition_indeces = []
        species.interclump_interclump_transition_indeces = []
        species.interclump_transition_frequencies = []
        species.interclump_transition_wavelengths = []
    else:
        species.clump_transitions = []
        species.clump_transition_indeces = []
        species.clump_interclump_transition_indeces = []
        species.clump_transition_frequencies = []
        species.clump_transition_wavelengths = []
  
    for transition in transitions:
  
        if transition in species.interclump_transitions and \
                transition in species.clump_transitions:
            pass
    
        elif (transition in constants.transitions) and  \
             (transition in constants.interclump_transitions):
            if interclump:
                species.interclump_transitions.append(transition)
                species.interclump_transition_indeces.append(
                        np.where(constants.interclump_transitions == transition))
            else:
                species.clump_transitions.append(transition)
                species.clump_transition_indeces.append(
                        np.where(constants.transitions == transition))
      
            if transition.split()[0] == '13CO':
                file = open(constants.MOLECULARPATH+'C13O16.INP').readlines()
      
            elif transition.split()[0] == 'H3O+':
                file = open(constants.MOLECULARPATH+'p-H3O+.INP').readlines()
      
            else:
                file = open(constants.MOLECULARPATH+'{}.INP'.format(
                    transition.split(' ')[0])).readlines()
      
            info = file[1].split()
            energies = []
            levels = []
            for i in range(int(info[0])):
                levels.append(file[2+i].split()[0])
                energies.append(float(file[2+i].split()[2]))
            energies = np.array(energies)
            levels = np.array(levels)
            transition_levels = file[1+int(info[0])+int(transition.split(' ')[1])].split()[:2]
            nINI = levels == transition_levels[0]
            nFIN = levels == transition_levels[1]
            if interclump:
                species.interclump_transition_frequencies.append((2.99792458*10**10
                    *(energies[nINI]-energies[nFIN]))[0])
                species.interclump_transition_wavelengths.append((1./100/(energies[nINI]
                    -energies[nFIN]))[0])
            else:
                species.clump_transition_frequencies.append((2.99792458*10**10
                    *(energies[nINI]-energies[nFIN]))[0])
                species.clump_transition_wavelengths.append((1./100/(energies[nINI]
                    -energies[nFIN]))[0])
    
        else:
            print(f'MODEL ERROR: Species transition {transition} not available. '
                  'Please use a different grid or select a different molecule.')
  
    return


def species_indeces(transitions=[], interclump=False):

    if isinstance(transitions, str):
        transitions = [transitions]
  
    if len(transitions) == 0:
        return
  
    elif transitions[0] == 'all' and interclump:
        indeces = np.arange(len(species.interclump_transitions))

    elif transitions[0] == 'all':
        indeces = np.arange(len(species.clump_transitions))

    else:
        indeces = []
        for transition in transitions:
            if interclump:
                i = np.where(np.asarray(species.interclump_transitions) == transition)[0]
            else:
                i = np.where(np.asarray(species.clump_transitions) == transition)[0]
            if len(i) == 0:
                transition = transition.split()
                if transition[0] == 'dust':
                    i = np.where(constants.dust_names == transition[1])[0]
                    if len(i) == 0:
                        print('Dust wavelength {} not found.'.format(transition))
                        continue
                    indeces.append(i[0])
                else:
                    print('Species transition {} not found.'.format(transition))
                    continue
            else:
                indeces.append(constants.wavelengths[constants.n_dust].size+i[0])
  
    return indeces
