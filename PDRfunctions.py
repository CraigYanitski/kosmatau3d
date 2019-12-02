# class functions
# collection of different functions and classes needed on PDR3D.py

import os
import numpy as np
from tkinter import *
import mpl_toolkits.mplot3d.axes3d as p3
import array
import datetime
import time
from copy import copy
import textwrap
from astropy.io import fits
# from visual import *
# from axis_angle_to_quaternion import axis_angle_to_quaternion
#from matplotlib.figure import Figure
import matplotlib.pyplot as plt #(plt)
# edit: Craig, 14.10.2019
import globals as gbl
#end edit
import fuv_absorption as fuva
import line_absorption_and_emission_velocity_dependent as lae  # script for radiative transfer
import rad_transfer_fixed_velocities as rt
from PDR_GUI_popUp import PDR_GUI  # GUI for maps

def directoryManager(directory, key, path): # add and update directory paths
  try:
    os.mkdir(path) # try to create directory
    directory[key]=path # add to directory
  except OSError as e:
    if e.errno == 17 and os.path.isdir(path):
      # errno == 17 : path exists
      # os.path.isdir returns true if path is en exisiting directory
      # print path,': directory exists'
      directory[key]=path # add to directory
    else:
      raise 
  return directory   

def progressBar(iteration, total, prefix='', suffix='', decimals=1, length=100, fill='#', printEnd='\r'):
  '''
  A progress routine found on StackOverflow.
  Implemented on 21.10.2019 by Craig.
  '''
  percent = ('{0:.' + str(decimals) + 'f}').format(100 * iteration/float(total))
  filledLength = int(length*iteration//total)
  bar = fill*filledLength+' '*(length-filledLength)
  print('\r{} {}it|{}| {}% {}'.format(prefix, iteration, bar, percent, suffix), end=printEnd)
  if iteration==total: print()
  return

class Ensemble:
  def __init__(self, number, metal = 100, 
                rho_ens_clumps = 10**5,  \
                rho_ens_inter = 10**4, \
                Mens_clumps = 10, \
                Mens_inter = 10,\
                FUV = 1, \
                Ml_cl = 0.001, \
                Mu_cl = 1, \
                Ml_inter = 0.001, \
                Mu_inter = 1,\
                B = 10., \
                velocity = np.array([0,0,0], float), \
                v_dispersion = 0., \
                inten = 0,  inten_tot = 0, \
                tau = 0.,  tau_tot = 0):
    self.number = number
    # for calculation of intensity
    self.metal = metal
    self.rho_ens_clumps = rho_ens_clumps
    self.rho_ens_inter = rho_ens_inter
    self.Mens_clumps = Mens_clumps
    self.Mens_cl_sum = None
    self.Mens_inter = Mens_inter
    self.Mens_int_sum = None
    self.FUV = FUV
    self.Ml_cl = Ml_cl
    self.Mu_cl = Mu_cl
    self.Ml_inter = Ml_inter
    self.Mu_inter = Mu_inter
    self.B = B
    # other parameters
    self.velocity = velocity
    # systematic velocity
    self.Mens = None
    self.Mens_sum = None
    self.rho_ens = None
    self.Mu = None
    self.Ml = None
    self.Afuv = None # averaged fuv extinction
    self.Afuv_tot = None # total averaged fuv extinction
    # species-specific parameters 
    self.inten = inten
    self.inten_tot = inten_tot
    self.tau = tau # define global absorption coefficient. 
    # value defines AREa BELOW gauss curve, in [cm^-1*km/s]!!!
    self.tau_tot = tau_tot

def ensemble_parameter():
  # calculate upper and lower mass threshold for ensemble interval
  #from PDR import _globals
  gbl._globals['compound']['number'] = []
  for i in range(gbl._globals['compound']['npoints']):
    B = gbl._globals['compound']['ens'][i].B
    alpha = gbl._globals['constants']['alpha']
    gamma = gbl._globals['constants']['gamma'] 
    MuLog = np.log10(gbl._globals['compound']['ens'][int(i)].Mu) * B
    MlLog = np.log10(gbl._globals['compound']['ens'][int(i)].Ml) * B
    # rho_ens = gbl._globals['compound']['ens'][int(i)].rho_ens
    Masses = np.arange(MlLog, MuLog + B, B)
    numberInt = len(Masses) 
    ml = np.zeros(numberInt) 
    mu = np.zeros(numberInt) # mass intervals
    j = 0
    for m in Masses:
      ml[j] = 10**(m/10.) * (1 - B**(1 - alpha))**(-1./(1 - alpha))
      mu[j] = ml[j] * B
      j = j+1
    # C = (4*np.pi/3.*(2 - alpha)/(1 + 3./gamma - alpha)*\
      # ((mu[numberInt - 1]**(1 + 3./gamma-alpha) - \
      # ml[0]**(1 + 3./gamma - alpha))/\
      # (mu[numberInt - 1]**(2 - alpha) - ml[0]**(2 - alpha))) * \
      # rho_ens * gbl._globals['constants']['M_H']???/\????
      # gbl._globals['constants']['M_sol'])**(gamma/3.)         
        # calculate contant C (Cubick et al) for the ensemble
        # C in M_s cm^-gamma
    # gbl._globals['compound']['ens'][i].C = C
    gbl._globals['compound']['ens'][i].ml = ml[0]
    gbl._globals['compound']['ens'][i].mu = mu[numberInt - 1]
  print('function ensemble_parameter calculated mu, ml for first ensemble: ',\
        gbl._globals['compound']['ens'][0].mu, '\n',\
        gbl._globals['compound']['ens'][0].ml)
  # print 'and C in M_s*cm^-gamma for first ensemble',\
          # gbl._globals['compound']['ens'][0].C
  # pause = input()


def find_transition_number_and_frequency():
  # find transition number as defined in
  # SetUpOrionBarModelEnvironment.v.1.1.nb
  #from PDR import _globals  
  gbl._globals['compound']['number'] = []
  gbl._globals['compound']['frequency'] = []
  for i in range(gbl._globals['compound']['nspe']):
    # for each species 
    number = -1
    species = gbl._globals['compound']['species'][i]
    transition = gbl._globals['compound']['transition'][i]
    if species == 'C+': number = 1
    elif species == 'C':  
      if transition > 3 or transition < 1:
        import sys    
        sys.exit('transition does not exist...exiting...')
      else: number = 1 + transition
    elif species == 'O':
      if transition > 3 or transition < 1:
        import sys    
        sys.exit('transition does not exist...exiting...')
      else: number = 4 + transition
    elif species == 'CO': 
      if transition > 49 or transition < 1:
        import sys    
        sys.exit('transition does not exist...exiting...')
      else: number = 7 + transition
    elif species == '13CO': 
      if transition > 49 or transition < 1:
        import sys    
        sys.exit('transition does not exist...exiting...')
      else: number = 56 + transition
    elif species == '13C+': number = 106
    elif species == '13C':
      if transition > 3 or transition < 1:
        import sys    
        sys.exit('transition does not exist...exiting...')
      else: number = 106 + transition
    elif species == 'HCO+':
      if transition > 15 or transition < 1:
        import sys    
        sys.exit('transition does not exist...exiting...')
      else: number = 109 + transition
    elif species == 'H13CO+':
      if transition > 30 or transition < 1:
        import sys    
        sys.exit('transition does not exist...exiting...')
      else: number = 124 + transition
    elif species == 'H3O+':
      if transition > 17 or transition < 1:
        import sys    
        sys.exit('transition does not exist...exiting...')
      else: number = 154 + transition
    elif species == 'C18O':
      if transition > 49 or transition < 1:
        import sys    
        sys.exit('transition does not exist...exiting...')
      else: number = 171 + transition
    elif species == 'Dust':
      if transition > 51 or transition < 1:
        import sys
        sys.exit('PDR_functions:this Dust transition is not implemented...exiting...')
      else: number = 191 + transition
#    elif species == 'CH+': 
#        if transition > 10 or transition < 1:
#            import sys    
#            sys.exit('transition does not exist...exiting...')
#        else: number = 171 + transition
    if number == -1:
      print('species: ', species)
      print('transition: ', transition)
      import sys    
      sys.exit('species/transition combination does not exist...exiting...')            
    gbl._globals['compound']['number'].append(number) 
    with open(gbl.INPUTPATH+'frequencies.dat','r') as f:
      for i, line in enumerate(f):
        if i == number:
          try: freq = float(line.split()[3]) 
          # read frequency from frequency.dat    
          except IndexError:
            import sys
            sys.exit('Problem :( frequency information needs to be added to file frequencies.dat ...exiting...')          
          gbl._globals['compound']['frequency'].append(freq)
          #print 'hier'
          
        #elif i > number:
          #print 'break'
          #break
    #gbl._globals['compound']['frequency'].append(freq)
  if len(gbl._globals['compound']['transition']) != len( gbl._globals['compound']['frequency']):
    print('frequencies found ' , len( gbl._globals['compound']['frequency']))
    print('frequencies seeked ', len(gbl._globals['compound']['transition']))
    import sys
    sys.exit('Problem :( Frequency information does not match. Maybe needs to be added to file frequencies.dat ...exiting...')
  
    # print 'number', gbl._globals['compound']['number']
    # print 'frequency',  gbl._globals['compound']['frequency'] 
    # pause = input('..ok?..remove again..in PDRfunctions..')
    return gbl._globals 
      #gbl._globals['compound']['frequency'].append(freq)
  # print 'number', gbl._globals['compound']['number']
  # print 'frequency',  gbl._globals['compound']['frequency'] 
  # pause = input('..ok?..remove again..in PDRfunctions..')
  return gbl._globals

def runKOSMAt():
  RUN = True
  while RUN:
    print('\nPlease select an option...')
    print('  (1) Setup model')
    print('  (2) Initialise velocities')
    print('  (3) Calculate FUV absorption')
    print('  (4) Calculate line absorption and emission')
    print('  (5) Calculate radiative transfer')
    print('  (6) Show plots')
    print('  (7) Execute all')
    INPUT = input('\nSelection: ')
    if INPUT=='1': setupModel()
    elif INPUT=='2': initVelocities()
    elif INPUT=='3': calculateFUV()
    elif INPUT=='4': lineAE()
    elif INPUT=='5': radTransfer()
    elif INPUT=='6':
      PLOTS = True
      while PLOTS:
        print('\nWhich plot would you like to see?')
        print('  (1) Total mass')
        print('  (2) Clump mass')
        print('  (3) Interclump mass')
        print('  (4) Total Intensity')
        print('  (5) Clump Intensity')
        print('  (6) Interclump intensity')
        print('  (7) Density')
        print('  (8) Z velocity')
        print('  (9) Total velocity')
        print('  (10) Dispersion')
        print('  (11) Total dispersion')
        print('  (12) FUV distribution')
        PLOT = input('Selection: ')
        if PLOT=='1':gbl.plotFlags['Total mass'] = True
        elif PLOT=='2': gbl.plotFlags['Clump mass'] = True
        elif PLOT=='3': gbl.plotFlags['Interclump mass'] = True
        elif PLOT=='4': gbl.plotFlags['Total intensity'] = True
        elif PLOT=='5': gbl.plotFlags['Clump intensity'] = True
        elif PLOT=='6': gbl.plotFlags['Interclump intensity'] = True
        elif PLOT=='7': gbl.plotFlags['Density'] = True
        elif PLOT=='8': gbl.plotFlags['Z velocity'] = True
        elif PLOT=='9': gbl.plotFlags['Total velocity'] = True
        elif PLOT=='10': gbl.plotFlags['Dispersion'] = True
        elif PLOT=='11': gbl.plotFlags['Total dispersion'] = True
        elif PLOT=='12': gbl.plotFlags['FUV distribution'] = True
        else: break
        try: showPlots()
        except: continue
        gbl.plotFlags['Total mass'] = gbl.plotFlags['Clump mass'] = gbl.plotFlags['Interclump mass'] = \
          gbl.plotFlags['Total Intensity'] = gbl.plotFlags['Clump Intensity'] = gbl.plotFlags['Interclump Intensity'] = \
          gbl.plotFlags['Density'] = gbl.plotFlags['Z velocity'] = gbl.plotFlags['Total velocity'] = \
          gbl.plotFlags['Dispersion'] = gbl.plotFlags['Total dispersion'] = gbl.plotFlags['FUV distribution'] = \
          gbl.plotFlags['Clump mass'] = gbl.plotFlags['Z velocity'] = gbl.plotFlags['Total velocity'] = \
          gbl.plotFlags['Interclump mass'] = gbl.plotFlags['Total dispersion'] = gbl.plotFlags['FUV distribution'] = \
          gbl.plotFlags['Clump intensity'] = gbl.plotFlags['Z velocity'] = gbl.plotFlags['Total velocity'] = \
          gbl.plotFlags['Interclump intensity'] = gbl.plotFlags['Total dispersion'] = gbl.plotFlags['FUV distribution'] = \
          gbl.plotFlags['Z velocity'] = gbl.plotFlags['Z velocity'] = gbl.plotFlags['Total velocity'] = \
          gbl.plotFlags['Total velocity'] = gbl.plotFlags['Total dispersion'] = gbl.plotFlags['FUV distribution'] = \
          gbl.plotFlags['FUV distribution'] = gbl.plotFlags['Z velocity'] = gbl.plotFlags['Total velocity'] = \
          gbl.plotFlags['Total dispersion'] = gbl.plotFlags['Total dispersion'] = gbl.plotFlags['FUV distribution'] = False
    elif INPUT=='7':
      setupModel()
      initVelocities()
      calculateFUV()
      lineAE()
      #radTransfer()
      #showPlots()
    else: RUN = False
  print('\nAuf Wiedersehen !!\n')
  return

def globalInit():
  return

def setupModel(timed=True):
  '''
  This is a separate function to import and setup the model.
  Created on 17.10.2019 by Craig.
  '''
  t_init = time.time()
  # time mark: start of Simulation
  t_start = datetime.datetime.utcnow()
  if gbl._globals['verbose']: print('time mark: starting 3D PDR simulation at ', t_start)
  # start time in absulute seconds after 1970
  gbl._globals['runtimes']['start'] = time.time()
  # list for times of line_absortption of each voxel
  gbl._globals['runtimes']['line_absorption'] = []
  ##########################################################
  ### import compound including user defined parameters
  t = datetime.datetime.utcnow()  
  if gbl._globals['verbose']: print('Importing compound information...elapsed time: ',\
    (t-t_start).seconds,' seconds, ',\
    (t-t_start).microseconds, ' microseconds')
  __import__(gbl._globals['Files']['compound_filename'])
  #pause = input('model gut durch')
  __import__(gbl._globals['Files']['parameter_filename'])
  #pause = input('transitions and observer durch')
  print('\nModel set up successfully.')
  if timed: print('  in', time.time()-t_init, 's')
  return

def initVelocities(timed=True):
  '''
  This is a function to separate this part of the code.
  Created on 17.10.2019 by Craig.
  '''
  t_init = time.time()
  #print '......test.....', gbl._globals['compound']['r_disk']
  #########################################################################
  gbl.velocities = np.linspace(gbl._globals['compound']['vel'] - gbl._globals['compound']['d_vel'], \
                               gbl._globals['compound']['vel'] + gbl._globals['compound']['d_vel'], \
                               gbl._globals['compound']['nstep'])
  # array with sampling velocities
  
  ###########################################################################
  if (np.linalg.norm(gbl._globals['compound']['los_start'] \
      - gbl._globals['compound']['los_end']) \
      <= float(gbl._globals['compound']['dgrid'])/2.):
      sys.exit('line of sight to short for integration...exiting...')
  
  # interrupt program if start- and endpoint of the line of 
  # sight are closer together then dgrid/2
  if not 'losoffsets' in gbl._globals['compound']:
    gbl._globals['compound']['losoffsets'] = np.array([0, 0])   
  # set offsets to [0,0] if not defined otherwise
  find_transition_number_and_frequency()
  if gbl._globals['verbose']: print('numbers related to species and transition: ', gbl._globals['compound']['number'])
  if gbl._globals['verbose']: print('related frequency: ', gbl._globals['compound']['frequency'])
  print('Velocities initialised.')
  if timed: print('  in', time.time()-t_init, 's')
  #pause1 = input('Ist de selected species korrekt?')
  gbl._globals['compound']['ranges'] = {}
  gbl._globals['compound']['ranges']['npoints'] = range(gbl._globals['compound']['npoints'])
  gbl._globals['compound']['ranges']['nspe'] = range(gbl._globals['compound']['nspe'])
  # gbl._globals['compound']['ranges']['npoints'] = range(gbl._globals['compound']['nspec'])
  # gbl._globals['compound']['ranges']['npoints'] = range(gbl._globals['compound']['nspec'])
  # gbl._globals['compound']['ranges']['npoints'] = range(gbl._globals['compound']['nspec'])
  # gbl._globals['compound']['ranges']['npoints'] = range(gbl._globals['compound']['nspec'])
  # gbl._globals['compound']['ranges']['npoints'] = range(gbl._globals['compound']['nspec'])
  # gbl._globals['compound']['ranges']['npoints'] = range(gbl._globals['compound']['nspec'])
  return

def calculateFUV(timed=True):
  '''
  A function to separate this part of the code.
  Created on 17.10.2019 by Craig.
  Modified for optimisation on 18.10.2019 by Craig.
  '''
  t_init = time.time()
  npoints = range(gbl._globals['compound']['npoints'])
  # Positive velocities indicate velocities into the direction of the los...
  # Determine averaged FUV absorption for each pixel
  #first only clumps:
  if gbl._globals['Flags']['reuse_FUV']:
      # import fuv field strength from table
      i = 0
      with open(gbl.KOSMAPATH+'temp/FUV.dat') as fuv:
        for line in fuv:
          gbl._globals['compound']['ens'][i].FUV = float(line.split()[0])
          i = i + 1
  else:
    # recalculate incident FUV flux for each voxel 
    for i in npoints:
      gbl._globals['compound']['ens'][i].Mens = \
        gbl._globals['compound']['ens'][i].Mens_clumps
      gbl._globals['compound']['ens'][i].rho_ens = \
        gbl._globals['compound']['ens'][i].rho_ens_clumps
      gbl._globals['compound']['ens'][i].Mu = \
        gbl._globals['compound']['ens'][i].Mu_cl
      gbl._globals['compound']['ens'][i].Ml = \
        gbl._globals['compound']['ens'][i].Ml_cl
    gbl._globals['statistic']['Nmin'] = gbl._globals['statistic']['Nmin_clumps']  
    gbl._globals['constants']['alpha'] = gbl._globals['constants']['alpha_cl']
    print('\nClump FUV calculation.')
    ##
    fuva.fuvAbsorption()
    ##
    for i in npoints:
      gbl._globals['compound']['ens'][i].Afuv_tot = gbl._globals['compound']['ens'][i].Afuv
    if gbl._globals['verbose']: print('clump fuv ext. pixel 1: ', gbl._globals['compound']['ens'][1].Afuv_tot)
    # pause = input('calculation of pixel averaged fuv extinction by clumps finished')
    # SECOND: contibution of inter-clump medium
    if sum(gbl._globals['compound']['ens'][i].Mens_inter for i in \
          range(gbl._globals['compound']['npoints'])) == 0: 
      if gbl._globals['verbose']: print('no interclump medium found for fuv absorption...')
      # pause = input('ok?')
    else: 
      for i in npoints:
        gbl._globals['compound']['ens'][i].Mens = \
          copy(gbl._globals['compound']['ens'][i].Mens_inter)
        gbl._globals['compound']['ens'][i].rho_ens = \
          copy(gbl._globals['compound']['ens'][i].rho_ens_inter)
        gbl._globals['compound']['ens'][i].Mu = \
          copy(gbl._globals['compound']['ens'][i].Mu_inter)
        gbl._globals['compound']['ens'][i].Ml = \
          copy(gbl._globals['compound']['ens'][i].Ml_inter)
      gbl._globals['statistic']['Nmin'] = gbl._globals['statistic']['Nmin_inter']  
      gbl._globals['constants']['alpha'] = gbl._globals['constants']['alpha_intercl']
      print('\nInterclump FUV calculation.')
      ##
      fuva.fuvAbsorption()
      ##
      gbl._globals['constants']['alpha'] = None
      # print 'epsilon inter-clump', gbl._globals['compound']['ens'][1].epsilon
      if gbl._globals['verbose']: print('fuv ext. by interclump pixel 1: ', \
                        gbl._globals['compound']['ens'][0].Afuv)
      emissivity_tot = []
      for i in npoints:
        gbl._globals['compound']['ens'][i].Afuv_tot = \
          gbl._globals['compound']['ens'][i].Afuv_tot + \
          gbl._globals['compound']['ens'][i].Afuv
        gbl._globals['compound']['ens'][i].Afuv = None
      if gbl._globals['verbose']: print('total fuv ext. voxel 1: ', gbl._globals['compound']['ens'][0].Afuv_tot)
      # pause = input('calculation of pixel averaged fuv extinction by interclump medium finished')
    __import__(gbl._globals['Files']['fuv_filename'])
    # calculate fuv field strength for each pixel
    if gbl._globals['verbose']: print('FUV for first ensemble: ', gbl._globals['compound']['ens'][1].FUV)
    print('writing fuv field strength for each pixel to file ./temp/FUV.dat...')
    try: os.remove(gbl.KOSMAPATH+'temp/FUV.dat')
    except OSError: pass
    # remove old fuv file
    with open(gbl.KOSMAPATH+'temp/FUV.dat','w') as t:
      # write fuv field strength for each pixel to file. enables reuse        
      for i in npoints:
        t.write(str(gbl._globals['compound']['ens'][i].FUV) + ' ')
        t.write(str(gbl._globals['compound']['coordinates'][i]) + '\n')
  print('\nFUV calculated successfully.')
  if timed: print('  in', time.time()-t_init, 's')
  #print 'fuv pixel 0_m2_0:', gbl._globals['compound']['ens'][2222].FUV
  # pause = input('FUV done')
  return

###########################################################################
# column density weighted FUV map
# import columnDensityWeightedFUV

def lineAE(timed=True):
  '''
  A function to separate this part of the code.
  Created on 17.10.2019 by Craig.
  Modified (partially) for optimisation on 18.10.2019 by Craig.
  '''
  t_init = time.time()    #a varriable to time this section
  npoints = range(gbl._globals['compound']['npoints'])
  nspe = range(gbl._globals['compound']['nspe'])
  # calculate line opacities and emissivities
  if gbl._globals['Flags']['reuse_kappa_epsilon']: 
  # read optical depth and intensity in from file. for each pixel and each species
  # one array with tau(vobs) and inten(vobs)
    with open(gbl.KOSMAPATH+'temp/tau_line_tot.dat', 'r') as t:
      for pix in npoints:
        li = []
        for species in nspe:
          l = t.readline().split()[2 : len(gbl.velocities) + 2]
          li.append(l)
          # print 'li', li
        gbl._globals['compound']['ens'][pix].tau_tot =  np.array(li, float)
    with open(gbl.KOSMAPATH+'temp/inten_line_tot.dat', 'r') as t:
      for pix in npoints:
        li = []
        for species in nspe:
          l = t.readline().split()[2 : len(gbl.velocities) + 2]
          li.append(l)   
          print('species: ', species)
          print('l: ', l)          
        gbl._globals['compound']['ens'][pix].inten_tot =  np.array(li, float)
    if gbl._globals['verbose']: print("gbl._globals['compound']['ens'][0].tau_tot: ", gbl._globals['compound']['ens'][0].tau_tot)
    # pause = input('.check else')
  
  else:
    t = datetime.datetime.utcnow()  
    if gbl._globals['verbose']: print('Calculating opacities...: ',\
          (t-t_start).seconds,' seconds, ',\
          (t-t_start).microseconds, ' microseconds')
    #pause = input('pause opacities and emissivities')
    for i in npoints:
      gbl._globals['compound']['ens'][(i)].Mens = \
        copy(gbl._globals['compound']['ens'][i].Mens_clumps)
      gbl._globals['compound']['ens'][(i)].rho_ens = \
        copy(gbl._globals['compound']['ens'][i].rho_ens_clumps)
      gbl._globals['compound']['ens'][(i)].Mu = \
        copy(gbl._globals['compound']['ens'][i].Mu_cl)
      gbl._globals['compound']['ens'][(i)].Ml = \
        copy(gbl._globals['compound']['ens'][i].Ml_cl)
    gbl._globals['statistic']['Nmin'] = gbl._globals['statistic']['Nmin_clumps']
    gbl._globals['constants']['alpha'] = copy(gbl._globals['constants']['alpha_cl'])
    gbl._globals['sigma']['sigma_j'] = copy(gbl._globals['sigma']['sigma_cl_j'])
    gbl._globals['sigma']['sigma_ens_j'] = copy(gbl._globals['sigma']['sigma_ens_cl_j'])
    # Module line_absorption
    temp = time.time()
    print('\nClump intensity calculation.')
    ##
    lae.lineAbsorptionEmission()
    ##
    gbl._globals['runtimes']['module_line_absorption'] = time.time()-temp
    for i in npoints:
      gbl._globals['compound']['ens'][(i)].Mens_cl_sum = \
        copy(gbl._globals['compound']['ens'][(i)].Mens_sum)
      gbl._globals['compound']['ens'][(i)].inten_tot = \
        copy(gbl._globals['compound']['ens'][(i)].inten)
      gbl._globals['compound']['ens'][(i)].tau_tot = \
        copy(gbl._globals['compound']['ens'][(i)].tau)
      print('Clump: {}\nMass: {}\nFUV: {}\nIntensity: {}\nOptical Depth: {}'.format(i, gbl._globals['compound']['ens'][(i)].Mens_cl_sum, gbl._globals['compound']['ens'][i].FUV, (gbl._globals['compound']['ens'][(i)].inten), (gbl._globals['compound']['ens'][(i)].tau)))
    # Write to data file (for dense clumps)
    input()
    with open(gbl.KOSMAPATH+'temp/tau_line_clumps.dat', 'w') as t:
      # write tau_tot to file for dense clumps
      # textwrap ensures that the arrays are not wraped
      for pix in npoints:
        for species in nspe:          
          t.write(textwrap.fill(str(pix) + ' ' + str(species) + ' ' + str(gbl._globals['compound']['ens'][pix].tau[species,:])[1:-1], width = gbl._globals['compound']['nstep'] *20))
          t.write('\n') 
    with open(gbl.KOSMAPATH+'temp/inten_line_clumps.dat', 'w') as t:
      # write inten_tot to file for dense clumps
      # textwrap ensures that the arrays are not wraped
      for pix in npoints:
        for species in nspe:                   
          t.write(textwrap.fill(str(pix) + ' ' + str(species) + ' ' + str(gbl._globals['compound']['ens'][(pix)].inten[species,:])[1:-1], width = gbl._globals['compound']['nstep'] *20))
          t.write('\n') 
    if gbl._globals['verbose']: print('epsilon cl ' , gbl._globals['compound']['ens'][0].inten_tot)
    # pause = input('...')
    # Interclump
    if sum(gbl._globals['compound']['ens'][(i)].Mens_inter for i in npoints) == 0: 
      print('no interclump medium found for interclump ensemble {}...'.format(i))
      #pause = input('ok?')
    else:
      # case with interclump medium
      for i in npoints:
        gbl._globals['compound']['ens'][(i)].Mens = \
          copy(gbl._globals['compound']['ens'][(i)].Mens_inter)
        gbl._globals['compound']['ens'][(i)].rho_ens = \
          copy(gbl._globals['compound']['ens'][(i)].rho_ens_inter)
        gbl._globals['compound']['ens'][(i)].Mu = \
          copy(gbl._globals['compound']['ens'][(i)].Mu_inter)
        gbl._globals['compound']['ens'][(i)].Ml = \
          copy(gbl._globals['compound']['ens'][(i)].Ml_inter)
      # Change FUV for interclumb medium here
      gbl._globals['compound']['ens'][(i)].FUV = gbl._globals['compound']['interclump_FUV']
        #copy(gbl._globals['compound']['ens'][(i)].FUV)
      gbl._globals['statistic']['Nmin'] = gbl._globals['statistic']['Nmin_inter']  
      gbl._globals['constants']['alpha'] = gbl._globals['constants']['alpha_intercl']
      gbl._globals['sigma']['sigma_j'] = copy(gbl._globals['sigma']['sigma_inter_j'])
      gbl._globals['sigma']['sigma_ens_j'] = copy(gbl._globals['sigma']['sigma_ens_inter_j'])
      print('\nInterclump intensity calculation.')
       ##
      lae.lineAbsorptionEmission() 
      ##
      gbl._globals['constants']['alpha'] = None
      for i in npoints:
        print('Interclump: {}\nMass: {}\nFUV: {}\nIntensity: {}\nOptical Depth: {}'.format(i, gbl._globals['compound']['ens'][(i)].Mens_cl_sum, gbl._globals['compound']['ens'][i].FUV, (gbl._globals['compound']['ens'][(i)].inten), (gbl._globals['compound']['ens'][(i)].tau)))
        gbl._globals['compound']['ens'][(i)].Mens = None
        gbl._globals['compound']['ens'][(i)].rho_ens = None
        gbl._globals['compound']['ens'][(i)].Mu = None
        gbl._globals['compound']['ens'][(i)].Ml = None
      # Write to file for interclump medium
      with open(gbl.KOSMAPATH+'temp/tau_line_interclump.dat', 'w') as t:
        # write tau_tot to file for dense clumps
        # textwrap ensures that the arrays are not wraped
        for pix in npoints:
          for species in nspe:                   
            t.write(textwrap.fill(str(pix) + ' ' + str(species) + ' ' + str(gbl._globals['compound']['ens'][(pix)].tau[species,:])[1:-1], width = gbl._globals['compound']['nstep'] *20))
            t.write('\n')
      with open(gbl.KOSMAPATH+'temp/inten_line_interclump.dat', 'w') as t:
        # write inten_tot to file for dense clumps
        # textwrap ensures that the arrays are not wraped
        for pix in npoints:
          for species in nspe:                   
            t.write(textwrap.fill(str(pix) + ' ' + str(species) + ' ' + str(gbl._globals['compound']['ens'][pix].inten[species,:])[1:-1], width = gbl._globals['compound']['nstep'] *20))
            t.write('\n')   
      for i in npoints:
        gbl._globals['compound']['ens'][(i)].Mens_inter_sum = \
          copy(gbl._globals['compound']['ens'][(i)].Mens_sum)
        gbl._globals['compound']['ens'][(i)].inten_tot = \
          copy(gbl._globals['compound']['ens'][(i)].inten_tot + gbl._globals['compound']['ens'][(i)].inten)
        gbl._globals['compound']['ens'][(i)].tau_tot = \
          copy(gbl._globals['compound']['ens'][(i)].tau_tot + gbl._globals['compound']['ens'][(i)].tau)
      if gbl._globals['verbose']: print('tau.tot: ' , gbl._globals['compound']['ens'][0].tau_tot)
      # pause = input('...')
      for i in npoints:
        gbl._globals['compound']['ens'][(i)].epsilon = None
        gbl._globals['compound']['ens'][(i)].kappa = None
    # Write to data file
    with open(gbl.KOSMAPATH+'temp/tau_line_tot.dat', 'w') as t:
      # write tau_tot to file
      # textwrap ensures that the arrays are not wraped
      for pix in npoints:
        for species in nspe:           
          t.write(textwrap.fill(str(pix) + ' ' + str(species) + ' ' + str(gbl._globals['compound']['ens'][(pix)].tau_tot[species,:])[1:-1], width = gbl._globals['compound']['nstep'] * 20)) 
          t.write('\n')
    with open(gbl.KOSMAPATH+'temp/inten_line_tot.dat', 'w') as t:
      # write inten_tot to file
      # textwrap ensures that the arrays are not wraped
      for pix in npoints:
        for species in nspe:                   
          t.write(textwrap.fill(str(pix) + ' ' + str(species) + ' ' + str(gbl._globals['compound']['ens'][(pix)].inten_tot[species,:])[1:-1], width = gbl._globals['compound']['nstep'] *20))
          t.write('\n')     
        # pause = input('...ensemble averaged emissivity and optical depth have been written into folder ./temp/... hit enter to go on!')
  print('\nLine absorption and emission calculated successfully.')
  if timed: print('  in', time.time()-t_init, 's')
  if gbl._globals['verbose']:
    for ens in gbl._globals['compound']['ens']:
      print('clump', ens.Mens_clumps, ens.Mens_cl_sum, '-- interclump', ens.Mens_inter, ens.Mens_inter_sum, '-- last run', ens.Mens_sum)
  return



#inten_tot = []
#tau_tot = []
#for i in range(gbl._globals['compound']['npoints']):
#    inten_tot.append(gbl._globals['compound']['ens'][i].inten_tot)
#    tau_tot.append(gbl._globals['compound']['ens'][i].tau_tot)

#print 'inten tot 0', gbl._globals['compound']['ens'][0].inten_tot
#print 'tau tot 0', gbl._globals['compound']['ens'][0].tau_tot
# pause = input('...ok?????...')
  
def radTransfer(timed=True):
  '''
  A function to separate the writing of the data files.
  Created on 17.10.2019 by Craig.
  '''
  t_init = time.time()    #a variable to time this section
  # Initialise the spectral cube
  x_length = gbl._globals['compound']['mapSize']['x']
  y_length = gbl._globals['compound']['mapSize']['y']
  n_spe    = gbl._globals['compound']['nspe']
  n_step   = gbl._globals['compound']['nstep']
  x_offset = gbl._globals['compound']['offsets']['x']
  y_offset = gbl._globals['compound']['offsets']['y']
  scale = gbl._globals['compound']['pixelsize']
  spec_cube = np.zeros(  (n_step,y_length,x_length) , dtype=np.float32 )
  hdu = fits.PrimaryHDU(spec_cube) # sets data as hdu-list
  # Define some ranges
  nspe = range(n_spe)
  nstep = range(n_step)
  xlength = range(x_length)
  ylength = range(y_length)
  # Write a general header
  hdu.header['NAME'] = "spectral cube of galactic Disk"
  hdu.header['CRPIX1'] = 1   #ref Element/Pixel
  hdu.header['CTYPE1'] = 'x-dim in pc' #Name 
  hdu.header['CRVAL1'] = ((-x_length+1)/2.)*scale #Reference Value in kpc
  hdu.header['CDELT1'] = 1. *scale #stepsize per counting index in kpc
  #
  hdu.header['CRPIX2'] = 1
  hdu.header['CTYPE2'] = "y-dim in pc"
  hdu.header['CRVAL2'] = ((-y_length+1)/2.)*scale 
  hdu.header['CDELT2'] = 1. *scale
  #
  hdu.header['CRPIX3'] = 1
  hdu.header['CTYPE3'] = "spectra in km/s"
  hdu.header['CRVAL3'] = (-gbl._globals['compound']['d_vel'] + gbl._globals['compound']['vel'])
  hdu.header['CDELT3'] = 2.*gbl._globals['compound']['d_vel']/(gbl._globals['compound']['nstep']-1.)
  #
  hdu.header['BUNIT'] = 'K'
  hdu.header['BZERO'] = 0
  hdu.header['BSCALE'] = 1
  #
  gbl._globals['runtimes']['start_rad_transfer'] = time.time()
  # Filling cube with data
  for sp in nspe: #for all Spezies
    print(sp+1, ' species out of: ', gbl._globals['compound']['nspe'])
    for x in xlength:   #for all x-Pos
      #print x+1, 'position out of:', x_length
      for y in ylength:   # for all y-Pos
        print('Species, x, y:', sp, x, y)
        yval = rt.rad_transfer(offset = [x - x_offset,y - y_offset], species = sp) #rad transfer for real pos
        for i in nstep:
          spec_cube[((i))][y][x] = yval[(i)]
    run = (time.time()-gbl._globals['runtimes']['start_rad_transfer'])
    if gbl._globals['verbose']: print('running ', run, 's of total estimated ',  run/((sp+1.)/gbl._globals['compound']['nspe']), 's')
    spec_name = gbl._globals['compound']['species'][sp] + '_' + \
      str(gbl._globals['compound']['transition'][sp]) + '-' + \
      str(gbl._globals['compound']['transition'][sp] - 1)\
      + str(gbl._globals['namesuffix'])  #+'_complement' #for other half of disk
    if gbl._globals['verbose']:print('spectra_name: ', spec_name)
    hdu.header['SPECIES']= spec_name
    hdu.header['FREQUENZ'] = gbl._globals['compound']['frequency'][sp]
    hdu.header['FREQUNIT'] = 'GHz'
    hdu.header['LAMBDA'] = gbl._globals['constants']['c']/gbl._globals['compound']['frequency'][sp]/10**3
    hdu.header['LAMBUNIT'] = 'micron'
    hdu.header['KONFAC1'] = 2*gbl._globals['constants']['kb'] *(gbl._globals['compound']['frequency'][sp]*10**9)**3 /(gbl._globals['constants']['c'])**3*1000*10**9
    hdu.header['KONFAC1e'] = 'multiply for K -> nW/m^2/sr' #integrated not nW/m^2/sr/Hz
    # Clearing old file
    #open(gbl.KOSMAPATH+'fits_file/data/'+ spec_name + '.fits', 'w').close()
    # Writing cube to file
    hdu.writeto(gbl.KOSMAPATH+'fits_file/data/'+ spec_name + '.fits', overwrite=True) 
  gbl._globals['runtimes']['end_rad_transfer'] = time.time()
  # Write log with runtimes
  # Anzahl Massen
  Mu = gbl._globals['compound']['ens'][0].Mu_cl
  Ml = gbl._globals['compound']['ens'][0].Ml_cl
  if gbl._globals['verbose']: print('Mu ', Mu)
  if gbl._globals['verbose']: print('np.log10(Mu) ', np.log10(Mu))
  MuLog = np.log10(Mu) * 10  # 0.1 -> -10 while 10 stays 10
  MlLog = np.log10(Ml) * 10
  number_masspoints = len(np.arange(MlLog, MuLog + 10, 10))
  # parameter for runtimes
  start = gbl._globals['runtimes']['start']
  end_model = gbl._globals['runtimes']['end_model']
  end_programm = time.time()
  #print gbl._globals['runtimes']['line_absorption']
  # calculates time for each voxel, and total time
  line_total = 0
  for i in range(len(gbl._globals['runtimes']['line_absorption'])-1):
    delta = gbl._globals['runtimes']['line_absorption'][(i)+1] - gbl._globals['runtimes']['line_absorption'][(i)]
    #print 'delta', delta
    line_total = line_total + delta
  with open(gbl.KOSMAPATH+'log/'+str(datetime.datetime.utcnow()).replace(' ', '_')+'_times.dat','w') as save:
    print('writing runtime-log')
    save.write(str ( 'Voxels_with_mass: ' + str(gbl._globals['compound']['npoints']) + '\n' ) )
    save.write(str ( 'Number_v_bins: ' + str(gbl._globals['compound']['nstep']) + '\n' ) )
    save.write(str ( 'number_species: ' + str(gbl._globals['compound']['nspe']) + '\n' ) )
    save.write(str ( 'number_mass_classes: ' + str(number_masspoints) + '\n\n' ) )
    save.write(str ('Setup_time_model_[s]: ' + str(end_model-start) + '\n'  ) )
    save.write(str ('loop_time_line_emisions_[s]: ' + str(line_total) + '\n'  ) )
    save.write(str ('total_modul_time_line_emisions_[s]: ' + str(gbl._globals['runtimes']['module_line_absorption']) + '\n'  ) )
    save.write(str ('total_radtransfer_[s]: ' + str(gbl._globals['runtimes']['end_rad_transfer']-gbl._globals['runtimes']['start_rad_transfer']) + '\n'  ) )
    save.write(str ('rest_[s]: ' + str(end_programm-start - (gbl._globals['runtimes']['end_rad_transfer']-gbl._globals['runtimes']['start_rad_transfer']) \
                    - gbl._globals['runtimes']['module_line_absorption']  - (end_model-start) ) + '\n'  ) )
    save.write(str ('total_time_of_programm_[s]: ' + str(end_programm-start) + '\n'  ) )
  print('\nWriting of files completed successfully.')
  if timed: print('  in', time.time()-t_init, 's')
  return

def showPlots():
  '''
  A function to separate the plotting portion of the code. It uses the global
  plotting flags set via globals.plotFlags, with options:
    Total mass plot
    Clump mass plot
    Interclump mass plot
    Total Intensity plot
    Clump Intensity plot
    Interclump intensity plot
    Density plot
    Z velocity plot
    Total velocity plot
    Dispersion plot
    Total dispersion plot
    FUV distribution plot
  Created on 17.10.2019 by Craig.
  Modified for optimisation on 18.10.2019 by Craig.
  '''
  
  npoints = range(gbl._globals['compound']['npoints'])
  nspe = range(gbl._globals['compound']['nspe'])
  
  # Mass plot
  if gbl.plotFlags['Total mass']:
    Mplot = []
    for i in npoints:
      #Mplot.append(ens[i].Mens_clumps + ens[i].Mens_inter)
      Mplot.append(np.log10(gbl._globals['compound']['ens'][i].Mens_clumps + gbl._globals['compound']['ens'][i].Mens_inter))
      #print "total Voxel mass", ens[i].Mens_clumps , ens[i].Mens_inter
    #print Mplot
    fig = plot_geo(gbl._globals['compound']['abs_coordinates'], \
                      wei = Mplot, limits = gbl._globals['plots']['plotrange'], \
                      title = 'total mass distribution', cbarLabel = "log S_mass per voxel", \
                      labelx ='x [pc]', labely ='y [pc]', labelz='z [pc]')
    plt.show()
  
  """
  ##################################################### in SM/PC^2
  show_mass_plot = 1
  if show_mass_plot == 1:
    Mplot = []
    for i in range(gbl._globals['compound']['npoints']):  
      #Mplot.append(ens[i].Mens_clumps + ens[i].Mens_inter)
      Mplot.append(np.log10((ens[i].Mens_clumps + ens[i].Mens_inter)/gbl._globals['compound']['pixelsize']**2)    )
      print "total Voxel mass", ens[i].Mens_clumps/1000/1000. , ens[i].Mens_inter
    print max(Mplot)
    #pause = input('stoppp')
    fig2 = plot_geo(gbl._globals['compound']['abs_coordinates'], \
            wei = Mplot, limits = gbl._globals['plots']['plotrange'], \
            title = 'total mass distribution', cbarLabel = r"$log( \Sigma [M_\odot / pc^2])$", \
            labelx ='x [pc]', labely ='y [pc]', labelz='z [pc]')
    p.show()
  """
  
  # Clump mass plot
  if gbl.plotFlags['Clump mass']:
    Mplot = []
    for i in npoints:  
      Mplot.append(np.log10(gbl._globals['compound']['ens'][i].Mens_clumps))
    fig = plot_geo(gbl._globals['compound']['abs_coordinates'], \
                      wei = Mplot, limits = gbl._globals['plots']['plotrange'], \
                      title = 'inter mass distribution', cbarLabel = "log Smass of clumps per voxel", \
                      labelx ='x [kpc]', labely ='y [kpc]', labelz='z [kpc]')
    plt.show()
  
  # Interclump mass plot
  if gbl.plotFlags['Interclump mass']:
    Mplot = []
    for i in npoints:  
      Mplot.append(np.log10(gbl._globals['compound']['ens'][i].Mens_inter))
    fig = plot_geo(gbl._globals['compound']['abs_coordinates'], \
                      wei = Mplot, limits = gbl._globals['plots']['plotrange'], \
                      title = 'inter mass distribution', cbarLabel = "log S_mass intermedium per voxel", \
                      labelx ='x [kpc]', labely ='y [kpc]', labelz='z [kpc]')
    plt.show()
  
  # Density plot
  if gbl.plotFlags['Density']:
    Mplot = []
    for i in npoints:  
      Mplot.append(np.log10(gbl._globals['compound']['ens'][i].rho_ens_clumps + gbl._globals['compound']['ens'][i].rho_ens_inter))
    fig = plot_geo(gbl._globals['compound']['abs_coordinates'], \
                      wei = Mplot, limits = gbl._globals['plots']['plotrange'], \
                      title = 'density distribution', cbarLabel = "log density [/cm^3]", \
                      labelx ='x [kpc]', labely ='y [kpc]', labelz='z [kpc]')
    plt.show()
  
  # Z velocity plot
  if gbl.plotFlags['Z velocity']:
    Mplot = []
    for i in npoints:
      Mplot.append(gbl._globals['compound']['ens'][i].velocity[2])
    fig = plot_geo(gbl._globals['compound']['abs_coordinates'], \
                      wei = Mplot, limits = gbl._globals['plots']['plotrange'], \
                      title = 'v distribution', cbarLabel = "velocity of z-component [km/s]", \
                      labelx ='x [kpc]', labely ='y [kpc]', labelz='z [kpc]')
    plt.show()
  
  # Total velocity plot
  if gbl.plotFlags['Total velocity']:
    Mplot = []
    for i in npoints:
      v_total = (gbl._globals['compound']['ens'][i].velocity[0]**2 + gbl._globals['compound']['ens'][i].velocity[1]**2 + gbl._globals['compound']['ens'][i].velocity[2]**2 )**0.5
      Mplot.append(v_total)
    fig = plot_geo(gbl._globals['compound']['abs_coordinates'], \
                      wei = Mplot, limits = gbl._globals['plots']['plotrange'], \
                      title = 'v distribution', cbarLabel = "v [km/s]", \
                      labelx ='x [kpc]', labely ='y [kpc]', labelz='z [kpc]')
    plt.show()
  
  # Velocity dispersion plot
  if gbl.plotFlags['Dispersion']:
    Mplot = []
    for i in npoints:
      #Mplot.append(ens[i].v_dispersion)
      Mplot.append(np.log10(gbl._globals['compound']['ens'][i].v_dispersion))
    fig = plot_geo(gbl._globals['compound']['abs_coordinates'], \
                      wei = Mplot, limits = gbl._globals['plots']['plotrange'], \
                      title = 'systematic v dispersion', cbarLabel = "log systematic veloctiy-dispersion z-component [km/s]", \
                      labelx ='x [kpc]', labely ='y [kpc]', labelz='z [kpc]')
    plt.show()
  
  # Total velocity dispersion plot
  if gbl.plotFlags['Total dispersion']:
    Mplot = []
    for i in npoints:
      hilf = ( float(gbl._globals['sigma']['sigma_ens_cl_j'][0])**2 + gbl._globals['compound']['ens'][i].v_dispersion**2 )**0.5
      Mplot.append(np.log10(hilf))
    fig = plot_geo(gbl._globals['compound']['abs_coordinates'], \
                      wei = Mplot, limits = gbl._globals['plots']['plotrange'], \
                      title = 'systematic v dispersion', cbarLabel = "log total veloctiy-dispersion z-component [km/s]", \
                      labelx ='x [kpc]', labely ='y [kpc]', labelz='z [kpc]')
    plt.show()
  
  # Clump intensity plot
  if gbl.plotFlags['Clump intensity']:
    clumps_inten_peak_plot = [None] * gbl._globals['compound']['nspe']
    for spe in npoints:
      clumps_inten_peak_plot[(spe)] = []
    with open(gbl.KOSMAPATH+'temp/inten_line_clumps.dat', 'r') as t:
      for pix in npoints:
        li = []
        for species in nspe:
          l = t.readline().split()[2 : len(gbl.velocities) + 2]
          clumps_inten_peak_plot[species].append(l[velocity_plot]) 
    if gbl._globals['verbose']:
      print('clumps_inten_peak_plot[0]: ', np.array(clumps_inten_peak_plot[0],float))
      pause = input('...')
    # plot clump intensity(velocity_plot) for each species
    for spe in nspe:
       #label = gbl._globals['compound']['species'][spe] + \
       #      str(gbl._globals['compound']['transition'][spe]) + \
       #      ' line intensity; clump' + ' v=' + str(velocities[velocity_plot])
      label =  gbl._globals['compound']['species'][spe] + '_' + \
                str(gbl._globals['compound']['transition'][spe]) + \
                ' line centre intensity [K] (only clumps)'# + ' v=' + str(velocities[velocity_plot])
      fig = plot_geo(gbl._globals['compound']['coordinates'], \
                        wei = np.array(clumps_inten_peak_plot[(spe)], float), \
                        species = spe ,\
                        limits = gbl._globals['plots']['plotrange'],\
                        title =  'title' , cbarLabel= label)
    plt.show()
  
  # Interclump intensity plot
  if gbl.plotFlags['Interclump intensity']:
    interclump_inten_plot = [None] * gbl._globals['compound']['nspe']
    for spe in nspe:
      interclump_inten_plot[(spe)] = []
    with open(gbl.KOSMAPATH+'temp/inten_line_interclump.dat', 'r') as t:
      for pix in npoints:
        li = []
        for species in nspe:
          l = t.readline().split()[2 : len(gbl.velocities) + 2]
          interclump_inten_plot[species].append(l[velocity_plot]) 
    if gbl._globals['verbose']:
      print('interclump_inten_plot[0]: ', np.array(interclump_inten_plot[0],float))
      pause = input('...')
    # plot interclump intensity(velocity_plot) for each species
    for spe in nspe:
      #label = gbl._globals['compound']['species'][spe] + \
      #       str(gbl._globals['compound']['transition'][spe]) + \
      #       ' line intensity; interclump' + ' v=' + str(velocities[velocity_plot]) 
      label = gbl._globals['compound']['species'][spe] + '_' + \
              str(gbl._globals['compound']['transition'][spe]) + \
              ' line centre intensity [K] (only interclump)'# + ' v=' + str(velocities[velocity_plot])
      fig = plot_geo(gbl._globals['compound']['coordinates'], \
                        wei = np.array(interclump_inten_plot[(spe)], float), \
                        species = spe, \
                        limits = gbl._globals['plots']['plotrange'],\
                        title =  'title' , cbarLabel= label)
    plt.show()
  
  # Total intensity plot
  if gbl.plotFlags['Total intensity']:
    inten_tot_peak_plot = [None] * gbl._globals['compound']['nspe']
    for spe in nspe:
      inten_tot_peak_plot[(spe)] = []
    with open(gbl.KOSMAPATH+'temp/inten_line_tot.dat', 'r') as t:
      for pix in npoints:
        li = []
        for species in nspe:
          l = t.readline().split()[2 : len(gbl.velocities) + 2]
          inten_tot_peak_plot[species].append(l[velocity_plot]) 
    if gbl._globals['verbose']:
      print('inten_tot_peak_plot[0]: ', np.array(inten_tot_peak_plot[0],float))
      pause = input('...')
    # plot total intensity(velocity_plot) for each species
    for spe in nspe:
      label = gbl._globals['compound']['species'][spe] + '_' + \
              str(gbl._globals['compound']['transition'][spe]) + \
              ' line centre intensity (K)'# + ' v=' + str(velocities[velocity_plot])
      fig = plot_geo(gbl._globals['compound']['coordinates'], \
                        wei = np.array(inten_tot_peak_plot[(spe)], float), \
                        species = spe, \
                        limits = gbl._globals['plots']['plotrange'],\
                        title =  'title' , cbarLabel= label)
    plt.show()
  
  # Plot FUV distribution
  if gbl.plotFlags['FUV distribution']:
    fuvTab = []
    with open(gbl.KOSMAPATH+'temp/FUV.dat') as fuv:
      for line in fuv:
        fuvTab.append(float(line.split()[0]))
    fuvTabLog = np.log10(np.array(fuvTab)) * 10
    # print 'fuvTab', fuvTab
    fig = plot_geo(gbl._globals['compound']['abs_coordinates'], \
                      wei = fuvTabLog / 10., limits = gbl._globals['plots']['plotrange'], \
                      title = 'FUV distribution', cbarLabel="Log$_{10}$[FUV flux (Draine)]", \
                      labelx ='x [kpc]', labely ='y [kpc]', labelz='z [kpc]')
    plt.show()
  
  return

############################################################################
class easy_choice:
  # Tkinter GUI for PDR3D. # tells the program what to do:
  # 1 - plot the compound
  # 2 - start the GUI to display spectra and map
  # 3 - more advanced spectra plotting (in this case the offsets 
    # of different los need to be defined when setting up the compund)
  # 4 - write maps to file. possible for an array of different species and 
  # transitions as input
  def __init__(self, root):
    self.root = root
    root.title('3D PDR Simulation: please choose!!')
    button_GUI = Button(root, text="Start Interface", command = self.GUI)
    button_GUI.config(bd=5)
    button_GUI.pack()  
    button_maps_to_file = Button(root, text="write maps to file (possible for many species)", command = self.maps_to_file)
    button_maps_to_file.pack()
    button_compound = Button(root, text="Show Geometry", command = self.compound)
    button_compound.pack()
    button_specs = Button(root, text= "Map Spectra (one species, spectra at different positions, offsets needed!!)", command = self.spec)
    button_specs.pack()
    button_specs_species = Button(root, text= "Spectra at one position, for many species", command = self.spec_species)
    button_specs_species.pack()
  def compound(self): 
    #from PDR import _globals
    gbl._globals['plot'] = 1
    self.root.destroy()
  def spec(self): 
    #from PDR import _globals
    gbl._globals['plot'] = 2
    self.root.destroy() 
  def GUI(self):
    #from PDR import _globals
    gbl._globals['plot'] = 3
    self.root.destroy()
  def maps_to_file(self):
    #from PDR import _globals
    gbl._globals['plot'] = 4
    self.root.destroy()
  def spec_species(self):
    #from PDR import _globals
    gbl._globals['plot'] = 5
    self.root.destroy()
##################################################################
# function which writes ensemble parameters 'component' 
# to file ./data/parameters.dat
# prviois parameters will be overwritten
def writeparameters(component):
  #from PDR import _globals
  with open(gbl._globals['directories']['data'] + 'parameters.dat','w') as par:
    par.write('pixel_edge_length_in_pc ' + str(gbl._globals['compound']['dint']/gbl._globals['constants']['pc']) + '\n')
    par.write('species ' + str(gbl._globals['compound']['species']) + '\n')
    par.write('transition ' + str(gbl._globals['compound']['transition'])  + '\n')
    for i in range(gbl._globals['compound']['npoints']):
      par.write('No._of_ensemble:' + str(i) + ' ')
      par.write(str(component[i].metal) + ' ')
      par.write(str(component[i].density) + ' ')
      par.write(str(component[i].Mens) + ' ')
      par.write(str(component[i].FUV) + ' ')
      par.write(str(component[i].mlower) + ' ')
      par.write(str(component[i].mupper) + '\n')

#######################################################################

# P L O T T I N G   F U N C T I O N S

# plots array of Koordinates with matplotlib. different weights for different
# points can be given. weight of coordinates will be proportional to surface
# of plotted points
# 
#
# Input:
#     coord: coordinates, np.array with [x,y,z] values
#     losoffset: offset for n  different lines of sight.
#                offset x and y positions given as n*2 array
#
# from testing: the size of the points changes with s: surface proportinal to s
#                and radius like sqrt(s)
#
# def types of input, errors
# marker ='o', '^','d':shape of points does not work??

def plot_geo(coord, wei = None, labelx='x',labely='y',labelz='z',title='Title', limits = [[0,1],[0,1],[0,1]], los_start = [0,0,0], los_end = [0,0,0], losoffset = None, cbarLabel = 'cbar_Label'):
  # scale size to fill pixel (?)
  # print 'weight of ensembles propotional to surface of points in plots'
  if losoffset != None and losoffset.size == 2: 
  # case only one los given
      losoffset = losoffset.flatten()
  coord = np.array(coord)
  los_start = np.array(los_start)
  los_end = np.array(los_end)
  xvalues = coord.transpose()[0]
  yvalues = coord.transpose()[1]
  zvalues = coord.transpose()[2]
  # npoints = xvalues.size # number of points to be ploted
  # si = 470000
  # if wei != None:
#      print 'wei!', wei
  #     wei = np.array(wei, float)
  #     wei = wei/sum(wei) * float(si)
#        wei = wei**2
  # else: wei = (np.ones(npoints)*float(si)/npoints)    
  # wei = (np.ones(npoints)*float(si)/npoints)  
  fig = plt.figure(figsize=(4,4))
  ax = p3.Axes3D(fig)
  # ax.scatter(xvalues, yvalues, zvalues, marker='s')
  #ax.set_scale('log')
  #scat = ax.scatter(xvalues, yvalues, zvalues, cmap=plt.cm.cubehelix, c=wei, vmin=min(wei), vmax = max(wei), marker='s')
  scat = ax.scatter(xvalues, yvalues, zvalues, cmap=plt.cm.jet, c=wei, vmin=min(wei), vmax = max(wei), marker="s", s=20, alpha=None, linewidths=None, edgecolors='none' )
  cb = fig.colorbar(scat, shrink=0.75)
  cb.set_label(cbarLabel)
  # colorscale???
  # ax.scatter(xvalues, yvalues, zvalues, s=wei, c='r', marker='s')
  # scale site with weight
  ax.set_xlabel(labelx) # labels cause error? 
  ax.set_ylabel(labely)
  ax.set_zlabel(labelz)
#    ax.legend()
#    fig.suptitle(title, fontsize=16)# if you need a title ...
#   if startint != None:
#       if endint != None:
  #     a = min(limits[0][1],limits[1][1],limits[2][1])
  #    alpha = np.linspace(0,a,4)
    #   print a
    #  print alpha
#       sx = r0[0] + alpha * lineos[0]
#      sy = r0[1] + alpha * lineos[1]
  #     sz = r0[2] + alpha * lineos[2]
  if losoffset == None:
    sx = np.array([los_start[0], los_end[0]])
    sy = np.array([los_start[1], los_end[1]])
    sz = np.array([los_start[2], los_end[2]])   
    ax.plot(sx, sy, sz, zdir = 'z')
  elif losoffset.size == 2:
    sx = np.array([los_start[0] + losoffset[0], los_end[0] + losoffset[0]])
    sy = np.array([los_start[1] + losoffset[1], los_end[1] + losoffset[1]])
    sz = np.array([los_start[2], los_end[2]])   
    ax.plot(sx, sy, sz, zdir='z')    
  else:
    l = losoffset.size/2
#        sx = np.zeros(l)
#        sy = np.zeros(l)
#        sz = np.zeros(l)        
    for k in range(l):
      sx = np.array([los_start[0] + losoffset[k][0], los_end[0] + losoffset[k][0]])
      sy = np.array([los_start[1] + losoffset[k][1], los_end[1] + losoffset[k][1]])
      sz = np.array([los_start[2], los_end[2]])               
      ax.plot(sx, sy, sz, zdir='z')  
  # define plotrange    
  # ax.set_xlim3d(limits[0][0],limits[0][1])
  ax.set_ylim3d(limits[1][0],limits[1][1])
  ax.set_zlim3d(limits[2][0],limits[2][1])
  for a in ax.w_xaxis.get_ticklines()+ax.w_xaxis.get_ticklabels():
    a.set_visible(False) 
  plt.rcParams.update({'font.size': 30}) ##22
  fig.add_axes(ax)
  plt.axis('off')
#  p.show()
  return fig

def plot_geo2D(coord, wei = None, species = -1, labely='Y',labelz='Z',title='Title', \
               limits = [[0,1],[0,1]], los_start = [0,0,0], los_end = [0,0,0], \
               losoffset = None, cbarLabel = 'cbar_Label'):
  # scale size to fill pixel (?)
  # print 'weight of ensembles propotional to surface of points in plots'
  from mpl_toolkits.axes_grid1 import make_axes_locatable
  coord = np.array(coord)
  # xvalues = coord.transpose()[0]
  yvalues = coord.transpose()[1]
  zvalues = coord.transpose()[2]
  fig, ax = pyplot.subplots(figsize=(6, 6))
  #fig = Figure(figsize=(4,4))
  #ax = fig.add_subplot(111)
  im= ax.scatter(yvalues, zvalues, cmap=pyplot.cm.cubehelix, c=wei, vmin=min(wei), vmax = max(wei), marker='s');
  divider = make_axes_locatable(ax)
  cax = divider.append_axes("right", size="5%", pad=0.05)
  cb = fig.colorbar(im, cax=cax)
  cb.set_label(cbarLabel)
  ax.set_xlabel(labely) 
  ax.set_ylabel(labelz)
  fig.add_axes(ax)
  pyplot.savefig(gbl.KOSMAPATH+'pics/intensity_'+str(species)+'.eps', format='eps', dpi=600, bbox_inches='tight') 
#  p.show()
  return fig

def plot_spec(ISumarray, losoffset = np.array([0,0]), xr = np.array([-10, 10]), nstep = 50):
  # Print statements
  if gbl._globals['verbose']: print('losoffset: ' , losoffset)
  if gbl._globals['verbose']: print('velocity range: ', xr)
  losoffset = np.array(losoffset)
  elem = losoffset.size/2
  color = ('b', 'g', 'r', 'c', 'm', 'y', 'k')
  co =  -1 
  # counter for loops...
  # array with x-values:
  xval = np.linspace(xr[0], xr[1], nstep)
  # initialise array to contain y-values for different lines of sight
  if gbl._globals['verbose']: print('velocities: ', xval)
  yval = np.zeros([elem, nstep])
  # print "yval.size", yval.size
  # pause = input('..plt spec...')
  n_losoffset = range(int(losoffset.size/2))
  n_nstep = range(nstep)
  n_nxs = range(nxs)
  n_nys = range(nys)
  n_elem = range(elem)
  for k in n_losoffset:
    for n in n_nstep:
      if gbl._globals['verbose']: print('-----velocity-----', xval[n], ' km/s ------')
      yval[k][n] = ISumarray[k](xval[n])
      # yval[k][n] = ISumarray[0](12)
      if gbl._globals['verbose']: print('yval: ', yval[k][n])    
  # pause = input('...plot spec ..')      
  if losoffset.size == 2: 
  # case x and y dimension = 1 i.e only one spectrum
    f, ax = plt.subplots()
    ax.plot(xval, yval[0])
  #      ax.set_title('Offset [0,0]')
  else:    
      nx = copy(losoffset[:,0])
      ny = copy(losoffset[:,1])
      nx.sort()
      ny.sort()
      nx = np.unique(nx)
      ny = np.unique(ny)
      nxs = nx.size
      nys = ny.size
      if nxs == 1: # case y-dimension = 1
        f, axarr = plt.subplots(nys, sharex = True)
        for y in n_nys:
          for n in n_elem:
            if losoffset[n,1] == ny[y]: nn = n
          axarr[y].plot(xval, yval[nn])
          axarr[y].set_title('y-Offset'+str([ny[y]]))     
      elif nys == 1: # case y-dimension = 1
        f, axarr = plt.subplots(nxs, sharex = True)
        for x in range(nxs):
          co = co + 1
          for n in range(elem):
            if losoffset[n,0] == nx[x]: nn = n
          axarr[x].plot(xval, yval[nn], color = color[int(co)])
          axarr[x].set_title('x-Offset' + str([nx[x]]))              
      else: # case x and y direction more then 1-dimensional   
        f, axarr = plt.subplots(nxs, nys)
        for x in range(nxs):
          for y in range(nys):
            nn = -1
            for n in range(elem):
              if losoffset[n,0] == nx[x] and \
                losoffset[n,1] == ny[y]:
                nn = n
            if nn == -1: pass
            else:
              axarr[x,y].plot(xval, yval[nn])
              axarr[x,y].set_title('Offset' + str([nx[x], ny[y]]))
            if x == nxs - 1: axarr[x,y].set_xlabel('Velocity')
            if y == 0:
              axarr[x,y].set_ylabel('Intensity')    
          # Fine-tune figure; hide x ticks for top plots and y ticks for right plots
    #      for x in range(nxs-1):
    #          plt.setp([a.get_xticklabels() for a in axarr[x,:]], visible=False)
      #     for y in range(nys-1)+1:    
      #        plt.setp([a.get_yticklabels() for a in axarr[:,y]], visible=False)
  plt.show()

# D I S T R I B U T I O N S
from copy import copy
from operator import mul # multiplication
from functools import reduce
class binomial():
  '''class calculation binomial coefficients (choose) and function'''
  def __init__(self, n, p):
    self.n = n
    self.p = p
    return
  
  def comb(self, k):
    ''' calculate nCr - the binomial coefficient
    >>> comb(3,2)
    3
    >>> comb(9,4)
    126
    >>> comb(9,6)
    84
    >>> comb(20,14)
    38760
    '''
    if k > self.n-k:  # for smaller intermediate values 
                      # use (n choose k) = (n choose n-k)
      k = self.n-k
    return int(reduce( mul, range(int(self.n-k+1), int(self.n+1)), 1) /
               reduce( mul, range(1,int(k+1)), 1) )
  '''
  def choose(self, k):
      """
      A fast way to calculate binomial coefficients by 
      Andrew Dalke (contrib). (altered from original version).
      """
      print 'n', self.n
      print 'k', k
      nn = copy(self.n)
      if 0 <= k <= nn:
          ntok = 1
          ktok = 1
          for t in xrange(1, min(k, nn - k) + 1): # for runtime?
              ntok *= nn
              ktok *= t
              nn -= 1  # n=n-1
          print 'ntok', ntok
          print 'ktok', ktok
          return ntok // ktok
      else:
          return 0
  '''
  # I find that comb is more robust against 
  # large numbers compared to choose
  def binomfunc(self, k):
    # print 'comb', self.comb(k)
    #print (float(self.comb(k)) * self.p**k * (1-self.p)**(self.n-k))
    return float(self.comb(k)) * self.p**k * (1-self.p)**(self.n-k)

class gauss:
  def __init__(self, v0, sigma, area = 1):
    self.area = float(area) # area below the curve (integrated curve) 
    self.v0 = float(v0) # peak velocity
    self.sigma = float(sigma) # standard derivation
    return
  
  def gaussfunc(self, v):
    if self.sigma == 0:
      return 0
    else:
      import numpy as np   
      return self.area/(np.sqrt(2*np.pi)*self.sigma)*\
             np.exp(-.5 * ( (v - self.v0) / self.sigma)**2 )

#################################################################
class AllCombi():
  '''Input: array "lowerupper" of the form [[min_1, max_1], [min_2, max_2], ...] 
  with min_i and max_i being integers (max_i => min_i)
  array needs to have the form ([[1,2],[3,4],...]) even for one element:
  ([[1],[3],...])
  Function AllCombi calculates all possible combinations with one number from each interval,
  i.e. starting with [min_1, min_2, ...] and finishing with [max_1, max_2,...]
  '''
  def __init__(self, lowerupper):
    self.lowerupper = lowerupper
    self.list = [] # list of number of clumps in each interval
    self.numbercombi = 1 # for total number of combinations
    return

  def createArray(self):
    for line in self.lowerupper:
      l = np.arange(line[0], line[1] + 1)
      self.list.append(l)
      self.numbercombi = self.numbercombi * len(l) # total number of possible combinations
    self.numberintervals = len(self.list)
    self.combinations = np.zeros([self.numbercombi, self.numberintervals]) 
    # reserve memory for complete array
    #print 'numbercombi', self.numbercombi
    #print 'combinations', self.combinations
    for j in range(self.numbercombi): # each combination
      for i in range(self.numberintervals): # each column 
        self.combinations[j][i] = int(j/(self.product(i))) % self.list[i].size + self.list[i][0] # + self.list[i][0] for list wich do not start with 0
        #print 'j,i', j,i, 'self.product', self.product(i), 'list[i].size', self.list[i].size
    return
  
  def product(self, a):
    prod = 1
    while a < (self.numberintervals - 1):
      prod = prod * self.list[a+1].size
      #print 'self list', self.list[a+1]
      #print 'prod', prod
      a = a + 1
    return prod
  
  '''    
  def writeArray(self): 
      with open("all_combinations_static.dat","w") as a:
          for l in range(self.combinations.shape[0]):
              for m in range(self.combinations[l].size):
                  a.write(str(self.combinations[l][m]) + " ")
              a.write("\n")
  '''


# G E O M E T R Y   M A N I P U L A T I O N
def rot3D(grid, angle, axis):
  '''
  rotates array of 3D vectors around fixed axis by angle phi
  axis: rotation axis, 3D vector
  angle in rad, 0<angle<2pi
  grid: grid to be rotated, has to be array of 3D vectors

  comments...:
  idl program also works for vectors lager than 3d, but not necessary?
  necessary for array of angles? dont think so...
  '''
  length = 0
  axis = np.array(axis, float) #for higher precision
  grid = np.array(grid, float) #for higher precision
  angle = float(angle)
  if grid.shape[0] == 0:
    'grid is empty. Returning None'
    return
  if grid.shape[0] == 1:
    grid = grid.flatten()
  if gbl._globals['verbose']: print('grid ', grid)
  grid = np.array(grid)/float(1)  # for higher precision
  grid_rot = copy(grid) # array to store result
  leng = int( grid.size/3 )
  if leng == 1:
    grid_rot = np.array(rotate(grid, angle, axis))
    # rotate ensemble coordinates     
  else: 
    for i in range(leng):
      grid_rot[i] = np.array(rotate(grid[i], angle, axis))
  return grid_rot

#################################################################
#################################################################
# calculates the signum of an numpy.array or an scalar. 
# with sign[0]=sign[-0]=1

def signum(num):
  num = np.array(num)
  if num.size == 0:
    print('error in signum: no value given')
    return nan
  elif num.size == 1:
    if float(num) < 0: return -1
    elif float(num) == 0: return 0
    else: return 1    
  else:
    resultarray = numpy.zeros(narray.size)
    for i in numpy.arange(narray.size):
      if narray[i] == 0:
        resultarray[i] = 1
      else: resultarray[i] = narray[i]/abs(narray[i])
    return resultarray

#################################################################
#################################################################
# array_unique: find unique rows within 2D array
# WORKS ONLY for arrays having the shape (2,..)

def array_unique(ar):
  ar = np.array(ar)
  t = []
  for row in ar:
    t.append((row[0], row[1])) #
  print(t)
  t2 = list(set(t))
  print(t2)
  return np.array(t2)


##########################################################################
######################### poisson distribution #########################
import numpy as np
from math import factorial
class poisson():
  def __init__(self, la):
    self.la = float(la) # lambda = n * p
    return
  
  def poissonfunc(self, k):
    # print 'comb', self.comb(k)
    return self.la**k /factorial(k) * np.e**(-self.la)
