# read in tabulated functions Ereal and Eimag
#
#
#
import numpy as np
import math
from scipy.interpolate import griddata
#from PDR import _globals
# edit: Craig, 14.10.2019
import globals as gbl
#end edit

####################### for testing #######################
# pause = input('Etilde: test modus, deleting gbl._globals')
# gbl._globals = {}
# gbl._globals['Etilde'] = {}
# gbl._globals['Etilde']['switch'] = True
###########################################################

if gbl._globals['Etilde']['switch']:
  # Tables Eimag and Ereal havent been read in yet
  #data location stated in globals.py

  ErealTab = []
  EimagTab = []
  # lists for function values (to be read in)
  ErealX = []
  EimagX = []
  # list for x-values for Ereal and Eimag

  with open(gbl.INPUTPATH+'Ereal.dat','r') as e:
    for line in e:
      l = line.split()
      ErealX.append(float(l[0]))
      ErealTab.append(float(l[1]))

  with open(gbl.INPUTPATH+'Eimag.dat','r') as e:
    for line in e:
      l = line.split()
      EimagX.append(float(l[0]))
      EimagTab.append(float(l[1]))

  ErealTab = np.array(ErealTab)
  EimagTab = np.array(EimagTab)

  gbl._globals['Etilde']['ErealTab'] = ErealTab
  gbl._globals['Etilde']['EimagTab'] = EimagTab
  gbl._globals['Etilde']['ErealX'] = ErealX
  gbl._globals['Etilde']['EimagX'] = EimagX
  # store arrays to avoid multiple read-ins

  gbl._globals['Etilde']['switch'] = False

else:
  # Tabels have already been read in
  ErealTab = gbl._globals['Etilde']['ErealTab']
  EimagTab = gbl._globals['Etilde']['EimagTab']
  ErealX = gbl._globals['Etilde']['ErealX']
  EimagX = gbl._globals['Etilde']['EimagX']


###########################################################################
########################define functions ###########################################

def Ereal(x):
  # print 'Ereal: x', x
  if x.imag == 0: x = x.real
  # x should be a real number. remove imaginary party '0j' which
  # prevents ordering
  if x < 0.01:
    return 2*x/np.sqrt(np.pi)
  elif x > 8.0:
    return 1/(np.sqrt(np.pi) * x)
  else:
    return griddata(ErealX, ErealTab, np.array([x]) ,  method = 'linear')[0]

def Eimag(x):
  if x == abs(x)*1j:
    # maser case. treated in linear approximation
    x = abs(x)
    return 1 + 2*x/np.sqrt(np.pi)
  else:
    x = abs(x)
    # x needs to be real and positive, i.e. abs(a) or abs(b)
    if x < 0.01:
      return 1 - 2*x/np.sqrt(np.pi)
    elif x > 8.0:
      return 1/(np.sqrt(np.pi) * x)
    else:
      return griddata(EimagX, EimagTab, np.array([x]), method = 'linear')[0]
