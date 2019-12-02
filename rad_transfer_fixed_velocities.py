###
# calculates radiative transport along the line of sight
# Input: offset of line of sight respective to initial line of sight
# Output: function ISum: integrated intensity, velocity dependent
# Needs modlues and classes: makeList, array_equal, gauss
# Needs different parameters vom PDR3D

############################################################
import numpy as np
from copy import copy
import cmath
# edit: Craig, 14.10.2019
import globals as gbl
#end edit

def rad_transfer(offset = np.array([0,0]), species = 0):
  '''
  Adapted from the Bruckmann version.
  Modified for optimisation on 18.10.2019 by Craig.
  '''
  #from PDR import gbl._globals
  
  ds = gbl._globals['compound']['dint']

  velocities = np.linspace(gbl._globals['compound']['vel'] - gbl._globals['compound']['d_vel'], \
                           gbl._globals['compound']['vel'] + gbl._globals['compound']['d_vel'], \
                           gbl._globals['compound']['nstep'])
  # sample velocities

  # from PDRfunctions import gauss # , signum

  # round los end- and startpoint to numerical threshold
  los_start = np.around(np.array(gbl._globals['compound']['los_start']), \
              gbl._globals['threshold']['epsgrid'])
  los_end = np.around(np.array(gbl._globals['compound']['los_end']), \
            gbl._globals['threshold']['epsgrid'])

  # find z positions to start/end integration
  if los_start[2] >= 0: startz = int(np.ceil(los_start[2]))
  else: startz = int(np.floor(los_start[2]))
  if los_end[2] >= 0: endz = int(np.ceil(los_end[2]))
  else: endz = int(np.floor(los_end[2]))

  #print 'z coordinate to start integration:', startz
  #print 'z coordinate where integration ends:', endz

  steps = int(abs(endz - startz))
  # number of steps between voxels of los
  #print 'number of integration steps on line of sight:', steps
  
  # positions for integration steps on z axis
  if startz < endz:
    position = np.arange(steps + 1) + startz
  elif endz < startz:
    position = np.arange(steps + 1)[::-1] + endz    
  #print 'positions for integrating steps on z axis: ', position
  #print 'number of positions:', position.size

  # pause = input('..ok?..')

  # find numbers of pixels which lie on line of sight. Put these numbers 
  # into position_newpos (-1 if not existent)
  position_newpos = np.zeros(steps + 1, int)
  position_newpos[:] = -1
  # Define some ranges
  npoints = range(gbl._globals['compound']['npoints'])
  n_position = range(position.size)
  n_newpos = range(position_newpos.size)
  n_velocities = range(len(velocities))
  n_steps = range(steps)
  for i in npoints:
    if gbl._globals['compound']['coordinates'][i][0] == offset[0] and \
      gbl._globals['compound']['coordinates'][int(i)][1] == offset[1]: 
      # if x and y position lie on line of sight          
      for j in n_position:
        if gbl._globals['compound']['coordinates'][int(i)][2] == position[int(j)]:
        # if z position coincides with integration position[j]
        # ensemble number i goes to position_newpos[j]
        # (filter out esembles relevant for los integration
        # and put into correct order)
          position_newpos[int(j)] = int(i)
  # print 'position_newpos', position_newpos
  # print 'position_newpos.size', position_newpos.size
  # print 'steps', steps

  # read in tau(vobs) and intensity(vobs) for all positions in position_newpos
  # convert tau to absorption coefficient kappa and intensity to 
  # volume emission coefficient epsilon CHECK
  epsilon_vobs = [None] * position_newpos.size
  kappa_vobs = [None] * position_newpos.size
  
  for j in n_newpos:
    if position_newpos[j] == -1:
      epsilon_vobs[j] = np.zeros(len(velocities))
      kappa_vobs[j] = np.zeros(len(velocities))
    else:
      epsilon_vobs[j] = gbl._globals['compound']['ens'][position_newpos[j]].inten_tot[species] / float(ds)  
      kappa_vobs[j] = gbl._globals['compound']['ens'][position_newpos[j]].tau_tot[species] / float(ds)  

  # print epsilon_vobs[0][5]
  # print epsilon_vobs[21][5]

  # pause = input('..ok2?..')

  epsilon_vobs_step = [None] * (position_newpos.size - 1)
  kappa_vobs_step = [None] * (position_newpos.size - 1)
  for j in n_newpos[:-1]:
    #print 'j-wert', j    #17.9
    #print 'episi', epsilon_vobs[j]
    #print 'episi+1', epsilon_vobs[j+1] #17.9
    # print 'epsilon_vobs[j]\n', epsilon_vobs[j]
    # print 'epsilon_vobs[j+1]\n', epsilon_vobs[j+1]
    epsilon_vobs_step[j] = (epsilon_vobs[j + 1] - epsilon_vobs[j])/float(ds)
    kappa_vobs_step[j] = (kappa_vobs[j + 1] - kappa_vobs[j])/float(ds)
      

  """
  # create list of instances (absorption and emission) from class gauss
  # , for every position needed for los integration.
  # *_v indicates velocity dependence
  epsilon_v = []
  kappa_v = []

  # define functions of emissivity and absorption, velocity dependent:
  # therefore volume emissivity and absorption coefficient
  # (epsilon_0 and kappa_0, values from pararray2) are multiplied
  # with gaussian functions to obtain correct line shape. 
  # functions are normalised in a way that epsilon_0, kappa_0
  # are the INTEGRATED emissivity
  # or respectively absorption coefficient of the gaussian function
  for i in range(steps + 1):   
      epsilon_v.append(gauss(pararray2[1][i], pararray2[2][i], \
                        pararray2[0][i]))
      kappa_v.append(gauss(pararray2[1][i], pararray2[2][i], pararray2[3][i]))

      # print "epsilon_i(11.3)", epsilon_v[i].gaussfunc(11.3)
      # print "epsilon_i(10)", epsilon_v[i].gaussfunc(10)
      # print "epsilon_i(8)", epsilon_v[i].gaussfunc(8)
      # print "kappa_i(11.3)", kappa_v[i].gaussfunc(11.3)

  # pause = input()

  # the next part creates lists of the functions needed
  # during the los integration. this is done via summation along the
  # positions of the cartesian grid which are lying in the line of sight  

  # absorption coefficient at gridpoints
  def kappajv(j):
      def kappav(v):
          return float(kappa_v[j].gaussfunc(v)) 
      return kappav
  kappa2_v = []
  for j in range(steps + 1):
      kappa2_v.append(kappajv(j))

  # linear approximation for absorption coefficient
  def kappastepjv(j):
      def kappastepv(v):
          return (kappa_v[j+1].gaussfunc(v) - kappa_v[j].gaussfunc(v))\
              /float(gbl._globals['compound']['dint'])
      return kappastepv
  kappastep_v = []
  for j in range(steps):
      kappastep_v.append(kappastepjv(j))
      # print 'kappastep(11.3) ',kappastep_v[j](11.3) 
      # print 'kappastep',kappastep_v[j](10) 
  # emissivity at gridpoints
  def epsilonjv(j):
      def epsilonv(v):
          return float(epsilon_v[j].gaussfunc(v))
      return epsilonv
  epsilon2_v = []
  for j in range(steps+1):
      epsilon2_v.append(epsilonjv(j))

  # linear approximation of emissivity
  def epsilonstepjv(j):
      def epsilonstepv(v):
          return (epsilon_v[j+1].gaussfunc(v) - epsilon_v[j].gaussfunc(v))\
                  /float(gbl._globals['compound']['dint'])
      return epsilonstepv
  epsilonstep_v = []
  for j in range(steps):
      epsilonstep_v.append(epsilonstepjv(j))
      # print 'epsilonstep(11.3) ', epsilonstep_v[j](11.3) 
  """

  I_v = [] # list for integrated intesities (each entry for one 
  # sampling velocity)

  for vel in n_velocities:
    #print 'velocity number', vel
    #print 'velocity', velocities[vel]
    # pause = input('..vel in rad_trans..')
    I = [] # initialise list for accumulated intensity at
    # every gridpoint on the los
    Inext = gbl._globals['compound']['I_bg']
    I.append(Inext)
    # background intensity needed as starting value


    # integration along line of sight at velocity v
    for i in n_steps:
      k = kappa_vobs[int(i)][vel]
      kstep = kappa_vobs_step[int(i)][vel]
      ee = epsilon_vobs[int(i)][vel]
      eestep = epsilon_vobs_step[int(i)][vel]
      intensity = copy(Inext)  
      del Inext

      if kstep == 0 and abs(k * ds) < 10**(-10):
      # if abs(k * ds + kstep/2. * ds**2) < 10**(-8):
        # case no absorption
        # print 'no absorption'
        # print 'k:', k, 'kstep:', kstep
        # print 'ee', ee, 'estep', eestep
        Inext = ee * ds + 1/2. * eestep * ds**2. + intensity 
        # print 'Inext', Inext
        # pause = input('..no abs..')

      # elif k != 0 and kstep == 0 : 
      # elif (k > (1/float(p) * 1/2. * abs(kstep) * ds) or k != 0 and kstep == 0):
      elif k > 10**3. * abs(kstep) * ds:
        # case constant absorption 
        # print 'constant absorption'
        # print 'ds', ds
        # print 'k:',k,'kstep:', kstep
        # print 'ee', ee,'estep', eestep            
        Inext = np.exp(-k * ds) * (((ee * k + eestep*\
                    (k * ds - 1))/(k ** 2.)) * np.exp(k * ds) - \
                    ((ee * k - eestep)/(k**2.)) + intensity)
        # print 'Inext', Inext
        # pause = input('...constant...')

      else: 
      # case none-constant absorption (k1 != 0)
        from Etilde import Ereal, Eimag
        # print 'non-constant absorption'
        # print 'ds
        a = k/cmath.sqrt(2. * kstep)
        b = (k + kstep * ds)/cmath.sqrt(2. * kstep)
        print('Kappa, kappa step:', k, kstep)
        print('a, b:', a, b)
        # print 'k:', k,'kstep:', kstep 
        # print 'ee', ee,'estep', eestep
        # print 'a:', a,'b:', b
        # print 'ds:', ds
        # print 'k + kstep * ds', k + kstep * ds
        # print 'cmath.sqrt(2. * kstep)', cmath.sqrt(2. * kstep)
        # print 'complete', (k + kstep * ds)/cmath.sqrt(2. * kstep)
        if kstep > 0:       
          Inext = (eestep/kstep * (1 - np.exp(-k * ds - kstep/2.*ds**2.))\
                  -(ee * kstep - eestep * k)/kstep * \
                  cmath.sqrt(np.pi/(2.* abs(kstep)))* \
                  (np.exp(a**2. - b**2.) * Ereal(a) - Ereal(b)) + \
                  intensity * np.exp(-k * ds - kstep/2. * ds**2.)).real 
        elif kstep < 0:
          Inext = (eestep/kstep * (1 - np.exp(-k * ds - kstep/2. * ds**2.))\
                  -(ee * kstep - eestep * k)/kstep * \
                  cmath.sqrt(np.pi/(2.* abs(kstep)))* \
                  (np.exp(a**2. - b**2.) * Eimag(a) - Eimag(b)) + \
                  intensity * np.exp(-k * ds - kstep/2. * ds**2.)).real
        # print 'Inext', Inext
        # pause = input('...general...')

      I.append(Inext)
      #print 'Inext', Inext #Silence
    # print 'I', I
    # print 'I[-1]', I[-1]
    I_v.append(I[-1])
  print('Integrated radiative transfer:', I_v)
  input()

  #print 'integrated I for all velocities \n', I_v
  # pause = input('..rad_trans for one los done..')
  return I_v

