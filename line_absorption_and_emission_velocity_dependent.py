# calculate averaged line intensity and optical depth for each 
# ensemble. 
# Ensembles representing clump and interclump medium will be treated 
# seperately and summed up later  
##############################################################
##########################################################
import pdb
import numpy as np
import scipy as sc
import time
from copy import copy
from scipy.interpolate import griddata
# griddata: routine for interpolation on grid
from scipy.integrate import quad
import multiprocessing as mp    # enable parallel computations
from tqdm import tqdm
# general purpose integration
#from PDR import _globals
import sys
#sys.path.append('../')
import PDRfunctions as pf
# edit: Craig, 14.10.2019
import globals as gbl
#end edit

def lineAbsorptionEmission():
  """
  This is to temporarily separate this code into a function.
  Adapted from the Bruckmann version.
  Created on 15.10.2019 by Craig. Edit as soon as possible!!
  Modified for optimisation on 18.10.2019 by Craig.
  """
  
  # Initialise a global iterator
  ens_vel = []
  for i in gbl._globals['compound']['ranges']['npoints']:
    ens_vel.append(np.array(gbl._globals['compound']['ens'][i].velocity[2]))
  # velocity of each ensemble (z-component)
  gbl._globals['compound']['ranges']['ens_vel'] = ens_vel
  gbl._globals['constants']['pn_gauss'] = 5 
  gbl._globals['constants']['N_gauss'] = 1000 
  # binomial or poisson distribution will be replaced
  # by gauss if p*n>pn_gauss and N_j > N_gauss
  
  # number of different species to be simulated
  speciesnumber = gbl._globals['compound']['number']
  # number of species of interest * transition of interest
  gbl._globals['constants']['dust_numb'] = 192 #dust starts with 192 in frequency dat
  ################################################################################
  gbl._globals['compound']['ranges']['vbin'] = np.linspace(gbl._globals['compound']['vel'] - gbl._globals['compound']['d_vel'],\
                                                gbl._globals['compound']['vel'] + gbl._globals['compound']['d_vel'],\
                                                gbl._globals['compound']['nstep']) 
  # 'bin-velocities' 
  gbl._globals['compound']['ranges']['Dv'] = gbl._globals['compound']['ranges']['vbin'][1] - gbl._globals['compound']['ranges']['vbin'][0] 
  gbl._globals['compound']['ranges']['vobs'] = gbl._globals['compound']['ranges']['vbin']  # "observed" sampling velocities (for now: vbin=vobs) 
  gbl._globals['compound']['ranges']['number vbin'] = range(len(gbl._globals['compound']['ranges']['vbin']))
  gbl._globals['compound']['ranges']['number vobs'] = range(len(gbl._globals['compound']['ranges']['vobs']))
  # Read in grids
  rhoMassUV_tau = []
  negative_tau = np.nan
  with open(gbl.INPUTPATH+'tau_linecentre.dat') as ta:
    next(ta) # skips header
    for line in ta:
      rhoMassUV_tau.append(line.split()[0:3]) # table rho, mass, uv
      #rhoMassUV_tau.append(line.split()[0:1])
      #print rhoMassUV_tau
      #pause = input('rhoMassUV_tau')
  # Some more iteration ranges
  gbl._globals['compound']['ranges']['rhoMassUV_tau'] = np.array(rhoMassUV_tau).astype(np.int)
  gbl._globals['compound']['ranges']['number rhoMassUV_tau'] = range(len(gbl._globals['compound']['ranges']['rhoMassUV_tau']))
  
  clumps = gbl._globals['compound']['ens']
  number_ensembles = gbl._globals['compound']['npoints']
  alpha = gbl._globals['constants']['alpha']
  gamma = gbl._globals['constants']['gamma']
  nsigma = gbl._globals['constants']['nsigma']
  nspe = gbl._globals['compound']['nspe']
  npoints = range(gbl._globals['compound']['npoints'])
  n_nspe = range(nspe)
  rhoMassUV_tau = np.array(rhoMassUV_tau).astype(np.int)
  n_rhomass = range(len(rhoMassUV_tau))
  
  # Min/Max for Rho-Mass-UV
  #init with first Rho-Mass-UV element
  max_arr = [int(gbl._globals['compound']['ranges']['rhoMassUV_tau'][0][0]),int(gbl._globals['compound']['ranges']['rhoMassUV_tau'][0][1]),int(gbl._globals['compound']['ranges']['rhoMassUV_tau'][0][2])]
  min_arr = [int(gbl._globals['compound']['ranges']['rhoMassUV_tau'][0][0]),int(gbl._globals['compound']['ranges']['rhoMassUV_tau'][0][1]),int(gbl._globals['compound']['ranges']['rhoMassUV_tau'][0][2])]
  #find and replace with the absolute min and maxima
  for i in gbl._globals['compound']['ranges']['number rhoMassUV_tau']:
    if int(gbl._globals['compound']['ranges']['rhoMassUV_tau'][i][0]) > max_arr[0]:
      max_arr[0] = int(gbl._globals['compound']['ranges']['rhoMassUV_tau'][i][0])
    if int(gbl._globals['compound']['ranges']['rhoMassUV_tau'][i][0]) < min_arr[0]:
      min_arr[0] = int(gbl._globals['compound']['ranges']['rhoMassUV_tau'][i][0])
    if int(gbl._globals['compound']['ranges']['rhoMassUV_tau'][i][1]) > max_arr[1]:
      max_arr[1] = int(gbl._globals['compound']['ranges']['rhoMassUV_tau'][i][1])
    if int(gbl._globals['compound']['ranges']['rhoMassUV_tau'][i][0]) < min_arr[1]:
      min_arr[1] = int(gbl._globals['compound']['ranges']['rhoMassUV_tau'][i][1])
    if int(gbl._globals['compound']['ranges']['rhoMassUV_tau'][i][2]) > max_arr[2]:
      max_arr[2] = int(gbl._globals['compound']['ranges']['rhoMassUV_tau'][i][2])
    if int(gbl._globals['compound']['ranges']['rhoMassUV_tau'][i][2]) < min_arr[2]:
      min_arr[2] = int(gbl._globals['compound']['ranges']['rhoMassUV_tau'][i][2])
  gbl._globals['compound']['ranges']['min_arr'] = min_arr
  gbl._globals['compound']['ranges']['max_arr'] = max_arr
  if gbl._globals['verbose']: print('\nMaximum Rho,Mass,UV', max_arr, '\nMinimum Rho,Mass,UV: ', min_arr)
  #pause = input('rhoMassUV_tau read? press enter!...')
  ####################################################
  tauAv =  [None] * gbl._globals['compound']['nspe'] # .. tau array for each species
  for i in gbl._globals['compound']['ranges']['nspe']:  
    # loop over species
    # get tau for the right species
    tauAv[i] = [] #Taus of this species will be appended here
    ###########################################common grid loop
    with open(gbl.INPUTPATH+'tau_linecentre.dat') as ta:
      # import tau-profile grid.
      # line-centre (peak) values for tau.
      # FWHM may vary
      # print 'gbl._globals['compound']['number'][i]', gbl._globals['compound']['number'][i]
      next(ta) #skips header
      for j, line in enumerate(ta): #run through tau_linecentre
        tautemp = line.split()[(2 + gbl._globals['compound']['number'][i]):(3 + gbl._globals['compound']['number'][i])][0]
        tautemp = float(tautemp)
        #print tautemp
        if tautemp > 0: tauAv[i].append(tautemp)
        elif tautemp == 0: tauAv[i].append(1e-100)
        # tau must not be exactly zero as log(tau) is used
        else:
          tauAv[i].append(1e-100)
          if np.isnan(negative_tau): negative_tau = tautemp
          elif abs(tautemp) > abs(negative_tau):
          # find largest abs(tau) with a negative tau
            negative_tau = tautemp
  if gbl._globals['check_negative'] == 0:
    if not np.isnan(negative_tau):
      if gbl._globals['verbose']:
        print('\nFor some clumps negative opacities have been found. The negative opacity with the largest absolute value is \n')
        print(negative_tau, '\n')
      #pause = input("negative tau's will be set to 1e-100. calculation of ensemble averaged optical depth is not possible with the current set-up as grid interpolations are done in log-space. Do you want to go on?")
  tauAv = np.array(tauAv, float) #transform lists-> array
  tauAvLog = np.log10(tauAv) * 10 # is used for interpolation
  gbl._globals['compound']['ranges']['tauAv'] = tauAv
  gbl._globals['compound']['ranges']['log tauAv'] = tauAvLog
  #print 'tauAv', tauAv
  #print 'tauAvLog', tauAvLog
  #pause = input('...')
  # Read in intensities
  rhoMassUV_intensity = []
  with open(gbl.INPUTPATH+'Tb_linecentre.dat') as ep:
    next(ep) #skips header
    for line in ep:
      rhoMassUV_intensity.append(line.split()[0:3])
      # table: rho, mass, uv for intensities
  gbl._globals['compound']['ranges']['rhoMassUV_intensity'] = rhoMassUV_intensity
  intensityClGrid = [None] * gbl._globals['compound']['nspe']
  # Average the intensity for each clump
  for i in gbl._globals['compound']['ranges']['nspe']: 
    # import grid of line intensities for each species 
    # line intensities are line integrated
    intensityClGrid[i] = [] #creates list for each species
    #pause = input('intensitygrid loop')
    with open(gbl.INPUTPATH+'Tb_linecentre.dat') as ep:
      next(ep)
      for line in ep:
        intensity_temp = float(line.split()[(2 + gbl._globals['compound']['number'][i]):\
                              (3 + gbl._globals['compound']['number'][i])][0])                  
        if intensity_temp > 0: intensityClGrid[i].append(intensity_temp)
        if intensity_temp <= 0: intensityClGrid[i].append(1e-100)
        # must not be exactly zero as log is used
  intensityClGrid = np.array(intensityClGrid, float)
  intensityClGridLog = np.log10(intensityClGrid) * 10 # is used for interpolation
  gbl._globals['compound']['ranges']['intensity Cl Grid'] = intensityClGrid
  gbl._globals['compound']['ranges']['log intensity Cl Grid'] = intensityClGridLog
  # Calculate averaged opacity and emissivity for each pixel
  sigma_cl_j =  gbl._globals['sigma']['sigma_j']
  # intrinsic clump linewidth standard deviation)
  gbl._globals['compound']['iterator'] = []
  if gbl._globals['compound']['npoints']<0:
    print('Parallel')
    #pf.progressBar(gbl._globals['compound']['iterator'], gbl._globals['compound']['npoints'], length=50)
    with tqdm(total=gbl._globals['compound']['npoints'], desc='Progress') as gbl._globals['compound']['progress']:
      pool = mp.Pool(4)
      #print([pix for pix in gbl._globals['compound']['ranges']['npoints']])
      test=pool.starmap_async(I_calc, [(i, pix) for (i, pix) in enumerate(gbl._globals['compound']['ranges']['npoints'])]).get()
      pool.close()
      print('FUV')
      for result in test:
        print(gbl._globals['compound']['ens'][result[0]].FUV)
        gbl._globals['compound']['ens'][result[0]].tau = result[1]
        gbl._globals['compound']['ens'][result[0]].inten = result[2]
        gbl._globals['compound']['ens'][result[0]].Mens_sum = result[3]
    #input()
  else:
    print('Series')
    with tqdm(gbl._globals['compound']['npoints'], ascii=True, desc='Progress') as gbl._globals['compound']['progress']:
      for pix in gbl._globals['compound']['ranges']['npoints']:
        result = I_calc(pix, pix)
        gbl._globals['compound']['ens'][result[0]].tau = result[1]
        gbl._globals['compound']['ens'][result[0]].inten = result[2]
        gbl._globals['compound']['ens'][result[0]].Mens_sum = result[3]
  # pause = input('calculation of averaged opacities and emissivities finished')
  return
  
# P A R A L L E L I S E
def I_calc(idx, pix, debug=False):
  '''
  This is my attempt to parallelise this section. Don't judge me...
  Created on 22.10.2019 by Craig.
  '''
#for pix in gbl._globals['compound']['ranges']['npoints']: #loop over pixels
  #print(pix, gbl._globals['compound']['ranges']['npoints'], gbl._globals['compound']['ranges']['nspe'])
  # saves time of each voxel in list, used later for runtime log
  #clumps = gbl._globals['compound']['ens']
  #number_ensembles = gbl._globals['compound']['gbl._globals['compound']['ranges']['npoints']']
  #alpha = gbl._globals['constants']['alpha']
  #gamma = gbl._globals['constants']['gamma']
  #nsigma = gbl._globals['constants']['nsigma']
  #nspe = gbl._globals['compound']['nspe']

  gbl._globals['runtimes']['line_absorption'].append(time.time())
  run = time.time()- gbl._globals['runtimes']['line_absorption'][0]
  if gbl._globals['verbose']:
    print('pixel number ', pix, 'of', gbl._globals['compound']['npoints'])
    print('running ', run, 's of an estimated ', run/((pix+1.)/gbl._globals['compound']['npoints']), 's') 
  sigma_cl_j =  gbl._globals['sigma']['sigma_j']
  sigma_ens_j = ( (np.array(gbl._globals['sigma']['sigma_ens_j']))**2 + gbl._globals['compound']['ens'][pix].v_dispersion**2 )**0.5
  # sum of inner clump turbulence + ensemble v_dispersion
  Mens = gbl._globals['compound']['ens'][pix].Mens # ensemble mass in solar masses
  totalSUM = 0
  if Mens == 0:
    if gbl._globals['verbose']: print('Mens = 0; no mass contained in pixel ', pix, \
                      ', setting line optical depth and intensity to zero for this pixel')
    tau = np.zeros([gbl._globals['compound']['nspe'], len(vobs)])
    inten = np.zeros([gbl._globals['compound']['nspe'], len(vobs)])
  else:
    Mu = gbl._globals['compound']['ens'][pix].Mu    # highest clump mass
    Ml = gbl._globals['compound']['ens'][pix].Ml    # lowest clump mass
    #print(Mu, Ml)
    MuLog = np.log10(Mu) * 10  # 0.1 -> -10 while 10 stays 10
    MlLog = np.log10(Ml) * 10
    #print MuLog
    #pause = input('Mulog')
    if MlLog < gbl._globals['compound']['ranges']['min_arr'][1] or MuLog > gbl._globals['compound']['ranges']['max_arr'][1]:
      print('WARNING: Outside mass grid')
      pause = input('Warning read? press enter!...')
    MLogArray = np.arange(MlLog, MuLog + 10, 10)
    #print(MLogArray)
    uv_log = np.log10(gbl._globals['compound']['ens'][pix].FUV) * 10
    rho_ens = gbl._globals['compound']['ens'][pix].rho_ens
    # rho_ens in cm^-3; needed to calculate surface density
    # of different clumps
    interpolationPoints = []
    number_masspoints = len(MLogArray)
    n_masspoints = range(number_masspoints) #another range
    number = []
    for M_log in MLogArray:
      n_s_binned = ((10 ** (M_log / 10.)) ** (1 - 3 / gbl._globals['constants']['gamma']) * \
                  sum((10 ** (M_log2 / 10.)) ** (1 + 3 / gbl._globals['constants']['gamma'] - gbl._globals['constants']['alpha'])\
                  for M_log2 in MLogArray)) / \
                  (sum((10 ** (M_log2 / 10.0)) ** (2 - gbl._globals['constants']['alpha']) for M_log2\
                  in MLogArray)) * rho_ens / 1.91
      # clump surface density in cm^-3 for discrete set-up
      # factor 1.91 is necessary to convert averaged clump density 
      # density into surface density
      #print '...clump surface density', n_s_binned
      #raw = input ('density OKOK?')
      if n_s_binned > 10**(gbl._globals['compound']['ranges']['max_arr'][0]/10) or n_s_binned < 10**(gbl._globals['compound']['ranges']['min_arr'][0]/10):
      # CHECK: needs to be adapted if grid is changed- automatisiert - done
        print('Surface density: ', n_s_binned)
        print('At voxel ', pix)
        sys.exit('WARNING: surface density lies outside grid ...exiting...')
      n_s_binned_log = np.log10(n_s_binned) * 10
      if gbl._globals['verbose']: print('density', n_s_binned)
      interpolationPoints.append([n_s_binned_log, M_log, uv_log])
      number.append((Mens * (10 ** (M_log / 10.0)) ** (1 - gbl._globals['constants']['alpha'])) / \
                    (sum((10 ** (M_log2 / 10.0)) ** (2 - gbl._globals['constants']['alpha']) for \
                    M_log2 in MLogArray)))
    #print 'interpolatiogbl._globals['compound']['ranges']['npoints'] Rho,Cl-Mass,UV', interpolatiogbl._globals['compound']['ranges']['npoints']
    #print 'number of clumps for each clump mass', number
    #raw = input (' OKOK?')
    if debug: print('\n', interpolationPoints, '\n')
    tauCl =  [None] * gbl._globals['compound']['nspe']
    intensityCl =  [None] * gbl._globals['compound']['nspe']
    for sp in gbl._globals['compound']['ranges']['nspe']:  
    # interpolate opacity and intensity for each species   
      tauCl[sp] = []
      intensityCl[sp] = []
      # print 'tauAvLog[sp]', tauAvLog[sp]
      tauCl[sp].append(10 ** (griddata(gbl._globals['compound']['ranges']['rhoMassUV_tau'], gbl._globals['compound']['ranges']['log tauAv'][sp], \
                        interpolationPoints, method='linear')/10.))
      # interpolated line peak opacity for each clump 
      #print 'tauCl',tauCl[sp]
      intensityCl[sp].append(10 ** (griddata(gbl._globals['compound']['ranges']['rhoMassUV_intensity'], \
                              gbl._globals['compound']['ranges']['log intensity Cl Grid'][sp], interpolationPoints, method='linear')/10.))          
      # interpolated, line averaged intensity in Kkm/s
    if debug:
      print(intensityCl)
      print(tauCl)
    RclTab = []
    for i in n_masspoints:
      # calculate clump radii
      if gbl._globals['verbose']: print('log density', interpolationPoints[i][0])
      Rcl = ((3 / (4.0 * np.pi) * (10 ** (MLogArray[i] / 10.0) * \
            gbl._globals['constants']['M_sol']) / (10 ** (interpolationPoints[i][0] / \
            10.0) * gbl._globals['constants']['M_H'] * 1.91)) ** (1 / 3.0)) / \
            gbl._globals['constants']['pc']
      #print 'Rcl', Rcl
      RclTab.append(Rcl)
    if gbl._globals['verbose']: print('Radius:', RclTab)
    #print 'interpolated tau:\n ', tauCl
    # print tauCl
    # print 'tauCl[0][0]', tauCl[0][0] #list build as: tauCl[species][0][mass]
    # pause = input('tauCl ok?')
        ############### Dust loop for converting Flux[Jansky] to Intensity[K]
        ## 1 dust exception
        ## cloud fluxes were normalized to distances of 10 parsec
        ## goal to get markus Flux[jansky] to I[Kelvin] like all other treatment
    for i in n_masspoints: #all clumps need individual correction
      for sp in gbl._globals['compound']['ranges']['nspe']:
        if gbl._globals['compound']['number'][sp] >= gbl._globals['constants']['dust_numb']: #means species is dust
          #print 'clump size', RclTab[i]
          jansky_to_WmmHz = 10**-26 #Jansky -> W/m^2/Hz
          intensity_area_factor = np.pi * (np.arcsin(RclTab[i]/10.) )**2 #normalized for clump in 10PC distance
          intensityCl[sp][0][i] = intensityCl[sp][0][i] / intensity_area_factor * jansky_to_WmmHz #jansky_to_W/m^2/Hz
          rayleigh = 2 * gbl._globals['constants']['kb'] / gbl._globals['constants']['c']**2 * (gbl._globals['compound']['frequency'][sp]*10**9)**2
          intensityCl[sp][0][i] = intensityCl[sp][0][i] / rayleigh  #von I[J/s/m^2/Hz/sr] -> I[Kelvin] (2k_b*f^2/c^2)
          #print 'intensity TAB CL', intensityCl[sp][0][i]
    #print gbl._globals['compound']['species']
    #print 'intensity TAB CL', intensityCl
    #pause = input('intensity correct?')
    #print 'interpolated intensity:\n ',
    #print intensityCl
    #print 'calculated radii (RclTab)',
    #print RclTab
    #pause = input('tau and intensity ok?')
    if np.any(np.isnan(tauCl)):
      print('Variable interpolationPoints: ', interpolationPoints)
      sys.exit('found "nan" in tauCl...exiting...')
    if np.any(np.isnan(intensityCl)):
      sys.exit('found "nan" in intensityCl...exiting...')
    # S T A T I S T I C S
    # v-dependence
    # Calculate mass fractions for v-bins
    for i in range(1, len(sigma_ens_j)):
      # check whether all sigma_ens_j are equal for all j 
      # (in this case some of the following steps can be simplified)
      if sigma_ens_j[i] == sigma_ens_j[0]: sigma_j_allequal = True
      else: 
        sigma_j_allequal = False
        break
    DeltaN_ji = np.zeros([number_masspoints, len(gbl._globals['compound']['ranges']['number vbin'])])
    # scale each element of DeltaN_ji 
    #print 'number of clumps for each clump mass', number
    for ma in n_masspoints:
      for v in gbl._globals['compound']['ranges']['number vbin']:
        if abs((gbl._globals['compound']['ranges']['vbin'][v]-gbl._globals['compound']['ranges']['ens_vel'][pix])/sigma_ens_j[ma]) >= 4.7534243088229 : 
          #exception catch for beeing further away from gauss than 1 million portion
          DeltaN_ji[ma,v] = 0  
        else:
          DeltaN_ji[ma, v] = number[ma]/(np.sqrt(2 * np.pi) * sigma_ens_j[ma]) * np.exp(-1/2.*((gbl._globals['compound']['ranges']['vbin'][v]-gbl._globals['compound']['ranges']['ens_vel'][pix])/(sigma_ens_j[ma]))**2) * gbl._globals['compound']['ranges']['Dv']
          if np.isnan(DeltaN_ji[ma,v]) == True:
            print('Voxel number: ', pix)
            pause = input('NaN in DealtaN_ji found!!')
          #print'\n'
          #print 'vbin[v]', vbin[v]
          #print 'mass in % in v-bin', DeltaN_ji[ma,v] / number[ma] *100.
          #print 'clumb number total', number[ma]
          #print 'bin v, ensemble v',  vbin[v], ens_vel[pix]
          #print 'v-difference to cental v', vbin[v]-ens_vel[pix]
          #print 'difference in number of sigmatot', (vbin[v]-ens_vel[pix])/sigma_ens_j[ma]
            # exception warning of disperion so high that too much mass is outside the velocitiy range
            # is important for example for galactic center or selected v-range to small
      totalbinnedmass = sum(DeltaN_ji[ma,:])
      #print 'Anteil abgedeckter Masse [%]', totalbinnedmass/ number[ma] *100
      if totalbinnedmass / number[ma] <= 0.95: #less than 95%
        print('Achtung! Anteil abgedeckter Masse [%] nur:', totalbinnedmass/number[ma]*100)
        print('increase v-range or reduce v-dispersion')
        #pause = input('Proceed?')  
    #print(DeltaN_ji)
        ############### dust with mass
        ## 2 dust difference
        ## for dust all mass is at each velocity at the same time
        ## only first v would be enough calculated, rest could be filled
    #print 'gbl._globals['compound']['number']er', gbl._globals['compound']['number'][sp]
    #if gbl._globals['compound']['number'][sp] >= dust_numb: #means species is dust
    DeltaN_ji_dust = np.zeros([number_masspoints, len(gbl._globals['compound']['ranges']['number vbin'])])
    for mass in n_masspoints:
      for v in gbl._globals['compound']['ranges']['number vbin']:
        DeltaN_ji_dust[mass][v] = number[mass] #writes full masses of this clump type to first velocity bins
        #dust has total mass at every velocity only first will be needed for calculations
    #print 'dust-masses', DeltaN_ji_dust
    #number_v_dust = copy(DeltaN_ji_dust[:,v]) #scaling can be done normaly
    #print 'number_v_dust', number_v_dust
    #print('DeltaN_ji', DeltaN_ji)
    #pause = input('...DeltaN_ji OK?')
    tau_v = []
    inten_v = []       
    for v in gbl._globals['compound']['ranges']['number vbin']:
      #print('velocity bin', gbl._globals['compound']['ranges']['vbin'][v])
      #print 'MA',ma
      #print 'V', v
      if DeltaN_ji[number_masspoints-1,v] == 0: #no masses at this velocity
        # dust get always activated as long as there is any mass at any v
        tau_v_i_Av = np.zeros((gbl._globals['compound']['nspe'], len(gbl._globals['compound']['ranges']['number vbin'])))                
        I_v_i_Av = np.zeros((gbl._globals['compound']['nspe'], len(gbl._globals['compound']['ranges']['number vbin'])))
        #print len(tau_v_i_Av), type(tau_v_i_Av)
        tau_v.append(tau_v_i_Av)
        inten_v.append(I_v_i_Av)
        #print 'no mass at this speed, skippping v-bin loop'
      else:
        erster_durchlauf=1        
        #print ("\033c") #clear the terminal
        # loop over v-bins 
        #print 'scaling velocity interval number: ', v
        #print 'pixel number', pix
        # pause = input('...')
        Lbox = gbl._globals['compound']['dint'] / gbl._globals['constants']['pc']
        number_v = copy(DeltaN_ji[:, v]) #selects only masses at this velocity
        #print number_v
        #pause = input('number_v?')
        #if number_v[(number_masspoints-1)] == 0 : 
        #  print 'no mass with this v found'
        #print 'number_v', number_v
        #print 'Lbox, before scaling', Lbox
        #print 'number of clumps per mass interval, before scaling', number_v
        if gbl._globals['binomial_poisson']:
          #  Use BINOMIAL distribution
          Lbox = Lbox / np.sqrt(float(number_v[(number_masspoints - 1)]))
          number_v = number_v / float(number_v[(number_masspoints - 1)])
          #print 'number of clumps, scaled to N_nM=1',
          #print(number_v)
          Lbox = Lbox * np.sqrt(float(gbl._globals['statistic']['Nmin']))
          number_v = number_v * float(gbl._globals['statistic']['Nmin'])
          #print 'number of clumps, scaled to Nmin',
          #print number_v #describes number of clumps local-box
          #print 'Lbox', Lbox
          #print 'RclTab', RclTab
          if np.pi * max(RclTab) ** 2 / Lbox ** 2 >= 1:
              # pixel smaller than clump, p>1. increase pixel size
              scale = np.ceil(np.pi * max(RclTab) ** 2 / Lbox ** 2)
              Lbox = Lbox * np.sqrt(scale)
              #print 'Lbox to small...increasing size' , 'by scale', scale
              number_v = number_v * scale
          number_v = np.around(number_v) # round to integer value 
          # pause = input('Binomial ok?')
        else:
          # Use Poisson distribution
          scale = 1./float(number_v[number_masspoints - 1]) * 100
          # scale to high numbers of clumps, here:
          # 100 high mass clumps
          Lbox = Lbox * np.sqrt(scale)
          number_v = number_v * scale
          #print 'Lbox', Lbox
          #print 'RclTab', RclTab
          if np.pi * max(RclTab) ** 2 / Lbox ** 2 >= 1:
            #print 'Lbox to small...increasing size'
            scale = np.ceil(np.pi * max(RclTab) ** 2 / Lbox ** 2) * 100.
            #print 'scale',scale
            Lbox = Lbox * np.sqrt(scale)
            number_v = number_v * scale   
          number_v = np.around(number_v) # round to integer value
        # Scaling Finished
        #print 'number of clumps per mass interval (tau), round', number_v
        #print 'pmax', np.pi * RclTab[number_masspoints - 1] ** 2 / Lbox ** 2
        #pause = input('after scale...')
        pTab = []
        probabilityTab = []
        expectedValTab = []
        standardDeriTab = []
        lowerupperTab = []
        for ma in n_masspoints:
          pTab.append(np.pi * RclTab[ma] ** 2 / Lbox ** 2) #A_clump/A_box
          expectedValTab.append(number_v[ma] * pTab[ma])
          standardDeriTab.append(np.sqrt(number_v[ma] * pTab[ma] * (1 - pTab[ma])))
          lower = max(0, np.floor(expectedValTab[ma] - gbl._globals['constants']['nsigma'] * standardDeriTab[ma]))
          upper = min(number_v[ma], np.ceil(expectedValTab[ma] + gbl._globals['constants']['nsigma'] * standardDeriTab[ma]))
          # print 'lower', lower # untere Grenze der Kollisionen inerhalb der nsigma(3)
          # print 'upper', upper # obere Grenze der Kollisionen
          lowerupperTab.append([lower, upper])
          if gbl._globals['binomial_poisson']:
            if (expectedValTab[ma] > gbl._globals['constants']['pn_gauss'] and number_v[ma] > gbl._globals['constants']['N_gauss']):
              # use gauss!
              g = pf.gauss(expectedValTab[ma], standardDeriTab[ma])
              probabilityTab.append(g.gaussfunc)
              pause = input('gauss!!...')
            else:
              # use binomial 
              #print 'binomfunction input n and p',number_v[ma], pTab[ma]
              b = pf.binomial(number_v[ma], pTab[ma]) # n and p for binominal 
              # binominalfunctions if called can give chances for k hits (for certrain n and p) 
              probabilityTab.append(b.binomfunc)
              #print ' probability for a k=0',probabilityTab[ma](0) #k=0
              #print ' probability for a k=1',probabilityTab[ma](1)
          else:
            if (expectedValTab[ma] > gbl._globals['constants']['pn_gauss'] and number_v[ma] > gbl._globals['constants']['N_gauss']):
              # use gauss
              g = pf.gauss(expectedValTab[ma], standardDeriTab[ma])
              probabilityTab.append(g.gaussfunc)
              pause = input('gauss!!...')
            else:
              # use poisson
              po = poisson(expectedValTab[ma])
              probabilityTab.append(po.poissonfunc)
        #print(probabilityTab[0](0), probabilityTab[0](1), probabilityTab[1](0), probabilityTab[1](1))
        #print 'pTab', pTab
        #print 'expectedValTab', expectedValTab
        #print 'standardDeriTab', standardDeriTab
        #print 'lowerupperTab', lowerupperTab     
        #pause = input('statistic ok?')
        # Create combination tables only for first run, (not for first velocity bin)
        if erster_durchlauf == 1:
          erster_durchlauf = 0 #sets it to false to creates combination only one time 
          a = pf.AllCombi(lowerupperTab) #creates AllCombiclass with min and max Collisions
          a.createArray() # creates combination table  with length of all combinations, and for each masssize on other axes
          combis = a.combinations #takes combination table
          #print 'combis' , combis
          #pause = input('all combis ok?')
        else:
          if not sigma_j_allequal:
            a = pf.AllCombi(lowerupperTab)
            a.createArray()
            combis = a.combinations
          else: pass
        n_combis = range(len(combis))
        #pause = input('combinations created ok?')
        tau_v_i = []
        inten_v_i = []
        # initialise arrays for "final" results (for this velocity)
        # for i in range(nspe): 
        #     tau_v_i[i] = np.zeros([nspe, len(vobs)])
        #     eps_v_i[i] = np.zeros([nspe, len(vobs)])
        p = 0
        totalSUM = copy((combis[-1]*10**(MLogArray/10)).sum())
        #pdb.set_trace()
        for c in n_combis:
        # loop over possible combinations of clumps
          #print('combis:', combis)
          #input()
          tau_x_i = np.zeros([gbl._globals['compound']['nspe'], len(gbl._globals['compound']['ranges']['vobs'])])
          inten_x_i = np.zeros([gbl._globals['compound']['nspe'], len(gbl._globals['compound']['ranges']['vobs'])])
          # that means tautot[species, vobs]
          ptot = 1
          for ma in n_masspoints:
            # loop over mass intervals
            for sp in gbl._globals['compound']['ranges']['nspe']:
              #pdb.set_trace()
                            
                            
                            ##### for dust this v loop is not nessesary 
              for vo in gbl._globals['compound']['ranges']['number vobs']:
                  #print('vo:', gbl._globals['compound']['ranges']['vobs'][vo])
                  #print 'v', v
                  #print 'sp', sp
                  #print 'vo', vo
                  #print 'vbin[v]', vbin[v]
                  #print 'vobs[vo]', vobs[vo]
                  #print 'sigma_cl_j[ma]', sigma_cl_j[ma]
                  #print np.exp(-1/2. * ((vbin[v]-vobs[vo])/(sigma_cl_j[ma]))**2)                           
                  #pause = input('...')
                                ##############################################################################################
                  # version with integral. takes more computing time
                  """
                  tau_x_i[sp, vo] = tau_x_i[sp, vo] + tauCl[sp][0][ma] * \
                                  combis[c][ma] * \
                                  quad(lambda x: np.exp(-1/2. * ((x-vobs[vo])/(float(sigma_cl_j[ma])))**2), vbin[v]-Dv/2., vbin[v]+Dv/2.)[0]
                  inten_x_i[sp, vo] = inten_x_i[sp, vo] + intensityCl[sp][0][ma] * \
                                  combis[c][ma] * \
                                  quad(lambda x: np.exp(-1/2. * ((x-vobs[vo])/(float(sigma_cl_j[ma])))**2), vbin[v]-Dv/2., vbin[v]+Dv/2.)[0]
                  """
                                ##################################################################################
                  # version without integral.
                  if gbl._globals['compound']['number'][sp] < gbl._globals['constants']['dust_numb']: #means species is non dust  
                    tau_x_i[sp, vo] = tau_x_i[sp, vo] + tauCl[sp][0][ma] * combis[c][ma] * np.exp(-1/2. * ((gbl._globals['compound']['ranges']['vbin'][v]-gbl._globals['compound']['ranges']['vobs'][vo])/(sigma_cl_j[ma]))**2)
                    inten_x_i[sp, vo] = inten_x_i[sp, vo] + intensityCl[sp][0][ma] * combis[c][ma] * np.exp(-1/2. * ((gbl._globals['compound']['ranges']['vbin'][v]-gbl._globals['compound']['ranges']['vobs'][vo])/(sigma_cl_j[ma]))**2)
                                ##################################################################################
                              
                  # if dust
                  if gbl._globals['compound']['number'][sp] >= gbl._globals['constants']['dust_numb']: #means species is dust
                                    ##################################################################################
                    tau_x_i[sp, vo] = tau_x_i[sp, vo] + tauCl[sp][0][ma] * combis[c][ma] #* \
                                  #np.exp(-1/2. * ((vbin[v]-vobs[vo])/(sigma_cl_j[ma]))**2)
                    inten_x_i[sp, vo] = inten_x_i[sp, vo] + intensityCl[sp][0][ma] * combis[c][ma] #* \
                                  #np.exp(-1/2. * ((vbin[v]-vobs[vo])/(sigma_cl_j[ma]))**2) relative mass is always 1
                                    ########### no mass fraction nessesary
                    #tau_x_i_dust[sp, vo] =  tauCl[sp][0][ma] * combis[c][ma] * 1 # 
                    #inten_x_i_dust[sp,vo] = intensityCl[sp][0][ma] * combis[c][ma] * 1
                    for v_fill in range(len(tau_x_i[sp])):
                      #print v_fill
                      tau_x_i[sp, v_fill] = tau_x_i[sp, 0] # fills all array fields of dust with same value
                      inten_x_i[sp, v_fill] = inten_x_i[sp, 0] # ...
                    #input('pause dust')
                    break #one time in v is enough for dust
                                    ##################################################################################
              ptot = ptot * float(probabilityTab[ma](int(combis[c][ma])))
              #pdb.set_trace()
          p = p + ptot
          Itemp = copy(inten_x_i)
          Ttemp = copy(tau_x_i)
          if debug: print(ma, combis[c][ma], ptot)
          if len(probabilityTab)==1: ptot = float(probabilityTab[0](combis[c][0]))
          elif len(probabilityTab)==2: ptot = float(probabilityTab[0](int(combis[c][0]))*probabilityTab[1](int(combis[c][1])))
          elif len(probabilityTab)==3: ptot = float(probabilityTab[0](int(combis[c][0]))*probabilityTab[1](int(combis[c][1]))*probabilityTab[2](int(combis[c][2])))
          #print(combis[c], ptot)
          if debug: print(ma, combis[c][ma], probabilityTab[ma](int(combis[c][ma])))
          inten_v_i.append([ptot, Itemp])
          tau_v_i.append([ptot, Ttemp])
          if debug:
            print('\nIntensity\n', inten_v_i[-1])
            print('\nTau\n', tau_v_i[-1])
            input()
          # for each combination of clumps: save probability p_x_i and related tau and
          # emissivity for each vobs and for each species
          if debug:
            print('velocity {}\ncombination {}\np {}\nI 13CO, I C+, I CO\n{}\n{}\n{}'.format(gbl._globals['compound']['ranges']['vbin'][v], combis[c], ptot, inten_x_i[0], inten_x_i[1], inten_x_i[2]))
            print('velocity {}\ncombination {}\np {}\ntau 13CO, tau C+, tau CO\n{}\n{}\n{}'.format(gbl._globals['compound']['ranges']['vbin'][v], combis[c], ptot, tau_x_i[0], tau_x_i[1], tau_x_i[2]))
            #print(combis[c], ptot, np.amax(inten_x_i[0]), np.amax(inten_x_i[1]), np.amax(inten_x_i[2]))
            #print(combis[c], ptot, np.amax(tau_x_i[0]), np.amax(tau_x_i[1]), np.amax(tau_x_i[2]))
            input()
        #print 'p =',
        #print p
        #print 'to obtain exact results p needs to be close to 1. If p is too small increase nsigma'
        # determine ensemble averaged tau and intensity for each velocity bin
        #print 'tau_v_i_total', tau_v_i[:][:]
        #print 'tau_v_i[:][0]', tau_v_i[:][0] , ' tau_v_i[:][1]', tau_v_i[:][1]
        #input(inten_v_i)
        I_v_i_Av = sum(inten_v_i[i][0] * inten_v_i[i][1] for i in range(len(inten_v_i)))
        tau_v_i_Av = -np.log(sum(tau_v_i[i][0] * np.exp(-tau_v_i[i][1] ) for i in range(len(tau_v_i)))) 
        if debug:
          for tau in tau_v_i:
            print(tau[0])
          print(I_v_i_Av)
          print(tau_v_i_Av)
          input()
        #print 'tau_v_i_Av', tau_v_i_Av
        #print len(tau_v_i_Av), type(tau_v_i_Av)
        #print 'I_v_i_Av', I_v_i_Av  
        #pause = input('vi_Av ok?')
        tau_v.append(tau_v_i_Av)
        inten_v.append(I_v_i_Av)
    tau   = sum(tau_v[i] for i in range(len(tau_v)))
    inten = sum(inten_v[i] for i in range(len(inten_v)))
    if debug:
      print('Emission x_i (single masspoint at one observing velocity)')
      print(np.shape(inten_x_i))
      print(np.shape(tau_x_i))
      print('Emission v_i Av (ensemble-averaged at one observing velocity)')
      print(np.shape(I_v_i_Av))
      print(np.shape(tau_v_i_Av))
      print('Emission v_i (number of combinations, masspoints)')
      print(np.shape(inten_v_i))
      print(np.shape(tau_v_i))
      print('Emission v (ensemble-averaged in voxel at all contributing velocities)')
      print(np.shape(inten_v))
      print(np.shape(tau_v))
      print('Emission (ensemble-averaged in voxel at one observing velocity, summed over contributions)')
      print(np.shape(inten))
      print(np.shape(tau))
      input()
    #gbl._globals['compound']['ens'][pix].tau = tau
    #gbl._globals['compound']['ens'][pix].inten = inten 
    # tau[species, vobs]
    # Intensity[species, vobs]
    # with the current input these are line centre values!
    #print 'ensemble averaged tau[species, vobs] \n', tau
    #print 'ensemble averaged Intensity[species, vobs] \n', inten 
    #pause = input('tau and inten Av ok?')
  gbl._globals['runtimes']['line_absorption'][-1] = time.time()-gbl._globals['runtimes']['line_absorption'][-1]
  gbl._globals['compound']['iterator'].append(pix)
  gbl._globals['compound']['progress'].update(1)
  #pf.progressBar(gbl._globals['compound']['iterator'], gbl._globals['compound']['npoints'], prefix='Progress:', suffix='complete', length=50) 
  #print((pix+1)/gbl._globals['compound']['npoints'])
  #print(c, standardDeriTab)
  #gbl._globals['compound']['ens'][pix].tau = tau
  #gbl._globals['compound']['ens'][pix].inten = inten 
  return (idx, tau, inten, totalSUM)

