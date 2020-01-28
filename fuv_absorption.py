# Calculate the averaged fuv extinction for each pixel 
##############################################################################
import numpy as np
from copy import copy
from scipy.interpolate import griddata # for grid interpolation
#from PDR import _globals
# edit: Craig, 14.10.2019
import globals as gbl
#end edit
import sys
sys.path.append(gbl.KOSMAPATH+'../')
import PDRfunctions as pf
from tqdm import tqdm

def fuvAbsorption():
  '''
  This is a temporary function to separate this stage.
  Adapted from the Bruckmann version.
  Created on 17.10.2019 by Craig.
  Modified for optimisation on 18.10.2019 by Craig.
  '''
  pn_gauss = 5 
  N_gauss = 1000 
  # binomial or poisson distribution will be replaced
  # by gauss if p*n>pn_gauss and N_j > N_gauss
  
  
  clumps = gbl._globals['compound']['ens']
  clumplist_A = range(gbl._globals['compound']['npoints']) 
  number_ensembles = len(clumplist_A) 
  alpha = gbl._globals['constants']['alpha']
  gamma = gbl._globals['constants']['gamma']
  nsigma = gbl._globals['constants']['nsigma']
  # import clump averaged fuv absorption, depending on clump density and mass
  RhoMassAFUV = []
  with open(gbl.INPUTPATH+'RhoMassAFUV.dat') as ta:
    for line in ta:
      RhoMassAFUV.append(line.split()) 
      # table: surface density, mass, fuv_absorption A
  RhoMassAFUV = np.array(RhoMassAFUV, float)
  # ns, mass, A not depending on FUV
  # pause = input('fuv extinction 1...')
  with tqdm(number_ensembles, ascii=True, desc='Progress') as progress:
    # P A R A L L E L I S E
    for iunique in clumplist_A:
      # for each pixel ...
      if gbl._globals['verbose']: print('ensemble number (for fuv absorption): ', iunique)
      Mens = clumps[iunique].Mens # ensemble mass in solar masses
      if Mens == 0:
        if gbl._globals['verbose']: print('no mass found in pixel ', iunique, ', setting AFUV = 0')
        AFUVEnsAv = 0
    ##
      else:
        # print 'Mens', Mens # ensemble mass in pixel i (in M_sol)
        rho_ens = clumps[iunique].rho_ens
        # rho_ens in cm^-3
        Mu = gbl._globals['compound']['ens'][iunique].Mu 
        # lagest discrete clump mass (NOT interval edge)
        Ml = gbl._globals['compound']['ens'][iunique].Ml 
        # lowest clump mass
        if gbl._globals['verbose']: print('Mu: ', Mu, '\nMl: ', Ml)
        # pause = input('?...')
        MuLog = np.log10(Mu) * 10
        MlLog = np.log10(Ml) * 10
        MLogArray = np.arange(MlLog, MuLog + 10, 10)
        # calculate constant A (Cubick et al) for the ensemble
        # A = Mens * (2 - alpha)/(mu**(2 - alpha) - ml**(2 - alpha))
        interpolationPoints = []
        # array to collect gridpoints at which clump opacities are needed
        numberInt = len(MLogArray) 
        # number of mass intervals
        if gbl._globals['verbose']: print('MLogArray', MLogArray)
        if gbl._globals['verbose']: print('numberInt', numberInt)
        # pause = input('bis hier..')
        number = [] # array to store number of clumps per mass interval  
        for M_log in MLogArray: 
          # print 'Mlog: ', M_log
          # pause = input('mass?')
          # Rcl_H_cm = (10**(M_log/10.)/\
          # gbl._globals['compound']['ens'][iunique].C)**(1./gamma) 
          # print 'calculated H-radius in cm: ', Rcl_H_cm
          # pause = input('H-radius')
          # radius of clump (for hydrogen not for species of interest) in cm
          # print 10**(m_log/10.)
          # print gbl._globals['constants']['M_sol']/gbl._globals['constants']['M_H']
          # print 10**(m_log/10.)*gbl._globals['constants']['M_sol']/\
          #      gbl._globals['constants']['M_H']
          # print '4/3 pi r^3', (4/3. * np.pi*Rcl_H_cm**3)
          n_s_binned = ((10**(M_log/10))**(1 - 3/gamma) * \
                      sum((10**(M_log2/10.))**(1 + 3/gamma - alpha) \
                      for M_log2 in np.arange(MlLog, MuLog + 10, 10)))/\
                      (sum((10**(M_log2/10.))**(2 - alpha) for M_log2 in \
                      np.arange(MlLog, MuLog + 10, 10)))*rho_ens/1.91
          # clump surface dennsity in cm^-3 for binned set-up
          # Correction factor 1.91 is necessary to convert averaged clump 
          # density into surface density
          # print 'averaged surface density cm**-3', n_s_binned  
          # n_s = 10**(M_log/10.)*gbl._globals['constants']['M_sol']/\
          #  gbl._globals['constants']['M_H']/(4/3.*np.pi*Rcl_H_cm**3)/1.91 
          # non-binned clump SURFACE density in cm^-3
          # pause = input('rho reasonable?')
          if n_s_binned > 10**7 or n_s_binned < 10**3: #ad automatic check
            if gbl._globals['verbose']: print('surface density: ', n_s_binned)
            sys.exit('surface density lies outside grid ...exiting...')
          n_s_binned_log = np.log10(n_s_binned) * 10
          interpolationPoints.append([n_s_binned_log, M_log])
          number.append((Mens * (10**(M_log/10.))**(1 - alpha))/\
                        (sum((10**(M_log2/10.))**(2 - alpha) for M_log2 in \
                        np.arange(MlLog, MuLog + 10, 10)))) 
          # pause = input()
          # number of clumps per mass interval following mass-size relation
    ##
        if gbl._globals['verbose']: print('interpolationPoints: ', interpolationPoints)
        #print 'number', number
        AFUVCl = 10**(griddata(RhoMassAFUV[:,:2], \
                (np.log10(RhoMassAFUV[:,2]) * 10), \
                interpolationPoints, method = 'linear')/10.)
        if gbl._globals['verbose']: print('AFUVCl', AFUVCl)
        # pause = input('A fuv interpolation ok?...')
        RclTab = []
        # print 'pc', gbl._globals['constants']['pc']
        for i in range(numberInt):
          # calculate Rcl via masses and densities of single clumps
          # print i
          if gbl._globals['verbose']: print('MLog: ', MLogArray[int(i)])
          if gbl._globals['verbose']: print('suface density log: ', interpolationPoints[int(i)][0])
          Rcl = ((3/(4.* np.pi)*(10**(MLogArray[int(i)]/10.) * \
                gbl._globals['constants']['M_sol'])/\
                ((10**(interpolationPoints[int(i)][0]/10.) * \
                gbl._globals['constants']['M_H'])*1.91))**(1./3.))/\
                gbl._globals['constants']['pc'] 
          RclTab.append(Rcl)
        if gbl._globals['verbose']: print('RclTab: ', RclTab)
        # interpolate single clump opacities, emissivities and radii
        if gbl._globals['verbose']: print('interpolated AFUV:\n ', AFUVCl)
        # print 'Rcl:\n ', RclTab
        if gbl._globals['verbose']: print('ensemble no: ', iunique)
        # pause = input('interpolation done ... ')

        Lbox = gbl._globals['compound']['dint']/gbl._globals['constants']['pc']
        # pixel edge length in pc
        number_tau = copy(np.array(number)) 
        # number of clumps per mass interval
        # print 'number of clumps per mass interval (tau)', number_tau    
        # pause = input('scaling...')
        # scale number of clumps per mass interval and Lbox in a way that 
        # Ncl/Lbox^2 is kept constant to fulfill the following conditions (for binomial):
        # 1. a minimum of 1 clump in each (the highes) mass interval
        # 2. as few clumps as possible for fast computing
        # 3. interger number of clumps in highest mass interval
        #    to surpress (large) rounding errors
        # 4. p < 1 i.e. np.pi * max(RclTab)**2/Lbox**2 < 1
        if gbl._globals['verbose']: print('numberInt: ', numberInt)
        if gbl._globals['verbose']: print('number_tau: ', number_tau)
        if gbl._globals['verbose']: print('Lbox: ', Lbox)
        if gbl._globals['binomial_poisson']:
              ### use BINOMIAL distribution ###########################
          Lbox = Lbox/np.sqrt(float(number_tau[numberInt-1]))
          number_tau = number_tau/float(number_tau[numberInt-1]) 
          # print 'number of clumps, scaled to N_nM=1',
          # print number_tau

          Lbox = Lbox * np.sqrt(float(gbl._globals['statistic']['Nmin']))
          number_tau = number_tau * float(gbl._globals['statistic']['Nmin'])
          # print 'number of clumps, scaled to Nmin',
          # print number_tau
          # pause = input('..ok?..')

          if np.pi * max(RclTab)**2/Lbox**2 < 1:
            #print 'scaling finished...'
            pass
          elif np.pi * max(RclTab)**2/Lbox**2 >= 1:
            # print 'some clumps are lager then the box, scaling ...'
            scale = np.ceil(np.pi * max(RclTab)**2/Lbox**2)
            # print 'scale', scale
            Lbox = Lbox * np.sqrt(scale) 
            number_tau = number_tau * scale
        else: # use Poisson distribution ###############################
          scale = 1./float(number_tau[numberInt-1]) * 100
          # scale to high numbers of clumps, here:
          # 100 high mass clumps
          Lbox = Lbox * np.sqrt(scale)
          number_tau = np.around(number_tau * scale)
          #print('Lbox: ', Lbox)
          #print('RclTab: ', RclTab)
          if np.pi * max(RclTab) ** 2 / Lbox ** 2 >= 1:
            print('Lbox to small...increasing size')
            scale = np.ceil(np.pi * max(RclTab) ** 2 / Lbox ** 2) * 100.
            Lbox = Lbox * np.sqrt(scale)
            number_tau = number_tau * scale      
          ####################################################################
        pTab = []
        probabilityTab = [] 
        expectedValTab = [] # expected values 
        standardDeriTab = [] # related standard derivations
        lowerupperTab = [] 
        # numbers of clumps accounted for in further calculation. 
        # for computing speed.
        for i in range(numberInt):
          # print 'i', i
          pTab.append(np.pi * RclTab[int(i)]**2/Lbox**2)
          expectedValTab.append(round(number_tau[int(i)]) * pTab[int(i)])
          standardDeriTab.append(np.sqrt(round(number_tau[int(i)]) * pTab[int(i)] * \
                                (1 - pTab[int(i)])))
          lower = max(0, np.floor(expectedValTab[int(i)] - nsigma *\
                      standardDeriTab[int(i)]))
          upper = min(round(number_tau[int(i)]), np.ceil(expectedValTab[int(i)] + \
                      nsigma * standardDeriTab[int(i)]))
          lowerupperTab.append([lower, upper])  
          if  gbl._globals['binomial_poisson']:
            if (expectedValTab[int(i)] > pn_gauss and round(number_tau[int(i)]) > N_gauss):
              # use gauss
              g = pf.gauss(expectedValTab[int(i)], standardDeriTab[int(i)])
              probabilityTab.append(g.gaussfunc)
              pause = input('gauss!!...')
            else:
              # use binomial 
              b = pf.binomial(round(number_tau[int(i)]), pTab[int(i)])
              probabilityTab.append(b.binomfunc)
          else:
            if (expectedValTab[int(i)] > pn_gauss and round(number_tau[int(i)]) > N_gauss):
              # use gauss
              g = pf.gauss(expectedValTab[int(i)], standardDeriTab[int(i)])
              probabilityTab.append(g.gaussfunc)
              pause = input('gauss!!...')
            else:
              # use poisson
              po = pf.poisson(expectedValTab[int(i)])
              probabilityTab.append(po.poissonfunc)

        # print 'pTab', pTab
        if gbl._globals['verbose']: print('upperlower: ', lowerupperTab)
        # pause = input('lowerupperTab...')
        a = pf.AllCombi(lowerupperTab)
        a.createArray()
        combis = a.combinations
        # print 'combinations', a.combinations
        AFUVNoBin = []; # table: FUV extinction A against probability
        p = 0;
        # pause = input("after allcombi")
        for c in range(len(combis)):
          AFUVtot = 0
          ptot = 1
          for inter in range(numberInt):
            AFUVtot +=  AFUVCl[inter] * combis[c][inter]
            ptot *= float(probabilityTab[inter]\
                          (int(combis[c][inter])))
            if gbl._globals['verbose']: print('ptot for inter {}: {}'.format(inter, probabilityTab[inter](int(combis[c][inter]))))
            # print 'proTab', probabilityTab[inter](int(combis[c][inter]))
            # pause = input()
          p = p + ptot
          if p>1: input('Error in probability')
          AFUVNoBin.append([ptot, AFUVtot])
          # pause = input()
        if gbl._globals['verbose']: print('p = ', p)
        # pause = input()
        if gbl._globals['verbose']:
          print('to obtain exact results p needs to be close to 1. Otherwise increase nsigma')
          print('FUV p:', AFUVNoBin, '\n', np.array(AFUVNoBin)[:,0].sum(), '\n')
        AFUVEnsAv = (-np.log(sum((np.array(AFUVNoBin)[:,0])*\
                    (np.exp(-np.array(AFUVNoBin)[:,1])))))
        if gbl._globals['verbose']: print('AFUVEnsAv: ', AFUVEnsAv)
        # pause = input('A?')
    
      gbl._globals['compound']['ens'][iunique].Afuv = AFUVEnsAv
      # store ensemble averaged FUV extinction
      progress.update(1)
  
  with open(gbl.KOSMAPATH+'temp/fuvExtinction.dat', 'w') as t:
  # write ensemble averaged opacity (A)
      for i in range(number_ensembles):
          t.write(str(gbl._globals['compound']['ens'][int(i)].Afuv) + '\n')
