class Ensemble(object):
  '''
  This class owes itself largely to the work  done by Silke Andree-Labsch and
  Christoph Bruckmann. It contains a number of combinations of fractal mass
  configurations, and therefore contains the ensemble-averaged intensity and
  optical depth.
  It can be either a clump or an interclump ensemble.
  '''
  # PRIVATE
  def __init__(self, clumpType):
    self.__clumpType = clumpType     #type of mass in ensemble (clump or interclump medium)
    self.__combinations = []    #list of combinations
    self.__intensity = 0        #ensemble-averaged intensity for each velocity
    self.__opticalDepth = 0     #ensemble-averaged optical depth for each velocity
    self.__FUV = 0
    self.__mass = 0     #ensemble-averaged mass
    self.__radius = 0   #ensemble-averaged radius
    if clumpType=='clump':
      self.__massLimitUpper = constants.getUpperClumpMass()   #maximum mass in ensemble
      self.__massLimitLower = constants.getLowerClumpMass()   #minimum mass in ensemble
    elif clumpType=='interclump':
      self.__massLimitUpper = constants.getUpperInterclumpMass()   #maximum mass in ensemble
      self.__massLimitLower = constants.getLowerInterclumpMass()   #minimum mass in ensemble
    self.__masspoints = []
    self.__masspointNumberRange = []
    self.__masspointRadii = []
    self.__interpolationPoints = []
    return
  def __setMass(self, mass):
    self.__mass = mass
    return

  # PUBLIC
  def setMass(self, mass):
    self.__setMass(mass)
    return
  def calculateMasspoints(self):
    self.__masspoints = np.arange(self.__massLimitLower, self.__massLimitUpper, 10)
    return
  def calculateSurfaceDensity(self):
    return np.log10(self.__surfaceDensity)
  def getInterpolationPoints(self):

  def getMasspoints(self):
    return self.__masspoints
  def calculateRadii(self):
    self.__surfaceDensity = 10**(self.__masspoints/10.*(1-3/gbl._globals['constants']['gamma']))*sum(10**(self.__masspoints/10.**(1+3/gbl._globals['constants']['gamma']-gbl._globals['constants']['alpha']))) / \
                            sum(10**(self.__masspoints/10.0*(2-gbl._globals['constants']['alpha'])))*rho_ens/1.91
    if (self.__surfaceDensity>10**(gbl._globals['compound']['ranges']['max_arr'][0]/10)).any() or (self.__surfaceDensity<10**(gbl._globals['compound']['ranges']['min_arr'][0]/10)).any(): sys.exit('WARNING: surface density outside of KOSMA-tau grid!\n\n')
    self.__interpolationPoints = np.stack((self.__surfaceDensity, self.__masspoints, np.full(self.__masspoints.size, self.__FUV)))
    self.__masspointRadii = ((3/(4.0*np.pi)*(10**(self.__masspoints/10.0)*gbl._globals['constants']['M_sol'])/(10**(self.__interpolationPoints/10.0)*gbl._globals['constants']['M_H']*1.91))**(1./3.))/gbl._globals['constants']['pc']
    return
  def getCombinations(self):
    for i,mass in self.getMasspoints:
      number = ((Mens*(10**(M_log/10.0))**(1-gbl._globals['constants']['alpha']))/(sum((10**(M_log2/10.0))**(2-gbl._globals['constants']['alpha']) for M_log2 in MLogArray)))
      pTab = (np.pi*self.radius**2/constants.pixelWidth**2)
      expectedValTab = (number_v*pTab)
      standardDeriTab = (np.sqrt(number_v*pTab*(1-pTab)))
      lower = max(0, np.floor(expectedValTab-gbl._globals['constants']['nsigma']*standardDeriTab))
      upper = min(number_v, np.ceil(expectedValTab+gbl._globals['constants']['nsigma']*standardDeriTab))
      self.__masspointNumberRange.append([lower, upper])
    combinations = self.calculateCombinations(self.__masspointNumberRange)
    for combination in combinations: self.__combinations.append(Combination(combination))
    return
  def calculate(self):
    self.calculateMasspoints()
    self.getCombinations()
    for combination in self.__combinations:
      combination.calculateEmission()
      result = combination.getScaledCombinationEmission()
      self.__intensity += result[0]*result[1]
      self.__opticalDepth += result[0]*np.exp(-result[2])
      self.__FUV += result[0]*np.exp(-result[3])
    return
  def getEnsembleEmission(self):
    return (self.__intensity,self.__opticalDepth,self.__FUV)

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