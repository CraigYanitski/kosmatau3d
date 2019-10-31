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
    self.__combinationObjects = []    #list of instances of Combination()
    self.__probability = []
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
    self.__deltaNji = []
    return
  def __setMass(self, mass):
    '''Set the mass.'''
    self.__mass = mass
    return

  # PUBLIC
  def setMass(self, mass):
    '''Set the mass.'''
    self.__setMass(mass)
    return
  def calculateMasspoints(self):
    '''This is a function to get the clump masses in this ensemble. It will soon be depreciated
       as I will change the reading of the KOSMA-tau files to convert from the fortran 10*log(values).'''
    self.__masspoints = np.arange(self.__massLimitLower, self.__massLimitUpper, 10)
    return
  def calculateMasspointDensity(self):
    '''I am not entirely certain this is needed...'''
    return np.log10(self.__masspointDensity)
  def getInterpolationPoints(self):
    '''This might be removed. I'm not certain the interpolation points need to be accessed out-of-class.'''
    return
  def getMasspoints(self):
    '''This is pretty self-explanatory. Return the list of clump masses when called...'''
    return self.__masspoints
  def calculateRadii(self):
    '''This function calculates the interpolation points necessary for reading the KOSMA-tau files.'''
    self.__masspointDensity = self.__masspoints**(1-3/gbl._globals['constants']['gamma'])*sum(self.__masspoints**(1+3/gbl._globals['constants']['gamma']-gbl._globals['constants']['alpha'])) / \
                            sum(self.__masspoints**(2-gbl._globals['constants']['alpha']))*self.__density/1.91
    if (self.__masspointDensity>gbl._globals['compound']['ranges']['max_arr'][0]).any() or (self.__masspointDensity<gbl._globals['compound']['ranges']['min_arr'][0]).any(): sys.exit('WARNING: surface density outside of KOSMA-tau grid!\n\n')
    self.__interpolationPoints = np.stack((self.__masspointDensity, self.__masspoints, np.full(self.__masspoints.size, self.__FUV)))
    self.__masspointRadii = ((3/(4.0*np.pi)*(self.__masspoints*gbl._globals['constants']['M_sol'])/(self.__interpolationPoints[0]*gbl._globals['constants']['M_H']*1.91))**(1./3.))/gbl._globals['constants']['pc']
    return
  def calculateCombinations(self):
    '''This function calculates all of the different combinations of clump masses that may be in a line-of-sight.
      It is basically the essence of the probabilistic approach used to create the superposition of clumps.'''
    list = []
    for line in self.__masspointNumberRange:
      l = np.arange(line[0], line[1] + 1)
      list.append(l)
    self.__combinations = np.zeros([self.numbercombi, self.numberintervals]) 
    for j in range(np.array(self.__masspointNumberRange).prod(1).sum()): # each combination
      for i in range(len(self.__masspointNumberRange)): # each column 
        prod = 1
        while a < (len(self.numberintervals) - 1):
          prod = prod * self.list[a+1].size
          a = a + 1
        self.__combinations[j][i] = int(j/(prod)) % self.list[i].size + self.list[i][0]
    return
  def getCombinations(self):
    '''This function removes all of the unnecessary degenerate looping during this calculation.
       Of course it is possible because of the wonders of nump.ndarray(). . .'''
    self.__deltaNji = (self.__mass*(self.__masspoints)**(1-gbl._globals['constants']['alpha'])) / sum((self.__masspoints)**(2-gbl._globals['constants']['alpha']))
    pTab = (np.pi*self.__masspointRadii**2/constants.pixelWidth**2)
    expectedValTab = (self.__deltaNji*pTab)
    standardDeriTab = (np.sqrt(self.__deltaNji*pTab*(1-pTab)))
    lower = np.maximum(np.zeros(np.size(self.__masspoints)), np.floor(expectedValTab-gbl._globals['constants']['nsigma']*standardDeriTab))
    upper = np.minimum(self.__deltaNji, np.ceil(expectedValTab+gbl._globals['constants']['nsigma']*standardDeriTab))
    self.__masspointNumberRange = np.array([lower, upper]).T
    self.calculateCombinations()
    if self.flagBinomialPoisson:
      if np.any(expectedValTab>gbl._globals['constants']['pn_gauss'] and number>gbl._globals['constants']['N_gauss']):
        # use gauss!
        g = Gauss(expectedValTab[ma], standardDeriTab[ma])
        self.__probability.append(g.gaussfunc)
        pause = input('gauss!!...')
      else:
        # use binomial 
        b = Binomial(number_v[ma], pTab[ma]) # n and p for binominal 
        self.__probability.append(b.binomfunc)
    else:
      if np.any(expectedValTab>gbl._globals['constants']['pn_gauss'] and number>gbl._globals['constants']['N_gauss']):
        # use gauss
        g = Gauss(expectedValTab[ma], standardDeriTab[ma])
        self.__probability.append(g.gaussfunc)
        pause = input('gauss!!...')
      else:
        # use poisson
        po = Poisson(expectedValTab[ma])
        self.__probability.append(po.poissonfunc)
    return
  def calculate(self):
    '''Maybe <<PARALLELISE>> this??

       This is a function to cycle through the Combination instances to create a large numpy.ndarray,
       which is used to calculate the final sums needed for the voxel.'''
    self.calculateMasspoints()
    self.getCombinations()
    combinations = []
    result = []
    for combination in self.__combinations:
      self.__combinationObjects.append(Combination(combination, self.__masspoints, probability=self.__probability))
      self.__combinationObjects[-1].calculateEmission()
      result.append(self.__combinationObjects[-1].getScaledCombinationEmission())
    result = np.array(result)
    self.__intensity = result.sum(3)[0]
    self.__opticalDepth = result.sum(3)[1]
    self.__FUV = result.sum(3)[2]
    return
  def getEnsembleEmission(self):
    '''This returns the ensemble emission...nothing more.'''
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

class poisson():
  def __init__(self, la):
    self.la = float(la) # lambda = n * p
    return
  
  def poissonfunc(self, k):
    # print 'comb', self.comb(k)
    return self.la**k /factorial(k) * np.e**(-self.la)
