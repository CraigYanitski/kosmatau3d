import numpy as np
from math import factorial
from copy import copy
from operator import mul # multiplication
from functools import reduce

class poisson():
  def __init__(self, la):
    self.la = float(la) # lambda = n * p
    return
  
  def poissonfunc(self, k):
    # print 'comb', self.comb(k)
    return self.la**k /factorial(k) * np.e**(-self.la)

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
      return self.area/(np.sqrt(2*np.pi)*self.sigma)*\
             np.exp(-.5 * ( (v - self.v0) / self.sigma)**2 )

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
