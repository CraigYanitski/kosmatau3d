import numpy as np
import functools as ft
from operator import mul
class Binomial():
  '''class calculation binomial coefficients (choose) and function'''
  '''This is a classs taken directly from the work of Silke Andree-Labsch and Christoph Bruckmann.'''
  def __init__(self, n, p):
    self.n = np.array(n, dtype=np.int)
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
    k = np.array(k, dtype=np.int)
    i = k>self.n-k          # for smaller intermediate values 
    k[i] = self.n[i]-k[i]   # use (n choose k) = (n choose n-k)
    if np.any(k>0):
      range1 = range(self.n[0]-k[0]+1,self.n[0]+1)
      range2 = range(1,k[0]+1)
      range3 = range(self.n[1]-k[1]+1,self.n[1]+1)
      range4 = range(1,k[1]+1)
      return np.array([ft.reduce(mul, range1)/ft.reduce(mul, range2),ft.reduce(mul, range3)/ft.reduce(mul, range4)], dtype=np.float)
    else: return 0
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
    return self.comb(k) * self.p**k * (1-self.p)**(self.n-k)