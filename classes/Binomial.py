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
    i = (k>self.n-k).any(1)          # for smaller intermediate values
    #iComb = i.any(1)
    #print(k[iComb], self.n[i], i)
    if k[i].size: k[i] = self.n[i]-k[i]   # use (n choose k) = (n choose n-k)
    probability = np.zeros([self.n[0].size, 2])
    for i in range(self.n[0].size):
      if np.any(k>0):
        range1 = range(int(self.n[0][i]-k[0]+1),int(self.n[0][i]+1))
        range2 = range(1,int(k[0]+1))
        range3 = range(int(self.n[1][i]-k[1]+1),int(self.n[1][i]+1))
        range4 = range(1,int(k[1]+1))
        #print(range1, range2, range3, range4)
        probability[i] = np.array([ft.reduce(mul, range1, 1)/ft.reduce(mul, range2, 1),ft.reduce(mul, range3, 1)/ft.reduce(mul, range4, 1)], dtype=np.float)
      else: probability[i] = 0
    #print(probability)
    return probability
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
    k = np.array(k, dtype=np.int)
    #print (float(self.comb(k)) * self.p**k * (1-self.p)**(self.n-k))
    return self.comb(k).T * self.p.T**k * (1-self.p.T)**(self.n-k)