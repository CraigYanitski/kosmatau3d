import numpy as np
import functools as ft
from operator import mul
class Binomial():
  '''class calculation binomial coefficients (choose) and function'''
  '''This is a classs taken directly from the work of Silke Andree-Labsch and Christoph Bruckmann.'''
  def __init__(self, n, p, debug=False):
    self.debug = debug
    self.n = np.array(n, dtype=np.int).T
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
    i = (k>self.n-k).any(0)          # for smaller intermediate values
    #iComb = i.any(1)
    if self.debug:
      print('\nCombination, n, p, i:\n', k, '\n', self.n, '\n', self.p, '\n', i)
      input()
    if k[:,i].size: k[:,i] = self.n[:,i]-k[:,i]   # use (n choose k) = (n choose n-k)
    comb = np.zeros([self.n[:,0].size, k.size])
    for i in range(self.n[:,0].size):
      if np.any(k>0):
        range1 = range(int(self.n[:,0][i]-k[:,0]+1),int(self.n[:,0][i]+1))
        range2 = range(1,int(k[:,0]+1))
        range3 = range(int(self.n[:,1][i]-k[:,1]+1),int(self.n[:,1][i]+1))
        range4 = range(1,int(k[:,1]+1))
        if self.debug:
          print('\nRanges:\n', range1, '\n', range2, '\n', range3, '\n', range4)
          input()
        comb[i] = np.array([ft.reduce(mul, range1, 1)/ft.reduce(mul, range2, 1),ft.reduce(mul, range3, 1)/ft.reduce(mul, range4, 1)], dtype=np.float)
      else: comb[i] = 0
    if self.debug:
      print('\nComb:\n', comb)
      input()
    return comb
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
    probability = self.comb(k) * self.p**k * (1-self.p)**(self.n-k)
    if self.debug: print('\nProbability\n', probability)
    return np.array([probability.mean(0)])