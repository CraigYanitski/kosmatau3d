import math
import numpy as np
import scipy.special
class Poisson():
  '''This is a classs taken directly from the work of Silke Andree-Labsch and Christoph Bruckmann.'''
  def __init__(self, la):
    self.la = la # lambda = n * p
    return
  
  def poissonfunc(self, k, verbose=False):
    # print 'comb', self.comb(k)
    factorial = scipy.special.factorial(k)
    if debug:
      print('\nself.la, k, factorial:\n', self.la, k, factorial)
      input()
    return (self.la**k/factorial*np.e**(-self.la))