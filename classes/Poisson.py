import math
import numpy as np
class Poisson():
  '''This is a classs taken directly from the work of Silke Andree-Labsch and Christoph Bruckmann.'''
  def __init__(self, la):
    self.la = la # lambda = n * p
    return
  
  def poissonfunc(self, k):
    # print 'comb', self.comb(k)
    return self.la**k /math.factorial(k) * np.e**(-self.la)