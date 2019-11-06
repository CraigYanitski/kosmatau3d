class Gauss:
  '''This is a classs taken directly from the work of Silke Andree-Labsch and Christoph Bruckmann.'''
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