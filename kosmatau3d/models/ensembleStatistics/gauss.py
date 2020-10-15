import numpy as np
class Gauss:
  '''This is a classs taken directly from the work of Silke Andree-Labsch and Christoph Bruckmann.'''
  def __init__(self, v0, sigma, area=1, debug=False, verbose=False):
    if verbose: print(v0, sigma)
    self.debug = debug
    self.area = np.array(area) # area below the curve (integrated curve) 
    self.v0 = np.array(v0) # peak velocity
    self.sigma = np.array(sigma) # standard deviation
    self.sigmaIndeces = self.sigma==0
    # print(self.v0,self.sigma,self.sigmaIndeces)
    return
  
  def gaussfunc(self, v, verbose=False):
    #if np.any(np.array(self.sigma)==0):
    #  return 0
    #else:
    power = np.array(-.5*((v-self.v0)/self.sigma)**2)
    probability = self.area/(np.sqrt(2*np.pi)*self.sigma)*np.exp(power)
    # print(power,probability,self.sigmaIndeces)
    if self.sigmaIndeces.ndim==1:
        probability.flatten()[self.sigmaIndeces] = 0
    else:
        probability[self.sigmaIndeces] = 0
    probability = probability.T
    if self.debug: input('\nProbability:\n{}'.format(probability))
    if verbose:
        print(v, probability)
        input()
    return probability.flatten()
