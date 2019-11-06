class Poisson():
  '''This is a classs taken directly from the work of Silke Andree-Labsch and Christoph Bruckmann.'''
  def __init__(self, la):
    self.la = float(la) # lambda = n * p
    return
  
  def poissonfunc(self, k):
    # print 'comb', self.comb(k)
    return self.la**k /factorial(k) * np.e**(-self.la)