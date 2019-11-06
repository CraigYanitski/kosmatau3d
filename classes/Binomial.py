class Binomial():
  '''class calculation binomial coefficients (choose) and function'''
  '''This is a classs taken directly from the work of Silke Andree-Labsch and Christoph Bruckmann.'''
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