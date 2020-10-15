import numpy as np
import functools as ft
from operator import mul
class Binomial():
  '''class calculation binomial coefficients (choose) and function'''
  '''This is a classs taken directly from the work of Silke Andree-Labsch and Christoph Bruckmann.'''
  def __init__(self, n, p, debug=False):
    self.debug = debug
    self.n = n.reshape(1, n.size)
    self.nIndeces = np.unique(np.where(self.n>0)[0])   #this is used to isolate where in the velocity array there is a masspoint
    if self.debug:
      print(n)
      print(self.nIndeces, self.n[self.nIndeces])
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
    i = np.where(self.n[self.nIndeces]-k<k)         # for smaller intermediate values
    if self.debug: input(i)
    #iComb = i.any(1)
    k = np.full((self.n[:,0].size, k.size), k)
    if self.debug:
      input('\nCombination, n, p, i:\n{}\n{}\n{}\n{}\n'.format(k, self.n, self.p, i))
    if k[self.nIndeces][i].size: k[self.nIndeces][i] = self.n[self.nIndeces][i]-k[self.nIndeces][i]   # use (n choose k) = (n choose n-k)
    comb = np.zeros([self.n[:,0].size, k[0,:].size])
    '''
    # This is the working version I had for two masspoints.
    for i in range(self.n[:,0].size):
      if np.any(k>0):
        if self.debug:
          input('\n({}, {})\n'.format(int(self.n[:,0][i]-k[:,0]+1),int(self.n[:,0][i]+1)))
        range1 = range(int(self.n[:,0][i]-k[:,0]+1),int(self.n[:,0][i]+1))
        range2 = range(1,int(k[:,0]+1))
        range3 = range(int(self.n[:,1][i]-k[:,1]+1),int(self.n[:,1][i]+1))
        range4 = range(1,int(k[:,1]+1))
        if self.debug:
          print('\nRanges:\n', range1, '\n', range2, '\n', range3, '\n', range4)
          input()
        comb[i] = np.array([ft.reduce(mul, range1, 1)/ft.reduce(mul, range2, 1),ft.reduce(mul, range3, 1)/ft.reduce(mul, range4, 1)], dtype=np.float)
      else: comb[i] = 0
    '''
    range1 = np.zeros([1, self.n[0,:].size])
    if self.debug:
      print('\n{}\n'.format(range1))
    range2 = np.zeros(self.n.size)
    for i in range(self.n[:,0].size): #loop over velocities
      # This needs to be fixed in the future...
      if self.n[0,:].size==1:
        if self.debug:
          print('\n({}, {})\n'.format(int(self.n[:,0][i]-k[:,0][i]+1),int(self.n[:,0][i]+1)))
        range1 = range(int(self.n[:,0][i]-k[:,0][i]+1),int(self.n[:,0][i]+1))
        range2 = range(1,int(k[:,0][i]+1))
        if self.debug:
          print('\nRanges:\n1:  {}\n2:  {}\n'.format(range1, range2))
        if self.debug:
          print('\nRanges:\n{}\n{}\n'.format(range1, range2))
        try:
          comb[i] = np.array([ft.reduce(mul, range1, 1)/ft.reduce(mul, range2, 1)], dtype=np.float)
        except OverflowError:
          print('range', range1, range2)
          input()
      elif self.n[0,:].size==2:
        if self.debug:
          print('\n({}, {})\n'.format(int(self.n[:,1][i]-k[:,1][i]+1),int(self.n[:,1][i]+1)))
        range1 = range(int(self.n[:,0][i]-k[:,0][i]+1),int(self.n[:,0][i]+1))
        range2 = range(1,int(k[:,0][i]+1))
        range3 = range(int(self.n[:,1][i]-k[:,1][i]+1),int(self.n[:,1][i]+1))
        range4 = range(1,int(k[:,1][i]+1))
        if self.debug:
          print('\nRanges:\n1:  {}\n2:  {}\n3:  {}\n4:  {}\n'.format(range1, range2, range3, range4))
        comb[i] = np.array([ft.reduce(mul, range1, 1)/ft.reduce(mul, range2, 1), ft.reduce(mul, range3, 1)/ft.reduce(mul, range4, 1)], dtype=np.float)
      elif self.n[0,:].size==3:
        if self.debug:
          print('\n({}, {})\n'.format(int(self.n[:,2][i]-k[:,2]+1),int(self.n[:,2][i]+1)))
        range1 = range(int(self.n[:,0][i]-k[:,0]+1),int(self.n[:,0][i]+1))
        range2 = range(1,int(k[:,0]+1))
        range3 = range(int(self.n[:,1][i]-k[:,1]+1),int(self.n[:,1][i]+1))
        range4 = range(1,int(k[:,1]+1))
        range5 = range(int(self.n[:,2][i]-k[:,2]+1),int(self.n[:,2][i]+1))
        range6 = range(1,int(k[:,2]+1))
        if self.debug:
          print('\nRanges:\n{}\n{}\n{}\n{}\n{}\n{}\n'.format(range1, range2, range3, range4, range5, range6))
        comb[i] = np.array([ft.reduce(mul, range1, 1)/ft.reduce(mul, range2, 1), ft.reduce(mul, range3, 1)/ft.reduce(mul, range4, 1), \
                            ft.reduce(mul, range5, 1)/ft.reduce(mul, range6, 1)], dtype=np.float)
      elif self.n[0,:].size==4:
        if self.debug:
          print('\n({}, {})\n'.format(int(self.n[:,3][i]-k[:,3]+1),int(self.n[:,3][i]+1)))
        range1 = range(int(self.n[:,0][i]-k[:,0]+1),int(self.n[:,0][i]+1))
        range2 = range(1,int(k[:,0]+1))
        range3 = range(int(self.n[:,1][i]-k[:,1]+1),int(self.n[:,1][i]+1))
        range4 = range(1,int(k[:,1]+1))
        range5 = range(int(self.n[:,2][i]-k[:,2]+1),int(self.n[:,2][i]+1))
        range6 = range(1,int(k[:,2]+1))
        range7 = range(int(self.n[:,0][i]-k[:,3]+1),int(self.n[:,3][i]+1))
        range8 = range(1,int(k[:,3]+3))
        if self.debug:
          print('\nRanges:\n{}\n{}\n{}\n{}\n{}\n{}\n{}\n{}\n'.format(range1, range2, range3, range4, range5, range6, range7, range8))
        comb[i] = np.array([ft.reduce(mul, range1, 1)/ft.reduce(mul, range2, 1), ft.reduce(mul, range3, 1)/ft.reduce(mul, range4, 1), \
                            ft.reduce(mul, range5, 1)/ft.reduce(mul, range6, 1), ft.reduce(mul, range7, 1)/ft.reduce(mul, range8, 1)], dtype=np.float)
      elif self.n[0,:].size==5:
        if self.debug:
          print('\n({}, {})\n'.format(int(self.n[:,4][i]-k[:,4]+1),int(self.n[:,4][i]+1)))
        range1 = range(int(self.n[:,0][i]-k[:,0]+1),int(self.n[:,0][i]+1))
        range2 = range(1,int(k[:,0]+1))
        range3 = range(int(self.n[:,1][i]-k[:,1]+1),int(self.n[:,1][i]+1))
        range4 = range(1,int(k[:,1]+1))
        range5 = range(int(self.n[:,2][i]-k[:,2]+1),int(self.n[:,2][i]+1))
        range6 = range(1,int(k[:,2]+1))
        range7 = range(int(self.n[:,3][i]-k[:,3]+1),int(self.n[:,3][i]+1))
        range8 = range(1,int(k[:,3]+1))
        range9 = range(int(self.n[:,4][i]-k[:,4]+1),int(self.n[:,4][i]+1))
        range10 = range(1,int(k[:,4]+1))
        if self.debug:
          print('\nRanges:\n{}\n{}\n{}\n{}\n{}\n{}\n{}\n{}\n{}\n{}\n'.format(range1, range2, range3, range4, range5, range6, range7, range8, range9, range10))
        comb[i] = np.array([ft.reduce(mul, range1, 1)/ft.reduce(mul, range2, 1), ft.reduce(mul, range3, 1)/ft.reduce(mul, range4, 1), \
                            ft.reduce(mul, range5, 1)/ft.reduce(mul, range6, 1), ft.reduce(mul, range7, 1)/ft.reduce(mul, range8, 1), \
                            ft.reduce(mul, range9, 1)/ft.reduce(mul, range10, 1)], dtype=np.float)
      elif self.n[0,:].size==6:
        if self.debug:
          print('\n({}, {})\n'.format(int(self.n[:,5][i]-k[:,5]+1),int(self.n[:,5][i]+1)))
        range1 = range(int(self.n[:,0][i]-k[:,0]+1),int(self.n[:,0][i]+1))
        range2 = range(1,int(k[:,0]+1))
        range3 = range(int(self.n[:,1][i]-k[:,1]+1),int(self.n[:,1][i]+1))
        range4 = range(1,int(k[:,1]+1))
        range5 = range(int(self.n[:,2][i]-k[:,2]+1),int(self.n[:,2][i]+1))
        range6 = range(1,int(k[:,2]+1))
        range7 = range(int(self.n[:,3][i]-k[:,3]+1),int(self.n[:,3][i]+1))
        range8 = range(1,int(k[:,3]+1))
        range9 = range(int(self.n[:,4][i]-k[:,4]+1),int(self.n[:,4][i]+1))
        range10 = range(1,int(k[:,4]+1))
        range11 = range(int(self.n[:,5][i]-k[:,5]+1),int(self.n[:,5][i]+1))
        range12 = range(1,int(k[:,5]+1))
        if self.debug:
          print('\nRanges:\n{}\n{}\n{}\n{}\n{}\n{}\n{}\n{}\n{}\n{}\n{}\n{}\n'.format(range1, range2, range3, range4, range5, range6, range7, range8, range9, range10, range11, range12))
        comb[i] = np.array([ft.reduce(mul, range1, 1)/ft.reduce(mul, range2, 1), ft.reduce(mul, range3, 1)/ft.reduce(mul, range4, 1), \
                            ft.reduce(mul, range5, 1)/ft.reduce(mul, range6, 1), ft.reduce(mul, range7, 1)/ft.reduce(mul, range8, 1), \
                            ft.reduce(mul, range9, 1)/ft.reduce(mul, range10, 1), ft.reduce(mul, range11, 1)/ft.reduce(mul, range12, 1)], dtype=np.float)
      elif self.n[0,:].size==7:
        if self.debug:
          print('\n({}, {})\n'.format(int(self.n[:,6][i]-k[:,6]+1),int(self.n[:,6][i]+1)))
        range1 = range(int(self.n[:,0][i]-k[:,0]+1),int(self.n[:,0][i]+1))
        range2 = range(1,int(k[:,0]+1))
        range3 = range(int(self.n[:,1][i]-k[:,1]+1),int(self.n[:,1][i]+1))
        range4 = range(1,int(k[:,1]+1))
        range5 = range(int(self.n[:,2][i]-k[:,2]+1),int(self.n[:,2][i]+1))
        range6 = range(1,int(k[:,2]+1))
        range7 = range(int(self.n[:,3][i]-k[:,3]+1),int(self.n[:,3][i]+1))
        range8 = range(1,int(k[:,3]+1))
        range9 = range(int(self.n[:,4][i]-k[:,4]+1),int(self.n[:,4][i]+1))
        range10 = range(1,int(k[:,4]+1))
        range11 = range(int(self.n[:,5][i]-k[:,5]+1),int(self.n[:,5][i]+1))
        range12 = range(1,int(k[:,5]+1))
        range13 = range(int(self.n[:,6][i]-k[:,6]+1),int(self.n[:,6][i]+1))
        range14 = range(1,int(k[:,6]+1))
        if self.debug:
          print('\nRanges:\n{}\n{}\n{}\n{}\n{}\n{}\n{}\n{}\n{}\n{}\n{}\n{}\n{}\n{}\n'.format(range1, range2, range3, range4, range5, range6, range7, range8, range9, range10, range11, range12, range13, range14))
        comb[i] = np.array([ft.reduce(mul, range1, 1)/ft.reduce(mul, range2, 1), ft.reduce(mul, range3, 1)/ft.reduce(mul, range4, 1), \
                            ft.reduce(mul, range5, 1)/ft.reduce(mul, range6, 1), ft.reduce(mul, range7, 1)/ft.reduce(mul, range8, 1), \
                            ft.reduce(mul, range9, 1)/ft.reduce(mul, range10, 1), ft.reduce(mul, range11, 1)/ft.reduce(mul, range12, 1), \
                            ft.reduce(mul, range13, 1)/ft.reduce(mul, range14, 1)], dtype=np.float)
      else: comb[i] = 0
    
    if self.debug:
      print('\nComb:\n{}\n'.format(comb))
    return comb
  def binomfunc(self, k):
    # print 'comb', self.comb(k)
    k = np.array(k, dtype=np.int)
    #print (float(self.comb(k)) * self.p**k * (1-self.p)**(self.n-k))
    probability = self.comb(k) * self.p.T**k * (1-self.p.T)**(self.n-k)
    probability[probability==1] = 0
    if self.debug: print('\nProbability\n{}'.format(probability))
    return probability[0]