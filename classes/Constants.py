class Constants():
  '''
  This is a class to contain all of the constants and interpolation
  functions used throughout the program. It is defined as a class in
  order to create the functions from the inputs in the Observations()
  class.

  The method of interpolation is passed as an argument when initialising
  this class. The acceptabled values are 'linear' and 'radial'.
  The default method is 'linear'.
  '''
  # PRIVATE
  def __init__(self, interpolation='linear', clumpMassRange=[10,10], interclumpMassRange=[0.1,0.1], alpha=1.84, gamma=2.31):
    self.__interpolation = interpolation
  	self.intensityInterpolation = []
  	self.tauInterpolation = []
    self.clumpMassUpper = clumpMassRange[0]
  	self.clumpMassLower = clumpMassRange[1]
    self.interclumpMassUpper = interclumpMassRange[0]
    self.interclumpMassLower = interclumpMassRange[1]
    self.alpha = alpha
    self.gamma = gamma
    self.massH = 1.008*1.6605*10**-27  #in [kg]
    self.massSol = 1.998892*10**30  #in [kg]
    self.c = 2.998*10**8  #in [m/s]
    self.kB = 1.3806*10**-23  #in [J/K]
    self.pc = 3.0856776*10**16  #in [m]
    return
  
  # PUBLIC
  def calculateInterpolation(self, species):
    nI,massI,uvI,I = obs.tbCenterline()
    nTau,massTau,uvTau,Tau = obs.tauCenterline()
    if self.__interplation=='linear':
      for index in species.indeces():
        rInterpI = sp.interpolate.LinearNDInterpolation(nI, massI, uvI, I[index])
        rInterpTau = sp.interpolate.LinearNDInterpolation(nTau, massTau, uvTau, Tau[index])
        self.__intensityInterpolation.append()
    elif self.__interplation=='radial':
      for index in species.indeces():
        rInterpI = sp.interpolate.Rbf(nI, massI, uvI, I[index])
        rInterpTau = sp.interpolate.Rbf(nTau, massTau, uvTau, Tau[index])
        self.__intensityInterpolation.append()
    else: sys.exit('<<ERROR>>: There is no such method as {} to interpolate the KOSMA-tau grid.\n\n \
                   Exitting...\n\n'.format(self.__interpolation))
  	return
  def interpolateTau(self, points, species):
    if len(species)>1:
      intensity = []
      for i in species.number:
        intensity.append(self.__intensityInterpolation(i))
    elif len(species):
      intensity = self.__intensityInterpolation(points)
    else: sys.exit('<<ERROR>>: a zero-length species array has been used to interpolate the intensity.\n\n')
    return np.array(intensity)
  def interpolateTau(self, points, species):
    if len(species)>1:
      tau = []
      for i in range(species):
        tau.append(self.__tauInterpolation(i))
    elif len(species):
      tau = self.__tauInterpolation(points)
    else: sys.exit('<<ERROR>>: a zero-length species array has been used to interpolate the optical depth.\n\n')
    return np.array(tau)
