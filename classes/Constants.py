class Constants():
  '''
  This is a class to contain all of the constants and interpolation
  tables used throughout the program. It is defined as a class in
  order to create the tables from observational inputs.
  '''
  # PRIVATE
  def __init__(self):
  	self.__intensityInterpolation = []
  	self.__tauInterpolation = []
    self.__clumpMassUpper = 10
  	self.__clumpMassLower = 10
    self.__interclumpMassUpper = 0.1
    self.__interclumpMassLower = 0.1
    return
  
  # PUBLIC
  def interpolateI(self, point, species):
    nI,massI,uvI,I = obs.tbCenterline()
  	return griddata([nI,massI,uvI], I[:,self.__speciesNumber[species]], point, method='linear')
  def interpolateTau(self, point, species):
    nTau,massTau,uvTau,Tau = obs.tauCenterline()
    return griddata([nTau,massTau,uvTau], Tau[:,self.__speciesNumber[species]], point, method='linear')
