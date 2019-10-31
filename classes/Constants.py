class Constants():
  '''
  This is a class to contain all of the constants used throughout the program.
  '''
  # PRIVATE
  def __init__(self, interpolation='linear', clumpMassRange=[10,10], interclumpMassRange=[0.1,0.1], alpha=1.84, gamma=2.31):
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
