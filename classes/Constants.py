class Constants():
  '''
  This is a class to contain all of the constants used throughout the program.
  '''
  #PRIVATE
  def __init__(self, clumpMassRange=[10,10], interclumpMassRange=[0.1,0.1], pixelWidth=1000, alpha=1.84, gamma=2.31, nsigma=3, nGauss=1000, pnGauss=5):
    self.clumpMassLimits = clumpMassRange
    self.interclumpMassLimits = interclumpMassRange
    self.alpha = alpha
    self.gamma = gamma
    self.nSigma = nsigma
    self.nGauss = nGauss
    self.pnGauss = pnGauss
    self.pixelWidth = pixelWidth
    self.massH = 1.008*1.6605*10**-27  #in [kg]
    self.massSolar = 1.998892*10**30  #in [kg]
    self.c = 2.998*10**8  #in [m/s]
    self.kB = 1.3806*10**-23  #in [J/K]
    self.pc = 3.0856776*10**16  #in [m]
    # Grid boundaries
    self.densityLimits = [3, 7]
    self.massLimits = [-3, 3]
    self.uvLimits = [0, 6]
    return
