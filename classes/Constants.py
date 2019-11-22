class Constants(object):
  '''
  This is a class to contain all of the constants used throughout the program.
  '''
  #PRIVATE
  def __init__(self, clumpMassRange=[1,3], interclumpMassRange=[-2,-1], pixelWidth=1000, alpha=1.84, gamma=2.31, nsigma=3, nGauss=1000, pnGauss=5, clumpDispersion=1.67/2.3548, ensembleDispersion=10./2.3548):
    self.clumpMassLimits = clumpMassRange
    self.interclumpMassLimits = interclumpMassRange
    self.alpha = alpha
    self.gamma = gamma
    self.nSigma = nsigma
    self.nGauss = nGauss
    self.pnGauss = pnGauss
    self.pixelWidth = pixelWidth  #in pc
    self.clumpDispersion = clumpDispersion
    self.ensembleDispersion = ensembleDispersion
    self.massH = 1.008*1.6605*10**-24  #in [g]
    self.massSolar = 1.998892*10**33  #in [g]
    self.c = 2.998*10**10  #in [cm/s]
    self.kB = 1.3806*10**-23  #in [J/K]
    self.pc = 3.0856776*10**18  #in [cm]
    # Grid boundaries
    self.densityLimits = [3, 7]
    self.massLimits = [-3, 3]
    self.uvLimits = [0, 6]
    # UV adjustment
    self.normUV = 2.89433*10**39
    self.globalUV = 10
    return
