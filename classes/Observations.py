import numpy as np
class Observations(object):
  '''
  This class will contain the input data needed to properly simulate the PDR. All of the information specific
  to the object being simulated should be in their own folder in INPUTPATH. The KOSMA-tau grid data is located
  in the folder 'grid'.
  '''
  # PRIVATE
  def __init__(self, directory):
    self.__INPUTPATH = '/home/yanitski/Desktop/KOSMA-tau^3/input/'
    self.__GRIDPATH = '/home/yanitski/Desktop/KOSMA-tau^3/input/grid/'
    if directory[-1]=='/': self.__directory = directory
    else: self.__directory = directory + '/'
    self.__setup()
    return
  def __setup(self):
    self.__massProfile = self.massProfile()
    self.__densityProfile =  self.densityProfile()
    self.__FUVfield = self.FUVfield()
    self.__rotationProfile = self.rotationProfile()
    self.__tauCenterline = self.tauCenterline()
    self.__tbCenterline = self.tbCenterline()
    self.__rhoMassAFUV = self.rhoMassAFUV()
    self.__speciesData = self.speciesData()

  # PUBLIC
  def clumpMassProfile(self, file='mass_profile.dat'):
    # Open file for the mass profile (clump) of the object
    clumpMass = np.genfromtxt(self.__INPUTPATH+self.__directory+file, names=['radius', 'h2 mass'])
    return (clumpMass['radius'],clumpMass['h2 mass'])
  def interclumpMassProfile(self, file='mass_profile_inter.dat'):
    # Open file for the mass profile (interclump) of the object
    interclumpMass = np.genfromtxt(self.__INPUTPATH+self.__directory+file, names=['radius', 'h2 mass'])
    return (interclumpMass['radius'],interclumpMass['h2 mass'])
  def densityProfile(self, file='density_profile.dat'):
    # Open file for the density profile of the object
    density = np.genfromtxt(self.__INPUTPATH+self.__directory+file, names=['radius', 'h2 surface density'])
    self.densityGrid = sp.interp1d(density[0], density[1])      #density interpolation
    return (density['radius'],density['h2 surface density'])
  def FUVfield(self, file='galactic_fuv.dat'):
    # Open file for the FUV profile of the object
    fuv = np.genfromtxt(self.__INPUTPATH+self.__directory+file, names=False)#['radius', 'energy density 912', 'energy density 1350', 'energy density 1500', \
                                                      #'energy density 1650', 'energy density 2000', 'energy density 2200', 'energy density 2500', \
                                                      #'energy density 2800', 'energy density 3650'])
    return (fuv[:,0],fuv[:,1],fuv[:,2:])
  def rotationProfile(self, file='rot_milki2018_14.dat'):
    # Open file for the rotation profile of the object
    rotation = np.genfromtxt(self.__INPUTPATH+self.__directory+file, names=['radius', 'rotation velocity'])
    return (rotation['radius'],rotation['rotation velocity'])
  def tauCenterline(self, file='tau_linecentre.dat'):
    # Open file for KOSMA-tau simulations of optical depths
    # FORMAT: n, M, UV, tau[molecules then dust]
    tau = np.genfromtxt(self.__GRIDPATH+file)
    return (tau[:,0],tau[:,1],tau[:,2],tau[:,3:])
  def tbCenterline(self, file='Tb_linecentre.dat'):
    # Open file for KOSMA-tau simulations of line intensities
    # FORMAT: n, M, UV, intensity[molecules then dust]
    tb = np.genfromtxt(self.__GRIDPATH+file)
    return (tb[:,0],tb[:,1],tb[:,2],tb[:,3:])
  def rhoMassAFUV(self, file='RhoMassAFUV.dat'):
    pma = np.genfromtxt(self.__GRIDPATH+file, names=['number density', 'mass', 'FUV extinction'])
    return (pma['number density'],pma['mass'],pma['FUV extinction'])
  def speciesData(self, file='frequencies.dat'):
    # Open file containing the transition frequencies of various elements
    frequencies = np.genfromtxt(self.__GRIDPATH+file, names=['number', 'species', 'transition', 'frequency'])
    return (frequencies['number'],frequencies['species'],frequencies['transition'],frequencies['frequency'])