import numpy as np
class Observations(object):
  '''
  This class will contain the input data needed to properly simulate the PDR. All of the information specific
  to the object being simulated should be in their own folder in INPUTPATH. The KOSMA-tau grid data is located
  in the folder 'grid'.
  '''
  # PRIVATE
  def __init__(self, directory='MilkyWay'):
    self.__INPUTPATH = '/home/yanitski/Desktop/KOSMA-tau^3/input/'
    self.__GRIDPATH = '/home/yanitski/Desktop/KOSMA-tau^3/grid/'
    if directory[-1]=='/': self.__directory = directory
    else: self.__directory = directory + '/'
    # self.clumpMassProfile = None
    # self.interclumpMassProfile = None
    # self.densityProfile =  None
    # self.FUVfield = None
    # self.rotationProfile = None
    # self.tauCenterline = None
    # self.tbCenterline = None
    # self.rhoMassAFUV = None
    # self.speciesData = None
    self.__initialise()
    return
  def __clumpMassProfile(self, file='mass_profile.dat'):
    # Open file for the mass profile (clump) of the object
    clumpMass = np.genfromtxt(self.__INPUTPATH+self.__directory+file, names=['radius', 'h2_mass'])
    return (clumpMass['radius'],clumpMass['h2_mass'])
  def __interclumpMassProfile(self, file='mass_profile_inter.dat'):
    # Open file for the mass profile (interclump) of the object
    interclumpMass = np.genfromtxt(self.__INPUTPATH+self.__directory+file, names=['radius', 'h2_mass'])
    return (interclumpMass['radius'],interclumpMass['h2_mass'])
  def __densityProfile(self, file='densities_clouds.dat'):
    # Open file for the density profile of the object
    density = np.genfromtxt(self.__INPUTPATH+self.__directory+file, names=['radius', 'h2_surface_density'])
    return (density['radius'],density['h2_surface_density'])
  def __FUVfield(self, file='galactic_FUV.dat'):
    # Open file for the FUV profile of the object
    fuv = np.loadtxt(self.__INPUTPATH+self.__directory+file)#, names=False)#['radius', 'energy density 912', 'energy density 1350', 'energy density 1500', \
                                                      #'energy density 1650', 'energy density 2000', 'energy density 2200', 'energy density 2500', \
                                                      #'energy density 2800', 'energy density 3650'])
    return (fuv[:,0],fuv[:,1],fuv[:,2:])
  def __rotationProfile(self, file='rot_milki2018_14.dat'):
    # Open file for the rotation profile of the object
    rotation = np.genfromtxt(self.__INPUTPATH+self.__directory+file, names=['radius', 'rotation_velocity'])
    return (rotation['radius'],rotation['rotation_velocity'])
  def __tauCenterline(self, file='tau_linecentre.dat'):
    # Open file for KOSMA-tau simulations of optical depths
    # FORMAT: n, M, UV, tau[molecules then dust]
    tau = np.genfromtxt(self.__GRIDPATH+file)
    return (tau[:,:3],tau[:,3:])
  def __tbCenterline(self, file='Tb_linecentre.dat'):
    # Open file for KOSMA-tau simulations of line intensities
    # FORMAT: n, M, UV, intensity[molecules then dust]
    tb = np.genfromtxt(self.__GRIDPATH+file)
    return (tb[:,:3],tb[:,3:])
  def __rhoMassAFUV(self, file='RhoMassAFUV.dat'):
    pma = np.genfromtxt(self.__GRIDPATH+file, names=['number_density', 'mass', 'FUV_extinction'])
    return (pma['number_density'],pma['mass'],pma['FUV_extinction'])
  def __speciesData(self, file='frequencies.dat'):
    # Open file containing the transition frequencies of various elements
    frequencies = np.genfromtxt(self.__GRIDPATH+file, names=['number', 'species', 'transition', 'frequency'], dtype="i8,U8,i8,f8", delimiter=',')
    return (frequencies['number'],frequencies['species'],frequencies['transition'],frequencies['frequency'])

  # PUBLIC
  def __initialise(self):
    self.clumpMassProfile = self.__clumpMassProfile()
    self.interclumpMassProfile = self.__interclumpMassProfile()
    self.densityProfile =  self.__densityProfile()
    self.FUVfield = self.__FUVfield()
    self.rotationProfile = self.__rotationProfile()
    self.tauCenterline = self.__tauCenterline()
    self.tbCenterline = self.__tbCenterline()
    self.rhoMassAFUV = self.__rhoMassAFUV()
    self.speciesData = self.__speciesData()
    return