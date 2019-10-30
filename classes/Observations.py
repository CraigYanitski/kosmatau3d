class Observations():
  '''
  This class will contain the input data needed to properly simulate the PDR.
  '''
  # PRIVATE
  def __init__(self, directory):
    if directory[-1]=='/': self.__INPUTPATH = '/home/yanitski/Desktop/KOSMA-tau^3/input/' + directory
    else: self.__INPUTPATH = '/home/yanitski/Desktop/KOSMA-tau^3/input/' + directory + '/'
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
  def massProfile(self, clumpFile='mass_profile.dat', interclumpFile='mass_profile_inter.dat'):
    # Open file for the mass profile (clump and interclump) of the object
    clumpMass = np.genfromtxt(self.__INPUTPATH+file, names=['radius', 'h2 mass'])
    interclumpMass = np.genfromtxt(self.__INPUTPATH+file, names=['radius', 'h2 mass'])
    self.clumpMassGrid = sp.interp1d(clumpMass[0], clumpMass[1])                  #clump mass interpolation
    self.interclumpMassGrid = sp.interp1d(interclumpMass[0], interclumpMass[1])   #interclump mass interpolation
    return ((interclumpMass['radius'],interclumpMass['h2 mass']), (clumpMass['radius'],clumpMass['h2 mass']))
  def densityProfile(self, file='density_profile.dat'):
    # Open file for the density profile of the object
    density = np.genfromtxt(self.__INPUTPATH+file, names=['radius', 'h2 surface density'])
    self.densityGrid = sp.interp1d(density[0], density[1])      #density interpolation
    return (density['radius'],density['h2 surface density'])
  def FUVfile(self, file='galactic_fuv.dat'):
    # Open file for the FUV profile of the object
    fuv = np.genfromtxt(self.__INPUTPATH+file, names=False)#['radius', 'energy density 912', 'energy density 1350', 'energy density 1500', \
                                                      #'energy density 1650', 'energy density 2000', 'energy density 2200', 'energy density 2500', \
                                                      #'energy density 2800', 'energy density 3650'])
    return (fuv[:,0],fuv[:,1],fuv[:,2:])
  def rotationProfile(self, file='rot_milki2018_14.dat'):
    # Open file for the rotation profile of the object
    rotation = np.genfromtxt(self.__INPUTPATH+file, names=['radius', 'rotation velocity'])
    self.rotationGrid = sp.interp1d(rotation[0], rotation[1])   #rotation velocity interpolation
    return (rotation['radius'],rotation['rotation velocity'])
  def tauCenterline(self, file='tau_linecentre.dat'):
    # Open file for KOSMA-tau simulations of optical depths
    # FORMAT: n, M, UV, tau[molecules then dust]
    tau = np.genfromtxt(self.__INPUTPATH+file)
    return (tau[:,0],tau[:,1],tau[:,2],tau[:,3:])
  def tbCenterline(self, file='Tb_linecentre.dat'):
    # Open file for KOSMA-tau simulations of line intensities
    # FORMAT: n, M, UV, intensity[molecules then dust]
    tb = np.genfromtxt(self.__INPUTPATH+file)
    return (tb[:,0],tb[:,1],tb[:,2],tb[:,3:])
  def rhoMassAFUV(self, file='RhoMassAFUV.dat'):
    pma = np.genfromtxt(self.__INPUTPATH+file, names=['number density', 'mass', 'FUV extinction'])
    return (pma['number densitty'],pma['mass'],pma['FUV extinction'])
  def speciesData(self, file='frequencies.dat'):
    # Open file containing the transition frequencies of various elements
    frequencies = np.genfromtxt(self.__INPUTPATH+file, names=['number', 'species', 'transition', 'frequency'])
    return (frequencies['name'],frequencies['species'],frequencies['transition'],frequencies['frequency'])