import os
import inspect
import numpy as np
class Observations(object):
  '''
  This class will contain the input data needed to properly simulate the PDR. All of the information specific
  to the object being simulated should be in their own folder in INPUTPATH. The KOSMA-tau grid data is located
  in the folder 'grid'.
  '''
  # PRIVATE
  def __init__(self, scale, directory='MilkyWay'):
    self.__scale = scale
    filename = inspect.getframeinfo(inspect.currentframe()).filename
    self.__KOSMAPATH = os.path.abspath(os.path.dirname(filename)+'/../')
    self.__INPUTPATH = self.__KOSMAPATH + '/input/'#'/home/craig/Desktop/Köln/kosma-tau^3-develop/kosma-tau-3d/input/'#
    self.__GRIDPATH = self.__KOSMAPATH + '/grid/'#'/home/craig/Desktop/Köln/kosma-tau^3-develop/kosma-tau-3d/grid/'#
    if directory[-1]=='/': self.__directory = directory
    else: self.__directory = directory + '/'
    self.__initialise()
    return
  def __clumpMassProfile(self, file='mass_profile.dat'):
    # Open file for the mass profile (clump) of the object (Msol/pc**2)
    clumpMass = np.genfromtxt(self.__INPUTPATH+self.__directory+file, names=['radius', 'h2_mass'])
    return (clumpMass['radius']*1000,clumpMass['h2_mass']*self.__scale**2)
  def __interclumpMassProfile(self, file='mass_profile_inter.dat'):
    # Open file for the mass profile (interclump) of the object (Msol/pc**2)
    interclumpMass = np.genfromtxt(self.__INPUTPATH+self.__directory+file, names=['radius', 'h2_mass'])
    return (interclumpMass['radius']*1000,interclumpMass['h2_mass']*self.__scale**2)
  def __densityProfile(self, file='densities_clouds.dat'):
    # Open file for the number density profile of the object (n/cm**3)
    density = np.genfromtxt(self.__INPUTPATH+self.__directory+file, names=['radius', 'h2_density'])
    return ((density['radius']*1000),density['h2_density'])
  def __FUVfield(self, file='galactic_FUV.dat'):
    '''Open file for the FUV profile of the object
       'radius', 'energy density 912', 'energy density 1350', 'energy density 1500', 'energy density 1650', 'energy density 2000', \
       'energy density 2200', 'energy density 2500', 'energy density 2800', 'energy density 3650'])'''
    fuv = np.loadtxt(self.__INPUTPATH+self.__directory+file)
    #r = np.arange(0,20000,50)
    #fuv = 50/4/np.pi/r**2
    return (fuv[:,0],fuv[:,1],fuv[:,2:])
    #return (r,fuv)
  def __rotationProfile(self, file='rot_milki2018_14.dat'):
    # Open file for the rotation profile of the object
    rotation = np.genfromtxt(self.__INPUTPATH+self.__directory+file)
    return (rotation[:,0]*1000., rotation[:,1:])
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
    pma = np.genfromtxt(self.__GRIDPATH+file)
    return (pma[:,:2],pma[:,2])
  def __speciesData(self, file='frequencies.dat'):
    # Open file containing the transition frequencies of various elements
    frequencies = np.genfromtxt(self.__GRIDPATH+file, names=['number', 'species', 'transition', 'frequency'], dtype="i8,U8,i8,f8", delimiter=',')
    return (frequencies['number'],frequencies['species'],frequencies['transition'],frequencies['frequency'])
  def __eTildeReal(self, file='Ereal.dat'):
    eReal = np.genfromtxt(self.__GRIDPATH+file, names=['x', 'Ereal'])
    return (eReal['x'],eReal['Ereal'])
  def __eTildeImaginary(self, file='Eimag.dat'):
    eImaginary = np.genfromtxt(self.__GRIDPATH+file, names=['x', 'Eimaginary'])
    return (eImaginary['x'],eImaginary['Eimaginary'])
  def __str__(self):
    if not len(speciesData):
      return 'There are no available transitions yet.'
    else:
      printout = 'Available transitions:'
      i = np.isfinite(self.speciesData[0])&np.isfinite(self.speciesData[1])&np.isfinite(self.speciesData[2])&np.isfinite(self.speciesData[3])
      transitions = self.speciesData[2]
      elements = self.speciesData[1]
      for i in range(len(elements)):
        printout += '\n  ->{} {}'.format(element[i], transitions[i])
      return printout

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
    self.eTildeReal = self.__eTildeReal()
    self.eTildeImaginary = self.__eTildeImaginary()
    return