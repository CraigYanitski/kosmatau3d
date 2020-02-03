import numpy as np
import constants
import observations

def clumpMassProfile(file='mass_profile.dat'):
  # Open file for the mass profile (clump) of the object (Msol/pc**2)
  clumpMass = np.genfromtxt(constants.INPUTPATH+constants.directory+file, names=['radius', 'h2_mass'])
  return (clumpMass['radius']*1000,clumpMass['h2_mass']*scale**2)

def interclumpMassProfile(file='mass_profile_inter.dat'):
  # Open file for the mass profile (interclump) of the object (Msol/pc**2)
  interclumpMass = np.genfromtxt(constants.INPUTPATH+constants.directory+file, names=['radius', 'h2_mass'])
  return (interclumpMass['radius']*1000,interclumpMass['h2_mass']*scale**2)

def densityProfile(file='densities_clouds.dat'):
  # Open file for the number density profile of the object (n/cm**3)
  density = np.genfromtxt(constants.INPUTPATH+constants.directory+file, names=['radius', 'h2_density'])
  return ((density['radius']*1000),density['h2_density'])

def FUVfield(file='galactic_FUV.dat'):
  '''Open file for the FUV profile of the object
     'radius', 'energy density 912', 'energy density 1350', 'energy density 1500', 'energy density 1650', 'energy density 2000', \
     'energy density 2200', 'energy density 2500', 'energy density 2800', 'energy density 3650'])'''
  fuv = np.loadtxt(constants.INPUTPATH+constants.directory+file)
  #r = np.arange(0,20000,50)
  #fuv = 50/4/np.pi/r**2
  return (fuv[:,0],fuv[:,1],fuv[:,2:])
  #return (r,fuv)

def rotationProfile(file='rot_milki2018_14.dat'):
  # Open file for the rotation profile of the object
  rotation = np.genfromtxt(constants.INPUTPATH+constants.directory+file)
  return (rotation[:,0]*1000., rotation[:,1:])

def tauCenterline(file='tau_linecentre.dat'):
  # Open file for KOSMA-tau simulations of optical depths
  # FORMAT: n, M, UV, tau[molecules then dust]
  tau = np.genfromtxt(constants.GRIDPATH+file)
  return (tau[:,:3],tau[:,3:])

def tbCenterline(file='Tb_linecentre.dat'):
  # Open file for KOSMA-tau simulations of line intensities
  # FORMAT: n, M, UV, intensity[molecules then dust]
  tb = np.genfromtxt(constants.GRIDPATH+file)
  return (tb[:,:3],tb[:,3:])

def rhoMassAFUV(file='RhoMassAFUV.dat'):
  pma = np.genfromtxt(constants.GRIDPATH+file)
  return (pma[:,:2],pma[:,2])

def speciesData(file='frequencies.dat'):
  # Open file containing the transition frequencies of various elements
  frequencies = np.genfromtxt(constants.GRIDPATH+file, names=['number', 'species', 'transition', 'frequency'], dtype="i8,U8,i8,f8", delimiter=',')
  return (frequencies['number'],frequencies['species'],frequencies['transition'],frequencies['frequency'])

def eTildeReal(file='Ereal.dat'):
  eReal = np.genfromtxt(constants.GRIDPATH+file, names=['x', 'Ereal'])
  return (eReal['x'],eReal['Ereal'])

def eTildeImaginary(file='Eimag.dat'):
  eImaginary = np.genfromtxt(constants.GRIDPATH+file, names=['x', 'Eimaginary'])
  return (eImaginary['x'],eImaginary['Eimaginary'])