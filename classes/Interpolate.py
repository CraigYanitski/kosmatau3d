import numpy as np
import sys
import scipy.interpolate as interpolate
import importlib as il
import Observations
class Interpolate(object):
  '''
  This is a class that can be used for the interpolation of the input data.
  It will contain functions to interpolate the intensity or optical depth
  for any species, given an index.

  The method of interpolation is passed as an argument when initialising
  this class. The acceptabled values are 'linear', 'cubic', and 'radial'.
  The default method is 'linear'. For the large intensity and optical depth
  grids, 'cubic' and 'radial' are the same.
  '''
  # PRIVATE
  def __init__(self, species, observations, directory='MilkyWay', interpolate='linear', verbose=False):
    self.__species = species
    self.__indeces = np.append(species[0].getFileIndeces(), species[1].getFileIndeces())
    # for element in species:
    #   for i in element.getFileIndeces(): self.__indeces.append(i)
    self.__observations = observations
    self.__interpolation = interpolate
    self.__verbose = verbose
    self.__intensityInterpolation,self.__tauInterpolation = self.__calculateGridInterpolation()
    self.__rotationInterpolation = self.__calculateRotationVelocity()
    self.__dispersionInterpolation = self.__calculateVelocityDispersion()
    self.__densityInterpolation = self.__calculateDensity()
    self.__clumpMassInterpolation = self.__clumpMassProfile()
    self.__interclumpMassInterpolation = self.__interclumpMassProfile()
    self.__FUVextinctionInterpolation = self.__interpolateFUVextinction()
    self.__FUVfieldInterpolation = self.__interpolateFUVfield()

  def __calculateGridInterpolation(self):
    np.seterr(divide='ignore')
    nmuvI,I = self.__observations.tbCenterline
    nmuvTau,Tau = self.__observations.tauCenterline
    intensityInterpolation = []
    tauInterpolation = []
    nmuvI /= 10.     #begin 'decoding' the grid for the interpolation
    nmuvTau /= 10.
    with np.errstate(divide='ignore', invalid='ignore'):
      logI = np.log10(I)    #'encode' the intensity of the grid for interpolation
      logTau = np.log10(Tau)
    if self.__interpolation=='linear':
      for index in self.__indeces:
        if self.__verbose: print('Creating intensity grid interpolation')
        rInterpI = interpolate.LinearNDInterpolator(nmuvI, logI[:,index-1])
        if self.__verbose: print('Creating tau grid interpolation')
        rInterpTau = interpolate.LinearNDInterpolator(nmuvTau, logTau[:,index-1])
        intensityInterpolation.append(rInterpI)
        tauInterpolation.append(rInterpTau)
      return intensityInterpolation,tauInterpolation
    elif self.__interpolation=='radial' or self.__interpolation=='cubic':
      for index in self.__indeces:
        if self.__verbose: print('Creating intensity grid interpolation')
        rInterpI = interpolate.Rbf(nmuvI[:,0], nmuvI[:,1], nmuvI[:,2], logI[:,index-1])
        if self.__verbose: print('Creating tau grid interpolation')
        rInterpTau = interpolate.Rbf(nmuvTau[:,0], nmuvTau[:,1], nmuvTau[:,2], logTau[:,index-1])
        intensityInterpolation.append(rInterpI)
        tauInterpolation.append(rInterpTau)
      return intensityInterpolation,tauInterpolation
    else: sys.exit('<<ERROR>>: There is no such method as {} to interpolate the KOSMA-tau grid.\n\nExitting...\n\n'.format(self.__interpolation))
  def __calculateRotationVelocity(self):
    if self.__verbose: print('Creating rotation velocity interpolation')
    rotation = self.__observations.rotationProfile 
    if self.__interpolation=='linear':
      return interpolate.interp1d(rotation[0], rotation[1][:,0], kind='linear')    #rotation velocity interpolation
    if self.__interpolation=='cubic' or self.__interpolation=='radial':
      return interpolate.interp1d(rotation[0], rotation[1][:,0], kind='cubic')    #rotation velocity interpolation
    else: sys.exit('<<ERROR>>: There is no such method as {} to interpolate the velocity profile.\n\nExitting...\n\n'.format(self.__interpolation))
  def __calculateVelocityDispersion(self):
    if self.__verbose: print('Creating velocity dispersion interpolation')
    rotation = self.__observations.rotationProfile 
    if self.__interpolation=='linear':
      return interpolate.interp1d(rotation[0], rotation[1][:,1], kind='linear')    #rotation velocity interpolation
    if self.__interpolation=='cubic' or self.__interpolation=='radial':
      return interpolate.interp1d(rotation[0], rotation[1][:,1], kind='cubic')    #rotation velocity interpolation
    else: sys.exit('<<ERROR>>: There is no such method as {} to interpolate the velocity profile.\n\nExitting...\n\n'.format(self.__interpolation))
  def __calculateDensity(self):
    if self.__verbose: print('Creating density interpolation')
    density = self.__observations.densityProfile
    if self.__interpolation=='linear':
      return interpolate.interp1d(density[0], density[1], kind='linear')      #density interpolation
    elif self.__interpolation=='cubic' or self.__interpolation=='radial':
      return interpolate.interp1d(density[0], density[1], kind='cubic')      #density interpolation
    else: sys.exit('<<ERROR>>: There is no such method as {} to interpolate the KOSMA-tau grid.\n\nExitting...\n\n'.format(self.__interpolation))
  def __clumpMassProfile(self):
    if self.__verbose: print('Creating clump mass interpolation')
    clumpMass = self.__observations.clumpMassProfile
    if self.__interpolation=='linear':
      return interpolate.interp1d(clumpMass[0], clumpMass[1], kind='linear')  #clump mass interpolation
    elif self.__interpolation=='cubic' or self.__interpolation=='radial':
      return interpolate.interp1d(clumpMass[0], clumpMass[1], kind='cubic')  #clump mass interpolation
    else: sys.exit('<<ERROR>>: There is no such method as {} to interpolate the KOSMA-tau grid.\n\nExitting...\n\n'.format(self.__interpolation))
  def __interclumpMassProfile(self):
    if self.__verbose: print('Creating interclump mass interpolation')
    interclumpMass = self.__observations.clumpMassProfile
    if self.__interpolation=='linear':
      return interpolate.interp1d(interclumpMass[0], interclumpMass[1], kind='linear')   #interclump mass interpolation
    elif self.__interpolation=='cubic' or self.__interpolation=='radial':
      return interpolate.interp1d(interclumpMass[0], interclumpMass[1], kind='cubic')   #interclump mass interpolation
    else: sys.exit('<<ERROR>>: There is no such method as {} to interpolate the KOSMA-tau grid.\n\nExitting...\n\n'.format(self.__interpolation))
  def __interpolateFUVextinction(self):
    if self.__verbose: print('Creating A_UV grid interpolation')
    rhomass,AUV = self.__observations.rhoMassAFUV
    rhomass /= 10.
    logAUV = np.log10(AUV)
    if self.__interpolation=='linear':
      return interpolate.LinearNDInterpolator(rhomass, logAUV)
    elif self.__interpolation=='cubic' or self.__interpolation=='radial':
      return interpolate.Rbf(rhomass, logAUV)
    else: sys.exit('<<ERROR>>: There is no such method as {} to interpolate the extinction in the KOSMA-tau grid.\n\nExitting...\n\n'.format(self.__interpolation))
  def __interpolateFUVfield(self):
    if self.__verbose: print('Creating FUV interpolation')
    fuv = self.__observations.FUVfield
    if self.__interpolation=='linear':
      return interpolate.interp1d(fuv[0], fuv[1], kind='linear')
    if self.__interpolation=='cubic' or self.__interpolation=='radial':
      return interpolate.interp1d(fuv[0], fuv[1], kind='cubic')
    else: sys.exit('<<ERROR>>: There is no such method as {} to interpolate the KOSMA-tau grid.\n\nExitting...\n\n'.format(self.__interpolation))
  def __str__(self):
    return 'Available Interpolations:\n -Clump intensity\n -Clump optical depth\n -Clump mass (galactic)\n -Clump density (galactic)\n -Voxel rotation (galactic)\n -UV extinction\n -FUV field (galactic)'

  # PUBLIC
  def reloadModules(self):
    il.reload(Observations)
    return
  def getObservations(self):
    return self.__observations
  def interpolateIntensity(self, points, speciesNumber, verbose=False):
    verbose = verbose or self.__verbose
    #points = np.log10(points)
    if len(speciesNumber):
      intensity = []
      for i in speciesNumber:
        if self.__interpolation=='linear': intensity.append(10**self.__intensityInterpolation[i](points))
        elif self.__interpolation=='radial' or self.__interpolation=='cubic': intensity.append(10**self.__intensityInterpolation[i](points[0], points[1], points[2]))
        if np.isnan(intensity[-1]) or intensity[-1]==0: intensity[-1] = 10**-100
      if verbose:
        print('Calculated the intensity for {} species.'.format(len(speciesNumber)))
    else:
      if verbose:
        print('There are no species of this type adding to the intensity.')
      intensity = 0
    return np.array(intensity)
  def interpolateTau(self, points, speciesNumber, verbose=False):
    verbose = verbose or self.__verbose
    #points = np.log10(points)
    if len(speciesNumber):
      tau = []
      for i in speciesNumber:
        if self.__interpolation=='linear': tau.append(10**self.__tauInterpolation[i](points))
        elif self.__interpolation=='radial' or self.__interpolation=='cubic': tau.append(10**self.__tauInterpolation[i](points[0], points[1], points[2]))
        if np.isnan(tau[-1]): tau[-1] = 10**-100
        elif tau[-1]<=0:
          temp = tau[-1]
          tau[-1] = 10**-100
          input('\n<<ERROR>> Negative opacity {} found.\n'.format(temp))
      if verbose:
        print('Calculated the optical depth for {} species.'.format(len(speciesNumber)))
    else:
      if verbose:
        print('There are no species adding to the optical depth.')
      tau = 0
    return np.array(tau)
  def interpolateRotationalVelocity(self, radius):
    return self.__rotationInterpolation(radius)
  def interpolateVelocityDispersion(self, radius):
    return self.__dispersionInterpolation(radius)
  def interpolateDensity(self, radius):
    density = self.__densityInterpolation(radius)
    if (density<0).any():
      input('<<ERROR>> density {} at radius {} pc!'.format(density, radius))
      sys.exit()
    return density
  def interpolateClumpMass(self, radius):
    mass=self.__clumpMassInterpolation(radius)
    if (mass<0).any():
      input('<<ERROR>> clump mass {} at radius {} pc!'.format(mass, radius))
      sys.exit()
    return mass
  def interpolateInterclumpMass(self, radius):
    mass=self.__interclumpMassInterpolation(radius)
    if (mass<0).any():
      input('<<ERROR>> interclump mass {} at radius {} pc!'.format(mass, radius))
      sys.exit()
    return mass
  def interpolateFUVextinction(self, density, mass):
    return 10**self.__FUVextinctionInterpolation(density, mass)
  def interpolateFUVfield(self, radius):
    return self.__FUVfieldInterpolation(radius)