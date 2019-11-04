import numpy as np
import sys
import scipy.interpolate as interpolate
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
  def __init__(self, species, observations, directory='MilkyWay', interpolate='linear'):
    self.__species = species
    self.__observations = observations
    self.__interpolation = interpolation
    self.__observations = Observations()
    self.__intensityInterpolation,self.__tauInterpolation = self.__calculateGridInterpolation()
    self.__rotationInterpolation = self.__calculaterotationVelocity()
    self.__densityInterpolation = self.__calculateDensity()
    self.__clumpMassInterpolation = self.__clumpMassProfile()
    self.__interclumpMassInterpolation = self.__interclumpMassProfile()
    self.__FUVextinctionInterpolation = self.__interpolateFUVextinction()
    self.__FUVfieldInterpolation = self.__interpolateFUVfield()

  def __calculateGridInterpolation(self):
    nI,massI,uvI,I = self.__observations.tbCenterline()
    nTau,massTau,uvTau,Tau = self.__observations.tauCenterline()
    intensityInterpolation = []
    tauInterpolation = []
    if self.__interplation=='linear':
      for index in self.__species.indeces():
        rInterpI = interpolate.LinearNDInterpolation(nI, massI, uvI, I[index])
        rInterpTau = interpolate.LinearNDInterpolation(nTau, massTau, uvTau, Tau[index])
        intensityInterpolation.append(rInterpI)
        tauInterpolation.append(rInterpTau)
      return intensityInterpolation,tauInterpolation
    elif self.__interplation=='radial' or self.__interpolation=='cubic':
      for index in self.__species.indeces():
        rInterpI = interpolate.Rbf(nI, massI, uvI, I[index])
        rInterpTau = interpolate.Rbf(nTau, massTau, uvTau, Tau[index])
        intensityInterpolation.append(rInterpI)
        tauInterpolation.append(rInterpTau)
      return intensityInterpolation,tauInterpolation
    else: sys.exit('<<ERROR>>: There is no such method as {} to interpolate the KOSMA-tau grid.\n\n \
                   Exitting...\n\n'.format(self.__interpolation))
  def __calculateRotationVelocity(self):
    rotation = self.__observations.rotationProfile() 
    if self.__interpolation=='linear':
      return sp.interp1d(rotation[0], rotation[1], kind='linear')    #rotation velocity interpolation
    if self.__interpolation=='cubic' or self.__interpolation=='radial':
      return sp.interp1d(rotation[0], rotation[1], kind='cubic')    #rotation velocity interpolation
    else: sys.exit('<<ERROR>>: There is no such method as {} to interpolate the KOSMA-tau grid.\n\n \
                   Exitting...\n\n'.format(self.__interpolation))
  def __calculateDensity(self):
    density = self.__observations.densityProfile()
    if self.__interpolation=='linear':
      return interpolate.interp1d(density[0], density[1], kind='linear')      #density interpolation
    elif self.__interpolation=='cubic' or self.__interpolation=='radial':
      return interpolate.interp1d(density[0], density[1], kind='cubic')      #density interpolation
    else: sys.exit('<<ERROR>>: There is no such method as {} to interpolate the KOSMA-tau grid.\n\n \
                   Exitting...\n\n'.format(self.__interpolation))
  def __clumpMassProfile(self):
    clumpmass = self.__observations.interclumpMassProfile()
    if self.__interpolation=='linear':
      return interpolate.interp1d(clumpMass[0], clumpMass[1], kind='cubic')  #clump mass interpolation
    elif self.__interpolation=='cubic' or self.__interpolation=='radial':
      return interpolate.interp1d(clumpMass[0], clumpMass[1], kind='cubic')  #clump mass interpolation
    else: sys.exit('<<ERROR>>: There is no such method as {} to interpolate the KOSMA-tau grid.\n\n \
                   Exitting...\n\n'.format(self.__interpolation))
  def __interclumpMassProfile(self):
    interclumpmass = self.__observations.clumpMassProfile()
    if self.__interpolation=='linear':
      return interpolate.interp1d(interclumpMass[0], interclumpMass[1], kind='linear')   #interclump mass interpolation
    elif self.__interpolation=='cubic' or self.__interpolation=='radial':
      return interpolate.interp1d(interclumpMass[0], interclumpMass[1], kind='cubic')   #interclump mass interpolation
    else: sys.exit('<<ERROR>>: There is no such method as {} to interpolate the KOSMA-tau grid.\n\n \
                   Exitting...\n\n'.format(self.__interpolation))
  def __interpolateFUVextinction(self):
    afuv = self.__observations.rhoMassAFUV()
    if self.__interpolation=='linear':
      return interpolate.interp2d(afuv[:2], afuv[2], kind='linear')
    elif self.__interpolation=='cubic' or self.__interpolation=='radial':
      return interpolate.interp2d(afuv[:2], afuv[2], kind='cubic')
    else: sys.exit('<<ERROR>>: There is no such method as {} to interpolate the extinction in the KOSMA-tau grid.\n\n \
                   Exitting...\n\n'.format(self.__interpolation))
  def __interpolateFUVfield(self):
    fuv = self.__observations.FUVfield()
    if self.__interpolation=='linear':
      return interpolate.interp1d(fuv[0], fuv[1], kind='linear')
    if self.__interpolation=='cubic' or self.__interpolation=='radial':
      return interpolate.interp1d(fuv[0], fuv[1], kind='cubic')
    else: sys.exit('<<ERROR>>: There is no such method as {} to interpolate the KOSMA-tau grid.\n\n \
                   Exitting...\n\n'.format(self.__interpolation))

  # PUBLIC
  def interpolateIntensity(self, points, species):
    if len(species)>1:
      intensity = []
      for i in species.number:
        intensity.append(self.__intensityInterpolation(i))
    elif len(species):
      intensity = self.__intensityInterpolation(points)
    else: sys.exit('<<ERROR>>: a zero-length species array has been used to interpolate the intensity.\n\n')
    return np.array(intensity)
  def interpolateTau(self, points, species):
    if len(species)>1:
      tau = []
      for i in species.number:
        tau.append(self.__tauInterpolation(i))
    elif len(species):
      tau = self.__tauInterpolation(points)
    else: sys.exit('<<ERROR>>: a zero-length species array has been used to interpolate the optical depth.\n\n')
    return np.array(tau)
  def interpolateRotation(self, radius):
    return self.__rotationInterpolation(radius)
  def interpolateDensity(self, radius):
    return self.__densityInterpolation(radius)
  def interpolateClumpMass(self, radius):
    return self.__clumpMassInterpolation(radius)
  def interpolateInterclumpMass(self, radius):
    return self.__interclumpMassInterpolation(radius)
  def interpolateFUVextinction(self, density, mass):
    return self.__FUVextinctionInterpolation([density, mass])
  def interpolateFUVfield(self, radius):
    return self.__FUVfieldInterpolation(radius)