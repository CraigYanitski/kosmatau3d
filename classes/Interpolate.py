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
  def __init__(self, species, interpolate='linear'):
    self.__listSpecies = species
    self.__interpolation = interpolation
    self.intensityInterpolation = []
    self.tauInterpolation = []
    self.__calculateGridInterpolation()
    self.__rotationInterpolation = self.__calculaterotationVelocity()
    self.__densityInterpolation = self.__calculateDensity()
    self.__clumpMassInterpolation = self.__clumpMassProfile()
    self.__interclumpMassInterpolation = self.__interclumpMassProfile()
    self.__FUVextinctionInterpolation = self.__interpolateFUVextinction()
    self.__FUVfieldInterpolation = self.__interpolateFUVfield()

  def __calculateGridInterpolation(self):
    nI,massI,uvI,I = obs.tbCenterline()
    nTau,massTau,uvTau,Tau = obs.tauCenterline()
    if self.__interplation=='linear':
      for index in self.__listSpecies.indeces():
        rInterpI = sp.interpolate.LinearNDInterpolation(nI, massI, uvI, I[index])
        rInterpTau = sp.interpolate.LinearNDInterpolation(nTau, massTau, uvTau, Tau[index])
        self.__intensityInterpolation.append()
    elif self.__interplation=='radial' or self.__interpolation=='cubic':
      for index in self.__listSpecies.indeces():
        rInterpI = sp.interpolate.Rbf(nI, massI, uvI, I[index])
        rInterpTau = sp.interpolate.Rbf(nTau, massTau, uvTau, Tau[index])
        self.__intensityInterpolation.append()
    else: sys.exit('<<ERROR>>: There is no such method as {} to interpolate the KOSMA-tau grid.\n\n \
                   Exitting...\n\n'.format(self.__interpolation))
    return
  def __calculateRotationVelocity(self):
    rotation = obs.rotationProfile() 
    if self.__interpolation=='linear':
      return sp.interp1d(rotation[0], rotation[1], kind='linear')    #rotation velocity interpolation
    if self.__interpolation=='cubic' or self.__interpolation=='radial':
      return sp.interp1d(rotation[0], rotation[1], kind='cubic')    #rotation velocity interpolation
    else: sys.exit('<<ERROR>>: There is no such method as {} to interpolate the KOSMA-tau grid.\n\n \
                   Exitting...\n\n'.format(self.__interpolation))
  def __calculateDensity(self):
    density = obs.densityProfile()
    if self.__interpolation=='linear':
      return sp.interp1d(density[0], density[1], kind='linear')      #density interpolation
    elif self.__interpolation=='cubic' or self.__interpolation=='radial':
      return sp.interp1d(density[0], density[1], kind='cubic')      #density interpolation
    else: sys.exit('<<ERROR>>: There is no such method as {} to interpolate the KOSMA-tau grid.\n\n \
                   Exitting...\n\n'.format(self.__interpolation))
  def __clumpMassProfile(self):
    clumpmass = obs.interclumpMassProfile()
    if self.__interpolation=='linear':
      return sp.interp1d(clumpMass[0], clumpMass[1], kind='cubic')  #clump mass interpolation
    elif self.__interpolation=='cubic' or self.__interpolation=='radial':
      return sp.interp1d(clumpMass[0], clumpMass[1], kind='cubic')  #clump mass interpolation
    else: sys.exit('<<ERROR>>: There is no such method as {} to interpolate the KOSMA-tau grid.\n\n \
                   Exitting...\n\n'.format(self.__interpolation))
  def __interclumpMassProfile(self):
    interclumpmass = obs.clumpMassProfile()
    if self.__interpolation=='linear':
      return sp.interp1d(interclumpMass[0], interclumpMass[1], kind='linear')   #interclump mass interpolation
    elif self.__interpolation=='cubic' or self.__interpolation=='radial':
      return sp.interp1d(interclumpMass[0], interclumpMass[1], kind='cubic')   #interclump mass interpolation
    else: sys.exit('<<ERROR>>: There is no such method as {} to interpolate the KOSMA-tau grid.\n\n \
                   Exitting...\n\n'.format(self.__interpolation))
  def __interpolateFUVextinction(self):
    afuv = obs.rhoMassAFUV()
    if self.__interpolation=='linear':
      return sp.interpolate.interp2d(afuv[:2], afuv[2], kind='linear')
    elif self.__interpolation=='cubic' or self.__interpolation=='radial':
      return sp.interpolate.interp2d(afuv[:2], afuv[2], kind='cubic')
    else: sys.exit('<<ERROR>>: There is no such method as {} to interpolate the extinction in the KOSMA-tau grid.\n\n \
                   Exitting...\n\n'.format(self.__interpolation))
  def __interpolateFUVfield(self):
    fuv = obs.FUVfield()
    fuvInterp = []
    if self.__interpolation=='linear':
      for i in range(len(fuv[1])):
        fuvInterp.append(sp.interpolate.interp1d(afuv[i]), kind='linear')
    if self.__interpolation=='cubic' or self.__interpolation=='radial':
      for i in range(len(fuv[1])):
        fuvInterp.append(sp.interpolate.interp1d(afuv[i]), kind='cubic')
    else: sys.exit('<<ERROR>>: There is no such method as {} to interpolate the KOSMA-tau grid.\n\n \
                   Exitting...\n\n'.format(self.__interpolation))
    return fuvInterp

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
      for i in range(species):
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