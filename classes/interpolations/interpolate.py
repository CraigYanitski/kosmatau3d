from numba import jit_module
import numpy as np
import sys
import scipy.interpolate as interpolate
import importlib as il

from numba.errors import NumbaWarning, NumbaDeprecationWarning, NumbaPendingDeprecationWarning
import warnings
warnings.simplefilter('ignore', category=NumbaWarning)
warnings.simplefilter('ignore', category=NumbaDeprecationWarning)
warnings.simplefilter('ignore', category=NumbaPendingDeprecationWarning)

import constants
import interpolations
import observations
import species

def initialise():
  interpolations.intensityInterpolation,interpolations.tauInterpolation = calculateGridInterpolation()
  interpolations.rotationInterpolation = calculateRotationVelocity()
  interpolations.dispersionInterpolation = calculateVelocityDispersion()
  interpolations.densityInterpolation = calculateDensity()
  interpolations.clumpMassInterpolation = clumpMassProfile()
  interpolations.interclumpMassInterpolation = interclumpMassProfile()
  interpolations.FUVextinctionInterpolation = interpolateFUVextinction()
  interpolations.FUVfieldInterpolation = interpolateFUVfield()

def calculateGridInterpolation(verbose=False):
  np.seterr(divide='ignore', invalid='ignore')
  indeces = np.append(species.molecules.getFileIndeces(), species.dust.getFileIndeces())
  nmuvI,I = observations.tbCenterline
  nmuvTau,Tau = observations.tauCenterline
  intensityInterpolation = []
  tauInterpolation = []
  nmuvI /= 10.     #begin 'decoding' the grid for the interpolation
  nmuvTau /= 10.
  #with np.errstate(divide='ignore', invalid='ignore'):
  logI = np.log10(I)    #'encode' the intensity of the grid for interpolation
  logTau = np.log10(Tau)
  if constants.interpolation=='linear':
    for index in indeces:
      if verbose: print('Creating intensity grid interpolation')
      rInterpI = interpolate.LinearNDInterpolator(nmuvI, logI[:,index-1])
      if verbose: print('Creating tau grid interpolation')
      rInterpTau = interpolate.LinearNDInterpolator(nmuvTau, logTau[:,index-1])
      intensityInterpolation.append(rInterpI)
      tauInterpolation.append(rInterpTau)
    return intensityInterpolation,tauInterpolation
  elif constants.interpolation=='radial' or constants.interpolation=='cubic':
    for index in indeces:
      if verbose: print('Creating intensity grid interpolation')
      rInterpI = interpolate.Rbf(nmuvI[:,0], nmuvI[:,1], nmuvI[:,2], logI[:,index-1])
      if verbose: print('Creating tau grid interpolation')
      rInterpTau = interpolate.Rbf(nmuvTau[:,0], nmuvTau[:,1], nmuvTau[:,2], logTau[:,index-1])
      intensityInterpolation.append(rInterpI)
      tauInterpolation.append(rInterpTau)
    return intensityInterpolation,tauInterpolation
  else: sys.exit('<<ERROR>>: There is no such method as {} to interpolate the KOSMA-tau grid.\n\nExitting...\n\n'.format(interpolation))

def calculateRotationVelocity(verbose=False):
  if verbose: print('Creating rotation velocity interpolation')
  rotation = observations.rotationProfile 
  if constants.interpolation=='linear':
    return interpolate.interp1d(rotation[0], rotation[1][:,0], kind='linear')    #rotation velocity interpolation
  elif constants.interpolation=='cubic' or constants.interpolation=='radial':
    return interpolate.interp1d(rotation[0], rotation[1][:,0], kind='cubic')    #rotation velocity interpolation
  else: sys.exit('<<ERROR>>: There is no such method as {} to interpolate the velocity profile.\n\nExitting...\n\n'.format(interpolation))

def calculateVelocityDispersion(verbose=False):
  if verbose: print('Creating velocity dispersion interpolation')
  rotation = observations.rotationProfile 
  if constants.interpolation=='linear':
    return interpolate.interp1d(rotation[0], rotation[1][:,1], kind='linear')    #rotation velocity interpolation
  elif constants.interpolation=='cubic' or constants.interpolation=='radial':
    return interpolate.interp1d(rotation[0], rotation[1][:,1], kind='cubic')    #rotation velocity interpolation
  else: sys.exit('<<ERROR>>: There is no such method as {} to interpolate the velocity profile.\n\nExitting...\n\n'.format(interpolation))

def calculateDensity(verbose=False):
  if verbose: print('Creating density interpolation')
  density = observations.densityProfile
  if constants.interpolation=='linear':
    return interpolate.interp1d(density[0], density[1], kind='linear')      #density interpolation
  elif constants.interpolation=='cubic' or constants.interpolation=='radial':
    return interpolate.interp1d(density[0], density[1], kind='cubic')      #density interpolation
  else: sys.exit('<<ERROR>>: There is no such method as {} to interpolate the KOSMA-tau grid.\n\nExitting...\n\n'.format(interpolation))

def clumpMassProfile(verbose=False):
  if verbose: print('Creating clump mass interpolation')
  clumpMass = observations.clumpMassProfile
  if constants.interpolation=='linear':
    return interpolate.interp1d(clumpMass[0], clumpMass[1], kind='linear')  #clump mass interpolation
  elif constants.interpolation=='cubic' or constants.interpolation=='radial':
    return interpolate.interp1d(clumpMass[0], clumpMass[1], kind='cubic')  #clump mass interpolation
  else: sys.exit('<<ERROR>>: There is no such method as {} to interpolate the KOSMA-tau grid.\n\nExitting...\n\n'.format(interpolation))

def interclumpMassProfile(verbose=False):
  if verbose: print('Creating interclump mass interpolation')
  interclumpMass = observations.clumpMassProfile
  if constants.interpolation=='linear':
    return interpolate.interp1d(interclumpMass[0], interclumpMass[1], kind='linear')   #interclump mass interpolation
  elif constants.interpolation=='cubic' or constants.interpolation=='radial':
    return interpolate.interp1d(interclumpMass[0], interclumpMass[1], kind='cubic')   #interclump mass interpolation
  else: sys.exit('<<ERROR>>: There is no such method as {} to interpolate the KOSMA-tau grid.\n\nExitting...\n\n'.format(interpolation))

def interpolateFUVextinction(verbose=False):
  if verbose: print('Creating A_UV grid interpolation')
  rhomass,AUV = observations.rhoMassAFUV
  rhomass /= 10.
  logAUV = np.log10(AUV)
  if constants.interpolation=='linear':
    return interpolate.LinearNDInterpolator(rhomass, logAUV)
  elif constants.interpolation=='cubic' or constants.interpolation=='radial':
    return interpolate.Rbf(rhomass, logAUV)
  else: sys.exit('<<ERROR>>: There is no such method as {} to interpolate the extinction in the KOSMA-tau grid.\n\nExitting...\n\n'.format(interpolation))

def interpolateFUVfield(verbose=False):
  if verbose: print('Creating FUV interpolation')
  fuv = observations.FUVfield
  if constants.interpolation=='linear':
    return interpolate.interp1d(fuv[0], fuv[1], kind='linear')
  elif constants.interpolation=='cubic' or constants.interpolation=='radial':
    return interpolate.interp1d(fuv[0], fuv[1], kind='cubic')
  else: sys.exit('<<ERROR>>: There is no such method as {} to interpolate the KOSMA-tau grid.\n\nExitting...\n\n'.format(interpolation))

jit_module(nopython=False)
