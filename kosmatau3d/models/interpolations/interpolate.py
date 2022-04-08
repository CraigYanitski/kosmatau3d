import sys
import importlib as il

import numpy as np
import scipy.interpolate as interpolate
import dill
import warnings

from copy import copy
from numba import jit_module
from numba.core.errors import NumbaWarning, NumbaDeprecationWarning, NumbaPendingDeprecationWarning

from kosmatau3d.models import constants
from kosmatau3d.models import interpolations
from kosmatau3d.models import observations
from kosmatau3d.models import species


warnings.simplefilter('ignore', category=NumbaWarning)
warnings.simplefilter('ignore', category=NumbaDeprecationWarning)
warnings.simplefilter('ignore', category=NumbaPendingDeprecationWarning)


def initialise_grid(dilled=False):
    calculateGridInterpolation(dilled=dilled)
    interpolations.FUVextinctionInterpolation = interpolateFUVextinction(dilled=dilled)
    interpolations.initialised = True
    return


def initialise_model(dilled=False):
    interpolations.rotationInterpolation = calculateRotationVelocity(dilled=dilled)
    interpolations.dispersionInterpolation = 1  # calculateVelocityDispersion(dilled=dilled)
    interpolations.densityInterpolation = calculateDensity(dilled=dilled)
    interpolations.clumpMassInterpolation = clumpMassProfile(dilled=dilled)
    interpolations.interclumpMassInterpolation = interclumpMassProfile(dilled=dilled)
    interpolations.FUVfieldInterpolation = interpolateFUVfield(dilled=dilled)
    # interpolations.eTildeReal = interpolateETildeReal()
    # interpolations.eTildeImaginary = interpolateETildeImaginary()
    return


# Grid interpolation


def calculateGridInterpolation(verbose=False, dilled=False):
    '''
    This interpolates the needed species from the KOSMA-tau emission grids. It automatically interpolates
    the dust continuum as well, even if there are no molecules specified. This is kept as a separate interpolation
    function than the molecular emissions since it is a continuum. The same grid files are used, however, so it is
    necessary to note the order of the emissions in the grids.
  
    There are 504 emission values for a given gridpoint, consisting of 171 molecular transitions and the emission
    at 333 wavelengths for the dust continuum. the indeces used for the dust is therefore [170:503], to refer to
    for the corresponding emission indeces.
    '''
    np.seterr(divide='ignore', invalid='ignore')
    if dilled:
        with open(constants.GRIDPATH + 'dilled/intensity_interpolation', 'rb') as file:
            intensityInterpolation = dill.load(file)
        with open(constants.GRIDPATH + 'dilled/tau_interpolation', 'rb') as file:
            tauInterpolation = dill.load(file)
        with open(constants.GRIDPATH + 'dilled/dust_intensity_interpolation', 'rb') as file:
            dustIntensityInterpolation = dill.load(file)
        with open(constants.GRIDPATH + 'dilled/dust_tau_interpolation', 'rb') as file:
            dustTauInterpolation = dill.load(file)
        interpolations.intensityInterpolation = []
        interpolations.tauInterpolation = []
        if len(species.moleculeIndeces) == len(intensityInterpolation):
            interpolations.intensityInterpolation = intensityInterpolation
            interpolations.tauInterpolation = tauInterpolation
        else:
            for index in species.moleculeIndeces:
                interpolations.intensityInterpolation.append(intensityInterpolation[index[0][0]])
                interpolations.tauInterpolation.append(tauInterpolation[index[0][0]])
        if constants.dust:
            if np.where(constants.nDust)[0].size == len(dustIntensityInterpolation):
                interpolations.dustIntensityInterpolation = dustIntensityInterpolation
                interpolations.dustTauInterpolation = dustTauInterpolation
            else:
                interpolations.dustIntensityInterpolation = []
                interpolations.dustTauInterpolation = []
                for i in np.where(constants.nDust)[0]:
                    interpolations.dustIntensityInterpolation.append(copy(dustIntensityInterpolation[i]))
                    interpolations.dustTauInterpolation.append(copy(dustTauInterpolation[i]))
        return #speciesIntensityInterpolation, speciesTauInterpolation

    # indeces = species.molecules.getFileIndeces()
    crnmuvI, I = observations.tbCenterline
    crnmuvTau, Tau = observations.tauCenterline
    interpolations.intensityInterpolation = []
    interpolations.tauInterpolation = []
    
    # Correct for the negative emission values (from Silke's code)
    I[I <= 0] = 1e-100
    Tau[Tau <= 0] = 1e-100
    
    # Begin 'encoding' the intensity of the grid for interpolation
    logI = np.log10(I)
    logTau = np.log10(Tau)
  
    if constants.logEncoded:
        # Fully 'encode' the emission grid for interpolation
        logI *= 10
        logTau *= 10
    
    else:
        # Begin 'decoding' the grid for interpolation
        crnmuvI /= 10.
        crnmuvTau /= 10.
  
    if constants.interpolation == 'linear':
        if constants.dust:
            interpolations.dustIntensityInterpolation = []
            interpolations.dustTauInterpolation = []
            for i in np.where(constants.nDust)[0]:
                interpolations.dustIntensityInterpolation.append(interpolate.LinearNDInterpolator(crnmuvI,
                                                                                                  logI[:, constants.moleculeNumber+i]))
                interpolations.dustTauInterpolation.append(interpolate.LinearNDInterpolator(crnmuvTau,
                                                                                            logTau[:, constants.moleculeNumber+i]))
        for index in species.moleculeIndeces:
            if verbose:
                print('Creating intensity grid interpolation')
            rInterpI = interpolate.LinearNDInterpolator(crnmuvI, logI[:, index])
            if verbose:
                print('Creating tau grid interpolation')
            rInterpTau = interpolate.LinearNDInterpolator(crnmuvTau, logTau[:, index])
            interpolations.intensityInterpolation.append(rInterpI)
            interpolations.tauInterpolation.append(rInterpTau)
        return #intensityInterpolation, tauInterpolation
  
    elif constants.interpolation == 'radial' or constants.interpolation == 'cubic':
        if constants.dust:
            interpolations.dustIntensityInterpolation = interpolate.Rbf(crnmuvI[:, 0], crnmuvI[:, 1], crnmuvI[:, 2], crnmuvI[:, 3],
                                                                        logI[:, constants.moleculeNumber:][:, constants.nDust])
            interpolations.dustTauInterpolation = interpolate.Rbf(crnmuvTau[:, 0], crnmuvTau[:, 1], crnmuvTau[:, 2], crnmuvTau[:, 3], logTau[:, constants.moleculeNumber:][:, constants.nDust])
        for index in species.moleculeIndeces:
            if verbose:
                print('Creating intensity grid interpolation')
            rInterpI = interpolate.Rbf(crnmuvI[:, 0], crnmuvI[:, 1], crnmuvI[:, 2], crnmuvI[:, 3], logI[:, index])
            if verbose:
                print('Creating tau grid interpolation')
            rInterpTau = interpolate.Rbf(crnmuvTau[:, 0], crnmuvTau[:, 1], crnmuvTau[:, 2], crnmuvTau[:, 3], logTau[:, index])
            interpolations.intensityInterpolation.append(rInterpI)
            interpolations.tauInterpolation.append(rInterpTau)
        return #intensityInterpolation, tauInterpolation
      
    else:
        sys.exit('<<ERROR>>: There is no such method as {} to interpolate the '.format(constants.interpolation) +
                 'KOSMA-tau grid.\n\nExitting...\n\n')


def interpolateFUVextinction(verbose=False, dilled=False):
    if dilled:
        with open(constants.GRIDPATH + 'dilled/fuv_extinction_interpolation', 'rb') as file:
            return dill.load(file)
    if verbose:
        print('Creating A_UV grid interpolation')
    rhomass,AUV = observations.rhoMassAFUV
    rhomass /= 10.
    logAUV = np.log10(AUV)
    if constants.interpolation == 'linear':
        return interpolate.LinearNDInterpolator(rhomass, logAUV)
    elif constants.interpolation == 'cubic' or constants.interpolation == 'radial':
        return interpolate.Rbf(rhomass, logAUV)
    else:
        sys.exit('<<ERROR>>: There is no such method as {} to interpolate '.format(constants.interpolation) +
                 'the extinction in the KOSMA-tau grid.\n\nExitting...\n\n')


# Model interpolation


def calculateRotationVelocity(verbose=False, dilled=False):
    if constants.directory != '':
        if dilled:
            with open(constants.INPUTPATH+constants.directory + 'dilled/rotation_interpolation', 'rb') as file:
                return dill.load(file)
        if verbose:
            print('Creating rotation velocity interpolation')
        rotation = observations.rotationProfile
        if constants.interpolation == 'linear':
            return interpolate.interp1d(rotation[0], rotation[1][:, 0], kind='linear')
        elif constants.interpolation == 'cubic' or constants.interpolation == 'radial':
            return interpolate.interp1d(rotation[0], rotation[1][:, 0], kind='cubic')
        else:
            sys.exit('<<ERROR>>: There is no such method as {} to interpolate '.format(constants.interpolation) +
                     'the velocity profile.\n\nExitting...\n\n')
    return


def calculateVelocityDispersion(verbose=False, dilled=False):
    if constants.directory != '':
        if dilled:
            with open(constants.INPUTPATH+constants.directory + 'dilled/dispersion_interpolation', 'rb') as file:
                return dill.load(file)
        if verbose:
            print('Creating velocity dispersion interpolation')
        rotation = observations.rotationProfile
        if constants.interpolation == 'linear':
            return interpolate.interp1d(rotation[0], rotation[1][:, 1], kind='linear')
        elif constants.interpolation == 'cubic' or constants.interpolation == 'radial':
            return interpolate.interp1d(rotation[0], rotation[1][:, 1], kind='cubic')
        else:
            sys.exit('<<ERROR>>: There is no such method as {} to interpolate '.format(constants.interpolation) +
                     'the velocity profile.\n\nExitting...\n\n')
    return


def calculateDensity(verbose=False, dilled=False):
    if constants.directory != '':
        if dilled:
            with open(constants.INPUTPATH+constants.directory * 'dilled/density_interpolation', 'rb') as file:
                return dill.load(file)
        if verbose:
            print('Creating density interpolation')
        density = observations.densityProfile
        if constants.interpolation == 'linear':
            return interpolate.interp1d(density[0], density[1], kind='linear')
        elif constants.interpolation == 'cubic' or constants.interpolation == 'radial':
            return interpolate.interp1d(density[0], density[1], kind='cubic')
        else:
            sys.exit('<<ERROR>>: There is no such method as {} to interpolate '.format(constants.interpolation) +
                     'the KOSMA-tau grid.\n\nExitting...\n\n')
    return


def clumpMassProfile(verbose=False, dilled=False):
    if constants.directory != '':
        if dilled:
            with open(constants.INPUTPATH+constants.directory + 'dilled/clump_mass_interpolation', 'rb') as file:
                return dill.load(file)
        if verbose:
            print('Creating clump mass interpolation')
        clumpMass = observations.clumpMassProfile
        if constants.interpolation == 'linear':
            return interpolate.interp1d(clumpMass[0], clumpMass[1], kind='linear')
        elif constants.interpolation == 'cubic' or constants.interpolation == 'radial':
            return interpolate.interp1d(clumpMass[0], clumpMass[1], kind='cubic')
        else:
            sys.exit('<<ERROR>>: There is no such method as {} to interpolate '.format(constants.interpolation) +
                     'the KOSMA-tau grid.\n\nExitting...\n\n')
    return


def interclumpMassProfile(verbose=False, dilled=True):
    if constants.directory != '':
        if dilled:
            with open(constants.INPUTPATH+constants.directory + 'dilled/interclump_mass_interpolation', 'rb') as file:
                return dill.load(file)
        if verbose:
            print('Creating interclump mass interpolation')
        interclumpMass = observations.interclumpMassProfile
        if constants.interpolation == 'linear':
            return interpolate.interp1d(interclumpMass[0], interclumpMass[1], kind='linear')
        elif constants.interpolation == 'cubic' or constants.interpolation == 'radial':
            return interpolate.interp1d(interclumpMass[0], interclumpMass[1], kind='cubic')
        else:
            sys.exit('<<ERROR>>: There is no such method as {} to interpolate '.format(constants.interpolation) +
                     'the KOSMA-tau grid.\n\nExitting...\n\n')
    return


def interpolateFUVfield(verbose=False, dilled=False):
    if constants.directory != '':
        if dilled:
            with open(constants.INPUTPATH+constants.directory + 'dilled/fuv_field_interpolation', 'rb') as file:
                return dill.load(file)
        if verbose:
            print('Creating FUV interpolation')
        fuv = observations.FUVfield
        lam = np.array([912, 1350, 1500, 1650, 2000, 2200, 2500, 2800, 3650])
        wav = np.linspace(912, 2066, num=1000)
        f = interpolate.interp1d(lam, fuv[2], axis=1)
        if constants.interpolation == 'linear':
            return interpolate.LinearNDInterpolator(fuv[0], np.trapz(f(wav), wav))
        elif constants.interpolation == 'cubic' or constants.interpolation == 'radial':
            return interpolate.Rbf(fuv[0][:, 0], fuv[0][:, 1], np.trapz(f(wav), wav))
        else:
            sys.exit('<<ERROR>>: There is no such method as {} to interpolate '.format(constants.interpolation) +
                     'the KOSMA-tau grid.\n\nExitting...\n\n')
    return


def interpolateETildeReal():
    return interpolate.interp1d(observations.eTildeReal[0], observations.eTildeReal[1], kind='linear')


def interpolateETildeImaginary():
    return interpolate.interp1d(observations.eTildeImaginary[0], observations.eTildeImaginary[1], kind='linear')


# jit_module(nopython=False)
