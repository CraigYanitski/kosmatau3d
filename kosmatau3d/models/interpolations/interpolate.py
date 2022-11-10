import sys
import importlib as il

import numpy as np
import pandas as pd
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
    calculate_grid_interpolation(dilled=dilled)
    interpolations.taufuv_interpolation = calculate_taufuv(dilled=dilled)
    interpolations.hi_column_density_interpolation = calculate_hi_column_density(dilled=dilled)
    interpolations.h2_column_density_interpolation = calculate_h2_column_density(dilled=dilled)
    interpolations.tg_interpolation = calculate_tg(dilled=dilled)
    interpolations.td_interpolation = calculate_td(dilled=dilled)
    interpolations.initialised = True
    return


def initialise_model(dilled=False, like_clumps=False, 
                     average_fuv=False, l_range=(912, 2066)):
    interpolations.galaxy_rotation_interpolation = calculate_galaxy_rotation(dilled=dilled)
    interpolations.velocity_dispersion_interpolation = 1  # calculate_velocity_dispersion(dilled=dilled)
    interpolations.number_density_interpolation = calculate_number_density(dilled=dilled)
    interpolations.h2_scale_height = calculate_h2_height()
    interpolations.h2_mass_interpolation = calculate_h2_mass(dilled=dilled)
    interpolations.hi_scale_height = calculate_hi_height(like_clumps=like_clumps)
    interpolations.hi_mass_interpolation = calculate_hi_mass(dilled=dilled,
                                                             like_clumps=like_clumps)
    interpolations.fuv_interpolation = calculate_fuv(l_range=l_range, 
                                                     average_fuv=average_fuv, 
                                                     dilled=dilled)
    # interpolations.e_tilde_real = calculate_e_tilde_real()
    # interpolations.e_tilde_imaginary = calculate_e_tilde_imaginary()
    return


# Grid interpolation


def calculate_grid_interpolation(verbose=False, dilled=False):
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
    if interpolations.initialised:
        del interpolations.species_intensity_interpolation
        del interpolations.species_tau_interpolation
        del interpolations.dust_intensity_interpolation
        del interpolations.dust_tau_interpolation
    interpolations.species_intensity_interpolation = []
    interpolations.species_tau_interpolation = []
    interpolations.dust_intensity_interpolation = []
    interpolations.dust_tau_interpolation = []
    if dilled:
        with open(constants.GRIDPATH + 'dilled/intensity_interpolation', 'rb') as file:
            species_intensity_interpolation = dill.load(file)
        with open(constants.GRIDPATH + 'dilled/tau_interpolation', 'rb') as file:
            species_tau_interpolation = dill.load(file)
        with open(constants.GRIDPATH + 'dilled/dust_intensity_interpolation', 'rb') as file:
            dust_intensity_interpolation = dill.load(file)
        with open(constants.GRIDPATH + 'dilled/dust_tau_interpolation', 'rb') as file:
            dust_tau_interpolation = dill.load(file)
        # interpolations.intensityInterpolation = []
        # interpolations.tauInterpolation = []
        if len(species.molecule_indeces) == len(species_intensity_interpolation):
            interpolations.species_intensity_interpolation = species_intensity_interpolation
            interpolations.species_tau_interpolation = species_tau_interpolation
        else:
            for index in species.molecule_indeces:
                interpolations.species_intensity_interpolation.append(
                        species_intensity_interpolation[index[0][0]])
                interpolations.species_tau_interpolation.append(species_tau_interpolation[index[0][0]])
        if constants.dust:
            if np.where(constants.n_dust)[0].size == len(dust_intensity_interpolation):
                interpolations.dust_intensity_interpolation = dust_intensity_interpolation
                interpolations.dust_tau_interpolation = dust_tau_interpolation
            else:
                interpolations.dust_intensity_interpolation = []
                interpolations.dust_tau_interpolation = []
                for i in np.where(constants.n_dust)[0]:
                    interpolations.dust_intensity_interpolation.append(
                            copy(dust_intensity_interpolation[i]))
                    interpolations.dust_tau_interpolation.append(copy(dust_tau_interpolation[i]))
        return #speciesIntensityInterpolation, speciesTauInterpolation

    # indeces = species.molecules.getFileIndeces()
    crnmuvI, I = observations.tb_centerline
    crnmuvTau, Tau = observations.tau_centerline
    # interpolations.intensityInterpolation = []
    # interpolations.tauInterpolation = []
    
    # Correct for the negative emission values (from Silke's code)
    I[I <= 0] = 1e-100
    Tau[Tau <= 0] = 1e-100
    
    # Begin 'encoding' the intensity of the grid for interpolation
    logI = np.log10(I)
    logTau = np.log10(Tau)
  
    if constants.log_encoded:
        # Fully 'encode' the emission grid for interpolation
        logI = 10*logI
        logTau = 10*logTau
    
    else:
        # Begin 'decoding' the grid for interpolation
        crnmuvI = crnmuvI/10.
        crnmuvTau = crmnuvTau/10.
  
    if constants.interpolation == 'linear':
        if constants.dust:
            interpolations.dust_intensity_interpolation = []
            interpolations.dust_tau_interpolation = []
            for i in np.where(constants.n_dust)[0]:
                interpolations.dust_intensity_interpolation.append(
                        interpolate.LinearNDInterpolator(crnmuvI,
                                                         logI[:, constants.molecule_number+i]))
                interpolations.dust_tau_interpolation.append(
                        interpolate.LinearNDInterpolator(crnmuvTau, 
                                                         logTau[:, constants.molecule_number+i]))
        for index in species.molecule_indeces:
            if verbose:
                print('Creating intensity grid interpolation')
            rInterpI = interpolate.LinearNDInterpolator(crnmuvI, logI[:, index])
            if verbose:
                print('Creating tau grid interpolation')
            rInterpTau = interpolate.LinearNDInterpolator(crnmuvTau, logTau[:, index])
            interpolations.species_intensity_interpolation.append(rInterpI)
            interpolations.species_tau_interpolation.append(rInterpTau)
        return #intensityInterpolation, tauInterpolation
  
    elif constants.interpolation == 'radial' or constants.interpolation == 'cubic':
        if constants.dust:
            interpolations.dust_intensity_interpolation = interpolate.Rbf(
                    crnmuvI[:, 0], crnmuvI[:, 1], crnmuvI[:, 2], crnmuvI[:, 3], 
                    logI[:, constants.moleculeNumber:][:, constants.n_dust])
            interpolations.dust_tau_interpolation = interpolate.Rbf(
                    crnmuvTau[:, 0], crnmuvTau[:, 1], crnmuvTau[:, 2], crnmuvTau[:, 3], 
                    logTau[:, constants.moleculeNumber:][:, constants.n_dust])
        for index in species.molecule_indeces:
            if verbose:
                print('Creating intensity grid interpolation')
            rInterpI = interpolate.Rbf(
                    crnmuvI[:, 0], crnmuvI[:, 1], crnmuvI[:, 2], crnmuvI[:, 3], logI[:, index])
            if verbose:
                print('Creating tau grid interpolation')
            rInterpTau = interpolate.Rbf(
                    crnmuvTau[:, 0], crnmuvTau[:, 1], crnmuvTau[:, 2], crnmuvTau[:, 3], logTau[:, index])
            interpolations.species_intensity_interpolation.append(rInterpI)
            interpolations.species_tau_interpolation.append(rInterpTau)
        return #intensityInterpolation, tauInterpolation
      
    else:
        sys.exit('<<ERROR>>: There is no such method as ' + 
                 '{} to interpolate the '.format(constants.interpolation) +
                 'KOSMA-tau grid.\n\nExitting...\n\n')


def calculate_taufuv(verbose=False, dilled=False):
    if dilled:
        with open(constants.GRIDPATH + 'dilled/taufuv_interpolation', 'rb') as file:
            return dill.load(file)
    if verbose:
        print('Creating A_UV grid interpolation')
    rhomass, taufuv = observations.rho_mass_taufuv
    rhomass = rhomass/10.
    log_taufuv = np.log10(taufuv)
    if constants.interpolation == 'linear':
        return interpolate.LinearNDInterpolator(rhomass, log_taufuv)
    elif constants.interpolation == 'cubic' or constants.interpolation == 'radial':
        return interpolate.Rbf(rhomass, log_taufuv)
    else:
        sys.exit('<<ERROR>>: There is no such method as ' + 
                 '{} to interpolate '.format(constants.interpolation) +
                 'the extinction in the KOSMA-tau grid.\n\nExitting...\n\n')


def calculate_hi_column_density(verbose=False, dilled=False):
    if dilled:
        with open(constants.GRIDPATH + 'dilled/hi_column_density_interpolation', 'rb') as file:
            return dill.load(file)
    if verbose:
        print('Creating N(H) grid interpolation')
    nmchi, N = observations.column_density
    nmchi = nmchi/10.
    log_N = np.log10(N)
    if constants.interpolation == 'linear':
        return interpolate.LinearNDInterpolator(nmchi.to_numpy(), log_N['H'])
    elif constants.interpolation == 'cubic' or constants.interpolation == 'radial':
        return interpolate.Rbf(nmchi, log_N['H'])
    else:
        sys.exit('<<ERROR>>: There is no such method as ' + 
                 '{} to interpolate '.format(constants.interpolation) +
                 'the H column density in the KOSMA-tau grid.\n\nExitting...\n\n')


def calculate_h2_column_density(verbose=False, dilled=False):
    if dilled:
        with open(constants.GRIDPATH + 'dilled/h2_column_density_interpolation', 'rb') as file:
            return dill.load(file)
    if verbose:
        print('Creating N(H2) grid interpolation')
    nmchi, N = observations.column_density
    nmchi = nmchi/10.
    log_N = np.log10(N)
    if constants.interpolation == 'linear':
        return interpolate.LinearNDInterpolator(nmchi.to_numpy(), log_N['H2'])
    elif constants.interpolation == 'cubic' or constants.interpolation == 'radial':
        return interpolate.Rbf(nmchi, log_N['H2'])
    else:
        sys.exit('<<ERROR>>: There is no such method as ' + 
                 '{} to interpolate '.format(constants.interpolation) +
                 'the H2 column density in the KOSMA-tau grid.\n\nExitting...\n\n')


def calculate_tg(verbose=False, dilled=False):
    if dilled:
        with open(constants.GRIDPATH + 'dilled/tg_interpolation', 'rb') as file:
            return dill.load(file)
    if verbose:
        print('Creating T_g grid interpolation')
    nmchi, T = observations.temperature
    nmchi = nmchi/10.
    log_T = np.log10(T)
    if constants.interpolation == 'linear':
        return interpolate.LinearNDInterpolator(nmchi.values, log_T.Tg.values)
    elif constants.interpolation == 'cubic' or constants.interpolation == 'radial':
        return interpolate.Rbf(nmchi.values, log_T.Tg.values)
    else:
        sys.exit('<<ERROR>>: There is no such method as ' + 
                 '{} to interpolate '.format(constants.interpolation) +
                 'the gas temperature in the KOSMA-tau grid.\n\nExitting...\n\n')


def calculate_td(verbose=False, dilled=False):
    if dilled:
        with open(constants.GRIDPATH + 'dilled/td_interpolation', 'rb') as file:
            return dill.load(file)
    if verbose:
        print('Creating T_d grid interpolation')
    nmchi, T = observations.temperature
    nmchi = nmchi/10.
    log_T = np.log10(T)
    if constants.interpolation == 'linear':
        return interpolate.LinearNDInterpolator(nmchi.values, log_T.Td.values)
    elif constants.interpolation == 'cubic' or constants.interpolation == 'radial':
        return interpolate.Rbf(nmchi.values, log_T.Td.values)
    else:
        sys.exit('<<ERROR>>: There is no such method as ' + 
                 '{} to interpolate '.format(constants.interpolation) +
                 'the dust temperature in the KOSMA-tau grid.\n\nExitting...\n\n')


# Model interpolation


def calculate_galaxy_rotation(verbose=False, dilled=False):
    if constants.directory != '':
        if dilled:
            with open(constants.INPUTPATH+constants.directory + 
                      'dilled/rotation_interpolation', 'rb') as file:
                return dill.load(file)
        if verbose:
            print('Creating rotation velocity interpolation')
        rotation = observations.galaxy_rotation_profile
        if constants.interpolation == 'linear':
            return interpolate.interp1d(rotation[0], rotation[1][:, 0], kind='linear')
        elif constants.interpolation == 'cubic' or constants.interpolation == 'radial':
            return interpolate.interp1d(rotation[0], rotation[1][:, 0], kind='cubic')
        else:
            sys.exit('<<ERROR>>: There is no such method as ' + 
                     '{} to interpolate '.format(constants.interpolation) +
                     'the velocity profile.\n\nExitting...\n\n')
    return


def calculate_velocity_dispersion(verbose=False, dilled=False):
    if constants.directory != '':
        if dilled:
            with open(constants.INPUTPATH+constants.directory + 
                      'dilled/dispersion_interpolation', 'rb') as file:
                return dill.load(file)
        if verbose:
            print('Creating velocity dispersion interpolation')
        rotation = observations.galaxy_rotation_profile
        if constants.interpolation == 'linear':
            return interpolate.interp1d(rotation[0], rotation[1][:, 1], kind='linear')
        elif constants.interpolation == 'cubic' or constants.interpolation == 'radial':
            return interpolate.interp1d(rotation[0], rotation[1][:, 1], kind='cubic')
        else:
            sys.exit('<<ERROR>>: There is no such method as ' + 
                     '{} to interpolate '.format(constants.interpolation) +
                     'the velocity profile.\n\nExitting...\n\n')
    return


def calculate_number_density(verbose=False, dilled=False):
    if constants.directory != '':
        if dilled:
            with open(constants.INPUTPATH+constants.directory + 
                      'dilled/density_interpolation', 'rb') as file:
                return dill.load(file)
        if verbose:
            print('Creating density interpolation')
        density = observations.number_density_profile
        if constants.interpolation == 'linear':
            return interpolate.interp1d(density[0], density[1], kind='linear')
        elif constants.interpolation == 'cubic' or constants.interpolation == 'radial':
            return interpolate.interp1d(density[0], density[1], kind='cubic')
        else:
            sys.exit('<<ERROR>>: There is no such method as ' + 
                     '{} to interpolate '.format(constants.interpolation) +
                     'the KOSMA-tau grid.\n\nExitting...\n\n')
    return


def calculate_h2_height():
        r_h, h_h2 = observations.h2_scale_height_profile
        return interpolate.interp1d(r_h.values, h_h2.values, kind='linear', 
                                      fill_value='extrapolate')


def calculate_h2_mass(verbose=False, dilled=False):
    if constants.directory != '':
        if dilled:
            with open(constants.INPUTPATH+constants.directory + 
                      'dilled/h2_mass_interpolation', 'rb') as file:
                return dill.load(file)
        if verbose:
            print('Creating clump mass interpolation')
        r_s, sigma_h2 = observations.h2_surface_mass_profile
        f_sigma = interpolate.interp1d(r_s.values, sigma_h2.values, kind='linear', 
                                       bounds_error=False, 
                                       fill_value=(sigma_h2.values[0], sigma_h2.values[-1]))
        r_new = np.linspace(0, constants.rgal*1.1, num=100)
        half_height = np.sqrt(2*np.log(2))*interpolations.h2_scale_height(r_new)
        rho_0 = f_sigma(r_new)/2/half_height
        i_partial = constants.voxel_size >= 2*half_height
        i_full = constants.voxel_size < 2*half_height
        n_z = np.ceil(half_height/constants.voxel_size*10).astype(int)
        r_mass = np.hstack(list(r_new for _ in range(n_z.max()+1)))
        z_mass = np.hstack(list(np.full_like(r_new, _*0.1*constants.voxel_size, dtype=float) 
                                for _ in range(n_z.max()+1)))
        m_mass = list(np.zeros_like(r_new, dtype=float) for _ in range(n_z.max()+1))
        for _ in range(n_z.max()+1):
            lim_l = (_*0.1-0.5)*constants.voxel_size
            lim_u = (_*0.1+0.5)*constants.voxel_size
            i_full = (lim_l>=-half_height) & (lim_u<=half_height)
            i_partial = (lim_l>=-half_height) & (lim_u>half_height)
            i_all = (lim_l<-half_height) & (lim_u>half_height)
            m_mass[_][i_full] = rho_0[i_full] * constants.voxel_size**3
            m_mass[_][i_partial] = rho_0[i_partial] * constants.voxel_size**2 * (half_height[i_partial]-lim_l)
            m_mass[_][i_all] = rho_0[i_all] * constants.voxel_size**2 * 2*half_height[i_all]
        interpolations.h2_mass_full = pd.DataFrame(data={'r':r_mass, 
                                                         'z':z_mass, 
                                                         'M':np.hstack(m_mass)})
        return interpolate.LinearNDInterpolator(np.array([interpolations.h2_mass_full.r, 
                                                          interpolations.h2_mass_full.z]).T,
                                                interpolations.h2_mass_full.M)
        # if constants.interpolation == 'linear':
        #     return interpolate.LinearNDInterpolator(np.array([clumpMass.R, clumpMass.Z]), 
        #                                             clumpMass.M_h2)
        # elif constants.interpolation == 'cubic' or constants.interpolation == 'radial':
        #     return interpolate.Rbf(clumpMass.R, clumpMass.Z, clumpMass.M_h2)
        # else:
        #     sys.exit('<<ERROR>>: There is no such method as ' + 
        #              '{} to interpolate '.format(constants.interpolation) +
        #              'the KOSMA-tau grid.\n\nExitting...\n\n')
    return


def calculate_hi_height(like_clumps=False):
    if like_clumps:
        r_h, h_hi = observations.h2_scale_height_profile
    else:
        r_h, h_hi = observations.hi_scale_height_profile
    return interpolate.interp1d(r_h.values, h_hi.values, kind='linear', 
                                  fill_value='extrapolate')


def calculate_hi_mass(like_clumps=False, verbose=False, dilled=True):
    if constants.directory != '':
        if dilled:
            with open(constants.INPUTPATH+constants.directory + 
                      'dilled/hi_mass_interpolation', 'rb') as file:
                return dill.load(file)
        if verbose:
            print('Creating interclump mass interpolation')
        if like_clumps:
            r_s, sigma_hi = observations.h2_surface_mass_profile
            #r_h, h_hi = observations.h2_scale_height_profile
        else:
            r_s, sigma_hi = observations.hi_surface_mass_profile
            #r_h, h_hi = observations.hi_scale_height_profile
        f_sigma = interpolate.interp1d(r_s.values, sigma_hi.values, kind='linear', 
                                       bounds_error=False, 
                                       fill_value=(sigma_hi.values[0], sigma_hi.values[-1]))
        #f_h = interpolate.interp1d(r_h.values, h_hi.values, kind='linear', 
        #                           fill_value='extrapolate')
        r_new = np.linspace(0, constants.rgal*1.1, num=100)
        half_height = np.sqrt(2*np.log(2))*interpolations.hi_scale_height(r_new)
        rho_0 = f_sigma(r_new)/2/half_height
        i_partial = constants.voxel_size >= 2*half_height
        i_full = constants.voxel_size < 2*half_height
        n_z = np.ceil(half_height/constants.voxel_size*10).astype(int)
        r_mass = np.hstack(list(r_new for _ in range(n_z.max()+1)))
        z_mass = np.hstack(list(np.full_like(r_new, _*0.1*constants.voxel_size, dtype=float) 
                                for _ in range(n_z.max()+1)))
        m_mass = list(np.zeros_like(r_new, dtype=float) for _ in range(n_z.max()+1))
        for _ in range(n_z.max()+1):
            lim_l = (_*0.1-0.5)*constants.voxel_size
            lim_u = (_*0.1+0.5)*constants.voxel_size
            i_full = (lim_l>=-half_height) & (lim_u<=half_height)
            i_partial = (lim_l>=-half_height) & (lim_u>half_height)
            i_all = (lim_l<-half_height) & (lim_u>half_height)
            m_mass[_][i_full] = rho_0[i_full] * constants.voxel_size**3
            m_mass[_][i_partial] = rho_0[i_partial] * constants.voxel_size**2 * (half_height[i_partial]-lim_l)
            m_mass[_][i_all] = rho_0[i_all] * constants.voxel_size**2 * 2*half_height[i_all]
        interpolations.hi_mass_full = pd.DataFrame(data={'r':r_mass, 
                                                         'z':z_mass, 
                                                         'M':np.hstack(m_mass)})
        return interpolate.LinearNDInterpolator(np.array([interpolations.hi_mass_full.r, 
                                                          interpolations.hi_mass_full.z]).T,
                                                interpolations.hi_mass_full.M)
        # if constants.interpolation == 'linear':
        #     return interpolate.LinearNDInterpolator(np.asarray([interclumpMass.R, interclumpMass.Z]), 
        #                                             interclumpMass.M_hi)
        # elif constants.interpolation == 'cubic' or constants.interpolation == 'radial':
        #     return interpolate.Rbf(interclumpMass.R, interclumpMass.Z, interclumpMass.M_hi)
        # else:
        #     sys.exit('<<ERROR>>: There is no such method as ' + 
        #              '{} to interpolate '.format(constants.interpolation) +
        #              'the KOSMA-tau grid.\n\nExitting...\n\n')
    return


def calculate_fuv(l_range=(912, 2066), average_fuv=False, 
                  verbose=False, dilled=False):
    if constants.directory != '':
        if dilled:
            with open(constants.INPUTPATH+constants.directory + 
                      'dilled/fuv_field_interpolation', 'rb') as file:
                return dill.load(file)
        if verbose:
            print('Creating FUV interpolation')
        fuv = observations.fuv_profile
        lam = np.array([912, 1350, 1500, 1650, 2000, 2200, 2500, 2800, 3650])
        wav = np.linspace(l_range[0], l_range[1], num=1000)
        i_gc = fuv[0][:, 0] < constants.fuv_r_gc
        fuv[1][i_gc] = fuv[1][i_gc] * constants.fuv_scale_gc
        fuv[2][i_gc, :] = fuv[2][i_gc, :] * constants.fuv_scale_gc
        f = interpolate.interp1d(lam, fuv[2], axis=1)
        if constants.interpolation == 'linear':
            if average_fuv:
                fuv_temp = copy(fuv[1])
                # fuv_temp[i_gc] = fuv_temp[i_gc] * constants.fuv_scale_gc
                return interpolate.LinearNDInterpolator(fuv[0], fuv_temp)
            else:
                fuv_temp = np.trapz(f(wav), wav)
                # fuv_temp[i_gc] = fuv_temp[i_gc] * constants.fuv_scale_gc
                return interpolate.LinearNDInterpolator(fuv[0], fuv_temp)
        elif constants.interpolation == 'cubic' or constants.interpolation == 'radial':
            fuv_temp = np.trapz(f(wav), wav)
            # fuv_temp[i_gc] = fuv_temp[i_gc] * constants.fuv_scale_gc
            return interpolate.Rbf(fuv[0][:, 0], fuv[0][:, 1], fuv_temp)
        else:
            sys.exit('<<ERROR>>: There is no such method as ' + 
                     '{} to interpolate '.format(constants.interpolation) +
                     'the KOSMA-tau grid.\n\nExitting...\n\n')
    return


def calculate_e_tilde_real():
    return interpolate.interp1d(observations.e_tilde_real[0], observations.e_tilde_real[1], 
                                kind='linear')


def calculate_e_tilde_imaginary():
    return interpolate.interp1d(observations.e_tilde_imaginary[0], observations.e_tilde_imaginary[1], 
                                kind='linear')


# jit_module(nopython=False)
