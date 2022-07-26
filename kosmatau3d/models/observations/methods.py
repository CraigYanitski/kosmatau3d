import warnings

import numpy as np
import pandas as pd

from numba import jit_module
from numba.core.errors import NumbaDeprecationWarning, NumbaPendingDeprecationWarning

from kosmatau3d.models import constants
from kosmatau3d.models import observations


warnings.simplefilter('ignore', category=NumbaDeprecationWarning)
warnings.simplefilter('ignore', category=NumbaPendingDeprecationWarning)


def initialise_grid(tau_grid_file='clump_tau_LineCenter.dat', 
                    tb_grid_file='clump_Tmb_LineCenter', 
                    tau_fuv_grid_file='RhoMassAFUV.dat'):
    # Grid
    tau_centerline(tau_grid_file)
    tb_centerline(tb_grid_file)
    rho_mass_taufuv(tau_fuv_grid_file)
    species_data()

    observations.grid_initialised = True
    
    return


def initialise_input(h2_surface_density_file='h2_surface-density.dat', 
                     hi_surface_density_file='hi_surface-density.dat', 
                     h2_scale_height_file='h2_scale-height.dat', 
                     hi_scale_height_file='hi_scale-height.dat', 
                     h_number_density_file='h_number-density.dat', 
                     fuv_file='galactic_FUV_complete.dat', 
                     velocity_file='rot_milki2018_14.dat'):
    # Model
    h2_surface_mass_profile(h2_surface_density_file)
    hi_surface_mass_profile(hi_surface_density_file)
    h2_scale_height_profile(h2_scale_height_file)
    hi_scale_height_profile(hi_scale_height_file)
    number_density_profile(h_number_density_file)
    fuv_profile(fuv_file)
    galaxy_rotation_profile(velocity_file)

    observations.model_initialised = True

    return


# G R I D


def tb_centerline(file='clump_Tmb_LineCenter.dat'):
    # Open file for KOSMA-tau simulations of line intensities
    # FORMAT: n, M, UV, intensity[molecules then dust]
    header = []
    with open(constants.GRIDPATH+file) as tb:
        header.append(tb.readline())
        header.append(tb.readline())
    molecules = header[1].split(': ')[1]
    species = []
    for molecule in molecules.split(', '):
        for transition in np.arange(1, int(molecule.split(' ')[1])+1):
            species.append('{} {}'.format(molecule.split(' ')[0], transition))
    constants.setupMolecules(np.array(species))
    tb = np.genfromtxt(constants.GRIDPATH+file)
    observations.tb_centerline = (tb[:, :4], tb[:, 4:])
    return


def tau_centerline(file='clump_tau_LineCenter.dat'):
    # Open file for KOSMA-tau simulations of optical depths
    # FORMAT: n, M, UV, tau[molecules then dust]
    tau = np.genfromtxt(constants.GRIDPATH+file)
    observations.tau_centerline = (tau[:, :4], tau[:, 4:])
    return


def rho_mass_taufuv(file='RhoMassAFUV.dat'):
    pma = np.genfromtxt(constants.GRIDPATH+file)
    observations.rho_mass_taufuv = (pma[:, :2], pma[:, 2])
    return


def species_data(file='frequencies.dat'):
    # Open file containing the transition frequencies of various elements
    frequencies = np.genfromtxt(constants.GRIDPATH+file, names=['number', 'species', 'transition', 'frequency'],
                                dtype="i8,U8,i8,f8", delimiter=',')
    observations.species_data = (frequencies['number'], frequencies['species'],
                                 frequencies['transition'], frequencies['frequency'])
    return


def init_rad_transfer():
    eTildeReal()
    eTildeImaginary()
    return


def e_tilde_real(file='Ereal.dat'):
    eReal = np.genfromtxt(constants.GRIDPATH+file, names=['x', 'Ereal'])
    observations.e_tilde_real = (eReal['x'], eReal['Ereal'])
    return


def e_tilde_imaginary(file='Eimag.dat'):
    eImaginary = np.genfromtxt(constants.GRIDPATH+file, names=['x', 'Eimaginary'])
    observations.e_tilde_imaginary = (eImaginary['x'], eImaginary['Eimaginary'])
    return


# M O D E L


def h2_surface_mass_profile(file='h2_surface-density.dat'):
    # Open file for the mass profile (clump) of the object (Msol/pc**2)
    if constants.directory != '':
        sigma_h2 = pd.read_csv(constants.INPUTPATH+constants.directory+file, 
                               delim_whitespace=True, skiprows=1)
        # if constants.voxel_size > 2*constants.scale_height_h2:
        #     r = np.hstack((clumpMass['r'] * 1000, clumpMass['r'] * 1000))
        #     z = np.hstack((np.zeros_like(clumpMass['r'] * 1000, dtype=float), 
        #                    np.full_like(clumpMass['r'] * 1000, constants.voxel_size, dtype=float)))
        #     m = np.hstack((clumpMass['n_h2']*2*constants.scale_height_h2,
        #                    clumpMass['n_h2']*0))
        # else:
        #     n_vox = int(np.ceil((constants.scale_height_h2-constants.voxel_size/2) 
        #                         / constants.voxel_size))
        #     r = np.hstack(list(clumpMass['r']*1000 for _ in range(n_vox+1)))
        #     z = np.hstack(list(np.full_like(clumpMass['r'], _*constants.voxel_size, dtype=float) 
        #                        for _ in range(n_vox+1)))
        #     m = np.hstack(list(clumpMass['n_h2']*constants.voxel_size if 
        #                        (_+0.5)*constants.voxel_size<constants.scale_height_h2 
        #                        else 0 if (_-0.5)*constants.voxel_size>constants.scale_height_h2 else 
        #                        clumpMass['n_h2']*(constants.scale_height_h2%constants.voxel_size) 
        #                        for _ in range(n_vox+1)))
        observations.h2_surface_mass_profile = (sigma_h2['r']*1000, sigma_h2['Sigma_h2'])
        # observations.h2_surface_mass_profile = pd.DataFrame(data={'R':r, 'Z':z, 
        #                                                     'M_h2':m*constants.voxel_size**2})
        #observations.clumpMassProfile = (clumpMass['radius']*1000, clumpMass['h2_mass_density'] *
        #                                 constants.voxel_size**2*constants.scale_height)
    return


def hi_surface_mass_profile(file='hi_surface-density.dat'):
    # Open file for the mass profile (interclump) of the object (Msol/pc**2)
    if constants.directory != '':
        sigma_hi = pd.read_csv(constants.INPUTPATH+constants.directory+file, 
                               delim_whitespace=True, skiprows=1)
        # if constants.voxel_size > 2*constants.scale_height_hi:
        #     r = np.hstack((interclumpMass['r'] * 1000, interclumpMass['r'] * 1000))
        #     z = np.hstack((np.zeros_like(interclumpMass['r'] * 1000, dtype=float), 
        #                    np.full_like(interclumpMass['r'] * 1000, constants.voxel_size, dtype=float)))
        #     m = np.hstack((interclumpMass['n_hi']*2*constants.scale_height_hi,
        #                    interclumpMass['n_hi']*0))
        # else:
        #     n_vox = int(np.ceil((constants.scale_height_hi-constants.voxel_size/2) 
        #                         / constants.voxel_size))
        #     r = np.hstack(list(interclumpMass['r']*1000 for _ in range(n_vox+1)))
        #     z = np.hstack(list(np.full_like(interclumpMass['r'], _*constants.voxel_size, dtype=float) 
        #                        for _ in range(n_vox+1)))
        #     m = np.hstack(list(interclumpMass['n_hi']*constants.voxel_size if 
        #                        (_+0.5)*constants.voxel_size<constants.scale_height_hi 
        #                        else 0 if (_-0.5)*constants.voxel_size>constants.scale_height_hi else 
        #                        interclumpMass['n_hi']*(constants.scale_height_hi%constants.voxel_size) 
        #                        for _ in range(n_vox+1)))
        observations.hi_surface_mass_profile = (sigma_hi['r']*1000, sigma_hi['Sigma_hi'])
        # observations.hi_surface_mass_profile = pd.DataFrame(data={'R':r, 'Z':z, 
        #                                                           'M_hi':m*constants.voxel_size**2})
        #observations.interclumpMassProfile = (interclumpMass['radius']*1000, 
        #                                      interclumpMass['hi_mass_density'] *
        #                                      constants.voxel_size**2*constants.scale_height)
    return


def h2_scale_height_profile(file='h2_scale-height.dat'):
    # Open file for the H2 scale height profile of the galaxy (pc).
    if constants.directory != '':
        height = pd.read_csv(constants.INPUTPATH+constants.directory+file, 
                             delim_whitespace=True, skiprows=1)
        observations.h2_scale_height_profile = (height['r']*1000, height['h_h2'])
    return


def hi_scale_height_profile(file='hi_scale-height.dat'):
    # Open file for the HI scale height profile of the galaxy (pc).
    if constants.directory != '':
        height = pd.read_csv(constants.INPUTPATH+constants.directory+file, 
                             delim_whitespace=True, skiprows=1)
        observations.hi_scale_height_profile = (height['r']*1000, height['h_hi'])
    return


def number_density_profile(file='number-densities.dat'):
    # Open file for the number density profile of the object (n/cm**3)
    if constants.directory != '':
        density = np.genfromtxt(constants.INPUTPATH+constants.directory+file, 
                                names=['radius', 'h2_density'])
        observations.number_density_profile = ((density['radius']*1000), density['h2_density'])
    return


def fuv_profile(file='galactic_FUV_complete.dat'):
    '''
    Open file for the FUV profile of the object
    'radius', 'energy density 912', 'energy density 1350', 'energy density 1500', 'energy density 1650',
    'energy density 2000', 'energy density 2200', 'energy density 2500', 'energy density 2800', 'energy density 3650'])
    '''
    if constants.directory != '':
        fuv = np.loadtxt(constants.INPUTPATH+constants.directory+file)
        # r = np.arange(0,20000,50)
        # fuv = 50/4/np.pi/r**2
        observations.fuv_profile = (fuv[:, :2], fuv[:, 2] / (constants.pc*100)**3/constants.u_draine0,
                                    fuv[:, 3:] / (constants.pc*100)**3/constants.u_draine0)
        # return (r,fuv)
    return


def galaxy_rotation_profile(file='rot_milki2018_14.dat'):
    '''
    Open file for the rotation profile of the object. There is `rot_milki2018_14.dat`, which is taken from
    Bhattacharjee et al. (2014) and Eilers et al. (2018), and there is `rot_curve.dat`, for which I cannot find a
    reference. Christoph's latest version of KOSMA-tau-3D seemed to be using the former data file, but I will test
    both. The latter file has higher rotational velocities, though it seems somewhat artificially created (multiple
    entries with the same value).
    '''
    if constants.directory != '':
        rotation = np.genfromtxt(constants.INPUTPATH+constants.directory+file)
        observations.galaxy_rotation_profile = (rotation[:, 0]*1000., rotation[:, 1:])
    return


# jit_module(nopython=False)
