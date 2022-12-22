import sys
import warnings
import multiprocessing
from pprint import pprint
from multiprocessing import Pool
from functools import partial
import importlib as il
from time import time
import copy as c
import cmath

import numpy as np
from numba import jit_module
import matplotlib.pyplot as plt
# import matplotlib as mpl
# mpl.use('Qt4Agg')
import scipy.interpolate as interpolate
import scipy.optimize as op
from scipy.special import erfi, erfc
from astropy.io import fits
from tqdm import tqdm

from kosmatau3d.models import (
    constants,
    species,
    interpolations,
    observations,
    radiativeTransfer as rt,
    )


def eTildeReal(file='Ereal.dat'):
    ereal = np.genfromtxt(constants.GRIDPATH+file, names=['x', 'Ereal'])
    return (ereal['x'], ereal['Ereal'])


def eTildeImaginary(file='Eimag.dat'):
    eimaginary = np.genfromtxt(constants.GRIDPATH+file, names=['x', 'Eimaginary'])
    return (eimaginary['x'], eimaginary['Eimaginary'])


def open_voxel_positions(directory):
    '''
    Open the file containing the Cartesian positions of all voxels.
    This will save the contents in a sub-module variable.

    :param directory: path to the model files
    :return:

    shape -> (n_voxels, 3)
    '''

    rt.voxel_positions = fits.open(directory + 'voxel_position.fits', mode='denywrite')
    return


def open_voxel_velocities(directory):
    '''
    Open the file containing the observed velocity of all voxels.
    This will save the contents in a sub-module variable.

    :param directory: path to the model files
    :return:

    shape -> (n_voxels, 1)
    '''

    rt.voxel_velocities = fits.open(directory+'voxel_velocity.fits', mode='denywrite')
    return


def open_voxel_filling_factors(directory):
    '''
    Open the file containing the voxel-filling factors of all voxels.
    This will save the contents in a sub-module variable.

    :param directory: path to the model files
    :return:

    shape -> (n_voxels, 1)
    '''

    voxel_filling_factors = fits.open(directory+'voxel-filling_factor.fits', 
                                         mode='denywrite')[0].data.mean(1)
    rt.voxel_filling_factors = np.ones_like(voxel_filling_factors)
    return


def open_voxel_emission(directory, hi=False, verbose=False):
    '''
    Open the file containing the emissivity and absorption of all voxels.
    This will save the contents in a sub-module variable.

    :param directory: path to the model files
    :return:

    Emissivities in K/pc; absorptions in 1/pc

    species file shape -> (n_voxels, n_v_obs, n_transitions)

    dust file shape -> (n_voxels, n_wavelengths)
    '''
    if hi:
        species = 'hi'
    else:
        species = 'species'
    dust = 'dust'
    rt.temp_species_emissivity = fits.open(directory+f'{species}_emissivity.fits', mode='denywrite')
    rt.temp_species_absorption = fits.open(directory+f'{species}_absorption.fits', mode='denywrite')
    rt.temp_dust_emissivity = fits.open(directory+f'{dust}_emissivity.fits', mode='denywrite')
    rt.temp_dust_absorption = fits.open(directory+f'{dust}_absorption.fits', mode='denywrite')

    if verbose:
        print(rt.temp_species_emissivity[0].data.mean(0).mean(0))
        print(rt.temp_dust_emissivity[0].data.mean(0))

    return


def calculateObservation(directory='', dim='xy', terminal=True, hi=False, vel_pool=False, plotV=False, 
                         slRange=[(-np.pi,np.pi), (-np.pi/2,np.pi/2)], nsl=[50,25],
                         pencil_beam=True, d_lon=None, d_lat=None, multiprocessing=0, 
                         debug=False, verbose=False):

    if debug:
        sl = [5, 5]

    rt.pencil_beam = pencil_beam
    if d_lon:
        rt.d_lon = d_lon
    else:
        rt.d_lon = abs(slRange[0][0]-slRange[0][1])/(nsl[0]-1)
    if d_lat:
        rt.d_lat = d_lat
    else:
        rt.d_lat = abs(slRange[1][0]-slRange[1][1])/(nsl[1]-1)
    
    # constants.velocityRange = np.linspace(-300, 300, 500)
    
    # print('Load data')

    rt.open_voxel_filling_factors(directory)
    rt.open_voxel_positions(directory)
    rt.open_voxel_velocities(directory)
    rt.open_voxel_emission(directory, hi=hi, verbose=verbose)
  
    nDust = rt.temp_dust_emissivity[0].shape[1]
    if nDust > 10:
        # default interpolation is along the last axis, which is what we need
        rt.interpDust = interpolate.interp1d(constants.wavelengths[:nDust],
                                             rt.temp_dust_emissivity[0].data,
                                             fill_value='extrapolate')
    
    constants.velocity_range = np.linspace(rt.temp_species_emissivity[0].header['CRVAL2'] -
                                           (rt.temp_species_emissivity[0].header['CRPIX2']) *
                                           rt.temp_species_emissivity[0].header['CDELT2'],
                                           rt.temp_species_emissivity[0].header['CRVAL2'] +
                                           (rt.temp_species_emissivity[0].header['CRPIX2']) *
                                           rt.temp_species_emissivity[0].header['CDELT2'],
                                           num=rt.temp_species_emissivity[0].header['NAXIS2'])
    
    # observations.methods.initialise_grid()
    # species.addMolecules(rt.tempSpeciesEmissivity[0].header['SPECIES'].split(', '))
    # print(rt.tempSpeciesEmissivity[0].header['SPECIES'].split(', '), '\n', species.moleculeWavelengths)
    
    # print('Data loaded :-)')
  
    xPositions, yPositions, zPositions = rt.voxel_positions[0].data[:, 0], \
        rt.voxel_positions[0].data[:, 1], \
        rt.voxel_positions[0].data[:, 2]
    r = np.sqrt((xPositions-constants.rgal_earth)**2 + yPositions**2)
  
    radGrid = np.sqrt((xPositions-constants.rgal_earth)**2 + yPositions**2 + zPositions**2)
    lonGrid = np.arctan2(yPositions, -(xPositions-constants.rgal_earth))
    rPolar = np.sqrt((xPositions-constants.rgal_earth)**2+yPositions**2)
    latGrid = np.arctan2(zPositions, rPolar)
  
    np.set_printoptions(threshold=100000)
    # print('\nLongitude\n', lonGrid, '\nLattitude\n', latGrid)
  
    # print('\nx\n', np.unique(xArray), '\ny\n', np.unique(yArray), '\nz\n', np.unique(zArray))
  
    if constants.from_earth:
        # For an observation from Earth, the data is modified to account for Earth's position at (8750, 0, 0) pc.
        #  The integrated intensity is then calculated in the y-z plane to account for different viewing angles across
        #  the galaxy. This can be post-processed to convert to galactic coordinates.
    
        # Define the boundaries separating the inner and outer disk
        xBoundary = (xPositions > 0) & (r > constants.rgal_earth)
        yBoundary = (yPositions < constants.rgal_earth) & (yPositions > -constants.rgal_earth)
    
        # np.set_printoptions(threshold=1000000)
    
        # Define sightlines calculated
        longrid = np.linspace(-np.pi, np.pi, num=nsl[0])
        latgrid = np.linspace(-np.pi/2, np.pi/2, num=nsl[1])
        # grid = np.meshgrid(lon, lat)
        # grid = np.array([grid[0].flatten(), grid[1].flatten()])
    
        # rt.vTqdm = tqdm(total=constants.velocityRange.size, desc='Observing velocity',
        #                                miniters=1, dynamic_ncols=True)
        # if terminal: rt.slTqdm = tqdm(total=longrid.size*latgrid.size, desc='Sightline',
        #                                miniters=1, dynamic_ncols=True)
    
        VintensityMapSpecies = []
        VintensityMapDust = []
        Vpositions = []
    
        rt.sightlines = np.zeros((longrid.size, latgrid.size))
    
        if debug:
            velocityRange = [0]
        else:
            velocityRange = constants.velocity_range
    
        result = multiprocessCalculation(slRange=slRange, nsl=nsl, multiprocessing=multiprocessing,
                                         vel_pool=vel_pool, dim=dim, debug=debug)
        
        result_positions = result[0]
        result_intensity_species = result[1][0]
        result_intensity_dust = result[2][0]
        result_optical_depth_species = result[1][1]
        result_optical_depth_dust = result[2][1]
        rt.sightlines = np.asarray(result[3]).max(0)
        vmin, vmax = result[4]
        
        # Save sightline lengths
        np.savetxt(directory+'/sightlines.csv', rt.sightlines, delimiter=',')
        
        # Convert to numpy arrays
        result_positions = np.asarray(result_positions)
        result_intensity_species = np.asarray(result_intensity_species)
        result_intensity_dust = np.asarray(result_intensity_dust)
        result_optical_depth_species = np.asarray(result_optical_depth_species)
        result_optical_depth_dust = np.asarray(result_optical_depth_dust)
    
        if verbose:
            rt.logger.info('Map position shape:')
            if result_positions.ndim > 1:
                rt.logger.info(result_positions.shape)
            else:
                for p in result_positions:
                    rt.logger.info(p.shape)
              
            rt.logger.info('Map intensity shapes:')
            rt.logger.info('Species')
            if result_intensity_species.ndim > 1:
                rt.logger.info(result_intensity_species.shape)
            else:
                for intensity in result_intensity_species:
                    rt.logger.info(intensity.shape)
            rt.logger.info('Dust')
            if result_intensity_dust.ndim > 1:
                rt.logger.info(result_intensity_dust.shape)
            else:
                for intensity in result_intensity_dust:
                    rt.logger.info(intensity.shape)
    
        # Setup the data to be saved in a FITS file. It will create an HDU list with position, species, and dust HDUs.
        if not debug:
    
            hdul_intensity = fits.HDUList()
            hdul_optical_depth = fits.HDUList()
            
            # Create HDUs for the map position and intensity and add the velocity in the headers
            position_hdu = fits.ImageHDU(result_positions)
            position_hdu.header['TYPE1'] = 'Angle'
            position_hdu.header['TYPE2'] = 'Position'
            position_hdu.header['DIREC'] = 'Radial'
            
            # print(VintensityMapSpecies.shape, np.shape(VintensityMapSpecies[0, 0, 0]))
            if verbose:
                print(np.shape(result_intensity_species), np.shape(result_intensity_dust))
            rt.intensity_species = result_intensity_species
            intensity_hdu_species = fits.ImageHDU(result_intensity_species)
            intensity_hdu_species = edit_hdu_header(intensity_hdu_species, sl_range=slRange, 
                                                    species=True, intensity=True,
                                                    verbose=verbose)
            
            # print(VintensityMapDust.shape, np.shape(VintensityMapDust[0, 0, 0]))
            rt.intensity_dust = result_intensity_dust
            intensity_hdu_dust = fits.ImageHDU(result_intensity_dust)
            intensity_hdu_dust = edit_hdu_header(intensity_hdu_dust, sl_range=slRange, 
                                                 dust=True, intensity=True,
                                                 verbose=verbose)
            
            # print(VintensityMapSpecies.shape, np.shape(VintensityMapSpecies[0, 0, 0]))
            if verbose:
                print(np.shape(result_optical_depth_species), np.shape(result_optical_depth_dust))
            rt.optical_depth_species = result_optical_depth_species
            optical_depth_hdu_species = fits.ImageHDU(result_optical_depth_species)
            optical_depth_hdu_species = edit_hdu_header(optical_depth_hdu_species, sl_range=slRange, 
                                                        species=True, optical_depth=True, 
                                                        verbose=verbose)
            
            # print(VintensityMapDust.shape, np.shape(VintensityMapDust[0, 0, 0]))
            rt.optical_depth_dust = result_optical_depth_dust
            optical_depth_hdu_dust= fits.ImageHDU(result_optical_depth_dust)
            optical_depth_hdu_dust = edit_hdu_header(optical_depth_hdu_dust, sl_range=slRange, 
                                                     dust=True, optical_depth=True,
                                                     verbose=verbose)
      
      
            # IntensityHDUSpecies.header['TYPE'] = 'Species transitions'
            # IntensityHDUSpecies.header['BUNIT'] = 'K'
            # IntensityHDUSpecies.header['CTYPE1'] = 'Wavelength'
            # IntensityHDUSpecies.header['CUNIT1'] = 'm'
            # IntensityHDUSpecies.header['CRVAL1'] = 'N/A'
            # IntensityHDUSpecies.header['CDELT1'] = 'N/A'
            # IntensityHDUSpecies.header['CRPIX1'] = 'N/A'
            # IntensityHDUSpecies.header['CTYPE2'] = 'GLON'
            # IntensityHDUSpecies.header['CUNIT2'] = 'rad'
            # if IntensityHDUSpecies.header['NAXIS2'] > 1:
            #     IntensityHDUSpecies.header['CRVAL2'] = (slRange[0][1]+slRange[0][0])/2.
            #     IntensityHDUSpecies.header['CDELT2'] = ((slRange[0][1]-slRange[0][0]) /
            #                                             (IntensityHDUSpecies.header['NAXIS2']-1))
            #     IntensityHDUSpecies.header['CRPIX2'] = (IntensityHDUSpecies.header['NAXIS2'])/2.
            # else:
            #     IntensityHDUSpecies.header['CRVAL2'] = slRange[0][0]
            #     IntensityHDUSpecies.header['CDELT2'] = 0
            #     IntensityHDUSpecies.header['CRPIX2'] = (IntensityHDUSpecies.header['NAXIS2'])
            # IntensityHDUSpecies.header['CTYPE3'] = 'GLAT'
            # IntensityHDUSpecies.header['CUNIT3'] = 'rad'
            # if IntensityHDUSpecies.header['NAXIS3'] > 1:
            #     IntensityHDUSpecies.header['CRVAL3'] = (slRange[1][1]+slRange[1][0])/2.
            #     IntensityHDUSpecies.header['CDELT3'] = ((slRange[1][1]-slRange[1][0]) /
            #                                             (IntensityHDUSpecies.header['NAXIS3']-1))
            #     IntensityHDUSpecies.header['CRPIX3'] = (IntensityHDUSpecies.header['NAXIS3'])/2.
            # else:
            #     IntensityHDUSpecies.header['CRVAL3'] = slRange[1][0]
            #     IntensityHDUSpecies.header['CDELT3'] = 0
            #     IntensityHDUSpecies.header['CRPIX3'] = (IntensityHDUSpecies.header['NAXIS3'])
            # IntensityHDUSpecies.header['CTYPE4'] = 'Velocity'
            # IntensityHDUSpecies.header['CUNIT4'] = 'km/s'
            # IntensityHDUSpecies.header['CRVAL4'] = (vmax+vmin)/2.
            # IntensityHDUSpecies.header['CDELT4'] = (vmax-vmin)/(IntensityHDUSpecies.header['NAXIS4']-1)
            # IntensityHDUSpecies.header['CRPIX4'] = (IntensityHDUSpecies.header['NAXIS4'])/2.
            # IntensityHDUSpecies.header['DIREC'] = 'Radial'
            # IntensityHDUSpecies.header['SPECIES'] = rt.tempSpeciesEmissivity[0].header['SPECIES']
      
            # IntensityHDUDust.header['TYPE'] = 'Dust continuum'
            # IntensityHDUDust.header['BUNIT'] = 'K'
            # IntensityHDUDust.header['CTYPE1'] = 'Wavelength'
            # IntensityHDUDust.header['CUNIT1'] = 'm'
            # IntensityHDUDust.header['CRVAL1'] = 'N/A'
            # IntensityHDUDust.header['CDELT1'] = 'N/A'
            # IntensityHDUDust.header['CRPIX1'] = 'N/A'
            # IntensityHDUDust.header['CTYPE2'] = 'GLON'
            # IntensityHDUDust.header['CUNIT2'] = 'rad'
            # if IntensityHDUDust.header['NAXIS2'] > 1:
            #     IntensityHDUDust.header['CRVAL2'] = (slRange[0][1]+slRange[0][0])/2
            #     IntensityHDUDust.header['CDELT2'] = ((slRange[0][1]-slRange[0][0]) /
            #                                          (IntensityHDUDust.header['NAXIS2']-1))
            #     IntensityHDUDust.header['CRPIX2'] = (IntensityHDUDust.header['NAXIS2'])/2.
            # else:
            #     IntensityHDUDust.header['CRVAL2'] = slRange[0][0]
            #     IntensityHDUDust.header['CDELT2'] = 0
            #     IntensityHDUDust.header['CRPIX2'] = (IntensityHDUDust.header['NAXIS2'])
            # IntensityHDUDust.header['CTYPE3'] = 'GLAT'
            # IntensityHDUDust.header['CUNIT3'] = 'rad'
            # if IntensityHDUDust.header['NAXIS3'] > 1:
            #     IntensityHDUDust.header['CRVAL3'] = (slRange[1][1]+slRange[1][0])/2
            #     IntensityHDUDust.header['CDELT3'] = ((slRange[1][1]-slRange[1][0]) /
            #                                          (IntensityHDUDust.header['NAXIS3']-1))
            #     IntensityHDUDust.header['CRPIX3'] = (IntensityHDUDust.header['NAXIS3'])/2.
            # else:
            #     IntensityHDUDust.header['CRVAL3'] = slRange[1][0]
            #     IntensityHDUDust.header['CDELT3'] = 0
            #     IntensityHDUDust.header['CRPIX3'] = (IntensityHDUDust.header['NAXIS3'])
            # IntensityHDUDust.header['DIREC'] = 'Radial'
            # IntensityHDUDust.header['DUST'] = rt.tempDustEmissivity[0].header['DUST']

            if hi:
                file_id = '_hi_'
            else:
                file_id = '_'
      
            hdul_intensity.append(c.deepcopy(position_hdu))
            hdul_intensity.append(intensity_hdu_species)
            hdul_intensity.append(intensity_hdu_dust)
            hdul_intensity.writeto(directory+f'/synthetic{file_id}intensity.fits', overwrite=True)
      
            hdul_optical_depth.append(c.deepcopy(position_hdu))
            hdul_optical_depth.append(optical_depth_hdu_species)
            hdul_optical_depth.append(optical_depth_hdu_dust)
            hdul_optical_depth.writeto(directory+f'/synthetic{file_id}optical_depth.fits', overwrite=True)
      
            print('Intensity map written successfully :-)')
    
        rt.voxel_positions.close()
        rt.voxel_velocities.close()
        rt.temp_species_emissivity.close()
        rt.temp_species_absorption.close()
        rt.temp_dust_emissivity.close()
        rt.temp_dust_absorption.close()
    
    return


def edit_hdu_header(hdu, sl_range=(0, 0), species=False, dust=False, 
                    intensity=False, optical_depth=False, verbose=False):
    if verbose:
        pprint(hdu.header)
    if species:
        hdu.header['TYPE'] = 'Species transitions'
    elif dust:
        hdu.header['TYPE'] = 'Dust continuum'
    if intensity:
        hdu.header['BUNIT'] = 'K'
    elif optical_depth:
        hdu.header['BUNIT'] = 'N/A'
    hdu.header['DIREC'] = 'Radial'
    hdu.header['CTYPE1'] = 'Wavelength'
    hdu.header['CUNIT1'] = 'm'
    hdu.header['CRVAL1'] = 'N/A'
    hdu.header['CDELT1'] = 'N/A'
    hdu.header['CRPIX1'] = 'N/A'
    hdu.header['CTYPE2'] = 'GLON'
    hdu.header['CUNIT2'] = 'rad'
    if hdu.header['NAXIS2'] > 1:
        hdu.header['CRVAL2'] = (sl_range[0][1]+sl_range[0][0])/2.
        hdu.header['CDELT2'] = ((sl_range[0][1]-sl_range[0][0]) / (hdu.header['NAXIS2']-1))
        hdu.header['CRPIX2'] = (hdu.header['NAXIS2'])/2.
    else:
        hdu.header['CRVAL2'] = sl_range[0][0]
        hdu.header['CDELT2'] = 0
        hdu.header['CRPIX2'] = (hdu.header['NAXIS2'])
    hdu.header['CTYPE3'] = 'GLAT'
    hdu.header['CUNIT3'] = 'rad'
    if hdu.header['NAXIS3'] > 1:
        hdu.header['CRVAL3'] = (sl_range[1][1]+sl_range[1][0])/2.
        hdu.header['CDELT3'] = ((sl_range[1][1]-sl_range[1][0]) / (hdu.header['NAXIS3']-1))
        hdu.header['CRPIX3'] = (hdu.header['NAXIS3'])/2.
    else:
        hdu.header['CRVAL3'] = sl_range[1][0]
        hdu.header['CDELT3'] = 0
        hdu.header['CRPIX3'] = (hdu.header['NAXIS3'])
    if species:
        hdu.header['CTYPE4'] = 'Velocity'
        hdu.header['CUNIT4'] = 'km/s'
        hdu.header['CRVAL4'] = rt.temp_species_emissivity[0].header['CRVAL2']#(vmax+vmin)/2.
        hdu.header['CDELT4'] = rt.temp_species_emissivity[0].header['CDELT2']#(vmax-vmin)/(hdu.header['NAXIS4']-1)
        hdu.header['CRPIX4'] = rt.temp_species_emissivity[0].header['CRPIX2']#(hdu.header['NAXIS4'])/2.
        hdu.header['SPECIES'] = rt.temp_species_emissivity[0].header['SPECIES']
    elif dust:
        hdu.header['DUST'] = rt.temp_dust_emissivity[0].header['DUST']
    return hdu


def multiprocessCalculation(slRange=[(-np.pi,np.pi), (-np.pi/2,np.pi/2)], nsl=[50,25], multiprocessing=0,
                            dim='spherical', vel_pool=False, debug=False):


    if vel_pool:
        v_positions = []
        v_intensity_map_species = []
        v_intensity_map_dust = []
        v_opticaldepth_map_species = []
        v_opticaldepth_map_dust = []
        sightlines = []
        args = (list(enumerate(constants.velocityRange)), constants.velocityRange.size)
        velChannel = partial(calculateVelocityChannel, slRange=slRange, nsl=nsl, dim=dim, debug=debug,
                             multiprocess=multiprocessing, vel_pool=vel_pool)
    else:
        lon = np.linspace(slRange[0][0], slRange[0][1], num=nsl[0])
        lat = np.linspace(slRange[1][0], slRange[1][1], num=nsl[1])
        longrid, latgrid = (arr.flatten() for arr in np.meshgrid(lon, lat))
        positions = []
        intensity_map_species = np.zeros((rt.temp_species_emissivity[0].shape[1], lat.size, lon.size,
                                          rt.temp_species_emissivity[0].shape[2]))
        optical_depth_map_species = np.zeros((rt.temp_species_absorption[0].shape[1], lat.size, lon.size,
                                              rt.temp_species_absorption[0].shape[2]))
        intensity_map_dust = np.zeros((lat.size, lon.size, rt.temp_dust_emissivity[0].shape[1]))
        optical_depth_map_dust = np.zeros((lat.size, lon.size, rt.temp_dust_absorption[0].shape[1]))
        sightlines = np.zeros((lon.size, lat.size))
        args = (list(enumerate(zip(longrid, latgrid))), longrid.size)
        calc_los = partial(calculate_sightline, slRange=slRange, nsl=nsl, dim=dim, debug=debug,
                           multiprocess=multiprocessing)
    
    t0 = time()
    
    if multiprocessing:
        pool = Pool(processes=multiprocessing)
        chunksize = max(int(args[1]/multiprocessing/100), 1)
        if vel_pool:
            intensity = pool.imap(velChannel, args[0], chunksize)
        else:
            intensity = pool.imap(calc_los, args[0], chunksize)
    else:
        intensity = []
        if vel_pool:
            vTqdm = tqdm(total=constants.velocity_range.size, desc='Observing velocity', 
                         miniters=1, dynamic_ncols=True)
            rt.slTqdm = tqdm(total=nsl[0]*nsl[1], desc='Sightline', 
                             miniters=1, dynamic_ncols=True)
            for iv in args[0]:
                intensity.append(velChannel(iv))
                #print(np.shape(intensity[-1][1][0]))
                #input()
                vTqdm.update()
        else:
            rt.slTqdm = tqdm(total=args[1], desc='Sightline', 
                             miniters=1, dynamic_ncols=True)
            for arg in args[0]:
                intensity.append(calc_los(arg))
                #print(np.shape(intensity[-1][1][0]))
                #input()
                # print(intensity[-1][0][0], '-->', intensity[-1][3])
                rt.slTqdm.update()
            print('\n\n')
        
    vmin = constants.velocity_range.max()
    vmax = constants.velocity_range.min()
    
    if multiprocessing:
        if vel_pool:
            vTqdm = tqdm(total=args[1], desc='Observing velocity', 
                         miniters=1, dynamic_ncols=True)
        else:
            rt.slTqdm = tqdm(total=args[1], desc='Sightline', 
                             miniters=1, dynamic_ncols=True)

    # Loop through all processes
    for n, i in enumerate(intensity):
        #   i.wait()
        #   print(i)
        # # This allows for discarding the velocity channels without an emission (assuming vel_pool is True)
        if vel_pool:
            if constants.velocityRange[n] < vmin:
                vmin = constants.velocityRange[n]
            if constants.velocityRange[n] > vmax:
                vmax = constants.velocityRange[n]

            if vel_pool:
                v_positions.append(i[0])
                v_intensity_map_species.append(i[1][0])
                v_opticaldepth_map_species.append(i[1][1])
                v_intensity_map_dust.append(i[2][0])
                v_opticaldepth_map_dust.append(i[2][1])
                sightlines.append(i[3])
        else:
            positions.append(i[0][1:])
            # print(positions)
            # print(i[1][0])
            i_lon = np.where(lon == i[0][1])[0][0]
            i_lat = np.where(lat == i[0][2])[0][0]
            intensity_map_species[:, i_lat, i_lon, :] = i[1][0]
            optical_depth_map_species[:, i_lat, i_lon, :] = i[1][1]
            intensity_map_dust[i_lat, i_lon, :] = i[2][0]
            optical_depth_map_dust[i_lat, i_lon, :] = i[2][1]
            sightlines[i_lon, i_lat] = i[3]

            # intensity[i] = None
          
        if multiprocessing:
            if vel_pool:
                vTqdm.update()
            else:
                rt.slTqdm.update()
       
    if vel_pool:
        vTqdm.close()
    else:
        rt.slTqdm.close()

    # print('\n\nTotal evaluation time for {} sightlines and {} velocity channels: {}\n\n'.format(nsl[0]*nsl[1], vNum,
    #                                                                                             time()-t0))

    if vel_pool:
        return (v_positions[0], (v_intensity_map_species, v_opticaldepth_map_species), 
                (v_intensity_map_dust[0], v_opticaldepth_map_dust[0]), sightlines, (vmin, vmax))
    else:
        return (positions, (intensity_map_species, optical_depth_map_species),
                (intensity_map_dust, optical_depth_map_dust), sightlines, (vmin, vmax))


def sightlength(x, l):
    return constants.rgal_earth**2 - constants.rgal**2 + x**2 - 2*constants.rgal_earth*x*np.cos(l)

# for i_vel,velocity in enumerate(constants.velocityRange):


def calculateVelocityChannel(ivelocity, slRange=[(-np.pi,np.pi), (-np.pi/2,np.pi/2)], nsl=[50,25],
                             dim='spherical', debug=False, multiprocess=0, vel_pool=True):

    # Convert the tuple to the desired index and velocity
    i_vel = ivelocity[0]
    
    # if multiprocess:
    #   t0 = time()
    #   print('\nCalculating velocity channel at {:.2f} km\s'.format(ivelocity[1]))
    
    # Setup the sightlines that are calculated
    longrid = np.linspace(slRange[0][0], slRange[0][1], num=nsl[0])
    latgrid = np.linspace(slRange[1][0], slRange[1][1], num=nsl[1])
    sightlines = np.zeros((longrid.size, latgrid.size))
    
    # print('lon/lat arrays created:', time()-t0)
  
    # Find the voxels that exist at the observing velocity
    nDust = rt.temp_dust_emissivity[0].shape[1]
    if nDust > 10:
        base = rt.interpDust(species.moleculeWavelengths)[:, 0]
    else:
        base = rt.temp_dust_emissivity[0].data[0, 0]
      
    rt.i_vox = ((rt.temp_species_emissivity[0].data[:, i_vel, :] > base).any(1) |
                               rt.temp_dust_emissivity[0].data[:, :].any(1))
    
    # print('Voxels selected (shape={}):'.format(i_vox[i_vox].shape), time()-t0)
  
    # Update velocity progress bar
    # rt.vTqdm.update()
    
    # Reset sightline progress bar
    if multiprocess == 0:
        rt.slTqdm.reset()
  
    # Initialise the intensities and map
    position = []
    intensity_map_species = []
    intensity_map_dust = []
    opticaldepth_map_species = []
    opticaldepth_map_dust = []
  
    # Get indeces of voxels at the correct observing velocity
    # iClumpV = np.where(iCV)
    # iInterclumpV = np.where(iIV)
  
    if rt.i_vox.any() is False:
        # print('\n\n', [], '\n\n')
        return 0, 0, 0, []  # sightlines
  
    # The voxel positions can be any of the voxels
    # rt.tempSpeciesEmissivity = tempSpeciesEmissivity#[iV,i,:] / constants.pc/100
    # rt.tempSpeciesAbsorption = tempSpeciesAbsorption#[iV,i,:] / constants.pc/100
    # rt.tempDustEmissivity = tempDustEmissivity#[iV,i,:] / constants.pc/100
    # rt.tempDustAbsorption = tempDustAbsorption#[iV,i,:] / constants.pc/100
    # rt.tempPosition = rt.voxelPositions[0].data[iV,:]
    # rt.tempClumpVelocity = clumpVelocity[iClumpV[0],:]
    # rt.tempInterclumpEmission = interclumpEmission[:,iInterclumpV[0],iInterclumpV[1],:]
    # rt.tempInterclumpPosition = voxelPositions[iInterclumpV[0],:]
    # rt.tempInterclumpVelocity = interclumpVelocity[iInterclumpV[0],:]
  
    for j, lat in enumerate(latgrid):
  
        # Initialise a list containing the intensities at each lattitude
        positionintensity_species = []
        positionintensity_dust = []
        positionopticaldepth_species = []
        positionopticaldepth_dust = []
    
        for i, lon in enumerate(longrid):
          
            # print('lon,lat before scipy call', lon, lat, ':', time()-t0)
      
            # Calculate sightline length
            Rslh = op.root_scalar(sightlength, args=lon, x0=constants.rG_gal_earth, x1=constants.rgal).root
            thetaC = np.arctan(constants.hd/Rslh)
            # if verbose:
            #   print('\n',thetaC,'\n')
            if abs(lat) < thetaC:
                Rsl = Rslh/np.cos(lat)
            else:
                Rsl = constants.hd/np.sin(abs(lat))
      
            # print('lon,lat after scipy call', lon, lat, ':', time()-t0)
            
            # Try to isolate voxels in LoS and work with all transitions, else do them one at a time
            try:
                result, vox = setLOS(lon=lon, lat=lat, i_vel=i_vel, dim=dim, debug=debug)
            except OSError:
                rt.logger.critical('OSError!!')
                sys.exit()
      
            position.append([lon, lat])
            
            if vox:
                normfactor = Rsl/constants.voxel_size/vox  # to normalise to disk shape...
                # if sightlines[i,j]<vox: sightlines[i,j] = normfactor*vox
                sightlines[i, j] = vox
            
            # Integrate along LoS
            if vox == 1:
                positionintensity_species.append(rt.intensity_species[0] * normfactor)
                positionintensity_dust.append(rt.intensity_dust * normfactor)
                positionopticaldepth_species.append(rt.opticaldepth_species[0] * normfactor)
                positionopticaldepth_dust.append(rt.opticaldepth_dust * normfactor)
            elif vox > 1:
                calculatert(scale=normfactor)
                positionintensity_species.append(rt.intensity_species[0])
                positionintensity_dust.append(rt.intensity_dust)
                positionopticaldepth_species.append(rt.opticaldepth_species[0] * normfactor)
                positionopticaldepth_dust.append(rt.opticaldepth_dust * normfactor)
                # intensity_species = []
                # intensity_dust = []
            else:
                positionintensity_species.append(np.zeros(rt.temp_species_emissivity[0].shape[-1]))
                positionintensity_dust.append(np.zeros(rt.temp_dust_emissivity[0].shape[-1]))
                positionopticaldepth_species.append(np.zeros(rt.temp_species_emissivity[0].shape[-1]))
                positionopticaldepth_dust.append(np.zeros(rt.temp_dust_emissivity[0].shape[-1]))
            
            if multiprocess == 0:
                rt.slTqdm.update()
            
            # if len(np.shape(positionintensity_species[-1])) > 1:
            #     rt.logger.error('Error {}'.format(np.shape(positionintensity_species[-1])))
            #     input()
        
        # Save intensities for each latitude
        intensity_map_species.append(positionintensity_species)
        intensity_map_dust.append(positionintensity_dust)
        opticaldepth_map_species.append(positionopticaldepth_species)
        opticaldepth_map_dust.append(positionopticaldepth_dust)
    
    # Save map for velocity channel
    # VintensityMapSpecies.append(intensityMapSpecies)
    # VintensityMapDust.append(intensityMapDust)
    # Vpositions.append(position)
  
    # if verbose:
    #   print('Evaluating {} km/s HDU'.format(velocity))
  
    # if plotV:
    #   fig = plt.figure()
    #   ax = fig.add_subplot(111, projection='mollweide')
    #   cb = ax.scatter(np.array(position)[:,0], np.array(position)[:,1],
    #                   c=np.array(intensityMapSpecies)[:,:,0,23].flatten(), s=64, marker='s')
    #   plt.ion()
    #   fig.colorbar(cb)
    #   ax.grid(True)
    #   ax.set_xticklabels([])
    #   ax.set_yticklabels([])
    #   plt.show(block=True)
    
    # print('\n\n', sightlines.shape, '\n\n')
    return (position, (intensity_map_species, opticaldepth_map_species), 
            (intensity_map_dust, opticaldepth_map_dust), sightlines)


def calculate_sightline(los, slRange=((-np.pi, np.pi), (-np.pi/2, np.pi/2)), nsl=(50, 25),
                        dim='spherical', debug=False, multiprocess=0):

    # # Convert the tuple to the desired index and velocity
    # i_vel = ivelocity[0]
    i_pos = los[0]
    lon = los[1][0]
    lat = los[1][1]
    position = (i_pos, lon, lat)

    # if multiprocess:
    #   t0 = time()
    #   print('\nCalculating velocity channel at {:.2f} km\s'.format(ivelocity[1]))

    # # Setup the sightlines that are calculated
    # longrid = np.linspace(slRange[0][0], slRange[0][1], num=nsl[0])
    # latgrid = np.linspace(slRange[1][0], slRange[1][1], num=nsl[1])
    # sightlines = np.zeros((longrid.size, latgrid.size))

    # print('lon/lat arrays created:', time()-t0)

    # # Find the voxels that exist at the observing velocity
    # nDust = rt.tempDustEmissivity[0].shape[2]
    # if nDust > 1:
    #     base = rt.interpDust(species.moleculeWavelengths)[:, :, 0]
    # else:
    #     base = rt.tempDustEmissivity[0].data[0, 0]

    rt.i_vox = rt.temp_dust_emissivity[0].data.any(1)

    # print('Voxels selected (shape={}):'.format(i_vox[i_vox].shape), time()-t0)

    # Update velocity progress bar
    # rt.vTqdm.update()

    # # Reset sightline progress bar
    # if multiprocess == 0:
    #     rt.slTqdm.reset()

    # # Initialise the intensities and map
    # position = []
    # intensityMapSpecies = []
    # intensityMapDust = []

    # Get indeces of voxels at the correct observing velocity
    # iClumpV = np.where(iCV)
    # iInterclumpV = np.where(iIV)

    if rt.i_vox.any() is False:
        # print('\n\n', [], '\n\n')
        return 0, (0, 0), (0, 0), 0  # sightlines

    # The voxel positions can be any of the voxels
    # rt.tempSpeciesEmissivity = tempSpeciesEmissivity#[iV,i,:] / constants.pc/100
    # rt.tempSpeciesAbsorption = tempSpeciesAbsorption#[iV,i,:] / constants.pc/100
    # rt.tempDustEmissivity = tempDustEmissivity#[iV,i,:] / constants.pc/100
    # rt.tempDustAbsorption = tempDustAbsorption#[iV,i,:] / constants.pc/100
    # rt.tempPosition = rt.voxelPositions[0].data[iV,:]
    # rt.tempClumpVelocity = clumpVelocity[iClumpV[0],:]
    # rt.tempInterclumpEmission = interclumpEmission[:,iInterclumpV[0],iInterclumpV[1],:]
    # rt.tempInterclumpPosition = voxelPositions[iInterclumpV[0],:]
    # rt.tempInterclumpVelocity = interclumpVelocity[iInterclumpV[0],:]

    # # Initialise a list containing the intensities at each lattitude
    # positionintensity_species = []
    # positionintensity_dust = []
    # positionopticaldepth_species = []
    # positionopticaldepth_dust = []

    # print('lon,lat before scipy call', lon, lat, ':', time()-t0)

    # Calculate sightline length
    Rslh = op.root_scalar(sightlength, args=lon, x0=constants.rgal_earth, x1=constants.rgal).root
    thetaC = np.arctan(constants.hd/Rslh)
    # if verbose:
    #   print('\n',thetaC,'\n')
    if abs(lat) < thetaC:
        Rsl = Rslh/np.cos(lat)
    else:
        Rsl = constants.hd/np.sin(abs(lat))

    # print('lon,lat after scipy call', lon, lat, ':', time()-t0)

    # Try to isolate voxels in LoS and work with all transitions, else do them one at a time
    try:
        result, vox = set_los(lon=lon, lat=lat, dim=dim, debug=debug)
    except OSError:
        rt.logger.critical('OSError!!')
        sys.exit()

    # position.append([lon, lat])

    if vox > 0:
        normfactor = Rsl/constants.voxel_size/vox  # to normalise to disk shape...
        # if sightlines[i,j]<vox: sightlines[i,j] = normfactor*vox

    # Integrate along LoS
    if vox == 1:
        # positionintensity_species.append(rt.intensity_species * normfactor)
        # positionintensity_dust.append(rt.intensity_dust * normfactor)
        # positionopticaldepth_species.append(rt.opticaldepth_species * normfactor)
        # positionopticaldepth_dust.append(rt.opticaldepth_dust * normfactor)
        result = ((rt.intensity_species * normfactor, rt.opticaldepth_species * normfactor),
                  (rt.intensity_dust * normfactor, rt.opticaldepth_dust * normfactor))
    elif vox > 1:
        calculatert(scale=normfactor)
        # positionintensity_species.append(rt.intensity_species)
        # positionintensity_dust.append(rt.intensity_dust)
        # positionopticaldepth_species.append(rt.k_species.sum(0) * normfactor*constants.voxel_size*100*constants.pc)
        # positionopticaldepth_dust.append(rt.kopticaldepth_dust.sum(0) * normfactor*constants.voxel_size*100*constants.pc)
        result = ((c.copy(rt.intensity_species), rt.opticaldepth_species * normfactor),
                  (c.copy(rt.intensity_dust), rt.opticaldepth_dust * normfactor))
        # intensity_species = []
        # intensity_dust = []
    else:
        # positionintensity_species.append(np.zeros(rt.tempSpeciesEmissivity[0].shape[1:]))
        # positionintensity_dust.append(np.zeros(rt.tempDustEmissivity[0].shape[1]))
        # positionopticaldepth_species.append(np.zeros(rt.tempSpeciesAbsorption[0].shape[1:]))
        # positionopticaldepth_dust.append(np.zeros(rt.tempDustAbsorption[0].shape[1]))
        result = ((np.zeros(rt.temp_species_emissivity[0].shape[1:]), np.zeros(rt.temp_species_absorption[0].shape[1:])),
                  (np.zeros(rt.temp_dust_emissivity[0].shape[1]), np.zeros(rt.temp_dust_absorption[0].shape[1])))

    #if len(np.shape(positionintensity_species[-1])) > 1:
    #    rt.logger.error('Error {}'.format(np.shape(positionintensity_species[-1])))
    #    input()

    # # Save intensities for each latitude
    # intensity_species.append(positionintensity_species)
    # intensity_dust.append(positionintensity_dust)

    # Save map for velocity channel
    # VintensityMapSpecies.append(intensityMapSpecies)
    # VintensityMapDust.append(intensityMapDust)
    # Vpositions.append(position)

    # if verbose:
    #   print('Evaluating {} km/s HDU'.format(velocity))

    # if plotV:
    #   fig = plt.figure()
    #   ax = fig.add_subplot(111, projection='mollweide')
    #   cb = ax.scatter(np.array(position)[:,0], np.array(position)[:,1],
    #                   c=np.array(intensityMapSpecies)[:,:,0,23].flatten(), s=64, marker='s')
    #   plt.ion()
    #   fig.colorbar(cb)
    #   ax.grid(True)
    #   ax.set_xticklabels([])
    #   ax.set_yticklabels([])
    #   plt.show(block=True)

    # print('\n\n', sightlines.shape, '\n\n')
    return (position, *result, vox)


def set_los(x=0, y=0, z=0, lon=0, lat=0, i_vox=[], i_vel=None, i_spe=None, i_dust=None,
           dim='xy', reverse=False, debug=False, verbose=False):
    '''
    The emission dimensions should be velocity x species x 2 x voxels.
    Axis 1 should be split for intensity and optical depth.
    The positions dimensions should be 3 x voxels.

    A function to find the voxels in a particular line-of-sight. This takes a position (specified in Cartesian or
    spherical coordinates) as input, and uses the previously-opened sub-module variables for the voxel positions.
    It the determines the voxels in the line-of-sight, orders them from farthest to closest, and saves the
    sub-module variables necessary for the RT calculation.

    :param x:
    :param y:
    :param z:
    :param lon:
    :param lat:
    :param i_vox:
    :param i_vel:
    :param i_spe:
    :param i_dust:
    :param dim:
    :param reverse:
    :param debug:
    :param verbose:
    :return:
        :param calculation_code: 0 for no voxels in the line-of sight, 1 for one voxel in the line-of-sight, and
            2 for more then one voxel in the line-of-sight. This is important for the radiative transfer calculation.
        :param n_vel: the number of the voxels in the line-of-sight. The original use (velocities) mostly defunct due
            to how the velocity dimension is handled, but rather it is used to get the number of voxels in the
            line-of-sight.
    '''
  
    scale = c.copy(constants.voxel_size)
  
    # #print(rt.tempClumpEmission.shape)
    #
    # # This block matches the voxel positions to add the ensembles belonging to the same voxel
    # nrows, ncols = rt.tempClumpPosition.shape
    # dtype={'names':['f{}'.format(i) for i in range(ncols)], \
    #        'formats':ncols * [rt.tempClumpPosition.dtype]}
    # common,iCommonClump,iCommonInterclump = np.intersect1d(rt.tempClumpPosition.view(dtype),
    #                                                        rt.tempInterclumpPosition.view(dtype), return_indices=True)
    #
    # # print(iCommonClump.max())
    #
    # # Intensity and optical depth have shape (voxel, wavelength)
    # gridIntensity = rt.tempInterclumpEmission[0,:,:]
    # gridIntensity[iCommonInterclump,:] = gridIntensity[iCommonInterclump,:] \
    #                                      + rt.tempClumpEmission[0,iCommonClump,:]# emission[0,:,:]#
    #
    # gridOpticalDepth = rt.tempInterclumpEmission[1,:,:]# emission[1,:,:]#
    # gridOpticalDepth[iCommonInterclump,:] += rt.tempClumpEmission[1,iCommonClump,:]
    
    # gridIntensity -> rt.tempSpeciesEmissivity/tempDustEmissivity,
    # gridOpticalDepth -> rt.tempSpeciesAbsorption/tempDustAbsorption
  
    # Specify the voxel positions relative to Earth
    xGrid, yGrid, zGrid = (rt.voxel_positions[0].data[:, 0],
                           rt.voxel_positions[0].data[:, 1],
                           rt.voxel_positions[0].data[:, 2],)
  
    if dim == 'spherical':
        # Set sightline position
        x1LoS = lon
        x2LoS = lat
    
        # Convert voxel positions to spherical
        radGrid = np.sqrt((xGrid-constants.rgal_earth)**2 + yGrid**2 + zGrid**2)
        lonGrid = np.arctan2(yGrid, -(xGrid-constants.rgal_earth))
        if lon < 0:
            lonGrid[lonGrid > 0] = lonGrid[lonGrid > 0] - 2*np.pi
        if lon > 0:
            lonGrid[lonGrid < 0] = lonGrid[lonGrid < 0] + 2*np.pi
        rPolar = np.sqrt((xGrid-constants.rgal_earth)**2+yGrid**2)
        latGrid = np.arctan2(zGrid, rPolar)
        if lat < 0:
            latGrid[latGrid > 0] = latGrid[latGrid > 0] - np.pi
        if lat > 0:
            latGrid[latGrid < 0] = latGrid[latGrid < 0] + np.pi
    
        # Choose the voxels in the sightline
        # adjustments for voxel orientation
        scaling = np.sqrt(1+rt.voxel_filling_factors)  # np.sqrt(2)
        width = scaling*constants.voxel_size*np.max([np.sin(np.abs(lonGrid-np.pi/4)),
                                                     np.sin(np.abs(lonGrid+np.pi/4))], axis=0)
        height = scaling*constants.voxel_size*np.max([np.sin(np.abs(latGrid-np.pi/4)),
                                                      np.sin(np.abs(latGrid+np.pi/4))], axis=0)
        # angular size of voxels (theoretically correct)
        if rt.pencil_beam:
            d_lon = 0.5*np.arctan(width/radGrid)
            d_lon[(lat > 1) | (lat < -1)] = 3.142  # check if voxel in los but not in selection criteria 
            d_lat = 0.5*np.arctan(height/radGrid)
            d_lat[(lat > 1) | (lat < -1)] = 1.571
        # beam width (observationally correct)
        else:
            d_lon = rt.d_lon
            d_lat = rt.d_lat
        i_los = np.where((abs(lonGrid-x1LoS) <= d_lon) & (abs(latGrid-x2LoS) <= d_lat) & rt.i_vox)[0]
        # i_los = np.where((abs(lonGrid-x1LoS) <= d_lon) & (abs(latGrid-x2LoS) <= d_lat))[0]
    
    elif 'disk' in dim:
        x1LoS = y
        x2LoS = z
        i_los = np.where((zGrid == z) & (yGrid == y))[0]
    
    elif ('x' in dim) and ('y' in dim):
        x1LoS = x
        x2LoS = y
        i_los = np.where((xGrid == x) & (yGrid == y))[0]
    
    elif ('x' in dim) and ('z' in dim):
        x1LoS = x
        x2LoS = z
        i_los = np.where((xGrid == x) & (zGrid == z))[0]
    
    elif ('z' in dim) and ('y' in dim):
        x1LoS = z
        x2LoS = y
        i_los = np.where((zGrid == z) & (yGrid == y))[0]
  
    else:
        rt.logger.error('\nPlease enter valid dimensions.\n')
        return False
    
    rt.logger.info(i_los)

    n_vel = rt.temp_species_emissivity[0].shape[1]
    n_spe = rt.temp_species_emissivity[0].shape[2]
    n_dust = rt.temp_dust_emissivity[0].shape[1]
    
    if i_los.size == 1:

        if not i_vel is None:
            ix_species = np.ix_(i_los, [i_vel], np.arange(n_spe))
            ix_dust = np.ix_(i_los, np.arange(n_dust))
        else:
            ix_species = np.ix_(i_los, np.arange(n_vel), np.arange(n_spe))
            ix_dust = np.ix_(i_los, np.arange(n_dust))
        
        eps_species = rt.temp_species_emissivity[0].data[ix_species][0, :, :] \
                      / rt.voxel_filling_factors[i_los]
        kap_species = rt.temp_species_emissivity[0].data[ix_species][0, :, :] \
                      / rt.voxel_filling_factors[i_los]
        eps_dust = rt.temp_dust_emissivity[0].data[ix_dust][0, :] / rt.voxel_filling_factors[i_los]
        kap_dust = rt.temp_dust_emissivity[0].data[ix_dust][0, :] / rt.voxel_filling_factors[i_los]

        #rt.intensity_species = rt.temp_species_emissivity[0].data[ix_species][0, :, :]
        rt.intensity_species = c.copy(eps_species)
        i_nan = kap_species < 1e-16
        rt.intensity_species = rt.intensity_species[~i_nan] / kap_species[~i_nan] \
                               * (1-np.exp(-kap_species[~i_nan]*scale))
        #rt.intensity_dust = rt.temp_dust_emissivity[0].data[ix_dust][0, :]
        rt.intensity_dust = c.copy(eps_dust)
        i_nan = kap_dust < 1e-16
        rt.intensity_dust= rt.intensity_dust[~i_nan] / kap_dust[~i_nan] \
                           * (1-np.exp(-kap_dust[~i_nan]*scale))
        rt.opticaldepth_species = scale * kap_species
        rt.opticaldepth_dust = scale * kap_dust
        
        # print('\n\n', (rt.tempSpeciesEmissivity[0].data[iLOS,i_vel,:]>0).any(), '\n\n')
    
        vel = rt.voxel_velocities[0].data[i_los]
        
        return 1, vel.size  # ,[intensity_species,intensity_dust]
    
    elif i_los.size > 1:
      
        if 'spherical' == dim:
            x3LoS = radGrid[i_los]
        elif not 'x' in dim:
            x3LoS = xGrid[i_los]
        elif not 'y' in dim:
            x3LoS = yGrid[i_los]
        elif not 'z' in dim:
            x3LoS = zGrid[i_los]
          
        if reverse:
            # order the voxels in the line-of-sight according to increasing radial distance,
            # compute RT from observer to edge of model
            i = np.argsort(x3LoS)
            i_los_ordered = i_los[i]
        else:
            # order the voxels in the line-of-sight according to decreasing radial distance,
            # compute RT from edge of model to observer
            i = np.argsort(x3LoS)[::-1]
            i_los_ordered = i_los[i]
        
        # print('\n\n', iLoS, i, iLoS_ordered, i_vel,
        #       (rt.tempSpeciesEmissivity[0].data[iLoS_ordered,i_vel,:]>0).any(), '\n\n')
        # print(rt.i_vox)
        # print(rt.x3LoS)
        # input()

        if not i_vel is None:
            ix_species = np.ix_(i_los_ordered, [i_vel], np.arange(n_spe))
            ix_dust = np.ix_(i_los_ordered, np.arange(n_dust))
        else:
            ix_species = np.ix_(i_los_ordered, np.arange(n_vel), np.arange(n_spe))
            ix_dust = np.ix_(i_los_ordered, np.arange(n_dust))
        
        # Species transitions
        #voxel_filling = rt.voxel_filling_factors[i_los_ordered][:-1]
        rt.e_species = c.copy(rt.temp_species_emissivity[0].data[ix_species][:-1] 
                / rt.voxel_filling_factors[i_los_ordered][:-1, np.newaxis, np.newaxis])
        rt.de_species = (rt.temp_species_emissivity[0].data[ix_species][1:]
                         - rt.temp_species_emissivity[0].data[ix_species][:-1]) / (scale)
        rt.k_species = c.copy(rt.temp_species_absorption[0].data[ix_species][:-1]
                / rt.voxel_filling_factors[i_los_ordered][:-1, np.newaxis, np.newaxis])
        rt.dk_species = (rt.temp_species_absorption[0].data[ix_species][1:]
                         - rt.temp_species_absorption[0].data[ix_species][:-1]) / (scale)
        eps0_species = rt.e_species[0]
        kap0_species = rt.k_species[0]
        rt.i0_species = c.copy(eps0_species)
        i_nan = kap0_species < 1e-16
        rt.i0_species[~i_nan] = eps0_species[~i_nan] / kap0_species[~i_nan] \
                                * (1-np.exp(-kap0_species[~i_nan])*scale)
        rt.opticaldepth_species = scale * rt.temp_species_absorption[0].data[ix_species].sum(0)
            
        # Dust continuum
        rt.e_dust = c.copy(rt.temp_dust_emissivity[0].data[ix_dust][:-1]
                / rt.voxel_filling_factors[i_los_ordered][:-1, np.newaxis])
        rt.de_dust = (rt.temp_dust_emissivity[0].data[ix_dust][1:]
                      - rt.temp_dust_emissivity[0].data[ix_dust][:-1]) / (scale)
        rt.k_dust = c.copy(rt.temp_dust_absorption[0].data[ix_dust][:-1]
                / rt.voxel_filling_factors[i_los_ordered][:-1, np.newaxis])
        rt.dk_dust = (rt.temp_dust_absorption[0].data[ix_dust][1:]
                      - rt.temp_dust_absorption[0].data[ix_dust][:-1]) / (scale)
        eps0_dust = rt.e_dust[0]
        kap0_dust = rt.k_dust[0]
        rt.i0_dust = c.copy(eps0_dust)
        i_nan = kap0_dust < 1e-16
        rt.i0_dust[~i_nan] = eps0_dust[~i_nan] / kap0_dust[~i_nan] \
                             * (1-np.exp(-kap0_dust[~i_nan]*scale))
        rt.opticaldepth_dust = scale * rt.temp_dust_absorption[0].data[ix_dust].sum(0)
    
    else:
        return 0, 0  # ,[]
  
    vel = rt.voxel_velocities[0].data[i_los]
    
    rt.logger.info(f'voxels: {i}')
    
    return 2, vel.size
    # ,(e_species,de_species,k_species,dk_species,
    #   e_dust,de_dust,k_dust,dk_dust)


def calculatert(scale=1, background_species_intensity=None, background_dust_intensity=None, 
                species=True, dust=True, verbose=False, test=True):

    # Modify data and initialise observed intensity array
    # if not dust:
    #   rt.e_species = rt.e_species.reshape((-1,1))
    #   rt.de_species = rt.de_species.reshape((-1,1))
    #   rt.k_species = rt.k_species.reshape((-1,1))
    #   rt.dk_species = rt.dk_species.reshape((-1,1))
    #   intensity_species = np.array([backgroundI])
    #   nSteps = rt.de_species.shape[0]
    # elif not species:
    #   rt.e_dust = rt.e_dust.reshape((-1,1))
    #   rt.de_dust = rt.de_dust.reshape((-1,1))
    #   rt.k_dust = rt.k_dust.reshape((-1,1))
    #   rt.dk_dust = rt.dk_dust.reshape((-1,1))
    #   intensity_dust = np.array([backgroundI])
    #   nSteps = rt.de_dust.shape[0]
    # else:
    
    if isinstance(background_species_intensity, (float, int)):
        rt.intensity_species = np.full(rt.k_species.shape[1:], background_species_intensity)
    else:
        #rt.intensity_species = c.copy(rt.i0_species)
        rt.intensity_species = np.full(rt.k_species.shape[1:], 0.)
    if isinstance(background_dust_intensity, (float, int)):
        rt.intensity_dust = np.full(rt.k_dust.shape[1], background_dust_intensity)
    else:
        #rt.intensity_dust = c.copy(rt.i0_dust)
        rt.intensity_dust = np.full(rt.k_dust.shape[1], 0.)
    
    n_steps = rt.e_species.shape[0]
    
    # Adjust according to the average voxel depth
    # if species:
    #     rt.e_species /= scale
    #     rt.de_species /= scale**2
    #     rt.k_species /= scale
    #     rt.dk_species /= scale**2
    # if dust:
    #     rt.e_dust /= scale
    #     rt.de_dust /= scale**2
    #     rt.k_dust /= scale
    #     rt.dk_dust /= scale**2
  
    scale = scale*constants.voxel_size  # scale voxel size
    # print(scale)
    
    #np.set_printoptions(threshold=100000)
    warnings.filterwarnings('error')
    
    # # Boolean indeces to separate how the intensity is calculated
    # k0 = (rt.kappaStep==0)&(abs(rt.kappa[:-1]*constants.resolution)<10**-10)
    # kg = rt.kappa[:-1]>10**3*abs(rt.kappaStep)*constants.resolution
    # kE = ~(k0|kg)
    # kEg = ~(k0|kg)&(rt.kappaStep>0)
    # kEl = ~(k0|kg)&(rt.kappaStep<0)
    
    # Calculate the variables needed to utilise the E tilde tables
    with warnings.catch_warnings():
        warnings.simplefilter('ignore')
        if species:
            a_species = (rt.k_species / np.sqrt(2*np.abs(rt.dk_species)))
            b_species = ((rt.k_species + rt.dk_species*scale) / np.sqrt(2*np.abs(rt.dk_species)))
        if dust:
            a_dust = (rt.k_dust / np.sqrt(2*np.abs(rt.dk_dust)))
            b_dust = ((rt.k_dust + rt.dk_dust*scale) / np.sqrt(2*np.abs(rt.dk_dust)))
    
    rt.logger.info(rt.k_species.shape[0])
        
    warnings.filterwarnings('ignore')
    # if test:
    #     print(n_steps)
    
    for i in range(n_steps):
      
        # << Compute radiative transfer for the species transitions >>
        if species:
    
            # Determine which form of integration the species transitions require
            k0_species = ((rt.dk_species[i, :, :] == 0) &
                          (np.abs(rt.k_species[i, :, :]*scale) < 10**-10))
            kg_species = (~k0_species & 
                          (np.abs(rt.k_species[i, :, :])
                          > (10**3*np.abs(rt.dk_species[i, :, :])*scale)))
            keg_species = (~k0_species & ~kg_species & (rt.dk_species[i, :, :] > 0))
            kel_species = (~k0_species & ~kg_species & (rt.dk_species[i, :, :] < 0))

            ix_k0 = (i, *np.where(k0_species))
            ix_kg = (i, *np.where(kg_species))
            ix_keg = (i, *np.where(keg_species))
            ix_kel = (i, *np.where(kel_species))
            
            # Integrate using the required method
            if k0_species.any():
                rt.intensity_species[ix_k0[1:]] = (rt.e_species[ix_k0] * scale
                                                        + 0.5*rt.de_species[ix_k0] * scale**2.
                                                        + rt.intensity_species[ix_k0[1:]])
        
            if kg_species.any():
                rt.intensity_species[ix_kg[1:]] = ((rt.e_species[ix_kg]*rt.k_species[ix_kg]
                                                    + rt.de_species[ix_kg]
                                                    * (rt.k_species[ix_kg]*scale-1.))
                                                    / (rt.k_species[ix_kg]**2.)) \
                                                  - np.exp(-rt.k_species[ix_kg]*scale) \
                                                  * ((rt.e_species[ix_kg]*rt.k_species[ix_kg]
                                                     - rt.de_species[ix_kg])
                                                     / (rt.k_species[ix_kg]**2.)) \
                                                  + (rt.intensity_species[ix_kg[1:]]
                                                     * np.exp(-rt.k_species[ix_kg]*scale))
      
            if keg_species.any():
                rt.logger.info('\na, b:\n{}\n{}'.format(a_dust[i], b_dust[i]))
                a_error = erfi(a_species[ix_keg])
                b_error = erfi(b_species[ix_keg])
                rt.intensity_species[ix_keg[1:]] = (rt.de_species[ix_keg]/rt.dk_species[ix_keg]
                                                    * (1.-np.exp(-rt.k_species[ix_keg]*scale
                                                                 - 0.5*rt.dk_species[ix_keg]*scale**2.))
                                                    - (rt.e_species[ix_keg]*rt.dk_species[ix_keg]
                                                       - rt.k_species[ix_keg]*rt.de_species[ix_keg])
                                                    / rt.dk_species[ix_keg]
                                                    * np.sqrt(np.pi/2./np.abs(rt.dk_species[ix_keg]))
                                                    * np.exp(-b_species[ix_keg]**2.) * (a_error - b_error)
                                                    + rt.intensity_species[ix_keg[1:]]
                                                    * np.exp(-rt.k_species[ix_keg]*scale
                                                             - 0.5*rt.dk_species[ix_keg]*scale**2.))
      
            if kel_species.any():
                rt.logger.info('\na, b:\n{}\n{}'.format(a_dust[i], b_dust[i]))
                a_error = erfc(a_species[ix_kel])
                b_error = erfc(b_species[ix_kel])
                rt.intensity_species[ix_kel[1:]] = (rt.de_species[ix_kel]/rt.dk_species[ix_kel]
                                                    * (1.-np.exp(-rt.k_species[ix_kel]*scale
                                                                 - 0.5*rt.dk_species[ix_kel]*scale**2.))
                                                    - (rt.e_species[ix_kel]*rt.dk_species[ix_kel]
                                                       - rt.k_species[ix_kel]*rt.de_species[ix_kel])
                                                    / rt.dk_species[ix_kel]
                                                    * np.sqrt(np.pi/2./np.abs(rt.dk_species[ix_kel]))
                                                    * np.exp(b_species[ix_kel]**2.) * (a_error - b_error)
                                                    + rt.intensity_species[ix_kel[1:]]
                                                    * np.exp(-rt.k_species[ix_kel]*scale
                                                             - 0.5*rt.dk_species[ix_kel]*scale**2.))
        
        # << Compute radiative transfer for the dust continuum >>
        if dust:
    
            # Determine which form of integration the dust continuum requires
            k0_dust = ((rt.dk_dust[i, :] == 0) &
                       (abs(rt.k_dust[i, :]*scale) < 10**-10))
            kg_dust = (~k0_dust & (np.abs(rt.k_dust[i, :])
                                   > (10**3*abs(rt.dk_dust[i, :])*scale)))
            keg_dust = (~k0_dust & ~kg_dust & (rt.dk_dust[i, :] > 0))
            kel_dust = (~k0_dust & ~kg_dust & (rt.dk_dust[i, :] < 0))

            ix_k0 = (i, *np.where(k0_dust))
            ix_kg = (i, *np.where(kg_dust))
            ix_keg = (i, *np.where(keg_dust))
            ix_kel = (i, *np.where(kel_dust))
            
            # Integrate using the required method
            if k0_dust.any():
                rt.intensity_dust[ix_k0[1:]] = (rt.e_dust[ix_k0]*scale
                                                + 0.5*rt.de_dust[ix_k0]*scale**2.
                                                + rt.intensity_dust[ix_k0[1:]])
        
            if kg_dust.any():
                rt.intensity_dust[ix_kg[1:]] = ((rt.e_dust[ix_kg]*rt.k_dust[ix_kg]
                                                 + rt.de_dust[ix_kg]*(rt.k_dust[ix_kg]*scale-1.))
                                                / (rt.k_dust[ix_kg]**2.)) \
                                               - np.exp(-rt.k_dust[ix_kg]*scale) \
                                               * ((rt.e_dust[ix_kg]*rt.k_dust[ix_kg]
                                                   - rt.de_dust[ix_kg])
                                                  / (rt.k_dust[ix_kg]**2.)) \
                                                          + (rt.intensity_dust[ix_kg[1:]]
                                                  * np.exp(-rt.k_dust[ix_kg]*scale))
      
            if keg_dust.any():
                rt.logger.info('\na, b:\n{}\n{}'.format(a_dust[i], b_dust[i]))
                a_error = erfi(a_dust[ix_keg])
                b_error = erfi(b_dust[ix_keg])
                rt.intensity_dust[ix_keg[1:]] = (rt.de_dust[ix_keg]/rt.dk_dust[ix_keg]
                                                 * (1.-np.exp(-rt.k_dust[ix_keg]*scale
                                                              - 0.5*rt.dk_dust[ix_keg]*scale**2.))
                                                 - (rt.e_dust[ix_keg]*rt.dk_dust[ix_keg]
                                                    - rt.k_dust[ix_keg]*rt.de_dust[ix_keg])
                                                 / rt.dk_dust[ix_keg]
                                                 * np.sqrt(np.pi/2./np.abs(rt.dk_dust[ix_keg]))
                                                 * np.exp(-b_dust[ix_keg]**2.) * (a_error - b_error)
                                                 + rt.intensity_dust[ix_keg[1:]]
                                                 * np.exp(-rt.k_dust[ix_keg]*scale
                                                          - 0.5*rt.dk_dust[ix_keg]*scale**2.))
      
            if kel_dust.any():
                rt.logger.info('\na, b:\n{}\n{}'.format(a_dust[i], b_dust[i]))
                a_error = erfc(a_dust[ix_kel])
                b_error = erfc(b_dust[ix_kel])
                rt.intensity_dust[ix_kel[1:]] = (rt.de_dust[ix_kel]/rt.dk_dust[ix_kel]
                                                 * (1.-np.exp(-rt.k_dust[ix_kel]*scale
                                                              - 0.5*rt.dk_dust[ix_kel]*scale**2.))
                                                 - (rt.e_dust[ix_kel] * rt.dk_dust[ix_kel]
                                                    - rt.k_dust[ix_kel]*rt.de_dust[ix_kel])
                                                 / rt.dk_dust[ix_kel]
                                                 * np.sqrt(np.pi/2./np.abs(rt.dk_dust[ix_kel]))
                                                 * np.exp(b_dust[ix_kel]**2.) * (a_error - b_error)
                                                 + rt.intensity_dust[ix_kel[1:]]
                                                 * np.exp(-rt.k_dust[ix_kel]*scale
                                                          - 0.5*rt.dk_dust[ix_kel]*scale**2.))
    
    rt.logger.info('Species intensity shape: {}'.format(np.shape(rt.intensity_species)))
    rt.logger.info('Dust intensity shape: {}'.format(np.shape(rt.intensity_dust)))
        
    # if (rt.intensity_species > 10**10).any() or (rt.intensity_species < 0).any():
    #
    #     print('\n\nSome of the species have either suspiciously large or negative intensities...')
    #
    #     dir = r'c:\users\cyani\KOSMA-tau^3\tests\full model'
    #
    #     i = np.where((rt.intensity_species > 10**10) | (rt.intensity_species < 0))[0]
    #     print('The indices are:', i)
    #     print('intensity:', rt.intensity_species[i])
    #     print('\n')
    #     np.save(dir+r'\ds_species.npy', scale)
    #     np.save(dir+r'\I_species.npy', rt.intensity_species[i])
    #     np.save(dir+r'\e_species.npy', rt.e_species[:, i])
    #     np.save(dir+r'\de_species.npy', rt.de_species[:, i])
    #     np.save(dir+r'\k_species.npy', rt.k_species[:, i])
    #     np.save(dir+r'\dk_species.npy', rt.dk_species[:, i])
    #     np.save(dir+r'\a_species.npy', a_species[:, i])
    #     np.save(dir+r'\b_species.npy', b_species[:, i])
    #
    # if (rt.intensity_dust > 10**10).any() or (rt.intensity_dust < 0).any():
    #
    #     print('\n\nSome of the dust has either suspiciously large or negative intensities...')
    #
    #     dir = r'c:\users\cyani\KOSMA-tau^3\tests\full model'
    #
    #     i = np.where((rt.intensity_dust > 10**10) | (rt.intensity_dust < 0))[0]
    #     print('The indices are:', i)
    #     print('intensity:', rt.intensity_dust[i])
    #     print('\n')
    #     np.save(dir+r'\ds_dust.npy', scale)
    #     np.save(dir+r'\I_dust.npy', rt.intensity_dust[i])
    #     np.save(dir+r'\e_dust.npy', rt.e_dust[:, i])
    #     np.save(dir+r'\de_dust.npy', rt.de_dust[:, i])
    #     np.save(dir+r'\k_dust.npy', rt.k_dust[:, i])
    #     np.save(dir+r'\dk_dust.npy', rt.dk_dust[:, i])
    #     np.save(dir+r'\a_dust.npy', a_dust[:, i])
    #     np.save(dir+r'\b_dust.npy', b_dust[:, i])
  
    return  # (intensity_species, intensity_dust)


def e_real(x, verbose=False):
  
    rt.logger.info('E real input: {}'.format(x))
  
    e_r = np.zeros_like(x)
  
    il = x < 0.01
    ig = ~il & (x > 8.0)
    ib = ~(ig | il)
  
    if il.any():
  
        rt.logger.info('x less than grid')
    
        e_r[il] = (2*x[il]/np.sqrt(np.pi))
    
    if ig.any():
  
        rt.logger.info('x greater than grid')
    
        e_r[ig] = (1/(np.sqrt(np.pi) * x[ig]))
    
    if ib.any():
  
        rt.logger.info('x interpolated')
    
        e_r[ib] = rt.eTildeReal(x[ib])
  
    return e_r


def e_imag(x, verbose=False):
  
    rt.logger.info('E imaginary input:'.format(x))
  
    e_i = np.zeros_like(x)
  
    im = x > 0
    il = ~im & (np.abs(x) < 0.01)
    ig = ~im & (np.abs(x) > 8.0)
    ib = ~(ig | il | im)
  
    # Force x to be a positive real value
    # x = np.abs(x)
    
    if im.any():
  
        rt.logger.info('Maser case')
    
        # MASER case: treat with linear approximation
        e_i[im] = 1 + 2*np.abs(x[im])/np.sqrt(np.pi)
  
    if il.any():
  
        rt.logger.info('x less than grid')
        
        e_i[il] = (1 - 2*np.abs(x[il])/np.sqrt(np.pi))
  
    if ig.any():
      
        rt.logger.info('x greater than grid')
        
        e_i[ig] = (1/(np.sqrt(np.pi) * np.abs(x[ig])))
  
    if ib.any():
      
        rt.logger.info('x interpolated')
        
        e_i[ib] = rt.eTildeImaginary(np.abs(x[ib]))
  
    return e_i


if __name__ == '__main__':

    rt.logger.info('spawned process')
  
    directory = r'C:\Users\cyani\projects\pdr\KT3_history\MilkyWay\r250_cm1-1_d1_uv10'
  
    constants.velocityRange = np.linspace(-300, 300, 500)
  
    rt.voxel_positions = fits.open(directory+'/voxel_position.fits', mode='denywrite')
    rt.voxel_velocities = fits.open(directory+'/voxel_velocity.fits', mode='denywrite')
    rt.temp_species_emissivity = fits.open(directory+'/species_emissivity.fits', mode='denywrite')
    rt.temp_species_absorption = fits.open(directory+'/species_absorption.fits', mode='denywrite')
    rt.temp_dust_emissivity = fits.open(directory+'/dust_emissivity.fits', mode='denywrite')
    rt.temp_dust_absorption = fits.open(directory+'/dust_absorption.fits', mode='denywrite')
  
    multiprocessCalculation(slRange=[(-np.pi, np.pi), (-np.pi/2, np.pi/2)], nsl=[50, 25], dim='spherical',
                            multiprocessing=2)
    
    # multiprocessing.freeze_support()

# jit_module(nopython=False)
