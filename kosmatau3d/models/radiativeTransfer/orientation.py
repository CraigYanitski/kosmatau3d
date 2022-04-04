import sys
import warnings
import multiprocessing
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

    rt.voxelPositions = fits.open(directory + 'voxel_position.fits', mode='denywrite')
    return


def open_voxel_velocities(directory):
    '''
    Open the file containing the observed velocity of all voxels.
    This will save the contents in a sub-module variable.

    :param directory: path to the model files
    :return:

    shape -> (n_voxels, 1)
    '''

    rt.voxelVelocities = fits.open(directory+'voxel_velocity.fits', mode='denywrite')
    return


def open_voxel_emission(directory):
    '''
    Open the file containing the emissivity and absorption of all voxels.
    This will save the contents in a sub-module variable.

    :param directory: path to the model files
    :return:

    Emissivities in K/pc; absorptions in 1/pc

    species file shape -> (n_voxels, n_v_obs, n_transitions)

    dust file shape -> (n_voxels, n_wavelengths)
    '''
    rt.tempSpeciesEmissivity = fits.open(directory+'species_emissivity.fits', mode='denywrite')
    rt.tempSpeciesAbsorption = fits.open(directory+'species_absorption.fits', mode='denywrite')
    rt.tempDustEmissivity = fits.open(directory+'dust_emissivity.fits', mode='denywrite')
    rt.tempDustAbsorption = fits.open(directory+'dust_absorption.fits', mode='denywrite')
    return


def calculateObservation(directory='', dim='xy', slRange=[(-np.pi,np.pi), (-np.pi/2,np.pi/2)], nsl=[50,25],
                         terminal=True, plotV=False, multiprocessing=0, debug=False, verbose=False):

    if debug:
        sl = [5, 5]
    
    # constants.velocityRange = np.linspace(-300, 300, 500)
    
    # print('Load data')

    rt.open_voxel_positions(directory)
    rt.open_voxel_velocities(directory)
    rt.open_voxel_emission(directory)
  
    nDust = rt.tempDustEmissivity[0].shape[2]
    if nDust > 1:
        # default interpolation is along the last axis, which is what we need
        rt.interpDust = interpolate.interp1d(constants.wavelengths[:nDust],
                                                            rt.tempDustEmissivity[0].data[:, 0, :],
                                                            fill_value='extrapolate')
    
    constants.velocityRange = np.linspace(rt.tempSpeciesEmissivity[0].header['CRVAL2'] -
                                          (rt.tempSpeciesEmissivity[0].header['CRPIX2']) *
                                          rt.tempSpeciesEmissivity[0].header['CDELT2'],
                                          rt.tempSpeciesEmissivity[0].header['CRVAL2'] +
                                          (rt.tempSpeciesEmissivity[0].header['CRPIX2']) *
                                          rt.tempSpeciesEmissivity[0].header['CDELT2'],
                                          num=rt.tempSpeciesEmissivity[0].header['NAXIS2'])
    
    observations.methods.initialise()
    species.addMolecules(rt.tempSpeciesEmissivity[0].header['SPECIES'].split(', '))
    # print(rt.tempSpeciesEmissivity[0].header['SPECIES'].split(', '), '\n', species.moleculeWavelengths)
    
    # print('Data loaded :-)')
  
    xPositions, yPositions, zPositions = rt.voxelPositions[0].data[:, 0], \
        rt.voxelPositions[0].data[:, 1], \
        rt.voxelPositions[0].data[:, 2]
    r = np.sqrt((xPositions-constants.rGalEarth)**2 + yPositions**2)
  
    radGrid = np.sqrt((xPositions-constants.rGalEarth)**2 + yPositions**2 + zPositions**2)
    lonGrid = np.arctan2(yPositions, -(xPositions-constants.rGalEarth))
    rPolar = np.sqrt((xPositions-constants.rGalEarth)**2+yPositions**2)
    latGrid = np.arctan2(zPositions, rPolar)
  
    np.set_printoptions(threshold=100000)
    # print('\nLongitude\n', lonGrid, '\nLattitude\n', latGrid)
  
    # print('\nx\n', np.unique(xArray), '\ny\n', np.unique(yArray), '\nz\n', np.unique(zArray))
  
    if constants.fromEarth:
        # For an observation from Earth, the data is modified to account for Earth's position at (8750, 0, 0) pc.
        #  The integrated intensity is then calculated in the y-z plane to account for different viewing angles across
        #  the galaxy. This can be post-processed to convert to galactic coordinates.
    
        hdul = fits.HDUList()
    
        # Define the boundaries separating the inner and outer disk
        xBoundary = (xPositions > 0) & (r > constants.rGalEarth)
        yBoundary = (yPositions < constants.rGalEarth) & (yPositions > -constants.rGalEarth)
    
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
            velocityRange = constants.velocityRange
    
        result = multiprocessCalculation(slRange=slRange, nsl=nsl, multiprocessing=multiprocessing,
                                         dim=dim, debug=debug)
        
        Vpositions = result[0]
        VintensityMapSpecies = result[1]
        VintensityMapDust = result[2]
        rt.sightlines = np.asarray(result[3]).max(0)
        vmin, vmax = result[4]
        
        # Save sightline lengths
        np.savetxt(directory+'/sightlines.csv', rt.sightlines, delimiter=',')
        
        # Convert to numpy arrays
        Vpositions = np.array(Vpositions[0])
        VintensityMapSpecies = np.array(VintensityMapSpecies)
        VintensityMapDust = np.array(VintensityMapDust)
    
        if verbose:
            rt.logger.info('Map position shape:')
            if Vpositions.ndim > 1:
                rt.logger.info(Vpositions.shape)
            else:
                for p in VmapPositions:
                    rt.logger.info(p.shape)
              
            rt.logger.info('Map intensity shapes:')
            rt.logger.info('Species')
            if VintensityMapSpecies.ndim > 1:
                rt.logger.info(VintensityMapSpecies.shape)
            else:
                for intensity in intensityMapSpecies:
                    rt.logger.info(intensity.shape)
            rt.logger.info('Dust')
            if VintensityMapDust.ndim > 1:
                rt.logger.info(VintensityMapDust.shape)
            else:
                for intensity in intensityMapDust:
                    rt.logger.info(intensity.shape)
    
        # Setup the data to be saved in a FITS file. It will create an HDU list with position, species, and dust HDUs.
        if not debug:
            # Create HDUs for the map position and intensity and add the velocity in the headers
            PositionHDU = fits.ImageHDU(Vpositions)
            
            # print(VintensityMapSpecies.shape, np.shape(VintensityMapSpecies[0, 0, 0]))
            rt.intensity_species = VintensityMapSpecies
            IntensityHDUSpecies = fits.ImageHDU(VintensityMapSpecies)
            
            # print(VintensityMapDust.shape, np.shape(VintensityMapDust[0, 0, 0]))
            rt.intensity_dust = VintensityMapDust
            IntensityHDUDust = fits.ImageHDU(VintensityMapDust)
            # input('Integrated intensity shape: {}'.format(VintensityMap.shape))
      
            PositionHDU.header['TYPE1'] = 'Angle'
            PositionHDU.header['TYPE2'] = 'Position'
            PositionHDU.header['DIREC'] = 'Radial'
      
            IntensityHDUSpecies.header['TYPE'] = 'Species transitions'
            IntensityHDUSpecies.header['BUNIT'] = 'K'
            IntensityHDUSpecies.header['CTYPE1'] = 'Wavelength'
            IntensityHDUSpecies.header['CUNIT1'] = 'm'
            IntensityHDUSpecies.header['CRVAL1'] = 'N/A'
            IntensityHDUSpecies.header['CDELT1'] = 'N/A'
            IntensityHDUSpecies.header['CRPIX1'] = 'N/A'
            IntensityHDUSpecies.header['CTYPE2'] = 'GLON'
            IntensityHDUSpecies.header['CUNIT2'] = 'rad'
            if IntensityHDUSpecies.header['NAXIS2'] > 1:
                IntensityHDUSpecies.header['CRVAL2'] = (slRange[0][1]+slRange[0][0])/2.
                IntensityHDUSpecies.header['CDELT2'] = ((slRange[0][1]-slRange[0][0]) /
                                                        (IntensityHDUSpecies.header['NAXIS2']-1))
                IntensityHDUSpecies.header['CRPIX2'] = (IntensityHDUSpecies.header['NAXIS2'])/2.
            else:
                IntensityHDUSpecies.header['CRVAL2'] = slRange[0][0]
                IntensityHDUSpecies.header['CDELT2'] = 0
                IntensityHDUSpecies.header['CRPIX2'] = (IntensityHDUSpecies.header['NAXIS2'])
            IntensityHDUSpecies.header['CTYPE3'] = 'GLAT'
            IntensityHDUSpecies.header['CUNIT3'] = 'rad'
            if IntensityHDUSpecies.header['NAXIS3'] > 1:
                IntensityHDUSpecies.header['CRVAL3'] = (slRange[1][1]+slRange[1][0])/2.
                IntensityHDUSpecies.header['CDELT3'] = ((slRange[1][1]-slRange[1][0]) /
                                                        (IntensityHDUSpecies.header['NAXIS3']-1))
                IntensityHDUSpecies.header['CRPIX3'] = (IntensityHDUSpecies.header['NAXIS3'])/2.
            else:
                IntensityHDUSpecies.header['CRVAL3'] = slRange[1][0]
                IntensityHDUSpecies.header['CDELT3'] = 0
                IntensityHDUSpecies.header['CRPIX3'] = (IntensityHDUSpecies.header['NAXIS3'])
            IntensityHDUSpecies.header['CTYPE4'] = 'Velocity'
            IntensityHDUSpecies.header['CUNIT4'] = 'km/s'
            IntensityHDUSpecies.header['CRVAL4'] = (vmax+vmin)/2.
            IntensityHDUSpecies.header['CDELT4'] = (vmax-vmin)/(IntensityHDUSpecies.header['NAXIS4']-1)
            IntensityHDUSpecies.header['CRPIX4'] = (IntensityHDUSpecies.header['NAXIS4'])/2.
            IntensityHDUSpecies.header['DIREC'] = 'Radial'
            IntensityHDUSpecies.header['SPECIES'] = rt.tempSpeciesEmissivity[0].header['SPECIES']
      
            IntensityHDUDust.header['TYPE'] = 'Dust continuum'
            IntensityHDUDust.header['BUNIT'] = 'K'
            IntensityHDUDust.header['CTYPE1'] = 'Wavelength'
            IntensityHDUDust.header['CUNIT1'] = 'm'
            IntensityHDUDust.header['CRVAL1'] = 'N/A'
            IntensityHDUDust.header['CDELT1'] = 'N/A'
            IntensityHDUDust.header['CRPIX1'] = 'N/A'
            IntensityHDUDust.header['CTYPE2'] = 'GLON'
            IntensityHDUDust.header['CUNIT2'] = 'rad'
            if IntensityHDUDust.header['NAXIS2'] > 1:
                IntensityHDUDust.header['CRVAL2'] = (slRange[0][1]+slRange[0][0])/2
                IntensityHDUDust.header['CDELT2'] = ((slRange[0][1]-slRange[0][0]) /
                                                     (IntensityHDUDust.header['NAXIS2']-1))
                IntensityHDUDust.header['CRPIX2'] = (IntensityHDUDust.header['NAXIS2'])/2.
            else:
                IntensityHDUDust.header['CRVAL2'] = slRange[0][0]
                IntensityHDUDust.header['CDELT2'] = 0
                IntensityHDUDust.header['CRPIX2'] = (IntensityHDUDust.header['NAXIS2'])
            IntensityHDUDust.header['CTYPE3'] = 'GLAT'
            IntensityHDUDust.header['CUNIT3'] = 'rad'
            if IntensityHDUDust.header['NAXIS3'] > 1:
                IntensityHDUDust.header['CRVAL3'] = (slRange[1][1]+slRange[1][0])/2
                IntensityHDUDust.header['CDELT3'] = ((slRange[1][1]-slRange[1][0]) /
                                                     (IntensityHDUDust.header['NAXIS3']-1))
                IntensityHDUDust.header['CRPIX3'] = (IntensityHDUDust.header['NAXIS3'])/2.
            else:
                IntensityHDUDust.header['CRVAL3'] = slRange[1][0]
                IntensityHDUDust.header['CDELT3'] = 0
                IntensityHDUDust.header['CRPIX3'] = (IntensityHDUDust.header['NAXIS3'])
            IntensityHDUDust.header['CTYPE4'] = 'Velocity'
            IntensityHDUDust.header['CUNIT4'] = 'km/s'
            IntensityHDUDust.header['CRVAL4'] = (vmax+vmin)/2.
            IntensityHDUDust.header['CDELT4'] = (vmax-vmin)/(IntensityHDUDust.header['NAXIS4']-1)
            IntensityHDUDust.header['CRPIX4'] = (IntensityHDUDust.header['NAXIS4'])/2.
            IntensityHDUDust.header['DIREC'] = 'Radial'
            IntensityHDUDust.header['DUST'] = rt.tempDustEmissivity[0].header['DUST']
      
            hdul.append(PositionHDU)
            hdul.append(IntensityHDUSpecies)
            hdul.append(IntensityHDUDust)
      
            hdul.writeto(directory+'/channel_intensity.fits', overwrite=True)
      
            print('Intensity map written successfully :-)')
    
        rt.voxelPositions.close()
        rt.voxelVelocities.close()
        rt.tempSpeciesEmissivity.close()
        rt.tempSpeciesAbsorption.close()
        rt.tempDustEmissivity.close()
        rt.tempDustAbsorption.close()
    
    return


def multiprocessCalculation(slRange=[(-np.pi,np.pi), (-np.pi/2,np.pi/2)], nsl=[50,25], multiprocessing=0,
                            dim='spherical', vel_pool=False, debug=False):
  
    Vpositions = []
    VintensityMapSpecies = []
    VintensityMapDust = []
    sightlines = []

    if vel_pool:
        velChannel = partial(calculateVelocityChannel, slRange=slRange, nsl=nsl, dim=dim, debug=debug,
                             multiprocess=multiprocessing, vel_pool=vel_pool)
        vNum = constants.velocityRange.size
    else:
        lon = np.linspace(slRange[0][0], slRange[0][1], num=nsl[0])
        lat = np.linspace(slRange[1][0], slRange[1][1], num=nsl[1])
        longrid, latgrid = (arr.flatten() for arr in np.meshgrid(lon, lat))
        sightlines = np.zeros((lon.size, lat.size))
        args = list(zip(enumerate()))
        calc_los = partial(calculateSightline, slRange=slRange, nsl=nsl, dim=dim, debug=debug,
                           multiprocess=multiprocessing, vel_pool=vel_pool)
    
    t0 = time()
    
    if multiprocessing:
        pool = Pool(processes=multiprocessing)
        chunksize = max(int(vNum/multiprocessing/100), 1)
        if vel_pool:
            intensity = pool.imap(velChannel, list(enumerate(constants.velocityRange)), chunksize)
        else:
            intensity = pool.imap(calc_los, list())
    else:
        intensity = []
        if vel_pool:
            vTqdm = tqdm(total=constants.velocityRange.size, desc='Observing velocity', miniters=1, dynamic_ncols=True)
        rt.slTqdm = tqdm(total=nsl[0]*nsl[1], desc='Sightline', miniters=1, dynamic_ncols=True)
        for iv in enumerate(constants.velocityRange):
            intensity.append(velChannel(iv))
            vTqdm.update()
        
    vmin = constants.velocityRange.max()
    vmax = constants.velocityRange.min()
    
    if multiprocessing:
        if vel_pool:
            vTqdm = tqdm(total=constants.velocityRange.size, desc='Observing velocity', miniters=1, dynamic_ncols=True)
        else:
            rt.slTqdm = tqdm(total=lon.size, desc='Sightline', miniters=1, dynamic_ncols=True)
    
    for n, i in enumerate(intensity):
        #   i.wait()
        #   print(i)
        # This allows for discarding the velocity channels without an emission (assuming vel_pool is True)
        if len(i[3]):
            if constants.velocityRange[n] < vmin:
                vmin = constants.velocityRange[n]
            if constants.velocityRange[n] > vmax:
                vmax = constants.velocityRange[n]
        
            Vpositions.append(i[0])
            VintensityMapSpecies.append(i[1])
            VintensityMapDust.append(i[2])
            sightlines.append(i[3])
            
            # intensity[i] = None
          
        if multiprocessing:
            if vel_pool:
                vTqdm.update()
            else:
                rt.slTqdm.update()
        
    # print('\n\nTotal evaluation time for {} sightlines and {} velocity channels: {}\n\n'.format(nsl[0]*nsl[1], vNum,
    #                                                                                             time()-t0))
    
    return (Vpositions, VintensityMapSpecies, VintensityMapDust, sightlines, (vmin, vmax))


def sightlength(x, l):
    return constants.rGalEarth**2 - constants.rGal**2 + x**2 - 2*constants.rGalEarth*x*np.cos(l)

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
    nDust = rt.tempDustEmissivity[0].shape[2]
    if nDust > 1:
        base = rt.interpDust(species.moleculeWavelengths)[:, :, 0]
    else:
        base = rt.tempDustEmissivity[0].data[0, i_vel, 0]
      
    rt.i_vox = ((rt.tempSpeciesEmissivity[0].data[:, i_vel, :] > base).any(1) |
                               rt.tempDustEmissivity[0].data[:, i_vel, :].any(1))
    
    # print('Voxels selected (shape={}):'.format(i_vox[i_vox].shape), time()-t0)
  
    # Update velocity progress bar
    # rt.vTqdm.update()
    
    # Reset sightline progress bar
    if multiprocess == 0:
        rt.slTqdm.reset()
  
    # Initialise the intensities and map
    position = []
    intensityMapSpecies = []
    intensityMapDust = []
  
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
    
        for i, lon in enumerate(longrid):
          
            # print('lon,lat before scipy call', lon, lat, ':', time()-t0)
      
            # Calculate sightline length
            Rslh = op.root_scalar(sightlength, args=lon, x0=constants.rGalEarth, x1=constants.rGal).root
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
                positionintensity_species.append(rt.intensity_species * normfactor)
                positionintensity_dust.append(rt.intensity_dust * normfactor)
            elif vox > 1:
                calculatert(scale=normfactor)
                positionintensity_species.append(rt.intensity_species)
                positionintensity_dust.append(rt.intensity_dust)
                # intensity_species = []
                # intensity_dust = []
            else:
                positionintensity_species.append(np.zeros(rt.tempSpeciesEmissivity[0].shape[-1]))
                positionintensity_dust.append(np.zeros(rt.tempDustEmissivity[0].shape[-1]))
            
            if multiprocess == 0:
                rt.slTqdm.update()
            
            if len(np.shape(positionintensity_species[-1])) > 1:
                rt.logger.error('Error {}'.format(np.shape(positionintensity_species[-1])))
                input()
        
        # Save intensities for each latitude
        intensityMapSpecies.append(positionintensity_species)
        intensityMapDust.append(positionintensity_dust)
    
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
    return (position, intensityMapSpecies, intensityMapDust, sightlines)


def calculateSightline(los, slRange=[(-np.pi,np.pi), (-np.pi/2,np.pi/2)], nsl=[50,25],
                             dim='spherical', debug=False, multiprocess=0):

    # # Convert the tuple to the desired index and velocity
    # i_vel = ivelocity[0]
    i_lon = los[0][0]
    lon = los[0][1]
    i_lat = los[1][0]
    lat = los[1][1]

    # if multiprocess:
    #   t0 = time()
    #   print('\nCalculating velocity channel at {:.2f} km\s'.format(ivelocity[1]))

    # # Setup the sightlines that are calculated
    # longrid = np.linspace(slRange[0][0], slRange[0][1], num=nsl[0])
    # latgrid = np.linspace(slRange[1][0], slRange[1][1], num=nsl[1])
    # sightlines = np.zeros((longrid.size, latgrid.size))

    # print('lon/lat arrays created:', time()-t0)

    # Find the voxels that exist at the observing velocity
    nDust = rt.tempDustEmissivity[0].shape[2]
    if nDust > 1:
        base = rt.interpDust(species.moleculeWavelengths)[:, :, 0]
    else:
        base = rt.tempDustEmissivity[0].data[0, i_vel, 0]

    rt.i_vox = ((rt.tempSpeciesEmissivity[0].data[:, i_vel, :] > base).any(1) |
                               rt.tempDustEmissivity[0].data[:, i_vel, :].any(1))

    # print('Voxels selected (shape={}):'.format(i_vox[i_vox].shape), time()-t0)

    # Update velocity progress bar
    # rt.vTqdm.update()

    # Reset sightline progress bar
    if multiprocess == 0:
        rt.slTqdm.reset()

    # Initialise the intensities and map
    position = []
    intensityMapSpecies = []
    intensityMapDust = []

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

        for i, lon in enumerate(longrid):

            # print('lon,lat before scipy call', lon, lat, ':', time()-t0)

            # Calculate sightline length
            Rslh = op.root_scalar(sightlength, args=lon, x0=constants.rGalEarth, x1=constants.rGal).root
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
                positionintensity_species.append(rt.intensity_species * normfactor)
                positionintensity_dust.append(rt.intensity_dust * normfactor)
            elif vox > 1:
                calculatert(scale=normfactor)
                positionintensity_species.append(rt.intensity_species)
                positionintensity_dust.append(rt.intensity_dust)
                # intensity_species = []
                # intensity_dust = []
            else:
                positionintensity_species.append(np.zeros(rt.tempSpeciesEmissivity[0].shape[-1]))
                positionintensity_dust.append(np.zeros(rt.tempDustEmissivity[0].shape[-1]))

            if multiprocess == 0:
                rt.slTqdm.update()

            if len(np.shape(positionintensity_species[-1])) > 1:
                rt.logger.error('Error {}'.format(np.shape(positionintensity_species[-1])))
                input()

        # Save intensities for each latitude
        intensityMapSpecies.append(positionintensity_species)
        intensityMapDust.append(positionintensity_dust)

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
    return (position, intensityMapSpecies, intensityMapDust, sightlines)


def setLOS(x=0, y=0, z=0, lon=0, lat=0, i_vox=[], i_vel=0, i_spe=None, i_dust=None,
           dim='xy', reverse=True, debug=False, verbose=False):
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
  
    scale = constants.voxel_size*constants.pc*100   # pc should be in cm
  
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
    xGrid, yGrid, zGrid = (rt.voxelPositions[0].data[:, 0],
                           rt.voxelPositions[0].data[:, 1],
                           rt.voxelPositions[0].data[:, 2],
                           )
  
    if dim == 'spherical':
        # Set sightline position
        x1LoS = lon
        x2LoS = lat
    
        # Convert voxel positions to spherical
        radGrid = np.sqrt((xGrid-constants.rGalEarth)**2 + yGrid**2 + zGrid**2)
        lonGrid = np.arctan2(yGrid, -(xGrid-constants.rGalEarth))
        if lon < 0:
            lonGrid[lonGrid > 0] = lonGrid[lonGrid > 0] - 2*np.pi
        if lon > 0:
            lonGrid[lonGrid < 0] = lonGrid[lonGrid < 0] + 2*np.pi
        rPolar = np.sqrt((xGrid-constants.rGalEarth)**2+yGrid**2)
        latGrid = np.arctan2(zGrid, rPolar)
        if lat < 0:
            latGrid[latGrid > 0] = latGrid[latGrid > 0] - np.pi
        if lat > 0:
            latGrid[latGrid < 0] = latGrid[latGrid < 0] + np.pi
    
        # Choose the voxels in the sightline
        # adjustments for voxel orientation
        scaling = np.sqrt(2)  # 2#
        width = scaling*constants.voxel_size*np.max([np.sin(np.abs(lonGrid-np.pi/4)),
                                                     np.sin(np.abs(lonGrid+np.pi/4))], axis=0)
        height = scaling*constants.voxel_size*np.max([np.sin(np.abs(latGrid-np.pi/4)),
                                                      np.sin(np.abs(latGrid+np.pi/4))], axis=0)
        # angular size of voxels
        d_lon = 0.5*np.arctan(width/radGrid)
        d_lon[(lat > 1.0) | (lat < -1.0)] = np.pi
        d_lat = 0.5*np.arctan(height/radGrid)
        d_lat[(lat > 1.0) | (lat < -1.0)] = np.pi/2.
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
    
    if i_los.size == 1:
  
        rt.intensity_species = scale * rt.tempSpeciesEmissivity[0].data[i_los, i_vel, :][0, :] / constants.pc/100
        rt.intensity_dust = scale * rt.tempDustEmissivity[0].data[i_los, i_vel, :][0, :] / constants.pc/100
        
        # print('\n\n', (rt.tempSpeciesEmissivity[0].data[iLOS,i_vel,:]>0).any(), '\n\n')
    
        vel = rt.voxelVelocities[0].data[i_los]
        
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
            # order the voxels in the line-of-sight according to increasing radial distance
            i = np.argsort(x3LoS)
            iLoS_ordered = i_los[i]
        else:
            # order the voxels in the line-of-sight according to decreasing radial distance
            i = np.argsort(x3LoS)[::-1]
            iLoS_ordered = i_los[i]
        
        # print('\n\n', iLoS, i, iLoS_ordered, i_vel,
        #       (rt.tempSpeciesEmissivity[0].data[iLoS_ordered,i_vel,:]>0).any(), '\n\n')
        # print(rt.i_vox)
        # print(rt.x3LoS)
        # input()
        
        if (i_spe is None) & (i_dust is None):
            rt.e_species = c.copy(rt.tempSpeciesEmissivity[0].data[iLoS_ordered, i_vel, :]) / constants.pc/100
            rt.de_species = (rt.e_species[1:]-rt.e_species[:-1])/(scale)
            rt.k_species = c.copy(rt.tempSpeciesAbsorption[0].data[iLoS_ordered, i_vel, :]) / constants.pc/100
            rt.dk_species = (rt.k_species[1:]-rt.k_species[:-1])/(scale)
            
            rt.e_dust = c.copy(rt.tempDustEmissivity[0].data[iLoS_ordered, :]) / constants.pc/100
            rt.de_dust = (rt.e_dust[1:]-rt.e_dust[:-1])/(scale)
            rt.k_dust = c.copy(rt.tempDustAbsorption[0].data[iLoS_ordered, :]) / constants.pc/100
            rt.dk_dust = (rt.k_dust[1:]-rt.k_dust[:-1])/(scale)
    
    else:
        return 0, 0  # ,[]
  
    vel = rt.voxelVelocities[0].data[i_los]
    
    rt.logger.info('voxels:', i)
    
    return 2, vel.size
    # ,(e_species,de_species,k_species,dk_species,
    #   e_dust,de_dust,k_dust,dk_dust)


def calculatert(scale=1, background_intensity=0., species=True, dust=True, verbose=False, test=True):

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
    rt.intensity_species = np.full(rt.k_species.shape[1], background_intensity)
    rt.intensity_dust = np.full(rt.k_dust.shape[1], background_intensity)
    n_steps = rt.de_species.shape[0]
    
    # Adjust according to the average voxel depth
    if species:
        rt.e_species /= scale
        rt.de_species /= scale**2
        rt.k_species /= scale
        rt.dk_species /= scale**2
    if dust:
        rt.e_dust /= scale
        rt.de_dust /= scale**2
        rt.k_dust /= scale
        rt.dk_dust /= scale**2
  
    scale = scale*constants.voxel_size*constants.pc*100  # pc should be in cm
    # print(scale)
    
    np.set_printoptions(threshold=100000)
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
            a_species = (rt.k_species[:-1, :] / np.sqrt(2*np.abs(rt.dk_species)))
            b_species = ((rt.k_species[:-1, :] + rt.dk_species*scale) / np.sqrt(2*np.abs(rt.dk_species)))
        if dust:
            a_dust = (rt.k_dust[:-1, :] / np.sqrt(2*np.abs(rt.dk_dust)))
            b_dust = ((rt.k_dust[:-1, :] + rt.dk_dust*scale) / np.sqrt(2*np.abs(rt.dk_dust)))
    
    rt.logger.info(rt.k_species.shape[0])
        
    # if test:
    #     print(n_steps)
    
    for i in range(n_steps):
      
        # << Compute radiative transfer for the species transitions >>
        if species:
    
            # Determine which form of integration the species transitions require
            k0_species = ((rt.dk_species[i, :] == 0) &
                          (np.abs(rt.k_species[i, :]*scale) < 10**-10))
            kg_species = ~k0_species & (np.abs(rt.k_species[i, :])
                                        > (10**3*np.abs(rt.dk_species[i, :])*scale))
            keg_species = ~k0_species & ~kg_species & (rt.dk_species[i, :] > 0)
            kel_species = ~k0_species & ~kg_species & (rt.dk_species[i, :] < 0)
            
            # Integrate using the required method
            if k0_species.any():
                rt.intensity_species[k0_species] = (rt.e_species[i, k0_species] * scale
                                                    + 0.5*rt.de_species[i, k0_species] * scale**2.
                                                    + rt.intensity_species[k0_species])
        
            if kg_species.any():
                rt.intensity_species[kg_species] = ((rt.e_species[i, kg_species]*rt.k_species[i, kg_species]
                                                    + rt.de_species[i, kg_species]
                                                    * (rt.k_species[i, kg_species]*scale-1.))
                                                    / (rt.k_species[i, kg_species]**2.)) \
                                                  - np.exp(-rt.k_species[i, kg_species]*scale)\
                                                  * ((rt.e_species[i, kg_species]*rt.k_species[i, kg_species]
                                                      - rt.de_species[i, kg_species])
                                                     / (rt.k_species[i, kg_species]**2.)) \
                                                  + (rt.intensity_species[kg_species]
                                                     * np.exp(-rt.k_species[i, kg_species]*scale))
      
            if keg_species.any():
                rt.logger.info('\na, b:\n{}\n{}'.format(a_dust[i, :], b_dust[i, :]))
                a_error = erfi(a_species[i, keg_species])
                b_error = erfi(b_species[i, keg_species])
                rt.intensity_species[keg_species] = (rt.de_species[i, keg_species]/rt.dk_species[i, keg_species]
                                                     * (1.-np.exp(-rt.k_species[i, keg_species]*scale
                                                                  - 0.5*rt.dk_species[i, keg_species]*scale**2.))
                                                     - (rt.e_species[i, keg_species]*rt.dk_species[i, keg_species]
                                                        - rt.k_species[i, keg_species]*rt.de_species[i, keg_species])
                                                     / rt.dk_species[i, keg_species]
                                                     * np.sqrt(np.pi/2./np.abs(rt.dk_species[i, keg_species]))
                                                     * np.exp(-b_species[i, keg_species]**2.) * (a_error - b_error)
                                                     + rt.intensity_species[keg_species]
                                                     * np.exp(-rt.k_species[i, keg_species]*scale
                                                              - 0.5*rt.dk_species[i, keg_species]*scale**2.))
      
            if kel_species.any():
                rt.logger.info('\na, b:\n{}\n{}'.format(a_dust[i, :], b_dust[i, :]))
                a_error = erfc(a_species[i, :][kel_species])
                b_error = erfc(b_species[i, :][kel_species])
                rt.intensity_species[kel_species] = (rt.de_species[i, kel_species]/rt.dk_species[i, kel_species]
                                                     * (1.-np.exp(-rt.k_species[i, kel_species]*scale
                                                                  - 0.5*rt.dk_species[i, kel_species]*scale**2.))
                                                     - (rt.e_species[i, kel_species]*rt.dk_species[i, kel_species]
                                                        - rt.k_species[i, kel_species]*rt.de_species[i, kel_species])
                                                     / rt.dk_species[i, kel_species]
                                                     * np.sqrt(np.pi/2./np.abs(rt.dk_species[i, kel_species]))
                                                     * np.exp(b_species[i, kel_species]**2.) * (a_error - b_error)
                                                     + rt.intensity_species[kel_species]
                                                     * np.exp(-rt.k_species[i, kel_species]*scale
                                                              - 0.5*rt.dk_species[i, kel_species]*scale**2.))
        
        # << Compute radiative transfer for the dust continuum >>
        if dust:
    
            # Determine which form of integration the dust continuum requires
            k0_dust = ((rt.dk_dust[i, :] == 0) &
                       (abs(rt.k_dust[:-1, :][i, :]*scale) < 10**-10))
            kg_dust = ~k0_dust & (np.abs(rt.k_dust[:-1, :][i, :])
                                  > (10**3*abs(rt.dk_dust[i, :])*scale))
            keg_dust = ~k0_dust & ~kg_dust & (rt.dk_dust[i, :] > 0)
            kel_dust = ~k0_dust & ~kg_dust & (rt.dk_dust[i, :] < 0)
            
            # Integrate using the required method
            if k0_dust.any():
                rt.intensity_dust[k0_dust] = (rt.e_dust[i, k0_dust]*scale
                                              + 0.5*rt.de_dust[i, k0_dust]*scale**2.
                                              + rt.intensity_dust[k0_dust])
        
            if kg_dust.any():
                rt.intensity_dust[kg_dust] = ((rt.e_dust[i, kg_dust]*rt.k_dust[i, kg_dust]
                                               + rt.de_dust[i, kg_dust]*(rt.k_dust[i, kg_dust]*scale-1.))
                                              / (rt.k_dust[i, kg_dust]**2.)) \
                                             - np.exp(-rt.k_dust[i, kg_dust]*scale) \
                                             * ((rt.e_dust[i, kg_dust]*rt.k_dust[i, kg_dust]
                                                 - rt.de_dust[i, kg_dust])
                                                / (rt.k_dust[i, kg_dust]**2.)) \
                                             + (rt.intensity_dust[kg_dust]
                                                * np.exp(-rt.k_dust[i, kg_dust]*scale))
      
            if keg_dust.any():
                rt.logger.info('\na, b:\n{}\n{}'.format(a_dust[i, :], b_dust[i, :]))
                a_error = erfi(a_dust[i, keg_dust])
                b_error = erfi(b_dust[i, keg_dust])
                rt.intensity_dust[keg_dust] = (rt.de_dust[i, keg_dust]/rt.dk_dust[i, keg_dust]
                                               * (1.-np.exp(-rt.k_dust[i, keg_dust]*scale
                                                            - 0.5*rt.dk_dust[i, keg_dust]*scale**2.))
                                               - (rt.e_dust[i, keg_dust]*rt.dk_dust[i, keg_dust]
                                                  - rt.k_dust[i, keg_dust]*rt.de_dust[i, keg_dust])
                                               / rt.dk_dust[i, keg_dust]
                                               * np.sqrt(np.pi/2./np.abs(rt.dk_dust[i, keg_dust]))
                                               * np.exp(-b_dust[i, keg_dust]**2.) * (a_error - b_error)
                                               + rt.intensity_dust[keg_dust]
                                               * np.exp(-rt.k_dust[i, keg_dust]*scale
                                                        - 0.5*rt.dk_dust[i, keg_dust]*scale**2.))
      
            if kel_dust.any():
                rt.logger.info('\na, b:\n{}\n{}'.format(a_dust[i, :], b_dust[i, :]))
                a_error = erfc(a_dust[i, kel_dust])
                b_error = erfc(b_dust[i, kel_dust])
                rt.intensity_dust[kel_dust] = (rt.de_dust[i, kel_dust]/rt.dk_dust[i, kel_dust]
                                               * (1.-np.exp(-rt.k_dust[i, kel_dust]*scale
                                                            - 0.5*rt.dk_dust[i, kel_dust]*scale**2.))
                                               - (rt.e_dust[i, kel_dust] * rt.dk_dust[i, kel_dust]
                                                  - rt.k_dust[i, kel_dust]*rt.de_dust[i, kel_dust])
                                               / rt.dk_dust[i, kel_dust]
                                               * np.sqrt(np.pi/2./np.abs(rt.dk_dust[i][kel_dust]))
                                               * np.exp(b_dust[i, kel_dust]**2.) * (a_error - b_error)
                                               + rt.intensity_dust[kel_dust]
                                               * np.exp(-rt.k_dust[i, kel_dust]*scale
                                                        - 0.5*rt.dk_dust[i, kel_dust]*scale**2.))
    
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
  
    rt.voxelPositions = fits.open(directory+'/voxel_position.fits', mode='denywrite')
    rt.voxelVelocities = fits.open(directory+'/voxel_velocity.fits', mode='denywrite')
    rt.tempSpeciesEmissivity = fits.open(directory+'/species_emissivity.fits', mode='denywrite')
    rt.tempSpeciesAbsorption = fits.open(directory+'/species_absorption.fits', mode='denywrite')
    rt.tempDustEmissivity = fits.open(directory+'/dust_emissivity.fits', mode='denywrite')
    rt.tempDustAbsorption = fits.open(directory+'/dust_absorption.fits', mode='denywrite')
  
    multiprocessCalculation(slRange=[(-np.pi, np.pi), (-np.pi/2, np.pi/2)], nsl=[50, 25], dim='spherical',
                            multiprocessing=2)
    
    # multiprocessing.freeze_support()

# jit_module(nopython=False)
