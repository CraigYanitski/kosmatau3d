"""
A subpackage containing various plotting functions both for individual models
and :code:`kosmatau3d` in general.
The model plotting functions have mostly been depreciated by the methods in the
:code:`SyntheticModel` class.
"""

from pprint import pprint

import numpy as np
import os
from astropy.io import fits
import matplotlib.pyplot as plt
import matplotlib.colors as colors
from scipy.interpolate import interp2d

if "READTHEDOCS" not in os.environ:
    from .viewMap import *
from .diagram import *
from kosmatau3d.models import constants


def plotModel(plot='total intensity', ce=[], ie=[], grid=None,
              directory='/home/craig/projects/pdr/KOSMA-tau^3/history/MilkyWay/r1000.0_s3015/',
              species='13C+ 1', debug=False):
    """A depreciated function to plot the 3D distribution of a quantity."""
  
    allSpecies = ['13C 1', '13C 2', '13C 3', '13C+ 1', 'C 1', 'C 2', 'C+ 1', 'CO 1', 'CO 2', 'CO 3', 'CO 4', 'CO 5', 'CO 6', 'CO 7', '13CO 1', '13CO 2', '13CO 3', '13CO 4', '13CO 5', '13CO 6', '13CO 7', 'O 2', 'Dust 1', 'Dust 2', 'Dust 3', 'Dust 4', 'Dust 5', 'Dust 6', 'Dust 7', 'Dust 8', 'Dust 9', 'Dust 10', 'Dust 11', 'Dust 12', 'Dust 13', 'Dust 14', 'Dust 15', 'Dust 16', 'Dust 17', 'Dust 18', 'Dust 19', 'Dust 20', 'Dust 21', 'Dust 22', 'Dust 23', 'Dust 24', 'Dust 25', 'Dust 26', 'Dust 27', 'Dust 28', 'Dust 29', 'Dust 30', 'Dust 31', 'Dust 32', 'Dust 33', 'Dust 34', 'Dust 35', 'Dust 36', 'Dust 37', 'Dust 38', 'Dust 39', 'Dust 40']
  
    directory = constants.HISTORYPATH + constants.directory + directory
  
    if grid:
  
        voxelPositions = []
        mass = []
        density = []
        fuv = []
        FUVabsorption = []
        velocity = []
        clumpVelocity = []
        voxelSpeciesEmissivity = []
        voxelSpeciesAbsorption = []
        voxelSpeciesIntensity = []
        voxelDustEmissivity = []
        voxelDustAbsorption = []
        voxelDustIntensity = []
    
        for voxel in grid.allVoxels():
            voxelPositions.append(voxel.getPosition())
            mass.append(voxel.getEnsembleMass())
            density.append(voxel.getDensity())
            fuv.append(voxel.getFUV())
            FUVabsorption.append(voxel.getFUVabsorption())
            velocity.append(voxel.getVelocity())
            clumpVelocity.append(voxel.getClumpVelocity())
            voxelSpeciesEmissivity.append(voxel.getSpeciesEmissivity())
            voxelSpeciesAbsorption.append(voxel.getSpeciesAbsorption())
            voxelSpeciesIntensity.append(voxel.getSpeciesIntensity())
            voxelDustEmissivity.append(voxel.getDustEmissivity())
            voxelDustAbsorption.append(voxel.getDustAbsorption())
            voxelDustIntensity.append(voxel.getDustIntensity())
    
        voxelPositions = np.array(voxelPositions)
        mass = np.array(mass)[:, 0]
        density = np.array(density)[:, 0]
        fuv = np.array(fuv)[:, 0]
        FUVabsorption = np.array(FUVabsorption)
        velocity = np.array(velocity)
        clumpVelocity = np.array(clumpVelocity)
        voxelSpeciesEmissivity = np.array(voxelSpeciesEmissivity)
        voxelSpeciesAbsorption = np.array(voxelSpeciesAbsorption)
        voxelSpeciesIntensity = np.array(voxelSpeciesIntensity)
        voxelDustEmissivity = np.array(voxelDustEmissivity)
        voxelDustAbsorption = np.array(voxelDustAbsorption)
        voxelDustIntensity = np.array(voxelDustIntensity)
    
    elif len(ce) == 0 and len(ie) == 0:
  
        voxelPositions = fits.open(directory+'voxel_position.fits')[0].data
        fuv = fits.open(directory+r'voxel_fuv.fits')[0].data
        FUVabsorption = fits.open(directory+r'voxel_FUVabsorption.fits')[0].data
        velocity = fits.open(directory+r'voxel_velocity.fits')[0].data
        
        voxelSpeciesEmissivity = fits.open(directory+r'emissivity_clump_species.fits')[0].data
        voxelSpeciesAbsorption = fits.open(directory+r'absorption_clump_species.fits')[0].data
        voxelDustEmissivity = fits.open(directory+r'emissivity_clump_dust.fits')[0].data
        voxelDustAbsorption = fits.open(directory+r'absorption_clump_dust.fits')[0].data
    
        # maxClumpEmission = np.zeros(clumpEmission[:,:,:,0].shape)
        # maxInterclumpEmission = np.zeros(interclumpEmission[:,:,:,0].shape)
        # for emission in range(2):
        #   for voxel in range(clumpEmission[:,0,0,0].size):
        #     for species in range(clumpEmission[0,0,:,0].size):
        #       maxClumpEmission[voxel,emission,species] = norm.fit(clumpEmission[voxel,emission,species,:])[0]
        #       maxInterclumpEmission[voxel,emission,species] = norm.fit(interclumpEmission[voxel,emission,species,:])[0]
      
    else:
  
        voxelPositions = fits.open(directory+'voxel_position.fits')[0].data
        fuv = fits.open(directory+'voxel_fuv.fits')[0].data
        FUVabsorption = fits.open(directory+'voxel_FUVabsorption.fits')[0].data
        velocity = fits.open(directory+'voxel_velocity.fits')[0].data
    
        maxClumpEmission = ce
        maxInterclumpEmission = ie
  
    # positions = self.__grid.getVoxelPositions()
    limits = [voxelPositions.min(), voxelPositions.max()]
  
    if plot == 'total intensity':
        weights = (maxClumpEmission[:, 0, :]+maxInterclumpEmission[:, 0, :]).max(1)
        scale = r'$I \ (\chi)$'
        plotTitle = 'PDR Intensity within the Milky Way'
    
    elif plot == 'total optical depth':
        weights = (maxClumpEmission[:, 1, :]+maxInterclumpEmission[:, 1, :]).max(1)
        scale = r'$\tau$'
        plotTitle = 'PDR Optical Depth within the Milky Way'
    
    elif plot == 'species total intensity':
        i = allSpecies == species
        weights = voxelSpeciesEmissivity[:, :, i]*constants.voxel_size
        weights = weights/max(weights)
        scale = r'$I \ (\chi)$'
        plotTitle = species + ' intensity within the Milky Way'
    
    elif plot == 'species total optical depth':
        i = allSpecies == species
        weights = maxClumpEmission[:, 1, i]+maxInterclumpEmission[:, 1, i]
        weights = weights / weights.max()
        scale = r'$\tau$'
        plotTitle = species + ' optical depth within the Milky Way'
    
    elif plot == 'clump intensity':
        weights = (maxClumpEmission[:, 0, :]).sum(1)
        scale = r'$I \ (\chi)$'
        plotTitle = 'Clump intensity within the Milky Way'
    
    elif plot == 'clump optical depth':
        weights = (maxClumpEmission[:, 1, :]).sum(1)
        scale = r'$\tau$'
        plotTitle = 'Clump optical depth within the Milky Way'
    
    elif plot == 'interclump intensity':
        weights = (maxInterclumpEmission[:, 0, :]).sum(1)
        scale = r'$I \ (\chi)$'
        plotTitle = 'Interclump intensity within the Milky Way'
    
    elif plot == 'interclump optical depth':
        weights = (maxInterclumpEmission[:, 1, :]).sum(1)
        scale = r'$\tau$'
        plotTitle = 'Interclump optical depth within the Milky Way'
    
    elif plot == 'mass':
        weights = (mass)
        scale = r'$M_{vox} \ (M_\odot)$'
        plotTitle = 'Ensemble mass within the Milky Way'
    
    elif plot == 'density':
        weights = (density)
        scale = r'$n \ (cm^{-3})$'
        plotTitle = 'Density within the Milky Way'
    
    elif plot == 'FUV':
        weights = (fuv)
        scale = r'$FUV \ (\chi)$'
        plotTitle = 'FUV within the Milky Way'
    
    elif plot == 'FUV absorption':
        weights = (FUVabsorption)
        scale = r'$\tau_{FUV}$'
        plotTitle = r'$\tau_{FUV}$ within the Milky Way'
    
    elif plot == 'velocity':
        weights = (velocity)
        scale = r'$v_{rot} \ (\frac{km}{s})$'
        plotTitle = 'Rotational velocity within the Milky Way'
    
    fig = plt.figure()
    ax = fig.add_subplot(111, projection='3d')
    
    model = ax.scatter(voxelPositions[:, 0], voxelPositions[:, 1], voxelPositions[:, 2], c=weights.flatten(),
                       cmap=plt.cm.hot, marker='s', s=27, alpha=0.5, linewidths=0)
    
    ax.set_xlim(limits)
    ax.set_ylim(limits)
    ax.set_zlim(limits)
    cbar = plt.colorbar(model)
    
    ax.set_title(plotTitle)
    ax.set_xlabel('X (pc)')
    ax.set_ylabel('Y (pc)')
    ax.set_zlabel('Z (pc)')
    cbar.set_label(scale, rotation=0)
    
    plt.show()
    
    return maxClumpEmission, maxInterclumpEmission


def PVplot(directory='', file='/channel_intensity.fits', latRange=[-np.pi/9,np.pi/9], species=[], dust=[], newMap=[],
           save=False, verbose=False, title=''):
    """A depreciated function to plot the PV diagram for a model."""
  
    channelMap = fits.open(directory+file)
    
    if species:
        print('Species channel map')
        speciesMap = channelMap[1].data
    if dust:
        print('Dust channel map')
        dustMap = channelMap[2].data
  
    # pprint(channelMap[1].header)
    vel = np.linspace(channelMap[1].header['CRVAL4'] - (channelMap[1].header['CRPIX4']-.5) * channelMap[1].header['CDELT4'],
                      channelMap[1].header['CRVAL4'] + (channelMap[1].header['CRPIX4']-.5) * channelMap[1].header['CDELT4'],
                      num=channelMap[1].header['NAXIS4'])
    lat = np.linspace(channelMap[1].header['CRVAL3'] - (channelMap[1].header['CRPIX3']-.5) * channelMap[1].header['CDELT3'],
                      channelMap[1].header['CRVAL3'] + (channelMap[1].header['CRPIX3']-.5) * channelMap[1].header['CDELT3'],
                      num=channelMap[1].header['NAXIS3'])
    lon = np.linspace(channelMap[1].header['CRVAL2'] - (channelMap[1].header['CRPIX2']-.5) * channelMap[1].header['CDELT2'],
                      channelMap[1].header['CRVAL2'] + (channelMap[1].header['CRPIX2']-.5) * channelMap[1].header['CDELT2'],
                      num=channelMap[1].header['NAXIS2'])*180/np.pi
    
    velocity, longitude = np.meshgrid(vel, lon)
    if newMap:
        vel_new = np.linspace(vel.min(), vel.max(), newMap[1])
        lon_new = np.linspace(lon.min(), lon.max(), newMap[0])
        velocity_new, longitude_new = np.meshgrid(vel_new, lon_new)
    
    # print(velocity)
    # print(latitude)
    # print(longitude)
    i_min = np.floor(lat.size/2).astype(np.int)  # np.abs(lat-latRange[0]).argmin()
    i_max = np.floor(lat.size/2).astype(np.int)+1  # np.abs(lat-latRange[1]).argmin()
    print(i_max, i_min)
    
    if isinstance(species, str):
        species = [species]
    if isinstance(dust, str):
        dust = [dust]
  
    i_species = []
    i_dust = []
    
    for transition in species:
        print('Find species')
        allSpecies = channelMap[1].header['SPECIES'].split(', ')
        i_species.append(np.where(np.asarray(allSpecies) == transition)[0][0])
    for line in dust:
        print('Find dust')
        allDust = channelMap[2].header['DUST'].split(', ')
        i_dust.append(np.where(np.asarray(allDust) == line)[0][0])
    
    for i, i_transition in enumerate(i_species):
  
        fig, ax = plt.subplots(1, 1, figsize=(8, 5))
        intensityMap = speciesMap[:, i_min:i_max, :, i_transition]
        i_bg = np.where(speciesMap[:, i_min:i_max, :, i_transition] == channelMap[2].data[:, i_min:i_max, :, 0])
        intensityMap[i_bg] = np.nan
        intensityMap = intensityMap.sum(1)
        if verbose:
            print(np.shape(vel), np.shape(lon), np.shape(intensityMap))
    
        if verbose:
            print('intensityMap min/max:', np.nanmin(intensityMap), np.nanmax(intensityMap))
    
        if newMap:
            f = interp2d(vel, lon, intensityMap.T)
            intensityMap_new = f(vel_new, lon_new)
            intensityMap_new[intensityMap_new == 0] = np.nan
            cm = ax.pcolormesh(longitude_new, velocity_new, intensityMap_new, shading='nearest',
                               norm=colors.SymLogNorm(linthresh=0.001, vmin=np.nanmin(intensityMap),
                                                      vmax=np.nanmax(intensityMap)),
                               cmap='jet')
        else:
            intensityMap[intensityMap == 0] = np.nan
            cm = ax.pcolormesh(longitude, velocity, intensityMap.T, shading='nearest',
                               # norm=colors.SymLogNorm(linthresh=0.001, vmin=np.nanmin(intensityMap),
                               #                        vmax=np.nanmax(intensityMap)),
                               #                        cmap='jet')
                               norm=colors.SymLogNorm(linthresh=0.001, vmin=0.01, vmax=10), cmap='jet')
        cb = fig.colorbar(cm, ax=ax, fraction=0.02)
        # bounds = np.array([0.03, 0.05, 0.1, 0.5, 1, 6, 10])
        # cm = ax.pcolormesh(longitude*180/np.pi, velocity, intensityMap.T, shading='nearest',
        #                    norm=colors.BoundaryNorm(boundaries=bounds, ncolors=50), cmap='jet')
        # cb = fig.colorbar(cm, extend='max', ax=ax, fraction=0.02, ticks=[0.03, 0.5, 2.5, 6, 10])
        cb.ax.set_ylabel(r'Intensity $K$', fontsize=20, labelpad=20, rotation=270)
    
        ax.invert_xaxis()
        ax.set_xlabel(r'Longitude $\left( ^\circ \right)$', fontsize=20)
        ax.set_ylabel(r'Velocity $\left( \frac{km}{s} \right)$', fontsize=20)
        if title:
            ax.set_title(title, fontsize=26)
        else:
            ax.set_title(r'{} galactic position-velocity diagram'.format(species[i]), fontsize=26)
        if save:
            plt.savefig(directory+r'/plots/{}_pv_plot.png'.format(species[i].replace(' ', '_')))
            plt.close()
        else:
            plt.show(block=False)
  
    for i, i_line in enumerate(i_dust):
  
        fig, ax = plt.subplots(1, 1, figsize=(15, 10))
        intensityMap = dustMap[:, i_min:i_max, :, i_line].sum(1)
        if verbose:
            print(np.shape(vel), np.shape(lon), np.shape(intensityMap))
    
        if verbose:
            print('intensityMap min/max:', np.nanmin(intensityMap), np.nanmax(intensityMap))
        
        if newMap:
            f = interp2d(vel, lon, intensityMap.T)
            intensityMap_new = f(vel_new, lon_new)
            intensityMap_new[intensityMap_new == 0] = np.nan
            cm = ax.pcolormesh(longitude_new, velocity_new, intensityMap_new, shading='nearest',
                               norm=colors.SymLogNorm(linthresh=0.001, vmin=np.nanmin(intensityMap),
                                                      vmax=np.nanmax(intensityMap)),
                               cmap='jet')
        else:
            intensityMap[intensityMap == 0] = np.nan
            cm = ax.pcolormesh(longitude*180/np.pi, velocity, intensityMap.T,
                               norm=colors.SymLogNorm(linthresh=0.001, vmin=np.nanmin(intensityMap),
                                                      vmax=np.nanmax(intensityMap)),
                               cmap='jet')
        cb = fig.colorbar(cm, extend='max', ax=ax, fraction=0.02)
        cb.ax.set_ylabel(r'Intensity $log_{10} \left( K \right)$', fontsize=20, labelpad=20, rotation=270)
        
        ax.invert_xaxis()
        ax.set_xlabel(r'Longitude $\left( ^\circ \right)$', fontsize=20)
        ax.set_ylabel(r'Velocity $\left( \frac{km}{s} \right)$', fontsize=20)
        ax.set_title(r'Dust {} galactic position-velocity diagram'.format(dust[i]), fontsize=26)
        
        if save:
            plt.savefig(directory+r'/plots/{}_pv_plot.png'.format(dust[i].replace(' ', '_')))
            plt.close()
        else:
            plt.show(block=False)
    
    return
