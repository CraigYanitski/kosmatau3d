"""
This subsubpackage contains functions to regrid and resample observational data,
compute error, compare to grids of :code:`kosmatau3d` models, and plot the 
results.
"""

import cygrid
import numpy as np
from astropy.io import fits
from astropy.wcs import WCS
from astropy.convolution import Gaussian1DKernel, Gaussian2DKernel, convolve
from astropy.visualization.wcsaxes.frame import EllipticalFrame
from scipy.interpolate import interp1d, interp2d
from scipy.integrate import trapz
from scipy.io import readsav
from scipy.stats import binned_statistic
import matplotlib.pyplot as plt
import matplotlib.colors as colors
from matplotlib.ticker import StrMethodFormatter
from functools import lru_cache
from spectres import spectres
from pprint import pprint
from copy import copy, deepcopy
import os

from .observation import *
from kosmatau3d import models


sedigism_rms_window = {
    'G014_13CO21_Tmb_DR1.fits' : (0, 75),
    'G312_13CO21_Tmb_DR1.fits' : (-75, 0),
    'G326_13CO21_Tmb_DR1.fits' : (-125, 0),
    'G341_13CO21_Tmb_DR1.fits' : (-150, 0),
    'G000_13CO21_Tmb_DR1.fits' : (-175, 150),
    'G001_13CO21_Tmb_DR1.fits' : (-75, 150),
    'G002_13CO21_Tmb_DR1.fits' : (-75, 150),
    'G003_13CO21_Tmb_DR1.fits' : (-50, 150),
    'G004_13CO21_Tmb_DR1.fits' : (-50, 50),
    'G005_13CO21_Tmb_DR1.fits' : (-50, 50),
    'G006_13CO21_Tmb_DR1.fits' : (-50, 50),
    'G007_13CO21_Tmb_DR1.fits' : (-50, 75),
    'G008_13CO21_Tmb_DR1.fits' : (-50, 75),
    'G009_13CO21_Tmb_DR1.fits' : (-25, 75),
    'G010_13CO21_Tmb_DR1.fits' : (-25, 100),
    'G011_13CO21_Tmb_DR1.fits' : (-25, 100),
    'G012_13CO21_Tmb_DR1.fits' : (-25, 100),
    'G013_13CO21_Tmb_DR1.fits' : (-25, 75),
    'G015_13CO21_Tmb_DR1.fits' : (0, 75),
    'G016_13CO21_Tmb_DR1.fits' : (0, 75),
    'G017_13CO21_Tmb_DR1.fits' : (0, 100),
    'G030_13CO21_Tmb_DR1.fits' : (0, 150),
    'G301_13CO21_Tmb_DR1.fits' : (-75, 25),
    'G302_13CO21_Tmb_DR1.fits' : (-75, 25),
    'G303_13CO21_Tmb_DR1.fits' : (-75, 50),
    'G304_13CO21_Tmb_DR1.fits' : (-75, 50),
    'G305_13CO21_Tmb_DR1.fits' : (-75, 50),
    'G306_13CO21_Tmb_DR1.fits' : (-75, 0),
    'G307_13CO21_Tmb_DR1.fits' : (-75, 25),
    'G308_13CO21_Tmb_DR1.fits' : (-75, 25),
    'G309_13CO21_Tmb_DR1.fits' : (-75, 25),
    'G310_13CO21_Tmb_DR1.fits' : (-75, 25),
    'G311_13CO21_Tmb_DR1.fits' : (-75, 25),
    'G313_13CO21_Tmb_DR1.fits' : (-75, 25),
    'G314_13CO21_Tmb_DR1.fits' : (-75, 25),
    'G315_13CO21_Tmb_DR1.fits' : (-75, 25),
    'G316_13CO21_Tmb_DR1.fits' : (-75, 25),
    'G317_13CO21_Tmb_DR1.fits' : (-75, 25),
    'G318_13CO21_Tmb_DR1.fits' : (-75, 25),
    'G319_13CO21_Tmb_DR1.fits' : (-100, 25),
    'G320_13CO21_Tmb_DR1.fits' : (-100, 25),
    'G321_13CO21_Tmb_DR1.fits' : (-100, 25),
    'G322_13CO21_Tmb_DR1.fits' : (-100, 25),
    'G323_13CO21_Tmb_DR1.fits' : (-100, 25),
    'G324_13CO21_Tmb_DR1.fits' : (-100, 25),
    'G324_13CO21_Tmb_DR1.fits.1' : (-100, 25),
    'G325_13CO21_Tmb_DR1.fits' : (-100, 25),
    'G327_13CO21_Tmb_DR1.fits' : (-125, 0),
    'G328_13CO21_Tmb_DR1.fits' : (-125, 0),
    'G329_13CO21_Tmb_DR1.fits' : (-125, 0),
    'G330_13CO21_Tmb_DR1.fits' : (-125, 0),
    'G331_13CO21_Tmb_DR1.fits' : (-125, 0),
    'G332_13CO21_Tmb_DR1.fits' : (-150, 25),
    'G333_13CO21_Tmb_DR1.fits' : (-150, 25),
    'G334_13CO21_Tmb_DR1.fits' : (-150, 25),
    'G335_13CO21_Tmb_DR1.fits' : (-150, 25),
    'G336_13CO21_Tmb_DR1.fits' : (-150, 25),
    'G337_13CO21_Tmb_DR1.fits' : (-150, 25),
    'G338_13CO21_Tmb_DR1.fits' : (-150, 25),
    'G339_13CO21_Tmb_DR1.fits' : (-150, 25),
    'G340_13CO21_Tmb_DR1.fits' : (-150, 25),
    'G342_13CO21_Tmb_DR1.fits' : (-175, 25),
    'G343_13CO21_Tmb_DR1.fits' : (-175, 25),
    'G344_13CO21_Tmb_DR1.fits' : (-175, 25),
    'G345_13CO21_Tmb_DR1.fits' : (-175, 25),
    'G346_13CO21_Tmb_DR1.fits' : (-175, 25),
    'G347_13CO21_Tmb_DR1.fits' : (-150, 25),
    'G348_13CO21_Tmb_DR1.fits' : (-150, 25),
    'G349_13CO21_Tmb_DR1.fits' : (-150, 25),
    'G350_13CO21_Tmb_DR1.fits' : (-150, 25),
    'G351_13CO21_Tmb_DR1.fits' : (-125, 25),
    'G352_13CO21_Tmb_DR1.fits' : (-125, 50),
    'G353_13CO21_Tmb_DR1.fits' : (-150, 25),
    'G354_13CO21_Tmb_DR1.fits' : (-100, 150),
    'G355_13CO21_Tmb_DR1.fits' : (-100, 150),
    'G356_13CO21_Tmb_DR1.fits' : (-100, 150),
    'G357_13CO21_Tmb_DR1.fits' : (-100, 50),
    'G358_13CO21_Tmb_DR1.fits' : (-100, 50),
    'G359_13CO21_Tmb_DR1.fits' : (-175, 150)
}
thrumms_rms_window = {
    'dr3.s300.12co.fits' : (),
    'dr3.s306.12co.fits' : (),
    'dr3.s312.12co.fits' : (),
    'dr3.s318.12co.fits' : (),
    'dr3.s324.12co.fits' : (),
    'dr3.s330.12co.fits' : (),
    'dr3.s336.12co.fits' : (),
    'dr3.s342.12co.fits' : (),
    'dr3.s348.12co.fits' : (),
    'dr3.s354.12co.fits' : (),
    'dr4.s300.13co.fits' : (),
    'dr4.s306.13co.fits' : (),
    'dr4.s312.13co.fits' : (),
    'dr4.s318.13co.fits' : (),
    'dr4.s324.13co.fits' : (),
    'dr4.s330.13co.fits' : (),
    'dr4.s336.13co.fits' : (),
    'dr4.s342.13co.fits' : (),
    'dr4.s348.13co.fits' : (),
    'dr4.s354.13co.fits' : ()
}
sgps_rms_window = {
    'g010.hi.fits' : 1.8,
    'g015.hi.fits' : 1.4,
    'g258.hi.fits' : 1.4,
    'g268.hi.fits' : 1.4,
    'g278.hi.fits' : 1.3,
    'g288.hi.fits' : 1.6,
    'g298.hi.fits' : 1.4,
    'g308.hi.fits' : 1.5,
    'g318.hi.fits' : 1.6,
    'g328.hi.fits' : 1.7,
    'g338.hi.fits' : 1.9,
    'g348.hi.fits' : 1.9,
    'g353.hi.fits' : 2.6,
}
cobe_idl_linfrq = np.array([
    115.3, 230.5, 345.8, 424.8, 461.0, 492.2, 556.9, 576.3, 691.5, 808.1,
     1113,  1460,  2226,  1901,  2060,  2311,  2459,  2589, 921.8])
cobe_idl_transitions = np.array([
    'CO 1',       'CO 2', 'CO 3', 'CO 4',  'C 1', 'CO 5',
    'CO 6', 'CO 7 + C 2', 'C+ 1',  'O 1', 'CO 8'])


def determine_rms(hdul, mission='', file=''):
    '''
    determine RMS of selected surveys.
    '''
    import astrokit
    # print(mission, 'GOT C+', mission=='GOT C+')

    if mission == 'GOT C+':
        # Create velocity axis
        lon = astrokit.get_axis(1, hdul)
        lat = astrokit.get_axis(2, hdul)
        vel = astrokit.get_axis(3, hdul)

        # make an empty map
        map_size = np.zeros_like(hdul[0].data[:, :, 0])
        hdu = fits.PrimaryHDU(map_size)
        hdul_rms = fits.HDUList([hdu])
        hdul_rms[0].header = deepcopy(hdul[0].header)

        # remove 3d attributes from header
        for attribute in list(hdul_rms[0].header.keys()):
            if not (attribute == ''):
                if (attribute[-1] == '3'):
                    del hdul_rms[0].header[attribute]
                elif (attribute == 'NAXIS'):
                    hdul_rms[0].header['NAXIS'] = 2

        clean_noise = deepcopy(np.nan_to_num(hdul[0].data, nan=0))

        rms_mask=np.zeros_like(hdul[0].data)
        i_mask = (np.where((vel>=-150)&(vel<=50))[0].reshape(-1, 1, 1), np.where((lat<=90)&(lat>=-90))[0].reshape(1, -1, 1), np.where((lon<=-10))[0].reshape(1, 1, -1), )
        rms_mask[i_mask] = 1
        i_mask = (np.where((vel>=-100)&(vel<=150))[0].reshape(-1, 1, 1), np.where((lat<=90)&(lat>=-90))[0].reshape(1, -1, 1), np.where((lon>=-10)&(lon<=10))[0].reshape(1, 1, -1), )
        rms_mask[i_mask] = 1
        i_mask = (np.where((vel>=-50)&(vel<=200))[0].reshape(-1, 1, 1), np.where((lat<=90)&(lat>=-90))[0].reshape(1, -1, 1), np.where((lon>=10))[0].reshape(1, 1, -1), )
        rms_mask[i_mask] = 1

        clean_data = np.ma.masked_array(hdul[0].data, rms_mask)
        clean_noise[clean_data.mask] = np.nan
        rms = np.nanstd(clean_noise, axis=0)
        rms[rms == 0] = rms[rms != 0].mean()

        hdul_rms[0].data = deepcopy(rms)

    elif mission == 'COGAL':
        mean=True

        # Ensure data has the correct dimensions (in fits format: glon, glat, vel_lsr)
        cube_co_fix = astrokit.swap_cube_axis(hdul)

        # Create the velocity axis
        lon_co = astrokit.get_axis(1, cube_co_fix)
        lat_co = astrokit.get_axis(2, cube_co_fix)
        vel_co = astrokit.get_axis(3, cube_co_fix) / 1e3

        # make an empty map
        map_size = np.zeros_like(cube_co_fix[0].data[:, :, 0])
        hdu = fits.PrimaryHDU(map_size)
        hdul_rms = fits.HDUList([hdu])
        hdul_rms[0].header = deepcopy(cube_co_fix[0].header)

        # remove 3d attributes from header
        for attribute in list(hdul_rms[0].header.keys()):
            if not (attribute == ''):
                if (attribute[-1] == '3'):
                    del hdul_rms[0].header[attribute]
                elif (attribute == 'NAXIS'):
                    hdul_rms[0].header['NAXIS'] = 2

        if mean:
            med = np.nanmean(cube_co_fix[0].data, axis=0)
        #     print('mean', med, 'K')
        else:
            med = np.nanmedian(cube_co_fix[0].data, axis=0)
        #     print('median', med, 'K')
        # print('std for data < median:', cube_co_fix[0].data[cube_co_fix[0].data < med.reshape(1, *med.shape)].std(), 'K')

        clean_noise = deepcopy(np.nan_to_num(cube_co_fix[0].data, nan=0))

        rms_mask=np.zeros_like(cube_co_fix[0].data)
        # i_mask = raw_data > med
        # rms_mask[i_mask] = 1
        print('Shape', cube_co_fix[0].shape)
        i_mask = (np.where((vel_co>=-50)&(vel_co<=150))[0].reshape(-1, 1, 1), np.where((lat_co<=90)&(lat_co>=-90))[0].reshape(1, -1, 1), np.where((lon_co>=-180)&(lon_co<-70))[0].reshape(1, 1, -1), )
        rms_mask[i_mask] = 1
        i_mask = (np.where((vel_co>=-220)&(vel_co<=125))[0].reshape(-1, 1, 1), np.where((lat_co<=90)&(lat_co>=-90))[0].reshape(1, -1, 1), np.where((lon_co>=-70)&(lon_co<-7))[0].reshape(1, 1, -1), )
        rms_mask[i_mask] = 1
        i_mask = (np.where((vel_co>=-310)&(vel_co<=320))[0].reshape(-1, 1, 1), np.where((lat_co<=90)&(lat_co>=-90))[0].reshape(1, -1, 1), np.where((lon_co>=-7)&(lon_co<7))[0].reshape(1, 1, -1), )
        rms_mask[i_mask] = 1
        i_mask = (np.where((vel_co>=-80)&(vel_co<=200))[0].reshape(-1, 1, 1), np.where((lat_co<=90)&(lat_co>=-90))[0].reshape(1, -1, 1), np.where((lon_co>=7)&(lon_co<50))[0].reshape(1, 1, -1), )
        rms_mask[i_mask] = 1
        i_mask = (np.where((vel_co>=-60)&(vel_co<=75))[0].reshape(-1, 1, 1), np.where((lat_co<=90)&(lat_co>=-90))[0].reshape(1, -1, 1), np.where((lon_co>=50)&(lon_co<70))[0].reshape(1, 1, -1), )
        rms_mask[i_mask] = 1
        i_mask = (np.where((vel_co>=-90)&(vel_co<=25))[0].reshape(-1, 1, 1), np.where((lat_co<=90)&(lat_co>=-90))[0].reshape(1, -1, 1), np.where((lon_co>=70)&(lon_co<=180))[0].reshape(1, 1, -1), )
        rms_mask[i_mask] = 1

        clean_data = np.ma.masked_array(cube_co_fix[0].data, rms_mask)
        clean_noise[clean_data.mask] = np.nan
        rms = np.nanstd(clean_noise, axis=0)
        rms[rms == 0] = rms[rms != 0].mean()

        hdul_rms[0].data = deepcopy(rms)

    elif mission == 'SEDIGISM':
        # Create velocity axis
        vel = astrokit.get_axis(3, hdul)

        # Set velocity window of the emission
        window = sedigism_rms_window[file]

        # Calculate the rms noise
        hdul_rms = astrokit.rms_map(hdul, window=window, rms_range=None)

    else:
        print('{} not available for RMS calculation'.format(mission))
        return None

    return hdul_rms


def regrid_observations(path='/media/hpc_backup/yanitski/projects/' 
        + 'pdr/observational_data/MilkyWay/', 
        mission=None, nu_planck=713e9, target_header=None, target_kernel=None, 
        output_file='obs_regridded.fits', cal_error=0, sig_th=3):
    '''
    This function will regrid the specified mission data to the specified 
    target header. 
    This prepares the observational data for comparison to simulated data. 
    For now, the target header must

    :param path: The directory storing the mission folders. 
            Each mission folder should contain the relevant fits files.
    :param mission: The mission you wish to regrid. 
            The default functionality is to regrid all missions.
    :param target_header: The target header, which must have longitude on 
            axis 1 and latitude on axis 2. 
            An optional axis 3 may be used for the velocity. 
            If there is not a third dimension, the data will be 
            regridded with a spacing of 1 km/s.
    :param target_kernel: The kernel to use for regridding.
    :return: Float -- 0 for success; 1 for fail
    '''

    import healpy as hp

    if isinstance(mission, str):
        mission = [mission]

    if mission==None or mission[0]=='':
        if path[-1] != '/': path += '/'
        mission = os.listdir(path)
    elif mission[0] in os.listdir(path):
        print('Regridding {} survey'.format(mission[0]))
    else:
        print('Invalid mission... {} observations do not exist or ' 
                + 'are not downloaded'.format(mission[0]))
        return 1

    if target_header==None:
        print('Error: Please specify target_header to regrid the ' 
                + 'observations.')
        print('Regrid failed.')
        return 1

    if target_header['NAXIS'] == 3:
        target_vel = np.linspace(target_header['CRVAL3'] 
                - target_header['CDELT3'] * (target_header['CRPIX3'] - 1),
                target_header['CRVAL3'] + target_header['CDELT3'] 
                * (target_header['NAXIS3'] - target_header['CRPIX3']),
                num=target_header['NAXIS3'])
    else:
        print('Error: Please enter a target velocity to regrid the '
                + 'spectroscopic observations.')
        return 1
        target_vel = None

    min_vel = 0
    max_vel = 0

    CO1 = False
    CO2 = False
    iCO1 = False
    iCO2 = False
    hi = False
    dust = False
    cobe = False

    for survey in mission:

        # COBE data must be located in instrument-separated folders
        # if survey == 'COBE':
        #     survey = 'COBE-FIRAS'

        temp_header = copy(target_header)
        twod_header = copy(target_header)
        if 'CDELT3' in twod_header.keys():
            twod_header['NAXIS'] = 2
            del twod_header['NAXIS3']
            del twod_header['CTYPE3']
            del twod_header['CDELT3']
            del twod_header['CRVAL3']
            del twod_header['CRPIX3']

        print(survey)

        if survey[-1] == '/':
            survey[-1] = ''
        files = os.listdir(path + survey + '/')

        if os.path.exists(path + survey + '/regridded/temp/') == False:
            os.makedirs(path + survey + '/regridded/temp/')

        # Initialise cygrid gridders
        if survey == 'COGAL':
            co1_gridder = cygrid.WcsGrid(target_header)
            co1_gridder.set_kernel(*target_kernel)
            co1_gridder_err = cygrid.WcsGrid(twod_header)
            co1_gridder_err.set_kernel(*target_kernel)
            CO1 = True
        elif survey == 'Mopra':
            co1_gridder = cygrid.WcsGrid(target_header)
            co1_gridder.set_kernel(*target_kernel)
            co1_gridder_err = cygrid.WcsGrid(twod_header)
            co1_gridder_err.set_kernel(*target_kernel)
            CO1 = True
            ico1_gridder = cygrid.WcsGrid(target_header)
            ico1_gridder.set_kernel(*target_kernel)
            ico1_gridder_err = cygrid.WcsGrid(twod_header)
            ico1_gridder_err.set_kernel(*target_kernel)
            iCO1 = True
        elif survey == 'ThrUMMS':
            co1_gridder = cygrid.WcsGrid(target_header)
            co1_gridder.set_kernel(*target_kernel)
            co1_gridder_err = cygrid.WcsGrid(twod_header)
            co1_gridder_err.set_kernel(*target_kernel)
            CO1 = True
            ico1_gridder = cygrid.WcsGrid(target_header)
            ico1_gridder.set_kernel(*target_kernel)
            ico1_gridder_err = cygrid.WcsGrid(twod_header)
            ico1_gridder_err.set_kernel(*target_kernel)
            iCO1 = True
        elif survey == 'THOR':
            hi_gridder = cygrid.WcsGrid(target_header)
            hi_gridder.set_kernel(*target_kernel)
            hi_gridder_err = cygrid.WcsGrid(twod_header)
            hi_gridder_err.set_kernel(*target_kernel)
            hi = True
        elif survey == 'SEDIGISM':
            ico2_gridder = cygrid.WcsGrid(target_header)
            ico2_gridder.set_kernel(*target_kernel)
            if cal_error:
                ico2_gridder_err = cygrid.WcsGrid(target_header)
            else:
                ico2_gridder_err = cygrid.WcsGrid(twod_header)
            ico2_gridder_err.set_kernel(*target_kernel)
            iCO2 = True
        elif survey == 'Planck':
            dust = True
        elif survey == 'COBE-FIRAS' or survey == 'COBE-DIRBE':
            cobe = True
        else:
            print('Survey {} not available. Choose another.'.format(survey))
            continue

        for file in files:

            if file == 'regridded' or file == 'temp' or \
                    'RMS' in file or survey == 'HiGAL':
                continue

            # Grid data and RMS
            elif survey == 'COGAL' and 'interp' in file:
                print(file)

                # Specify transition
                transitions = ['CO 1']
                transition_indeces = [0]

                # Open file
                obs = fits.open(path + survey + '/' + file)

                # Get axes
                lon = np.linspace(obs[0].header['CRVAL2'] 
                        - obs[0].header['CDELT2'] 
                        * (obs[0].header['CRPIX2'] - 1),
                        obs[0].header['CRVAL2'] + obs[0].header['CDELT2'] 
                        * (obs[0].header['NAXIS2'] - obs[0].header['CRPIX2']),
                        num=obs[0].header['NAXIS2'])
                lat = np.linspace(obs[0].header['CRVAL3'] 
                        - obs[0].header['CDELT3'] 
                        * (obs[0].header['CRPIX3'] - 1),
                        obs[0].header['CRVAL3'] + obs[0].header['CDELT3'] 
                        * (obs[0].header['NAXIS3'] - obs[0].header['CRPIX3']),
                        num=obs[0].header['NAXIS3'])
                vel = np.linspace(obs[0].header['CRVAL1'] 
                        - obs[0].header['CDELT1'] 
                        * (obs[0].header['CRPIX1'] - 1),
                        obs[0].header['CRVAL1'] + obs[0].header['CDELT1'] 
                        * (obs[0].header['NAXIS1'] - obs[0].header['CRPIX1']),
                        num=obs[0].header['NAXIS1'])

                # Copy header
                temp_header['NAXIS'] = 3
                temp_header['CTYPE1'] = target_header['CTYPE1']
                temp_header['CTYPE2'] = target_header['CTYPE2']
                temp_header['NAXIS3'] = target_header['NAXIS3']
                temp_header['CTYPE3'] = target_header['CTYPE3']
                temp_header['CRVAL3'] = target_header['CRVAL3']
                temp_header['CDELT3'] = target_header['CDELT3']
                temp_header['CRPIX3'] = target_header['CRPIX3']

                obs_data = spectres(target_vel, vel, 
                        np.nan_to_num(obs[0].data, nan=0), fill=0, 
                        verbose=False)
                obs_data = np.nan_to_num(obs_data.reshape(-1, 
                    obs_data.shape[-1]), nan=0)
                obs_error = determine_rms(obs, 
                        mission=survey)[0].data.reshape(-1, 1)
                print(obs_error.shape)
                # obs_error = np.swapaxes(obs_error, 0, 2)
                # obs_error = np.swapaxes(obs_error, 0, 1)

                # Grid
                lon_mesh, lat_mesh = np.meshgrid(lon, lat)
                co1_gridder.grid(lon_mesh.flatten(), lat_mesh.flatten(), 
                        obs_data)
                co1_gridder_err.grid(lon_mesh.flatten(), lat_mesh.flatten(), 
                        obs_error[:, 0])

                if vel.min() < min_vel:
                    min_vel = vel.min()
                if vel.max() > max_vel:
                    max_vel = vel.max()
            elif survey == 'Mopra' and ('_Vfull.fits' in file) \
                    and (('12CO' in file) or ('13CO' in file)):
                print(file)

                transition = file.split('_')[0]

                # Specify transition
                if '12CO' in transition:
                    transitions = ['CO 1']
                elif '13CO' in transition:
                    transitions = ['13CO 1']
                transition_indeces = [0]

                # Open file
                obs = fits.open(path + survey + '/' + file)
                obs_error = np.nanmean(fits.open(path + survey + '/' 
                    + file.replace('_Vfull', '.sigma'))[0].data)

                # Get axes
                lon = np.linspace(obs[0].header['CRVAL1'] 
                        - obs[0].header['CDELT1'] 
                        * (obs[0].header['CRPIX1'] - 1),
                        obs[0].header['CRVAL1'] + obs[0].header['CDELT1'] 
                        * (obs[0].header['NAXIS1'] - obs[0].header['CRPIX1']),
                        num=obs[0].header['NAXIS1'])
                lat = np.linspace(obs[0].header['CRVAL2'] 
                        - obs[0].header['CDELT2'] 
                        * (obs[0].header['CRPIX2'] - 1),
                        obs[0].header['CRVAL2'] + obs[0].header['CDELT2'] 
                        * (obs[0].header['NAXIS2'] - obs[0].header['CRPIX2']),
                        num=obs[0].header['NAXIS2'])
                vel = np.linspace(obs[0].header['CRVAL3'] 
                        - obs[0].header['CDELT3'] 
                        * (obs[0].header['CRPIX3'] - 1),
                        obs[0].header['CRVAL3'] + obs[0].header['CDELT3'] 
                        * (obs[0].header['NAXIS3'] - obs[0].header['CRPIX3']),
                        num=obs[0].header['NAXIS3']) / 1e3

                # Copy header
                temp_header['NAXIS'] = 3
                temp_header['CTYPE1'] = target_header['CTYPE1']
                temp_header['CTYPE2'] = target_header['CTYPE2']
                temp_header['NAXIS3'] = target_header['NAXIS3']
                temp_header['CTYPE3'] = target_header['CTYPE3']
                temp_header['CRVAL3'] = target_header['CRVAL3']
                temp_header['CDELT3'] = target_header['CDELT3']
                temp_header['CRPIX3'] = target_header['CRPIX3']

                # Grid
                obs_data = np.swapaxes(obs[0].data, 0, 2)
                obs_data = np.swapaxes(obs_data, 0, 1)
                obs_data = spectres(target_vel, vel, 
                        np.nan_to_num(obs[0].data.T, nan=0), fill=0, 
                        verbose=False)
                obs_data = obs_data.reshape(-1, obs_data.shape[-1])
                obs_error = fits.open(path + survey + '/' 
                        + file.replace('_Vfull', '.sigma'))[0].data.reshape(-1)
                i_nan = np.isnan(obs_error)
                del obs

                # Grid
                lon_mesh, lat_mesh = np.meshgrid(lon, lat)
                if '12CO' in transition:
                    co1_gridder.grid(lon_mesh.flatten(), lat_mesh.flatten(), 
                            obs_data)
                    co1_gridder_err.grid(lon_mesh.flatten()[~i_nan.flatten()], 
                            lat_mesh.flatten()[~i_nan.flatten()],
                            obs_error[~i_nan])
                elif '13CO' in transition:
                    ico1_gridder.grid(lon_mesh.flatten(), lat_mesh.flatten(), 
                            obs_data)
                    ico1_gridder_err.grid(lon_mesh.flatten()[~i_nan.flatten()], 
                            lat_mesh.flatten()[~i_nan.flatten()],
                            obs_error[~i_nan])

                if vel.min() < min_vel:
                    min_vel = vel.min()
                if vel.max() > max_vel:
                    max_vel = vel.max()
            elif survey == 'ThrUMMS':
                print(file)

                transition = file.split('.')[2]

                # Specify transition
                if '12' in transition:
                    transitions = ['CO 1']
                elif '13' in transition:
                    transitions = ['13CO 1']
                else:
                    continue
                transition_indeces = [0]

                # Open file
                obs = fits.open(path + survey + '/' + file)

                # Get axes
                lon = np.linspace(obs[0].header['CRVAL1'] 
                        - obs[0].header['CDELT1'] 
                        * (obs[0].header['CRPIX1'] - 1),
                        obs[0].header['CRVAL1'] + obs[0].header['CDELT1'] 
                        * (obs[0].header['NAXIS1'] - obs[0].header['CRPIX1']),
                        num=obs[0].header['NAXIS1'])
                lat = np.linspace(obs[0].header['CRVAL2'] 
                        - obs[0].header['CDELT2'] 
                        * (obs[0].header['CRPIX2'] - 1),
                        obs[0].header['CRVAL2'] + obs[0].header['CDELT2'] 
                        * (obs[0].header['NAXIS2'] - obs[0].header['CRPIX2']),
                        num=obs[0].header['NAXIS2'])
                vel = np.linspace(obs[0].header['CRVAL3'] 
                        - obs[0].header['CDELT3'] 
                        * (obs[0].header['CRPIX3'] - 1),
                        obs[0].header['CRVAL3'] + obs[0].header['CDELT3'] 
                        * (obs[0].header['NAXIS3'] - obs[0].header['CRPIX3']),
                        num=obs[0].header['NAXIS3']) / 1e3

                # Copy header
                temp_header['NAXIS'] = 3
                temp_header['CTYPE1'] = target_header['CTYPE1']
                temp_header['CTYPE2'] = target_header['CTYPE2']
                temp_header['NAXIS3'] = target_header['NAXIS3']
                temp_header['CTYPE3'] = target_header['CTYPE3']
                temp_header['CRVAL3'] = target_header['CRVAL3']
                temp_header['CDELT3'] = target_header['CDELT3']
                temp_header['CRPIX3'] = target_header['CRPIX3']

                obs_data = np.swapaxes(obs[0].data, 0, 2)
                obs_data = np.swapaxes(obs_data, 0, 1)
                obs_data = spectres(target_vel, vel, 
                        np.nan_to_num(obs_data, nan=0), fill=0, verbose=False)
                obs_data = obs_data.reshape(-1, obs_data.shape[-1])
                if '12' in transition:
                    # from Barnes et al. (2015)
                    obs_error = np.full((obs_data.shape[0]), 1.3)
                elif '13' in transition:
                    # from Barnes et al. (2015)
                    obs_error = np.full((obs_data.shape[0]), 0.7)
                del obs

                # Grid
                lon_mesh, lat_mesh = np.meshgrid(lon, lat)
                if '12' in transition:
                    co1_gridder.grid(lon_mesh.flatten(), 
                            lat_mesh.flatten(), obs_data)
                    co1_gridder_err.grid(lon_mesh.flatten(), 
                            lat_mesh.flatten(), obs_error)
                elif '13' in transition:
                    ico1_gridder.grid(lon_mesh.flatten(), 
                            lat_mesh.flatten(), obs_data)
                    ico1_gridder_err.grid(lon_mesh.flatten(), 
                            lat_mesh.flatten(), obs_error)

                if vel.min() < min_vel:
                    min_vel = vel.min()
                if vel.max() > max_vel:
                    max_vel = vel.max()
            elif survey == 'THOR':
                if not ('.fits' in file):
                    continue
                print(file)

                # Specify transition
                transitions = ['HI']
                transition_indeces = [0]

                # Open file
                obs = fits.open(path + survey + '/' + file)

                # Get axes
                lon = np.linspace(obs[0].header['CRVAL1'] 
                        - obs[0].header['CDELT1'] 
                        * (obs[0].header['CRPIX1'] - 1),
                        obs[0].header['CRVAL1'] + obs[0].header['CDELT1'] 
                        * (obs[0].header['NAXIS1'] - obs[0].header['CRPIX1']),
                        num=obs[0].header['NAXIS1'])
                lat = np.linspace(obs[0].header['CRVAL2'] 
                        - obs[0].header['CDELT2'] 
                        * (obs[0].header['CRPIX2'] - 1),
                        obs[0].header['CRVAL2'] + obs[0].header['CDELT2'] 
                        * (obs[0].header['NAXIS2'] - obs[0].header['CRPIX2']),
                        num=obs[0].header['NAXIS2'])
                vel = np.linspace(obs[0].header['CRVAL3'] 
                        - obs[0].header['CDELT3'] 
                        * (obs[0].header['CRPIX3'] - 1),
                        obs[0].header['CRVAL3'] + obs[0].header['CDELT3'] 
                        * (obs[0].header['NAXIS3'] - obs[0].header['CRPIX3']),
                        num=obs[0].header['NAXIS3']) / 1e3

                # Copy header
                temp_header['NAXIS'] = 3
                temp_header['CTYPE1'] = target_header['CTYPE1']
                temp_header['CTYPE2'] = target_header['CTYPE2']
                temp_header['NAXIS3'] = target_header['NAXIS3']
                temp_header['CTYPE3'] = target_header['CTYPE3']
                temp_header['CRVAL3'] = target_header['CRVAL3']
                temp_header['CDELT3'] = target_header['CDELT3']
                temp_header['CRPIX3'] = target_header['CRPIX3']


                if 'THOR' in file:
                    obs_data = np.swapaxes(obs[0].data, 0, 2)
                    obs_data = (obs_data * 1.0492224297700237*0.21106**2
                            /obs[0].header['BMAJ']
                            /obs[0].header['BMAJ']) #Jy/beam->Tb
                elif '.pks.' in file:
                    obs_data = np.swapaxes(obs[0].data, 0, 2)
                elif 'CGPS' in file:
                    vel = vel[::-1]
                    obs_data = np.swapaxes(obs[0].data[0, ::-1, :, :], 0, 2)
                elif 'CAR' in file or 'SIN' in file:
                    obs_data = np.swapaxes(obs[0].data, 0, 2)
                else:
                    obs_data = np.swapaxes(obs[0].data[0], 0, 2)
                obs_data = np.swapaxes(obs_data, 0, 1)
                obs_data = spectres(target_vel, vel, np.nan_to_num(obs_data, 
                    nan=0), fill=0, verbose=False)
                obs_data = obs_data.reshape(-1, obs_data.shape[-1])
                # obs_error = determine_rms(obs, mission=survey, 
                #         file=file)[0].data.reshape(-1, 1)
                # i_nan = np.isnan(obs_error)
                # print(np.nanmean(obs_error), np.mean(obs_error[~i_nan]))
                # np.save('/home/yanitski/obs_error.npy', obs_error)
                # exit()
                del obs

                # Grid
                lon_mesh, lat_mesh = np.meshgrid(lon, lat)
                hi_gridder.grid(lon_mesh.flatten(), lat_mesh.flatten(), 
                        obs_data)
                hi_gridder_err.grid(lon_mesh.flatten(), lat_mesh.flatten(),
                                      np.full_like(lon_mesh.flatten(), 1.6))

                if vel.min() < min_vel:
                    min_vel = vel.min()
                if vel.max() > max_vel:
                    max_vel = vel.max()
            elif survey == 'SEDIGISM':
                print(file)

                # Specify transition
                transitions = ['13CO 2']
                transition_indeces = [0]

                # Open file
                obs = fits.open(path + survey + '/' + file)

                # Get axes
                lon = np.linspace(obs[0].header['CRVAL1'] 
                        - obs[0].header['CDELT1'] 
                        * (obs[0].header['CRPIX1'] - 1),
                        obs[0].header['CRVAL1'] + obs[0].header['CDELT1'] 
                        * (obs[0].header['NAXIS1'] - obs[0].header['CRPIX1']),
                        num=obs[0].header['NAXIS1'])
                lat = np.linspace(obs[0].header['CRVAL2'] 
                        - obs[0].header['CDELT2'] 
                        * (obs[0].header['CRPIX2'] - 1),
                        obs[0].header['CRVAL2'] + obs[0].header['CDELT2'] 
                        * (obs[0].header['NAXIS2'] - obs[0].header['CRPIX2']),
                        num=obs[0].header['NAXIS2'])
                vel = np.linspace(obs[0].header['CRVAL3'] 
                        - obs[0].header['CDELT3'] 
                        * (obs[0].header['CRPIX3'] - 1),
                        obs[0].header['CRVAL3'] + obs[0].header['CDELT3'] 
                        * (obs[0].header['NAXIS3'] - obs[0].header['CRPIX3']),
                        num=obs[0].header['NAXIS3']) / 1e3

                # Copy header
                temp_header['NAXIS'] = 3
                temp_header['CTYPE1'] = target_header['CTYPE1']
                temp_header['CTYPE2'] = target_header['CTYPE2']
                temp_header['NAXIS3'] = target_header['NAXIS3']
                temp_header['CTYPE3'] = target_header['CTYPE3']
                temp_header['CRVAL3'] = target_header['CRVAL3']
                temp_header['CDELT3'] = target_header['CDELT3']
                temp_header['CRPIX3'] = target_header['CRPIX3']


                obs_data = np.swapaxes(obs[0].data, 0, 2)
                obs_data = np.swapaxes(obs_data, 0, 1)
                obs_data = spectres(target_vel, vel, 
                        np.nan_to_num(obs_data, nan=0), fill=0, verbose=False)
                obs_data = obs_data.reshape(-1, obs_data.shape[-1])
                obs_rms_error = determine_rms(obs, mission=survey, 
                        file=file)[0].data.reshape(-1, 1)
                if cal_error:
                    idx_sig = obs_data > sig_th*obs_rms_error
                    obs_error = obs_rms_error
                    obs_error[idx_sig] = cal_error*obs_data[idx_sig] \
                            + obs_error[idx_sig]
                else:
                    obs_error = obs_rms_error
                i_nan = np.isnan(obs_error)
                # print(np.nanmean(obs_error), np.mean(obs_error[~i_nan]))
                # np.save('/home/yanitski/obs_error.npy', obs_error)
                # exit()
                del obs

                # Grid
                lon_mesh, lat_mesh = np.meshgrid(lon, lat)
                ico2_gridder.grid(lon_mesh.flatten(), lat_mesh.flatten(), 
                        obs_data)
                if cal_error:
                    ico2_gridder_err.grid(lon_mesh.flatten(), 
                            lat_mesh.flatten(), obs_error)
                else:
                    ico2_gridder_err.grid(lon_mesh.flatten()[~i_nan.flatten()], 
                            lat_mesh.flatten()[~i_nan.flatten()], 
                            obs_error[~i_nan])

                if vel.min() < min_vel:
                    min_vel = vel.min()
                if vel.max() > max_vel:
                    max_vel = vel.max()
            elif survey == 'Planck':
                print(file)

                # Specify observation type
                transitions = ['550um']

                # Open file and get data
                obs = fits.open(path + survey + '/' + file)

                if obs[1].header['ORDERING'] == 'NESTED':
                    nest = True
                else:
                    nest = False

                if 'commander' in file:
                    obs_tkj = obs[1].data['I_ML'] * 1e-6
                    obs_t = obs[1].data['TEMP_ML']
                    obs_beta = obs[1].data['BETA_ML']
                    obs_gamma = 6.646e-34/1.38e-23/obs_t
                    obs_data = obs_tkj 
                    # * (nu_planck/545e9)**(obs_beta+1) 
                    # * (np.exp(obs_gamma*545e9)-1)
                    # / (np.exp(obs_gamma*nu_planck)-1)
                    obs_error = obs[1].data['I_RMS'] * 1e-6
                elif 'GNILC' in file:
                    freq = int(file.split('F')[1].split('_')[0])
                    obs_data = obs[1].data['I'] * 32.56 / freq ** 2
                    obs_error = np.full_like(obs_data, obs_data.min())
                else:
                    print('File {} not available in survey {}'.format(
                        file, survey))
                    continue

                nside = obs[1].header['NSIDE']
                npix = hp.nside2npix(nside)

                # fix header
                if 'CDELT3' in temp_header.keys():
                    temp_header['NAXIS'] = 2
                    del temp_header['NAXIS3']
                    del temp_header['CTYPE3']
                    del temp_header['CDELT3']
                    del temp_header['CRVAL3']
                    del temp_header['CRPIX3']

                # Gridding parameters (split large files into multiple chuncks)
                chunck = 1000000
                n_chuncks = int(np.ceil(obs_data.size/chunck))
                lat_mesh, lon_mesh = np.degrees(hp.pix2ang(nside=nside, 
                    ipix=np.arange(npix), nest=nest))
                lat_mesh = 90 - lat_mesh

                # Grid
                dust_gridder = cygrid.WcsGrid(temp_header)
                dust_gridder.set_kernel(*target_kernel)
                for _ in range(n_chuncks):
                    dust_gridder.grid(lon_mesh[_*chunck:(_+1)*chunck],
                                      lat_mesh[_*chunck:(_+1)*chunck],
                                      obs_data[_*chunck:(_+1)*chunck])
                dust_gridder_err = cygrid.WcsGrid(temp_header)
                dust_gridder_err.set_kernel(*target_kernel)
                for _ in range(n_chuncks):
                    dust_gridder_err.grid(lon_mesh[_*chunck:(_+1)*chunck],
                                          lat_mesh[_*chunck:(_+1)*chunck],
                                          obs_error[_*chunck:(_+1)*chunck])

                print('save file')
                temp_header['TRANSL'] = transitions[0]
                temp_header['TRANSI'] = '0'
                grid_hdu = fits.PrimaryHDU(data=dust_gridder.get_datacube(), 
                        header=fits.Header(temp_header))
                grid_hdu_err = fits.PrimaryHDU(
                        data=dust_gridder_err.get_datacube(), 
                        header=fits.Header(twod_header))
                grid_hdu.writeto(path + survey + '/regridded/temp/planck_' 
                        + file.split('_')[2] + '_regridded.fits',
                        overwrite=True, output_verify='ignore')
                grid_hdu_err.writeto(path + survey + '/regridded/temp/planck_' 
                        + file.split('_')[2] + '_regridded_error.fits', 
                        overwrite=True, output_verify='ignore')
                del obs_data
                del obs_error
                del obs
            elif survey == 'COBE-FIRAS':
                print(file)

                # Specify transitions
                if 'HIGH' in file.split('.')[0]:
                    transitions = ['CO 6', 'C 2', 'H2O f1113', 'N+ 1', 
                            'H2O 2', 'C+ 1', 'O 1', 'Si', 'N+ 2', 'CH 2']   
                    #'CO 6', 'O 2'
                    transition_indeces = [0, 1, 2, 3, 4, 5, 6, 7, 8, 9] #0, 6
                elif 'HRES' in file.split('.')[0]:
                    transitions = ['CO 1', 'CO 2', 'CO 3', 'O2 13', 
                            'CO 4', 'C 1', 'H2O f557', 'CO 5']
                    transition_indeces = [0, 1, 2, 4, 5, 7]
                elif 'LOWF' in file.split('.')[0]:
                    transitions = ['CO 1', 'CO 2', 'CO 3', 
                            'CO 4', 'C 1', 'CO 5']
                    transition_indeces = [0, 1, 2, 4, 5, 7]

                # Open data and convert to brightness temperature
                obs = fits.open(path + survey + '/' + file)
                linfrq = np.array([obs[0].header[key] for key in \
                        obs[0].header.keys()
                                   if 'LINFRQ' in key])
                obs_data = (np.nan_to_num(obs[1].data['LINE_FLU'], nan=0)
                        * (2.9979**3) / (linfrq**3) / 2
                        / 1.38 * 10 ** 8) 
                #* np.pi**3*7**2/180**2  
                # corresponding beam size in sr
                obs_error = (np.nan_to_num(obs[1].data['LINE_FL2'], nan=0) 
                        * (2.9979**3) / (linfrq**3) / 2
                        / 1.38 * 10 ** 8) 
                #* np.pi**3*7**2/180**2  
                # corresponding beam size in sr

                # Get axes
                lon_mesh = obs[1].data['GAL_LON']
                lat_mesh = obs[1].data['GAL_LAT']

                # Fix header
                if 'CDELT3' in temp_header.keys():
                    # temp_header['NAXIS'] = 2
                    # del temp_header['NAXIS3']
                    del temp_header['CTYPE3']
                    del temp_header['CDELT3']
                    del temp_header['CRVAL3']
                    del temp_header['CRPIX3']
                temp_header['NAXIS'] = 3
                temp_header['NAXIS3'] = obs_data.shape[1]

                # Grid
                gridder = cygrid.WcsGrid(temp_header)
                gridder.set_kernel(*target_kernel)
                gridder.grid(lon_mesh, lat_mesh, obs_data)
                gridder_err = cygrid.WcsGrid(temp_header)
                gridder_err.set_kernel(*target_kernel)
                gridder_err.grid(lon_mesh, lat_mesh, obs_error)
                temp_header['TRANSL'] = ', '.join(transitions)
                temp_header['TRANSI'] = ', '.join('{}'.format(_) for 
                        _ in np.arange(len(transitions)))
                grid_hdu = fits.PrimaryHDU(data=gridder.get_datacube(), 
                        header=fits.Header(temp_header))
                grid_hdu_err = fits.PrimaryHDU(data=gridder_err.get_datacube(), 
                        header=fits.Header(temp_header))
                grid_hdu.writeto(path + survey + '/regridded/temp/' 
                        + file.replace('.FITS', '_regridded.fits'),
                        overwrite=True, output_verify='ignore')
                grid_hdu_err.writeto(path + survey + '/regridded/temp/' 
                        + file.replace('.FITS', '_regridded_error.fits'),
                        overwrite=True, output_verify='ignore')
            elif survey == 'COBE-DIRBE':
                if file == 'DIRBE_SKYMAP_INFO.FITS' or \
                        file == 'additional_files':
                    continue
                print(file)

                # Specify transitions
                transitions = ['240um']
                transition_indeces = [0]

                # Open data and convert to brightness temperature
                obs = fits.open(path + survey + '/' + file)
                pixcoord = fits.open(path + survey + 
                        '/additional_files/DIRBE_SKYMAP_INFO.FITS')
                linfrq = 2.9979e5 / 240  
                # convert the relevant DIRBE band to GHz
                obs_data = (np.nan_to_num([obs[1].data['Resid']], nan=0).T 
                        * (2.9979**2) / (linfrq**2) / 2 / 1.38 * 10)
                obs_error = (np.nan_to_num([obs[1].data['StdDev']], nan=0).T 
                        * (2.9979**2) / (linfrq**2) / 2 / 1.38 * 10) 

                # Get axes
                lon_mesh = pixcoord[1].data['GLON-CSC']
                lat_mesh = pixcoord[1].data['GLAT-CSC']

                # Fix header
                if 'CDELT3' in temp_header.keys():
                    # temp_header['NAXIS'] = 2
                    # del temp_header['NAXIS3']
                    temp_header['NAXIS'] = 2
                    del temp_header['NAXIS3']
                    del temp_header['CTYPE3']
                    del temp_header['CDELT3']
                    del temp_header['CRVAL3']
                    del temp_header['CRPIX3']

                # Grid
                print(obs_data.shape)
                gridder = cygrid.WcsGrid(temp_header)
                gridder.set_kernel(*target_kernel)
                gridder.grid(lon_mesh, lat_mesh, obs_data.flatten())
                gridder_err = cygrid.WcsGrid(temp_header)
                gridder_err.set_kernel(*target_kernel)
                gridder_err.grid(lon_mesh, lat_mesh, obs_error.flatten())
                temp_header['TRANSL'] = ', '.join(transitions)
                temp_header['TRANSI'] = ', '.join('{}'.format(_) for 
                        _ in np.arange(len(transitions)))
                grid_hdu = fits.PrimaryHDU(data=gridder.get_datacube(), 
                        header=fits.Header(temp_header))
                grid_hdu_err = fits.PrimaryHDU(data=gridder_err.get_datacube(), 
                        header=fits.Header(temp_header))
                grid_hdu.writeto(path + survey + '/regridded/temp/' 
                        + file.replace('.FITS', '_regridded.fits'),
                        overwrite=True, output_verify='ignore')
                grid_hdu_err.writeto(path + survey + '/regridded/temp/' 
                        + file.replace('.FITS', '_regridded_error.fits'),
                        overwrite=True, output_verify='ignore')
            else:
                # print('The specified survey {} is unavailable. 
                # Please add it to ´model selection´.'.format(survey))
                continue

        if CO1:
            temp_header['TRANSL'] = 'CO 1'
            temp_header['TRANSI'] = '0'
            grid_hdu = fits.PrimaryHDU(data=co1_gridder.get_datacube(), 
                    header=fits.Header(temp_header))
            grid_hdu_err = fits.PrimaryHDU(data=co1_gridder_err.get_datacube(), 
                    header=fits.Header(twod_header))
            grid_hdu.writeto(path + survey + '/regridded/temp/' +
                    'co1_test_regridded.fits', overwrite=True, 
                    output_verify='ignore')
            grid_hdu_err.writeto(path + survey + '/regridded/temp/' +
                    'co1_test_regridded_error.fits', overwrite=True, 
                    output_verify='ignore')
        if CO2:
            temp_header['TRANSL'] = 'CO 2'
            temp_header['TRANSI'] = '0'
            grid_hdu = fits.PrimaryHDU(data=co2_gridder.get_datacube(), 
                    header=fits.Header(temp_header))
            grid_hdu_err = fits.PrimaryHDU(data=co2_gridder_err.get_datacube(), 
                    header=fits.Header(twod_header))
            grid_hdu.writeto(path + survey + '/regridded/temp/' +
                    'co2_test_regridded.fits', overwrite=True, 
                    output_verify='ignore')
            grid_hdu_err.writeto(path + survey + '/regridded/temp/' +
                    'co2_test_regridded_error.fits', overwrite=True,
                    output_verify='ignore')
        if iCO1:
            temp_header['TRANSL'] = '13CO 1'
            temp_header['TRANSI'] = '0'
            grid_hdu = fits.PrimaryHDU(data=ico1_gridder.get_datacube(), 
                    header=fits.Header(temp_header))
            grid_hdu_err = fits.PrimaryHDU(
                    data=ico1_gridder_err.get_datacube(), 
                    header=fits.Header(twod_header))
            grid_hdu.writeto(path + survey + '/regridded/temp/' +
                    '13co1_test_regridded.fits', overwrite=True, 
                    output_verify='ignore')
            grid_hdu_err.writeto(path + survey + '/regridded/temp/' +
                    '13co1_test_regridded_error.fits', overwrite=True, 
                    output_verify='ignore')
        if iCO2:
            temp_header['TRANSL'] = '13CO 2'
            temp_header['TRANSI'] = '0'
            grid_hdu = fits.PrimaryHDU(data=ico2_gridder.get_datacube(), 
                    header=fits.Header(temp_header))
            grid_hdu_err = fits.PrimaryHDU(
                    data=ico2_gridder_err.get_datacube(), 
                    header=fits.Header(twod_header))
            grid_hdu.writeto(path + survey + '/regridded/temp/' +
                    '13co2_test_regridded.fits', overwrite=True, 
                    output_verify='ignore')
            if cal_error:
                grid_hdu_err.writeto(path + survey + '/regridded/temp/' + 
                        '13co2_test_regridded_complete_error.fits', 
                        overwrite=True, output_verify='ignore')
            else:
                grid_hdu_err.writeto(path + survey + '/regridded/temp/' + 
                        '13co2_test_regridded_error.fits', overwrite=True, 
                        output_verify='ignore')
        if hi:
            temp_header['TRANSL'] = 'HI'
            temp_header['TRANSI'] = '0'
            grid_hdu = fits.PrimaryHDU(data=hi_gridder.get_datacube(), 
                    header=fits.Header(temp_header))
            grid_hdu_err = fits.PrimaryHDU(data=hi_gridder_err.get_datacube(), 
                    header=fits.Header(twod_header))
            grid_hdu.writeto(path + survey + '/regridded/temp/' +
                    'hi_test_regridded.fits', overwrite=True, 
                    output_verify='ignore')
            grid_hdu_err.writeto(path + survey + '/regridded/temp/' +
                    'hi_test_regridded_error.fits', overwrite=True, 
                    output_verify='ignore')
        # if dust:
        #     temp_header['TRANSL'] = 'Dust'
        #     temp_header['TRANSI'] = '0'
        #     grid_hdu = fits.PrimaryHDU(data=dust_gridder.get_datacube(), 
        #             header=fits.Header(temp_header))
        #     grid_hdu_err = fits.PrimaryHDU(
        #             data=dust_gridder_err.get_datacube(), 
        #             header=fits.Header(twod_header))
        #     grid_hdu.writeto(path + survey 
        #             + '/regridded/temp/planck_dust_' + file 
        #             + '_regridded.fits',
        #             overwrite=True, output_verify='ignore')
        #     grid_hdu_err.writeto(path + survey + '/regridded/temp/' 
        #             + 'planck_dust_regridded_error.fits',
        #             overwrite=True, output_verify='ignore')

        CO1 = False
        CO2 = False
        iCO1 = False
        iCO2 = False
        dust = False
        cobe = False

    print('Regrid successfully completed.')

    return 0


def combine_regridded(path=None, regridded_path=None, target_header=None,
                      target_kernel=None, target_vel=None, output_file=None):
    '''
    Combine separate regridded files. Useful for large datasets.
    '''

    # if path==None or regridded_path==None or output_file==None\
    #         or isinstance(target_header, dict) or isinstance(target_vel, np.ndarray):
    #     print('use appropriate kwargs to combine regridded data')
    #     return

    files = os.listdir(path + regridded_path)

    if '/COBE/' in path or '/COBE-FIRAS/' in path:
        output_data = np.zeros((target_header['NAXIS2'], target_header['NAXIS1']))
    else:
        output_data = np.zeros((target_header['NAXIS3'], target_header['NAXIS2'], target_header['NAXIS1']))

    lon = np.linspace(target_header['CRVAL1'] - target_header['CDELT1'] * (target_header['CRPIX1'] - 1),
                      target_header['CRVAL1'] + target_header['CDELT1'] * (target_header['NAXIS1'] - target_header['CRPIX1']),
                      num=target_header['NAXIS1'])
    lat = np.linspace(target_header['CRVAL2'] - target_header['CDELT2'] * (target_header['CRPIX2'] - 1),
                      target_header['CRVAL2'] + target_header['CDELT2'] * (target_header['NAXIS2'] - target_header['CRPIX2']),
                      num=target_header['NAXIS2'])
    lon_mesh, lat_mesh = np.meshgrid(lon, lat)

    gridder = cygrid.WcsGrid(target_header)
    gridder.set_kernel(*target_kernel)

    for file in files:

        if file == output_file:
            continue

        if '/COBE/' in path or '/COBE-FIRAS/' in path:
            obs = fits.open(path + regridded_path + file)
            i = output_data == 0
            output_data[i] = output_data[i] + np.nan_to_num(obs[0].data)[i]
        else:
            obs = fits.open(path + regridded_path + file)
            # obs_data = np.reshape(obs[0].data, ())
            gridder.grid(lon_mesh.flatten(), lat_mesh.flatten(), np.nan_to_num(obs[0].data.reshape(-1, obs[0].shape[-1])))
            # i = output_data == 0
            # vel = np.linspace(obs[0].header['CRVAL3'] - obs[0].header['CDELT3'] * (obs[0].header['CRPIX3'] - 1),
            #                   obs[0].header['CRVAL3'] + obs[0].header['CDELT3'] * (
            #                               obs[0].header['NAXIS3'] - obs[0].header['CRPIX3']),
            #                   num=obs[0].header['NAXIS3'])
            # interp = interp1d(vel, np.nan_to_num(obs[0].data), axis=0, fill_value=0, bounds_error=False)
            # # output_data[i] = output_data[i] + interp(target_vel)[i]
            # output_data[i] = output_data[i] + np.nan_to_num(obs[0].data[i])

    output_data = gridder.get_datacube()

    # output = fits.PrimaryHDU(data=output_data, header=fits.Header(obs[0].header))
    output = fits.PrimaryHDU(data=output_data, header=fits.Header(obs[0].header))
    output.writeto(path + output_file, overwrite=True, output_verify='ignore')

    print('Saved as ' + output_file)

    return


def view_observation(path='/mnt/hpc_backup/yanitski/projects/pdr/observational_data/MilkyWay/',
                     mission='COGAL', transition='', regridded_path='/regridded/', filename='CO1_obs_regridded.fits',
                     plot='integrated', integrate_b=[], i_lat=None, list_observations=False, xlabel='', ylabel='',
                     clabel='', title='', fontsize=16, scale=1, logval=False, vmin=0, vmax=None, cmap='viridis',
                     xlim=None, ylim=None, save=False):
    '''
    This function will plot the specified (regridded) observation. By default it prints the integrated emission.

    :param path: Directory containing the observational data.
    :param mission: Name of the survey you wish to view.
    :param regridded_path: folder name for the regridded data. This should usually be '/regridded/', but some surveys
           might require '/regridded/temp/'
    :param filename: Directory containing the observational data.
    :param plot: The type of plot to create with this function.
           - 'integrated'
           - 'pv'
           By default it 'integrated'.
    :param title: The title of the produced plot.
    :param list_observations: Set this flag to True to print the available observation files.
    :return:
    '''

    if list_observations:
        print(path + mission + regridded_path)
        files = os.listdir(path + mission + regridded_path)
        print('   - '.join([file for file in files if file != 'temp']))
        return

    obs = fits.open(path + mission + regridded_path + filename)
    cygrid_data = obs[0].data
    cygrid_data[cygrid_data<0] = 0
    header = obs[0].header

    lon = np.linspace(header['CRVAL1'] - header['CDELT1'] * (header['CRPIX1'] - 1),
                      header['CRVAL1'] + header['CDELT1'] * (header['NAXIS1'] - header['CRPIX1']),
                      num=header['NAXIS1'])
    lat = np.linspace(header['CRVAL2'] - header['CDELT2'] * (header['CRPIX2'] - 1),
                      header['CRVAL2'] + header['CDELT2'] * (header['NAXIS2'] - header['CRPIX2']),
                      num=header['NAXIS2'])

    pprint(header)

    twod_header = copy(header)
    if (mission != 'COBE-FIRAS') and (mission != 'COBE-DIRBE') and (mission != 'Planck') and not ('error' in filename):
        print('spectroscopic data')
        if transition != header['TRANSL']:
            print('Line transition not in specified file. Please supply the correct line/file.')
            print('  - {}'.format(header['TRANSL']))
            return
        twod_header['NAXIS'] = 2
        del twod_header['NAXIS3']
        del twod_header['CTYPE3']
        del twod_header['CRPIX3']
        del twod_header['CRVAL3']
        del twod_header['CDELT3']
        vel = np.linspace(header['CRVAL3'] - header['CDELT3'] * (header['CRPIX3'] - 1),
                          header['CRVAL3'] + header['CDELT3'] * (header['NAXIS3'] - header['CRPIX3']),
                          num=header['NAXIS3'])
        cygrid_integrated_data = np.trapz(cygrid_data, vel, axis=0)
        cygrid_integrated_data[cygrid_integrated_data == 0] = np.nan
        cygrid_data[cygrid_data == 0] = np.nan
    elif 'error' in filename:
        print('error')
        if transition != header['TRANSL']:
            print('Line transition not in specified file. Please supply the correct line/file.')
            print('  - {}'.format(header['TRANSL']))
            return
        twod_header['NAXIS'] = 2
        del twod_header['CTYPE3']
        del twod_header['CRPIX3']
        del twod_header['CRVAL3']
        del twod_header['CDELT3']
        cygrid_integrated_data = copy(cygrid_data)
        cygrid_integrated_data[cygrid_integrated_data == 0] = np.nan
        pprint(twod_header)
    else:
        print('COBE or Planck')
        if twod_header['NAXIS'] == 3:
            twod_header['NAXIS'] = 2
            del twod_header['NAXIS3']
            # del twod_header['OBSERR']
            del twod_header['TRANSL']
            # del twod_header['TRANSI']
        transitions = np.asarray(header['TRANSL'].split(', '))
        # i_transitions = np.asarray(header['TRANSI'].split(', '))
        i_line = transitions == transition
        # print('\n', transitions, '\n', transition, i_line, int(lat.size/2), '\n')
        if not np.any(i_line):
            print('Line transition {} not in specified file. Please supply the correct line/file.'.format(transition))
            for line in transitions:
                print('  - {}'.format(line))
            return
        if mission == 'COBE-FIRAS':
            cygrid_integrated_data = cygrid_data[i_line, :, :][0]
        else:
            cygrid_integrated_data = cygrid_data
        pprint(twod_header)

    twod_wcs = WCS(twod_header)

    # for i_vel in [0, 720, 1440]:
    #     fig = plt.figure(figsize=(20, 20))
    #     ax = fig.add_subplot(111, projection=twod_wcs, frame_class=EllipticalFrame, slices=('x', 'y', 250))
    #     #     ax = fig.add_subplot(111)
    #     ax.imshow(cygrid_data[i_vel, :, :], vmin=0, vmax=10)
    #     plt.show()

    fig = plt.figure(figsize=(20, 10))
    if plot == 'integrated':
        print(np.nanmin(cygrid_integrated_data), np.nanmax(cygrid_integrated_data))
        if vmax == None:
            vmax = scale * cygrid_integrated_data.max()
        ax = fig.add_subplot(111, projection=twod_wcs, frame_class=EllipticalFrame)
        if logval:
            cm = ax.imshow(np.log10(cygrid_integrated_data), vmin=vmin, vmax=vmax, cmap=cmap)
        else:
            cm = ax.imshow(cygrid_integrated_data, vmin=vmin, vmax=vmax, cmap=cmap)
        # ax.set_aspect(np.abs(twod_wcs.wcs.cdelt[1] / twod_wcs.wcs.cdelt[0]))
        cb = fig.colorbar(cm, ax=ax, fraction=0.03, aspect=20)
        if clabel:
            cb.ax.set_ylabel(clabel, fontsize=fontsize)
        else:
            cb.ax.set_ylabel(r'Integrated Intensity ($K \frac{km}{s}$)', fontsize=fontsize)
    elif plot == 'pv' or plot == 'PV':
        print('Number of finite-valued sightlines:', cygrid_data.size-np.where(np.isnan(cygrid_data))[0].size)
        if type(i_lat) != int:
            i_lat = int(lat.size/2)
        if vmax == None:
            vmax = scale * cygrid_data.max()
        ax = fig.add_subplot(111)
        if integrate_b:
            print('integrate latitude')
            if 'error' in filename:
                print('min', np.nanmin(np.log10(np.nansum(cygrid_data[(integrate_b[0]-1):integrate_b[1], :], axis=0))),
                      '\nmax', np.nanmax(np.log10(np.nansum(cygrid_data[(integrate_b[0]-1):integrate_b[1], :], axis=0))),
                      '\nNaN?', (np.isnan(np.log10(np.nansum(cygrid_data[(integrate_b[0]-1):integrate_b[1], :], axis=0)))).any(),
                      np.isnan(np.log10(cygrid_data)).all())
                cm = ax.pcolormesh(lon, [1], np.nansum([cygrid_data[(integrate_b[0]-1):integrate_b[1], :]], axis=1),
                               vmin=vmin, vmax=vmax, shading='auto', cmap=cmap)
            else:
                print('min', np.nanmin(np.log10(np.nansum(cygrid_data[:, (integrate_b[0]-1):integrate_b[1], :], axis=1))),
                      '\nmax', np.nanmax(np.log10(np.nansum(cygrid_data[:, (integrate_b[0]-1):integrate_b[1], :], axis=1))),
                      '\nNaN?', (np.isnan(np.log10(np.nansum(cygrid_data[:, (integrate_b[0]-1):integrate_b[1], :], axis=1)))).any(),
                      np.isnan(np.log10(cygrid_data)).all())
                cm = ax.pcolormesh(lon, vel, np.log10(np.nansum(cygrid_data[:, (integrate_b[0]-1):integrate_b[1], :], axis=1)),
                               vmin=vmin, vmax=vmax, shading='auto', cmap=cmap)
        else:
            print('at latitude', lat[i_lat], 'deg')
            if 'error' in filename:
                print(np.nanmin(cygrid_data[int(lat.size/2), :]),
                      np.nanmax(cygrid_data[int(lat.size/2), :]),
                      (~np.isnan(cygrid_data[int(lat.size/2), :])).any())
                cm = ax.pcolormesh(lon, [1], cygrid_data[i_lat, :],
                               vmin=0, vmax=vmax, shading='auto', cmap=cmap)
            else:
                print(np.nanmin(cygrid_data[:, int(lat.size/2), :]),
                      np.nanmax(cygrid_data[:, int(lat.size/2), :]),
                      (~np.isnan(cygrid_data[:, int(lat.size/2), :])).any())
                cm = ax.pcolormesh(lon, vel, cygrid_data[:, i_lat, :],
                               vmin=0, vmax=vmax, shading='auto', cmap=cmap)
        if xlabel:
            ax.set_xlabel(xlabel, fontsize=fontsize)
        else:
            ax.set_xlabel(r'$\lambda \ \left( ^\circ \right)$', fontsize=fontsize)
        if ylabel:
            ax.set_ylabel(ylabel, fontsize=fontsize)
        else:
            ax.set_ylabel(r'$V_{LSR} \ \left( \frac{km}{s} \right)$', fontsize=fontsize)
        cb = fig.colorbar(cm, ax=ax, fraction=0.1, aspect=20)
        if clabel:
            cb.ax.set_ylabel(clabel, fontsize=fontsize)
        else:
            cb.ax.set_ylabel(r'Intensity ($K$)', fontsize=fontsize)
    elif plot == 'spectrum':
        if type(i_lat) != int:
            i_lat = int(lat.size/2)
        ax = fig.add_subplot(111)
        ax.step(lon, cygrid_integrated_data[i_lat, :].sum(0))
        ax.text(0.5, 0.01, '{} degrees'.format(lat[i_lat]), ha='center', va='bottom', transform=ax.transAxes,
                fontsize=fontsize)
        if xlabel:
            ax.set_xlabel(xlabel, fontsize=fontsize)
        else:
            ax.set_xlabel(r'$\lambda \ \left( ^\circ \right)$', fontsize=fontsize)
        if ylabel:
            ax.set_ylabel(ylabel, fontsize=fontsize)
        else:
            ax.set_ylabel(r'$T_{int} \ \left( K \frac{km}{s} \right)$', fontsize=fontsize)
    else:
        print('Please select a valid plotting method.')
        return

    if title:
        plt.title(title, fontsize=fontsize)
    else:
        plt.title(mission + ' ' + transition + ' ' + plot + ' plot', fontsize=fontsize)

    if xlim:
        ax.set_xlim(xlim)
    if ylim:
        ax.set_ylim(ylim)

    if save:
        plt.savefig(path + mission + regridded_path + transition.replace(' ', '') + '/' + filename.replace(".fits", "") + '_' + plot + '.png')
    else:
        plt.show()

    return


def error_correction(data, conf=''):
    '''
    Calculate configuration error
    '''
    if len(data.shape) == 0:
        return 0
    if data.shape[0] == 1 and data.ndim == 3:
        data_copy = np.nan_to_num(data[0], nan=0)
    else:
        data_copy = np.nan_to_num(data, nan=0)
    correction = np.zeros_like(data)

    if conf == 'axisymmetric':

        if data_copy.ndim == 1:
            for i in range(int(np.ceil(data_copy.shape[0]/2))):
                idx = np.ix_([i, -1-i])
                avg = np.mean(data_copy[idx])
                err = np.std(data_copy[idx]-avg, ddof=1)
                correction[idx] = err / np.sqrt(err.size)
        elif data_copy.ndim == 2:
            for i in range(int(np.ceil(data_copy.shape[0]/2))):
                for j in range(int(np.ceil(data_copy.shape[1]/2))):
                    idx = np.ix_([i, -1-i], [j, -1-j])
                    if np.isnan(data_copy[idx]).all():
                        continue
                    avg = np.nanmean(data_copy[idx])
                    err = np.std(data_copy[idx]-avg, ddof=1)
                    correction[idx] = err/np.sqrt(err.size)
        elif data_copy.ndim == 3:
            for i in range(int(np.ceil(data_copy.shape[0]/2))):
                for j in range(int(np.ceil(data_copy.shape[1]/2))):
                    for k in range(int(data_copy.shape[2])):
                        idx = np.ix_([i, -1-i], [j, -1-j], [k, k, -1-k, -1-k])
                        idx = (np.array([i, i, -1-i, -1-i]), np.array([j, -1-j, j, -1-j]), np.array([k, k, -1-k, -1-k]))
                        if np.isnan(data_copy[idx]).all():
                            continue
                        avg = np.nanmean(data_copy[idx])
                        err = np.std(data_copy[idx]-avg, ddof=1)
                        correction[idx] = err/np.sqrt(err.size)

    if data.shape != data_copy.shape:
        correction = np.asarray([correction])

    return correction



def model_selection_old(path='/mnt/hpc_backup/yanitski/projects/pdr/KT3_history/MilkyWay', missions=None, lat=None,
                    model_dir='', model_param=[[]], comp_type='pv', log_comp=True, cmap='gist_ncar',
                    PLOT=False, PRINT=False, debug=False):
    '''
    This function will compare the Milky Way models to the position-velocity observations depending on the
      dataset. It will utilise the sightlines of the dataset.

      :param path: directory containing the model folders.
      :param model_dir: directory format for each model in the grid.
      :param model_param: a list of lists containing the parameter space of the models.
      :param resolution: the voxel size (in pc) of the models used in the comparison; used in model folder name.
      :param missions: as this is setup for the Milky Way, the mission selection is
             - COBE-FIRAS
             - COBE-DIRBE
             - COGAL
             - Mopra
             - ThrUMMS
             - SEDIGISM
             - Planck
             - THOR (partial)
             - Hi-GAL (abandoned)
      :param cmap: color map for the plots.
      :param PLOT: flag to save plots.
      :param PRINT: flag to enable verbose output.
      :param debug: flag to run in debug mode.

      currently, this will look solely at the galactic plane.
    '''

    # Check that the missions are specified properly.
    if missions == '' or missions == None or missions == []:
        missions = os.listdir(path.replace('KT3_history', 'observational_data'))
    elif isinstance(missions, str):
        missions = [missions]
    elif isinstance(missions, list):
        # Use both COBE instruments if simply 'COBE' is specified
        if 'COBE' in missions:
            if 'COBE-FIRAS' in missions: missions.remove('COBE-FIRAS')
            if 'COBE-DIRBE' in missions: missions.remove('COBE-DIRBE')
            missions.append('COBE-FIRAS')
            missions.append('COBE-DIRBE')
    else:
        print('Please specify a list of missions to compare the models.')
        return

    if model_dir == '' or model_param == [[]]:
        print('Please specify both model_dir and model_param.')
        return

    model_params = np.meshgrid(*model_param)
    # model_params = zip(np.transpose([model_params[n].flatten() for n in range(len(model_param))]))
    model_params = np.transpose([model_params[n].flatten() for n in range(len(model_param))])

    if log_comp:
        comp = comp_type + '_logT'
    else:
        comp = comp_type

    obs_data = []
    lon = []
    vel = []

    for survey in missions:

        print('\n\n  {}'.format(survey))
        print('  ' + '-'*len(survey))

        if survey in ('HiGAL', 'GRS', 'THOR', 'CGPS', 'SGPS', 'VGPS'):
            continue

        if path[-1] != '/': path += '/'
        directory = path.replace('KT3_history', 'observational_data') + survey + '/regridded/temp/'
        # if mission == 'COBE-FIRAS':
        #     directory += 'temp/'
        files = os.listdir(directory)

        if len(files) == 0:
            print('no data is available.')
            exit()

        # Loop through the different models
        # for i_obs in range(len(obs_data)):
        for file in files:

            if not ('.fits' in file or '.idl' in file) or 'error' in file:
                continue

            print(f'\n  {file}')

            if '.fits' in file:
                obs = fits.open(directory + file)
                obs_error = fits.open(directory + file.replace('.fits', '_error.fits'))[0].data
            elif '.sav' in file or '.idl' in file:
                obs = readsav(directory + file)
            else:
                continue

            if survey == 'COBE-DIRBE' or survey == 'COBE-FIRAS':
                if file == 'craig.idl':

                    # This file requires a comparison to the galactic plane
                    # lat = 0

                    # these values are hard-coded since the data is not in the file
                    linfrq = np.array([115.3, 230.5, 345.8, 424.8, 461.0, 492.2, 556.9, 576.3, 691.5, 808.1,
                                       1113, 1460, 2226, 1901, 2060, 2311, 2459, 2589, 921.8])
                    transitions = np.array(['CO 1', 'CO 2', 'CO 3', 'CO 4', 'C 3', 'CO 5',
                                            'CO 6', 'CO 7 + C 1', 'C+ 1', 'O 1', 'CO 8'])
                    transition_indeces = np.array([0, 1, 2, 4, 5, 7, 8, 9, 13, 14, 18])

                    obs_data = obs['amplitude'] / (2.9979**3) * (linfrq**3) * 2 * 1.38 / 10 ** 8
                    obs_error = obs['sigma'] / (2.9979**3) * (linfrq**3) * 2 * 1.38 / 10 ** 8
                    lon_obs = obs['long']
                    lat_obs = np.array([0])
                    i_lat_obs_init = np.ones(1, dtype=bool)

                else:

                    obs_data = obs[0].data
                    if debug:
                        pprint(obs[0].header)
                    # obs_error = obs[0].header['OBSERR'].split()
                    transitions = obs[0].header['TRANSL'].split(', ')
                    transition_indeces = obs[0].header['TRANSI'].split(', ')
                    lon_obs = np.linspace(obs[0].header['CRVAL1'] - obs[0].header['CDELT1'] * (obs[0].header['CRPIX1'] - 1),
                                          obs[0].header['CRVAL1'] + obs[0].header['CDELT1'] * (
                                                      obs[0].header['NAXIS1'] - obs[0].header['CRPIX1']),
                                          num=obs[0].header['NAXIS1'])
                    lat_obs = np.linspace(obs[0].header['CRVAL2'] - obs[0].header['CDELT2'] * (obs[0].header['CRPIX2'] - 1),
                                          obs[0].header['CRVAL2'] + obs[0].header['CDELT2'] * (
                                                      obs[0].header['NAXIS2'] - obs[0].header['CRPIX2']),
                                          num=obs[0].header['NAXIS2'])
                    i_lat_obs_init = obs_data.any(0).any(1)
                vel = None

            elif survey == 'Planck':
                obs_data = obs[0].data
                if debug:
                    pprint(obs[0].header)
                # obs_error = obs[0].header['OBSERR'].split()
                transitions = obs[0].header['TRANSL'].split(', ')
                transition_indeces = obs[0].header['TRANSI'].split(', ')
                lon_obs = np.linspace(obs[0].header['CRVAL1'] - obs[0].header['CDELT1'] * (obs[0].header['CRPIX1'] - 1),
                                      obs[0].header['CRVAL1'] + obs[0].header['CDELT1'] * (
                                                  obs[0].header['NAXIS1'] - obs[0].header['CRPIX1']),
                                      num=obs[0].header['NAXIS1'])
                lat_obs = np.linspace(obs[0].header['CRVAL2'] - obs[0].header['CDELT2'] * (obs[0].header['CRPIX2'] - 1),
                                      obs[0].header['CRVAL2'] + obs[0].header['CDELT2'] * (
                                                  obs[0].header['NAXIS2'] - obs[0].header['CRPIX2']),
                                      num=obs[0].header['NAXIS2'])
                i_lat_obs_init = obs_data.any(1)
                vel = None
            else:
                obs_data = obs[0].data
                if debug:
                    pprint(obs[0].header)
                # obs_error = obs[0].header['OBSERR']
                transitions = obs[0].header['TRANSL'].split(', ')
                transition_indeces = obs[0].header['TRANSI'].split(', ')
                lon_obs = np.linspace(obs[0].header['CRVAL1'] - obs[0].header['CDELT1'] * (obs[0].header['CRPIX1'] - 1),
                                      obs[0].header['CRVAL1'] + obs[0].header['CDELT1'] * (
                                                  obs[0].header['NAXIS1'] - obs[0].header['CRPIX1']),
                                      num=obs[0].header['NAXIS1'])
                lat_obs = np.linspace(obs[0].header['CRVAL2'] - obs[0].header['CDELT2'] * (obs[0].header['CRPIX2'] - 1),
                                      obs[0].header['CRVAL2'] + obs[0].header['CDELT2'] * (
                                                  obs[0].header['NAXIS2'] - obs[0].header['CRPIX2']),
                                      num=obs[0].header['NAXIS2'])
                vel_obs = np.linspace(obs[0].header['CRVAL3'] - obs[0].header['CDELT3'] * (obs[0].header['CRPIX3'] - 1),
                                      obs[0].header['CRVAL3'] + obs[0].header['CDELT3'] * (
                                                  obs[0].header['NAXIS3'] - obs[0].header['CRPIX3']),
                                      num=obs[0].header['NAXIS3'])
                i_lat_obs_init = (obs_data.any(0)).any(1)

            for i, transition in enumerate(transitions):

                print('\n  fitting', transition)

                chi2_grid = []
                loglikelihood_grid = []
                params = []

                if (survey == 'COBE-DIRBE') or (survey == 'COBE-FIRAS'):
                    if '.idl' in file:
                        obs_error_conf = error_correction(obs_data[:, int(transition_indeces[i])],
                                                          conf='axisymmetric')
                        vmin = obs_data[:, int(transition_indeces[i])].min()
                        vmax = obs_data[:, int(transition_indeces[i])].max()
                    else:
                        obs_error_conf = error_correction(obs_data[int(transition_indeces[i]), :, :],
                                                          conf='axisymmetric')
                        vmin = obs_data[int(transition_indeces[i]), :, :].min()
                        vmax = obs_data[int(transition_indeces[i]), :, :].max()
                    # print(transition)
                    # print(vmin, vmax)
                else:
                    obs_error_conf = error_correction(obs_data, conf='axisymmetric')
                    vmin = np.nanmin(obs_data)
                    vmax = np.nanmax(obs_data)

                # print('Average error in observation: {}'.format(np.nanmean(obs_error)))
                # print('Correction due to axisymmetry: {}'.format(np.nanmean(obs_error_conf)))

                for param in model_params:

                    if PRINT:
                        print('  ' + model_dir.format(*param) + '    \r', end='')

                    # Check existance of model
                    dir_model = model_dir.format(*param) + 'synthetic_intensity.fits'
                        # '/r{}_cm{}-{}_d{}_uv{}/channel_intensity.fits'.format(resolution, fcl,
                        #                                                               ficl, fden, fuv)
                    if not os.path.isfile(path + dir_model):
                        print('  missing')
                        continue

                    # Open the model
                    params.append(param)
                    model = fits.open(path + dir_model)
                    # pprint(model[1].header)

                    # Locate the species transition (only one dust value is considered; constant background)
                    if survey in ('Planck', 'COBE-DIRBE'):
                        i_spec = np.array([0])
                    else:
                        i_spec = np.where(np.asarray(model[1].header['SPECIES'].split(', ')) == transition)[0]
                    if i_spec.size == 0:
                        print('  {} not in model.            '.format(transition))
                        break
                    # else:
                    #     i_spec = i_spec[0]

                    if os.path.isdir(path + 'fit_results/{}/'.format(survey)
                                     + '{}/{}/'.format(file.replace('.fits', ''), transition)) == False:
                        os.makedirs(path + 'fit_results/{}/'.format(survey)
                                    + '{}/{}/'.format(file.replace('.fits', ''), transition))

                    # Create arrays for the longitude and velocity axes
                    lat_model = np.linspace(
                        model[1].header['CRVAL3']
                        - model[1].header['CDELT3'] * (model[1].header['CRPIX3'] - 0.5),
                        model[1].header['CRVAL3']
                        + model[1].header['CDELT3'] * (model[1].header['NAXIS3'] - model[1].header['CRPIX3'] - 0.5),
                        num=model[1].header['NAXIS3']) * 180 / np.pi
                    lon_model = np.linspace(
                        model[1].header['CRVAL2']
                        - model[1].header['CDELT2'] * (model[1].header['CRPIX2'] - 0.5),
                        model[1].header['CRVAL2']
                        + model[1].header['CDELT2'] * (model[1].header['NAXIS2']-model[1].header['CRPIX2'] - 0.5),
                        num=model[1].header['NAXIS2']) * 180 / np.pi
                    # if not (mission == 'COGAL'):
                    #     lon_model[lon_model<0] += 360
                    vel_model = np.linspace(
                        model[1].header['CRVAL4']
                        - model[1].header['CDELT4'] * (model[1].header['CRPIX4'] - 0.5),
                        model[1].header['CRVAL4']
                        + model[1].header['CDELT4'] * (model[1].header['CRPIX4'] - 0.5),
                        num=model[1].header['NAXIS4'])
                    # breakpoint()

                    # Identify the common latitudes between the observations and the model
                    if ((isinstance(lat, int) or isinstance(lat, float)) and comp_type == 'pv' and
                        not (survey == 'COBE-FIRAS' or survey == 'COBE-DIRBE' or survey == 'Planck')):
                        lat_min = lat
                        lat_max = lat
                    else:
                        lat_min = lat_obs[i_lat_obs_init].min()
                        lat_max = lat_obs[i_lat_obs_init].max()

                    if survey == 'COBE-FIRAS' or survey == 'COBE-DIRBE' or survey == 'Planck':
                        i_lon_model = np.linspace(0, lon_model.size-1, num=lon_obs.size, dtype=int)
                        i_lat_model = (lat_model >= lat_min) \
                                      & (lat_model <= lat_max)
                        i_lat_obs = (lat_obs >= lat_model[i_lat_model].min()) \
                                    & (lat_obs <= lat_model[i_lat_model].max())
                    else:
                        i_lat_model = np.where((lat_model >= lat_min)
                                               & (lat_model <= lat_max))[0]
                        i_lat_obs = np.where((lat_obs >= lat_model[i_lat_model].min())
                                             & (lat_obs <= lat_model[i_lat_model].max()))[0]

                    # Interpolate at the observational axes
                    if (survey == 'COBE-DIRBE') or (survey == 'COBE-FIRAS'):
                        idx_t = np.ix_(np.arange(vel_model.size), i_lat_model, i_lon_model, i_spec)
                        # idx_d = np.ix_(np.arange(vel_model.size), i_lat_model, i_lon_model, [0])
                        idx_d = np.ix_(i_lat_model, i_lon_model, [0])
                        # print(i_lat_model.size, i_lon_model.size)
                        # print(vel_model.shape,
                        #       model[1].data[idx_t].shape,
                        #       model[2].data[idx_d].shape,
                        #       i_lat_obs.size)
                        model_data = model[1].data[idx_t] \
                                     - model[2].data[idx_d]
                        model_interp = trapz(model_data, vel_model, axis=model_data.shape.index(vel_model.size))
                        if '.idl' in file:
                            idx_obs = np.ix_(np.arange(obs_data.shape[0]), [int(transition_indeces[i])])
                            obs_data_final = obs_data[:, int(transition_indeces[i])].reshape(-1)#[:, int(transition_indeces[i])]
                            # print(obs_data_final.shape)
                            # print(obs_error[:, int(transition_indeces[i])].shape, obs_error_conf.shape)
                            # print(obs_error[idx_obs].shape, obs_error_conf.shape)
                            obs_error_final = np.sqrt(obs_error[:, int(transition_indeces[i])]**2
                                                      +obs_error_conf**2).reshape(-1)#[:, int(transition_indeces[i])]
                            # print(obs_error_final.shape)
                            model_interp = model_interp.reshape(-1)
                        else:
                            idx_obs = np.ix_([int(transition_indeces[i])], i_lat_obs, np.arange(obs_data.shape[2]))
                            obs_data_final = obs_data[idx_obs][0, :, :]#[int(transition_indeces[i]), i_lat_obs, :]
                            obs_error_final = np.sqrt(obs_error**2+obs_error_conf**2)[idx_obs][0, :, :]#[int(transition_indeces[i]), i_lat_obs, :]
                            model_interp = model_interp[:, :, 0]
                        # print(model_interp.shape, obs_data_final.shape)

                    elif survey == 'Planck':
                        # idx_d = np.ix_([0], i_lat_model, np.arange(lon_model.size), [0])
                        idx_d = np.ix_(i_lat_model, np.arange(lon_model.size), [0])
                        # model_interp = model[2].data[idx_d][0, :, :, 0]#[0, i_lat_model, :, i_spec]
                        model_interp = model[2].data[idx_d][:, :, 0]#[0, i_lat_model, :, i_spec]
                        idx_obs = np.ix_(i_lat_obs, np.arange(obs_data.shape[1]))
                        obs_data_final = obs_data[idx_obs]#[i_lat_obs, :]
                        obs_error_final = np.sqrt(obs_error**2+obs_error_conf**2)[idx_obs]#[i_lat_obs, :]
                        # print(model_interp.shape, obs_data_final.shape)

                    else:
                        # print(vel_model.shape,
                        #       model[1].data[:, i_lat_model, :, i_spec].shape,
                        #       model[2].data[:, i_lat_model, :, 0].shape,
                        #       obs_data[:, i_lat_obs, :].shape,
                        #       np.swapaxes(obs_data[:, i_lat_obs, :], 0, 1).shape)
                        # input()
                        idx_t = np.ix_(np.arange(vel_model.size), i_lat_model, np.arange(lon_model.size), i_spec)
                        # idx_d = np.ix_(np.arange(vel_model.size), i_lat_model, np.arange(lon_model.size), [0])
                        idx_d = np.ix_(i_lat_model, np.arange(lon_model.size), [0])
                        model_data = (model[1].data[idx_t]
                                      - model[2].data[idx_d])
                        idx_obs = np.ix_(np.arange(vel_obs.size), i_lat_obs, np.arange(lon_obs.size))
                        idx_map = np.ix_(i_lat_obs, np.arange(lon_obs.size))
                        obs_data_temp = obs_data[idx_obs]#np.swapaxes(obs_data[idx_obs], 0, 1)
                        if comp_type == 'integrated':
                            # This will give an incorrect result if either the latitude or longitude axes are the same
                            # size as the velocity axis
                            model_interp = np.trapz(model_data, vel_model, axis=model_data.shape.index(vel_model.size))
                            obs_data_final = np.trapz(obs_data_temp, vel_obs, axis=obs_data_temp.shape.index(vel_obs.size))
                            obs_error_final = len(vel_obs)*np.sqrt(obs_error**2+obs_error_conf**2)[idx_obs]
                        elif comp_type == 'pv':
                            model_interp = model_data[:, :, :, 0]#np.swapaxes(model_data, 0, 1)
                            obs_data_final = copy(obs_data_temp)
                            obs_error_final = np.sqrt(obs_error**2+obs_error_conf**2)[idx_obs]
                            # obs_error_final = np.swapaxes([obs_error[i_lat_obs, :]]*obs_data_temp.shape[1], 0, 1)
                        else:
                            print('ERROR >> Comparison type {} not valid; '.format(comp_type) +
                                  'please choose "pv" or "integrated"')
                            return 1

                        model_interp[model_interp == 0] = np.nan

                        if debug:
                            # print(vel)
                            # print(lon)
                            # print(vel_model)
                            # print(lon_model)
                            print(vmax, vmin)
                            print(obs_data.max(), obs_data.min())
                            print(vel.shape, lon.shape, obs_data.shape, model_interp.shape)

                            # lon_model_grid, vel_model_grid = np.meshgrid(lon_model, vel_model)
                            # plt.pcolormesh(lon_model_grid, vel_model_grid, model_data,
                            #                norm=colors.SymLogNorm(base=10, linthresh=0.1, vmin=0, vmax=vmax),
                            #                shading='auto', cmap=cmap)
                            # plt.show(block=True)
                            cb = plt.pcolormesh(lon, vel, model_interp,
                                                norm=colors.SymLogNorm(base=10, linthresh=0.1, vmin=0, vmax=vmax),
                                                shading='auto', cmap='viridis')
                            plt.colorbar(cb)
                            plt.show(block=True)
                            cb = plt.pcolormesh(lon, vel, obs_data,
                                                norm=colors.SymLogNorm(base=10, linthresh=0.1, vmin=vmin, vmax=vmax),
                                                shading='auto', cmap='viridis', alpha=0.25)
                            plt.colorbar(cb)
                            plt.show(block=True)
                            # i_obs = np.isnan(model_interp)
                            # obs_data[i_obs] = np.nan
                            kernel = Gaussian2DKernel(x_stddev=2)
                            obs_data_smoothed = convolve(obs_data, kernel)
                            vmin = obs_data_smoothed.min()
                            vmax = obs_data_smoothed.max()
                            cb = plt.pcolormesh(lon, vel, obs_data_smoothed,
                                                norm=colors.SymLogNorm(base=10, linthresh=0.1, vmin=vmin, vmax=vmax),
                                                shading='auto', cmap='viridis', alpha=0.25)
                            plt.colorbar(cb)
                            plt.show(block=True)

                            plt.contour(lon, vel, obs_data_smoothed, levels=[0.25*vmax, 0.5*vmax], colors='black')
                            plt.show(block=True)

                    if log_comp:
                        obs_error_final = obs_error_final / np.log(10) / obs_data_final
                        obs_data_final[obs_data_final <= 0] = 1e-100
                        obs_data_final = np.log10(obs_data_final)
                        model_interp[model_interp <= 0] = 1e-100
                        model_interp = np.log10(model_interp)

                    # Calculate the chi**2 and overall likelihood
                    chi2 = np.zeros_like(model_interp)
                    loglikelihood = np.zeros_like(model_interp)
                    # if ((survey == 'COBE') or (survey == 'COBE-FIRAS')):
                    if '.idl' in file:
                        i_signal = np.where(obs_data_final != 0)
                    else:
                        i_signal = np.where(obs_data_final != 0)
                        # print(model_interp.shape,
                        #       obs_data_final.shape,
                        #       obs_error_final.shape,
                        #       model_interp[i_signal].shape,
                        #       obs_data_final[i_signal].shape,
                        #       obs_error_final[i_signal].shape
                        #       )
                    # print(obs_data_final[i_signal].size - len(param))
                    # print((obs_data_final[i_signal] - model_interp[i_signal]).any())
                    # print(((obs_data_final[i_signal] - model_interp[i_signal])** 2 / \
                    #                  obs_error_final[i_signal] ** 2).max())
                    # print('chi2, i_signal, obs_data_final, model_interp, obs_error_final\n',
                    #       chi2.shape, np.shape(i_signal), obs_data_final.shape,
                    #       model_interp.shape, obs_error_final.shape)
                    chi2[i_signal] = (obs_data_final[i_signal] - model_interp[i_signal]) ** 2 / \
                                     obs_error_final[i_signal] ** 2
                    loglikelihood[i_signal] = -chi2[i_signal] / 2 \
                                              - 0.5 * np.log(2 * np.pi * obs_error_final[i_signal]**2)
                    # elif survey == 'Planck':
                    #     i_signal = np.where(obs_data_final > obs_error_final)
                    #     chi2[i_signal] = (obs_data_final[i_signal] - model_interp[i_signal]) ** 2 / \
                    #                      obs_error_final[i_signal] ** 2
                    #     loglikelihood[i_signal] = -chi2[i_signal] / 2 \
                    #                               - 0.5 * np.log10(np.sqrt(2 * np.pi) * obs_error_final[i_signal])
                    # else:
                    #     # i_signal = np.where(obs_data_temp >= obs_data_temp.min())
                    #     i_signal = np.where(obs_data_final > obs_error_final)
                    #     chi2[i_signal] = (obs_data_final[i_signal] - model_interp[i_signal]) ** 2 \
                    #                      / obs_error_final[i_signal] ** 2
                    #     loglikelihood[i_signal] = -chi2[i_signal] / 2 \
                    #                               - 0.5 * np.log10(np.sqrt(2 * np.pi) * obs_error_final[i_signal])

                    # print(model_interp)
                    # print(obs_data[~np.isnan(obs_data)])
                    # print(chi2)
                    chi2 = np.nan_to_num(chi2, nan=0)
                    loglikelihood = np.nan_to_num(loglikelihood, nan=0)
                    # likelihood[likelihood == 0] = 1e10
                    # input()
                    np.save(path + 'fit_results/{}/{}/{}/{}_chi2.npy'.format(survey, file.replace('.fits', ''),
                                                                             transition, dir_model.replace('/', '')
                                                                             .replace('channel_intensity.fits', '')
                                                                             + '_' + comp),
                            chi2)
                    np.save(path + 'fit_results/{}/{}/{}/{}_loglikelihood.npy'.format(survey, file.replace('.fits', ''),
                                                                                      transition,
                                                                                      dir_model.replace('/', '')
                                                                                      .replace('channel_intensity.fits',
                                                                                               '') + '_' + comp),
                            loglikelihood)

                    chi2_grid.append(chi2.sum())
                    loglikelihood_grid.append(loglikelihood.sum())
                    #                             print('  ', likelihood.min())

                    # print(obs_data.min(), obs_data.max())

                    if PLOT:
                        if (survey == 'COBE-DIRBE') or (survey == 'COBE-FIRAS'):
                            fig, ax = plt.subplots(1, 1, figsize=(7, 7))
                            fig2, ax2 = plt.subplots(1, 1, figsize=(7, 7))
                            cm = ax.scatter(lon_model, lat_model, c=np.asarray(model_interp), cmap=cmap,
                                            norm=colors.SymLogNorm(base=10, linthresh=0.1,
                                                                   vmin=0, vmax=vmax))
                            cm2 = ax2.scatter(lon_obs, lat_obs, c=obs_data[:, :, int(transition_indeces[i])],
                                              norm=colors.SymLogNorm(base=10, linthresh=0.1,
                                                                     vmin=0, vmax=vmax), cmap=cmap)
                            cb = fig.colorbar(cm, ax=ax, fraction=0.02)
                            cb2 = fig.colorbar(cm2, ax=ax2, fraction=0.02)
                            ax.set_xlabel(r'Longitude ($^\circ$)', fontsize=16)
                            ax.set_ylabel(r'Latitude ($^\circ$)', fontsize=16)
                            ax2.set_xlabel(r'Longitude ($^\circ$)', fontsize=16)
                            ax2.set_ylabel(r'Latitude ($^\circ$)', fontsize=16)
                        else:
                            # vel_grid, lon_grid = np.meshgrid(vel, lon)
                            # model_interp = f_model_interp(vel_grid, lon_grid)
                            # print(0.1*vmax)
                            i_nan = np.isnan(obs_data)
                            fig, ax = plt.subplots(1, 1, figsize=(7, 7))
                            fig2, ax2 = plt.subplots(1, 1, figsize=(7, 7))
                            if True:#mission == 'COGAL':
                                cm = ax.pcolormesh(lon_model, vel_model, model_interp,
                                                   norm=colors.SymLogNorm(base=10, linthresh=0.1,
                                                                          vmin=vmin, vmax=vmax),
                                                   shading='auto', cmap=cmap)
                                cm2 = ax2.pcolormesh(lon_obs, vel_obs, obs_data_temp,
                                                     norm=colors.SymLogNorm(base=10, linthresh=0.1,
                                                                            vmin=vmin, vmax=vmax),
                                                     shading='auto', cmap=cmap)
                                # extents = [lon.max(), lon.min(), vel.min(), vel.max()]
                                # cm = plt.imshow(model_interp, cmap=cmap, extent=extents, aspect='auto', origin='upper',
                                #      norm=colors.SymLogNorm(base=10, linthresh=0.1, vmin=0, vmax=vmax))
                                ax.contour(lon_obs, vel_obs, obs_data_temp, levels=[0.1*vmax], colors='xkcd:magenta')
                                ax.text(0.1, 0.9, '{:.3f}'.format(0.25*vmax), color='xkcd:magenta',
                                        fontsize=16, transform=ax.transAxes)
                            else:
                                cm = ax.pcolormesh(lon, vel, model_interp.T,
                                                   norm=colors.SymLogNorm(base=10, linthresh=0.1,
                                                                          vmin=0, vmax=vmax),
                                                   shading='gouraud', cmap=cmap)
                                ax.contour(lon, vel, obs_data, levels=[0.1*vmax], colors='black')
                            cb = fig.colorbar(cm, ax=ax, fraction=0.02)
                            cb2 = fig.colorbar(cm2, ax=ax2, fraction=0.02)
                            ax.set_xlabel(r'Longitude ($^\circ$)', fontsize=16)
                            ax.set_ylabel(r'Radial Velocity ($\frac{km}{s}$)', fontsize=16)
                            ax2.set_xlabel(r'Longitude ($^\circ$)', fontsize=16)
                            ax2.set_ylabel(r'Radial Velocity ($\frac{km}{s}$)', fontsize=16)
                            # ax.invert_xaxis()
                        fig.savefig(
                            path + 'fit_results/{}/{}/{}/synthetic_observation_{}.png'.format(
                                survey, file.replace('.fits', ''), transition, dir_model.split('/')[0]))
                        fig2.savefig(
                            path + 'fit_results/{}/{}/{}/synthetic_observation_smoothed.png'.format(
                                survey, file.replace('.fits', ''), transition))
                        plt.close('all')

                    model.close()
                    # del f_model_interp
                    del model_interp
                    del chi2
                    del loglikelihood
                    del lon_model
                    del lat_model
                    del vel_model

                if np.size(loglikelihood_grid) == 0:
                    continue
                # print(loglikelihood_grid)

                i_bestfit2 = np.where(loglikelihood_grid == np.max(loglikelihood_grid))[0][0]
                # print(i_bestfit2)
                print('\n    The best-fitting model for {} with transition {}\n'.format(file.replace('.fits', ''), transition) +
                      '  has parameters ' + model_dir.format(*params[i_bestfit2]))

                np.save(path + 'fit_results/{}/{}/{}/{}_chi2.npy'.format(survey, file.replace('.fits', ''),
                                                                         transition, comp),
                        chi2_grid)
                np.save(path + 'fit_results/{}/{}/{}/{}_loglikelihood.npy'.format(survey, file.replace('.fits', ''),
                                                                                  transition, comp),
                        loglikelihood_grid)
                np.save(path + 'fit_results/{}/{}/{}/{}_parameters.npy'.format(survey, file.replace('.fits', ''),
                                                                               transition, comp),
                        params)

            if '.fits' in file:
                obs.close()
            del obs_data
            del obs_error
            del lon_obs
            del lat_obs
            # if not vel is None:
            #     del vel

        try:
            i_bestfit1 = np.where(chi2_grid == np.min(chi2_grid))[0][0]
            print('  The best-fitting model has parameters' +
                  model_dir.format(*params[i_bestfit1]))
        except ValueError:
            pass
        # i_bestfit2 = np.where(loglikelihood_grid == np.max(loglikelihood_grid))[0][0]


    return


def model_selection_new(path='/mnt/yanitski_backup/yanitski/projects/pdr/KT3_history/MilkyWay', 
        missions=None, lat=None, f_idx=0, error_cal=0, model_dir='', 
        model_param=[[]], comp_type='pv', log_comp=True, cmap='gist_ncar', 
        spectra=True, PLOT=False, PRINT=False, debug=False):
    '''
    Like `model_selection_old`, but using `SyntheticModel()` and `Observation` instances.
    '''
    
    # Check that the missions are specified properly.
    if missions == '' or missions == None or missions == []:
        missions = os.listdir(path.replace('KT3_history', 'observational_data'))
    elif isinstance(missions, str):
        missions = [missions]
    elif isinstance(missions, list):
        # Use both COBE instruments if simply 'COBE' is specified
        if 'COBE' in missions:
            if 'COBE-FIRAS' in missions: missions.remove('COBE-FIRAS')
            if 'COBE-DIRBE' in missions: missions.remove('COBE-DIRBE')
            missions.append('COBE-FIRAS')
            missions.append('COBE-DIRBE')
    else:
        print('Please specify a list of missions to compare the models.')
        return
  
    if model_dir == '' or model_param == [[]]:
        print('Please specify both model_dir and model_param.')
        return

    if isinstance(f_idx, int):
        f_idx = len(missions) * [f_idx]
    elif isinstance(f_idx, (tuple, list)):
        if len(f_idx) == 1:
            f_idx = len(missions) * f_idx
        elif len(f_idx) != len(missions):
            print('A valid index must be specified for f_idx.')
            return
        else:
            pass
    else:
        print('An int or list(int) must be specified for f_idx.')
        return
  
    model_params = np.meshgrid(*model_param)
    # model_params = zip(np.transpose([model_params[n].flatten() for n in range(len(model_param))]))
    model_params = np.transpose([model_params[n].flatten() for n in range(len(model_param))])

    # pprint([model_dir.format(*m) for m in model_params])

    if log_comp:
        comp = comp_type + '_logT'
    else:
        comp = comp_type

    obs_data = []
    lon = []
    vel = []

    for s, survey in enumerate(missions):

        print('\n\n  {}'.format(survey))
        print('  ' + '-'*len(survey))

        if survey in ('HiGAL', 'GRS', 'CGPS', 'SGPS', 'VGPS'):
            continue

        if path[-1] != '/': path += '/'
        directory = path.replace('KT3_history', 'observational_data') + survey + '/regridded/temp/'
        obs = Observation(base_dir=path.replace('KT3_history', 'observational_data'))
        obs.load_survey(survey=survey)
        # if mission == 'COBE-FIRAS':
        #     directory += 'temp/'
        files = os.listdir(directory)

        if len(files) == 0:
            print('no data is available.')
            exit()

        # Loop through the different models
        # for i_obs in range(len(obs_data)):
        for f, file in enumerate(obs.files):

            if f < f_idx[s]:
                continue

            if not ('.fits' in file or '.idl' in file or '.csv' in file) or 'error' in file:
                continue

            print(file)

            transitions = obs.obs_iid[f]
            transition_indeces = obs.obs_i_iid[f]
            if comp_type == 'integrated':
                obs_data = obs.get_intensity(idx=f, integrated=True)
            else:
                obs_data = obs.get_intensity(idx=f, integrated=False)
            obs_data_temp = deepcopy(obs_data)
            obs_error = obs.obs_error_data[f]
            obs_error_conf = obs.obs_error_conf_data[f]
            lon_obs = obs.obs_lon[f]
            lon_obs[lon_obs>180] = lon_obs[lon_obs>180] - 360
            lat_obs = obs.obs_lat[f]
            vel_obs = obs.obs_vel[f]
            i_lat_obs_init = obs.get_obs_extent(idx=f, kind='index')[1]

            for i, transition in enumerate(transitions):

                print('\n  fitting', transition)

                chi2_grid = []
                loglikelihood_grid = []
                params = []

                # if (survey == 'COBE-DIRBE') or (survey == 'COBE-FIRAS'):
                #     if '.idl' in file:
                #         # obs_error_conf = error_correction(obs_data[:, int(transition_indeces[i])],
                #         #                                   conf='axisymmetric')
                #         vmin = obs_data[:, int(transition_indeces[i])].min()
                #         vmax = obs_data[:, int(transition_indeces[i])].max()
                #     else:
                #         # obs_error_conf = error_correction(obs_data[int(transition_indeces[i]), :, :],
                #         #                                   conf='axisymmetric')
                #         vmin = obs_data[int(transition_indeces[i]), :, :].min()
                #         vmax = obs_data[int(transition_indeces[i]), :, :].max()
                #     # print(transition)
                #     # print(vmin, vmax)
                # elif not (survey=='GOT_C+' and '.csv' in file):
                #     # obs_error_conf = error_correction(obs_data, conf='axisymmetric')
                #     vmin = np.nanmin(obs_data)
                #     vmax = np.nanmax(obs_data)
                # else:
                #     obs_error_conf = 0

                if PRINT:
                    print('    Average error in observation: {}'.format(np.nanmean(obs_error)))
                    print('    Correction due to axisymmetry: {}'.format(np.nanmean(obs_error_conf)))
                
                model = models.SyntheticModel(base_dir=path)

                for param in model_params:

                    if PRINT:
                        print('  ' + model_dir.format(*param) + '    \r', end='')

                    # Check existance of model
                    dir_model = model_dir.format(*param)#
                    if not os.path.isfile(path + dir_model + 'synthetic_intensity.fits'): 
                        print('No intensity')
                        continue

                    # Open the model
                    params.append(param)
                    model.load_model(directory=dir_model)

                    # Locate the species transition (only one dust value is considered; constant background)
                    # if transition in model.dust:
                    #     i_spec = np.where(model.dust == transition)[0]
                    # elif transition in model.species:
                    #     i_spec = np.where(model.species == transition)[0]
                    # elif transition == 'HI':
                    #     i_spec = np.zeros(1)
                    # else:
                    #     print(f'  {transition} not in model.  ')
                    #     break
                    # else:
                    #     i_spec = i_spec[0]

                    if os.path.isdir(path + f"fit_results/{survey}/{file.replace('.fits', '')}/{transition}/") == False:
                        os.makedirs(path + f"fit_results/{survey}/{file.replace('.fits', '')}/{transition}/")

                    # Create arrays for the longitude and velocity axes
                    lat_model = model.map_lat
                    lon_model = model.map_lon
                    vel_model = model.map_vel

                    # Identify the common latitudes between the observations and the model
                    if isinstance(lat, (int, float)):
                        lat_min = lat
                        lat_max = lat
                    else:
                        lat_min = lat_obs[i_lat_obs_init].min()
                        lat_max = lat_obs[i_lat_obs_init].max()
  
                    i_lon_model = np.arange(lon_model.size)
                    i_lat_model = np.where((lat_model >= lat_min)
                                            & (lat_model <= lat_max))[0]
                    # print(lat_obs, lat_model[i_lat_model])
                    i_lat_obs = np.where((lat_obs >= lat_model.min())
                                         & (lat_obs <= lat_model.max()))[0]

                    ix_err = np.ix_(i_lat_model, i_lon_model)
                    if survey == 'COBE-FIRAS':
                        bin_edge = np.array([*((lon_obs[:-1]+lon_obs[1:])/2), 182.5])
                        if ' + ' in transition:
                            model_intensity = model.get_species_intensity(transition=transition.split(' + '), include_dust=False, integrated=True).sum(0)
                        else:
                            model_intensity = model.get_species_intensity(transition=transition, include_dust=False, integrated=True)
                        lon_model_alt = deepcopy(lon_model)
                        lon_model_alt[lon_model_alt<-177.5] = lon_model_alt[lon_model_alt<-177.5] + 360
                        model_data_temp = binned_statistic(lon_model_alt, model_intensity, statistic='mean', bins=bin_edge)[0]
                        model_data = np.array([model_data_temp[:, -1], *model_data_temp.T]).T
                        obs_data = np.array(model_data_temp.shape[0]*[obs_data_temp[:, transition_indeces[i]]])
                        ix = np.ix_(i_lat_model, np.arange(lon_obs.size))
                        ix_obs = np.ix_(*(np.arange(_) for _ in obs_data.shape))
                        obs_error_final = np.array(model_data_temp.shape[0]*[np.sqrt(obs_error[:, transition_indeces[i]]**2+obs_error_conf[:, transition_indeces[i]]**2)])[ix_obs][ix]
                    elif survey == 'GOT_C+' and '.csv' in file:
                        obs_lon = np.unique(lon_obs)
                        obs_vel = np.unique(vel_obs)
                        # obs_lon[obs_lon>180] = obs_lon[obs_lon>180]-360
                        bin_mid = obs_lon
                        bin_edge = np.array([-180, *((bin_mid[:-1]+bin_mid[1:])/2), 180])
                        model_intensity = model.get_species_intensity(transition='C+ 1', include_dust=False, integrated=False)[:, 4, :]
                        model_data = binned_statistic(model.map_lon, model_intensity, statistic='mean', bins=bin_edge)[0]
                        obs_data = np.zeros_like(model_data)
                        obs_error_final = np.zeros_like(model_data)
                        for l, lon in enumerate(obs_lon):
                            obs_data[:, l] = obs_data_temp[lon_obs==lon]
                            obs_error_final[:, l] = np.sqrt(obs_error[lon_obs==lon].max()**2 + obs_error_conf[lon_obs==lon]**2)
                        ix = np.ix_(*(np.arange(_) for _ in obs_data.shape))
                        ix_obs = ix
                    elif transition == 'HI':
                        if comp_type == 'integrated':
                            ix = np.ix_(i_lat_model, i_lon_model)
                            ix_obs = np.ix_(i_lat_obs, i_lon_model)
                            model_data = model.get_hi_intensity(integrated=True)
                            obs_error_final = np.sqrt((np.sqrt(vel_model.size)*obs_error)**2+(np.sqrt((obs_error_conf**2).sum(axis=0))**2))[i_obs][ix_err]
                        elif comp_type == 'pv':
                            ix = np.ix_(np.arange(vel_model.size), i_lat_model, i_lon_model)
                            ix_obs = np.ix_(np.arange(vel_model.size), i_lat_obs, i_lon_model)
                            model_data = model.get_hi_intensity(integrated=False)
                            obs_error_final = np.sqrt(obs_error**2+obs_error_conf**2)[ix_obs][ix]
                    elif transition in model.species:
                        if comp_type == 'integrated':
                            ix = np.ix_(i_lat_model, i_lon_model)
                            ix_obs = np.ix_(i_lat_obs, i_lon_model)
                            model_data = model.get_species_intensity(transition=transition, include_dust=False, integrated=True)
                            obs_error_final = np.sqrt((np.sqrt(vel_model.size)*obs_error)**2+(np.sqrt((obs_error_conf**2).sum(axis=0))**2))[ix_obs][ix_err]
                        elif comp_type == 'pv':
                            ix = np.ix_(np.arange(vel_model.size), i_lat_model, i_lon_model)
                            ix_obs = np.ix_(np.arange(vel_model.size), i_lat_obs, i_lon_model)
                            model_data = model.get_species_intensity(transition=transition, include_dust=False, integrated=False)
                            obs_error_final = np.sqrt(obs_error**2
                                    +obs_error_conf**2
                                    +(obs_data*error_cal)**2)[ix_obs][ix]
                    elif transition in model.dust:
                        ix = np.ix_(i_lat_model, i_lon_model)
                        ix_obs = np.ix_(i_lat_obs, i_lon_model)
                        model_data = model.get_dust_intensity(wavelength=transition)
                        obs_error_final = np.sqrt(obs_error**2+obs_error_conf**2)[ix_obs][ix_err]

                    model_interp = model_data[ix]
                    obs_data_final = obs_data[ix_obs][ix]
  
                    model_interp[model_interp == 0] = np.nan

                    if np.isnan(obs_data_final).all():
                        print('Observation improperly indexed!')
                        break

                    # if ((survey == 'COBE') or (survey == 'COBE-FIRAS')):
                    if '.idl' in file:
                        i_signal = np.where(~((obs_data_final == 0) | np.isnan(obs_data_final)
                                              | np.isnan(model_interp)))
                    else:
                        i_signal = np.where(~((obs_data_final == 0) | np.isnan(obs_data_final)
                                              | np.isnan(model_interp)))
                    
                    if log_comp:
                        obs_error_final[i_signal] = obs_error_final[i_signal] / np.log(10) / obs_data_final[i_signal]
                        obs_data_final[obs_data_final <= 0] = 1e-100
                        obs_data_final[i_signal] = np.log10(obs_data_final[i_signal])
                        model_interp[model_interp <= 0] = 1e-100
                        model_interp[i_signal] = np.log10(model_interp[i_signal])

                    # Calculate the chi**2 and overall likelihood
                    chi2 = np.zeros_like(model_interp)
                    loglikelihood = np.zeros_like(model_interp)
                        
                    chi2[i_signal] = (obs_data_final[i_signal] - model_interp[i_signal]) ** 2 / \
                                    obs_error_final[i_signal] ** 2
                    loglikelihood[i_signal] = -chi2[i_signal] / 2 \
                                            - 0.5 * np.log(2 * np.pi * obs_error_final[i_signal]**2)
                    
                    chi2 = np.nan_to_num(chi2, nan=0)
                    loglikelihood = np.nan_to_num(loglikelihood, nan=0)

                    if ~chi2.any() or np.isnan(chi2).all():
                        print(chi2.any() or np.isnan(chi2).all())
                        print(dir_model)
                        print('Observation improperly indexed -- all NaN!')
                        return
                    
                    np.save(path + f"fit_results/{survey}/{file.replace('.fits', '')}/{transition}/"
                            f"{dir_model.replace('/', '').replace('channel_intensity.fits', '') + '_' + comp}_chi2.npy", chi2)
                    np.save(path + f"fit_results/{survey}/{file.replace('.fits', '')}/{transition}/"
                            f"{dir_model.replace('/', '').replace('channel_intensity.fits','') + '_' + comp}_loglikelihood.npy", loglikelihood)

                    chi2_grid.append(chi2.sum())
                    loglikelihood_grid.append(loglikelihood.sum())
  
                    # model.close()
                    # del f_model_interp
                    del model_interp
                    del chi2
                    del loglikelihood
                    del lon_model
                    del lat_model
                    del vel_model

                if np.size(loglikelihood_grid) == 0:
                    continue
                # print(loglikelihood_grid)

                i_bestfit2 = np.where(loglikelihood_grid == np.max(loglikelihood_grid))[0][0]
                # print(i_bestfit2)
                print(f"\n\n    The best-fitting model for {file.replace('.fits', '')} with transition {transition}"
                      f"  \n    has parameters {model_dir.format(*params[i_bestfit2])}")

                np.save(path + f"fit_results/{survey}/{file.replace('.fits', '')}/{transition}/{comp}_chi2.npy", chi2_grid)
                np.save(path + f"fit_results/{survey}/{file.replace('.fits', '')}/{transition}/{comp}_loglikelihood.npy", loglikelihood_grid)
                np.save(path + f"fit_results/{survey}/{file.replace('.fits', '')}/{transition}/{comp}_parameters.npy", params)

            # if '.fits' in file:
            #     obs.close()
            del obs_data
            del obs_error
            del lon_obs
            del lat_obs
  
        try:
            i_bestfit1 = np.where(chi2_grid == np.min(chi2_grid))[0][0]
            print(f'\n\n  The best-fitting model has parameters   {model_dir.format(*params[i_bestfit1])}')
        except ValueError:
            pass

    return


def line_ratio_comparison(obs_path='/mnt/hpc_backup/yanitski/projects/pdr/observational_data/MilkyWay/',
                          model_path='/mnt/hpc_backup/yanitski/projects/pdr/KT3_history/fit_results/MilkyWay/',
                          missions=[], transitions=[], survey_files=[], file_format='', model_param=[[]],
                          log_comp=True, lat=None, comp_type='integrated',
                          figsize=(), ylabel='', ylim=(), title='', notch=True, boxwidth=0.5,
                          label_rotation=30, label_fontsize=16, fontsize=20,
                          save_plot=False, output_file='', output_format='png',
                          debug=False, verbose=False, violin_width=1, violin_spacing=1, **kwargs):
    '''
    Plot comparison of line ratios.
    '''

    survey_paths = [[], []]  # list for paths to plotted data (in order observed, synthetic)
    survey_maps = [[], []]   # list for unprocessed (pre-calculation) data to plot (in order observed, synthetic)
    survey_i_lat_init = []   # list for initial observed lattitude indeces where there is an observation
    survey_i_vel_init = []   # list for initial observed velocity indeces where there is an observation
    survey_ix = [[], []]     # list of indeces used in comparison (in order observed, synthetic)
    survey_ratios = [[], []] # list of processed (post-calculation) data to plot (in order observed, synthetic)
    # list of labels for each ratio (in order observed, synthetic)
    survey_labels = [['/'.join(missions)], []]
    # list of available transitions in the data (in order observed, synthetic; used for debugging)
    survey_transitions = [[], []]
    survey_lons = [[], []]   # list of data longitude (in order observed, synthetic)
    survey_lats = [[], []]   # list of data lattitude (in order observed, synthetic)
    survey_vels = [[], []]   # list of data velocities (in order observed, synthetic; 0 for 2D maps)

    if ('COBE-FIRAS' in missions or 'COBE-DIRBE' in missions or 'Planck' in missions):
        comp_type = 'integrated'

    for i, survey in enumerate(missions):
        print(survey)
        survey_paths[0].append(obs_path + missions[i] + '/regridded/temp/' + survey_files[i])
        if '.fits' in survey_paths[0][i]:
            obs = fits.open(survey_paths[0][i])
            survey_transitions[0].append(np.asarray(obs[0].header['TRANSL'].split(', ')))
            i_trans_obs = survey_transitions[0][i] == transitions[i]
            if i_trans_obs.any() == False:
                print('Invalid transition {} for file {}'.format(transitions[i], survey_paths[i]))
                return
            survey_lons[0].append(np.linspace(obs[0].header['CRVAL1'] - obs[0].header['CDELT1'] * (
                                              obs[0].header['CRPIX1'] - 1),
                                              obs[0].header['CRVAL1'] + obs[0].header['CDELT1'] * (
                                              obs[0].header['NAXIS1'] - obs[0].header['CRPIX1']),
                                              num=obs[0].header['NAXIS1']))
            survey_lats[0].append(np.linspace(obs[0].header['CRVAL2'] - obs[0].header['CDELT2'] * (
                                              obs[0].header['CRPIX2'] - 1),
                                              obs[0].header['CRVAL2'] + obs[0].header['CDELT2'] * (
                                              obs[0].header['NAXIS2'] - obs[0].header['CRPIX2']),
                                              num=obs[0].header['NAXIS2']))
            if survey == 'COBE-FIRAS':
                survey_vels[0].append(np.asarray([0]))
                survey_maps[0].append(obs[0].data[i_trans_obs, :, :][0, :, :])
                survey_i_lat_init.append(survey_maps[0][i].any(1))
            elif survey == 'COBE-DIRBE':
                survey_vels[0].append(np.asarray([0]))
                survey_maps[0].append(obs[0].data[i_trans_obs, :, :][0, :, :])
                survey_i_lat_init.append(survey_maps[0][i].any(1))
            elif survey == 'Planck':
                survey_vels[0].append(np.asarray([0]))
                survey_maps[0].append(obs[0].data[:, :])
                survey_i_lat_init.append(survey_maps[0][i].any(1))
            elif ('COBE-FIRAS' in missions or 'Planck' in missions) or comp_type == 'integrated':
                survey_vels[0].append(np.linspace(obs[0].header['CRVAL3'] - obs[0].header['CDELT3'] * (
                                                  obs[0].header['CRPIX3'] - 1),
                                                  obs[0].header['CRVAL3'] + obs[0].header['CDELT3'] * (
                                                  obs[0].header['NAXIS3'] - obs[0].header['CRPIX3']),
                                                  num=obs[0].header['NAXIS3']))
                survey_maps[0].append(np.trapz(obs[0].data, survey_vels[0][i], axis=0))
                survey_i_lat_init.append(survey_maps[0][i].any(1))
            elif comp_type == 'pv':
                survey_vels[0].append(np.linspace(obs[0].header['CRVAL3'] - obs[0].header['CDELT3'] * (
                                                  obs[0].header['CRPIX3'] - 1),
                                                  obs[0].header['CRVAL3'] + obs[0].header['CDELT3'] * (
                                                  obs[0].header['NAXIS3'] - obs[0].header['CRPIX3']),
                                                  num=obs[0].header['NAXIS3']))
                survey_maps[0].append(obs[0].data)
                survey_i_lat_init.append(survey_maps[0][i].any(0).any(1))
                survey_i_vel_init.append(survey_maps[0][i].any(2).any(1))
            else:
                print('''Enter a valid string for `comp_type`: 'integrated' for a 2D comparison of maps '''
                      + '''and 'pv' for a 3D comparison of spectroscopic data.''')
                return
        elif '.idl' in survey_paths[0][i]:
            obs = readsav(survey_paths[0][i])
            survey_lons[0].append(obs['long'])
            survey_lats[0].append(np.asarray([0]))
            survey_transitions[0].append(transitions[i])
            i_trans_obs = cobe_idl_transitions == transitions[i]
            if i_trans_obs.any() == False:
                print()
                return
            survey_maps[0].append((obs['amplitude'][:, i_trans_obs] / (2.9979**3)
                                   * (cobe_idl_linfrq[i_trans_obs]**3) * 2 * 1.38 / 10 ** 8).reshape(0, -1))
            survey_i_lat_init.append(np.ones(1, dtype=bool))
        else:
            print('Invalid mission path {}.'.format(survey_paths[0][i]))
            return

    print(survey_maps[0][0].shape, survey_maps[0][1].shape)

    i_lon_obs = []
    i_lon_model = []
    i_lat_obs = []
    i_lat_model = []
    if ((isinstance(lat, int) or isinstance(lat, float)) and comp_type == 'pv' and
        not ('COBE-FIRAS' in missions or 'COBE-DIRBE' in missions or 'Planck' in missions)):
        lat_min = lat
        lat_max = lat
    else:
        lat_min = np.amin([survey_lats[0][i][survey_i_lat_init[i]] for i in range(2)])
        lat_max = np.amax([survey_lats[0][i][survey_i_lat_init[i]] for i in range(2)])
    if comp_type == 'pv':
        vel_min = np.amin([survey_vels[0][i][survey_i_vel_init[i]] for i in range(2)])
        vel_max = np.amax([survey_vels[0][i][survey_i_vel_init[i]] for i in range(2)])

    model_param_grid = np.meshgrid(*np.asarray(model_param, dtype=object))
    model_params = zip(*[model_param_grid[n].flatten() for n in range(len(model_param_grid))])

    for param in model_params:

        survey_labels[1].append(file_format.format(*param))
        model_dir = model_path + survey_labels[1][-1] + 'synthetic_intensity.fits'

        model = fits.open(model_dir)

        # Create arrays for the longitude and velocity axes of the synthetic observation if required
        if len(survey_lats[1]) == 0:
            survey_lats[1].append(
                np.linspace(model[1].header['CRVAL3']
                            - model[1].header['CDELT3'] * (model[1].header['CRPIX3'] - 0.5),
                            model[1].header['CRVAL3']
                            + model[1].header['CDELT3'] * (model[1].header['NAXIS3'] - model[1].header['CRPIX3'] - 0.5),
                            num=model[1].header['NAXIS3']) * 180 / np.pi)
            survey_lons[1].append(
                np.linspace(model[1].header['CRVAL2']
                            - model[1].header['CDELT2'] * (model[1].header['CRPIX2'] - 0.5),
                            model[1].header['CRVAL2']
                            + model[1].header['CDELT2'] * (model[1].header['NAXIS2'] - model[1].header['CRPIX2'] - 0.5),
                            num=model[1].header['NAXIS2']) * 180 / np.pi)
            # if not (mission == 'COGAL'):
            #     lon_model[lon_model<0] += 360
            survey_vels[1].append(
                np.linspace(model[1].header['CRVAL4']
                            - model[1].header['CDELT4'] * (model[1].header['CRPIX4'] - 0.5),
                            model[1].header['CRVAL4']
                            + model[1].header['CDELT4'] * (model[1].header['CRPIX4'] - 0.5),
                            num=model[1].header['NAXIS4']))

        survey_maps[1].append([])

        for i, transition in enumerate(transitions):
            if transition == 'Dust':
                survey_transitions[1].append(np.asarray(model[2].header['DUST'].split(', ')))
                # survey_maps[1][-1].append(deepcopy(model[2].data[0, :, :, 0]))
                survey_maps[1][-1].append(deepcopy(model[2].data[:, :, 0]))
            else:
                survey_transitions[1].append(np.asarray(model[1].header['SPECIES'].split(', ')))
                i_transition = np.where(survey_transitions[1][-1] == transition)[0]
                if i_transition.size == 1 and comp_type == 'integrated' \
                        and not ('COBE-FIRAS' in missions or 'Planck' in missions):
                    survey_maps[1][-1].append(np.trapz(model[1].data[:, :, :, i_transition][:, :, :, 0]
                                                       # - model[2].data[:, :, :, 0], survey_vels[1][0], axis=0))
                                                       - model[2].data[:, :, 0], survey_vels[1][0], axis=0))
                elif i_transition.size == 1:
                    survey_maps[1][-1].append(model[1].data[:, :, :, i_transition][:, :, :, 0]
                                              - model[2].data[:, :, 0])
                else:
                    print('Transition {} not found in model'.format(transition))
                    return

            if len(survey_ix[0]) < 2:
                i_lon_obs = np.linspace(0, survey_lons[0][i].size-1, num=survey_lons[0][i].size, dtype=int)
                i_lon_model = np.linspace(0, survey_lons[1][0].size-1, num=survey_lons[1][0].size, dtype=int)
                i_lat_model = np.where((survey_lats[1][0] >= lat_min)
                                       & (survey_lats[1][0] <= lat_max))[0]
                i_lat_obs = np.where((survey_lats[0][i] >= survey_lats[1][0][i_lat_model].min())
                                     & (survey_lats[0][i] <= survey_lats[1][0][i_lat_model].max()))[0]
                if comp_type == 'integrated': #compare 2D maps
                    survey_ix[0].append(np.ix_(i_lat_obs, i_lon_obs))
                    survey_ix[1].append(np.ix_(i_lat_model, i_lon_model))
                elif comp_type == 'pv': #compare 3D maps
                    i_vel_model = np.where((survey_vels[1][0] >= vel_min)
                                           & (survey_vels[1][0] <= vel_max))[0]
                    i_vel_obs = np.where((survey_vels[0][i] >= survey_vels[1][0][i_vel_model].min())
                                         & (survey_vels[0][i] <= survey_vels[1][0][i_vel_model].max()))[0]
                    survey_ix[0].append(np.ix_(i_vel_obs, i_lat_obs, i_lon_obs))
                    survey_ix[1].append(np.ix_(i_vel_model, i_lat_model, i_lon_model))

        if len(survey_ratios[0]) == 0:
            survey_ratios[0].append((survey_maps[0][0][survey_ix[0][0]]
                                     / survey_maps[0][1][survey_ix[0][1]]).flatten())
        # for i in range(len(survey_maps[1])):
        survey_ratios[1].append((survey_maps[1][-1][0][survey_ix[1][0]]
                                 / survey_maps[1][-1][1][survey_ix[1][1]]).flatten())

    labels = [*survey_labels[0], *survey_labels[1]]
    if len(figsize) == 0:
        figsize = (violin_spacing*(len(labels)), 7)
    fig, ax = plt.subplots(1, 1, figsize=(violin_spacing*len(labels), 7))
    violin_positions = np.arange(0, violin_spacing*len(labels), violin_spacing)
    if log_comp == True:
        comp = '_logT'
        i_nan_obs = (survey_ratios[0][0] <= 0) | np.isnan(survey_ratios[0][0])
        i_nan_model = (survey_ratios[1][0] <= 0) | np.isnan(survey_ratios[1][0])
        data = [np.log10(survey_ratios[0][0][~i_nan_obs]),
                *[np.log10(survey_ratios[1][i][~i_nan_model]) for i in range(len(survey_ratios[1]))]]
        ax.violinplot(data, positions=violin_positions, widths=violin_width)
        if ylabel == '' and comp_type == 'integrated':
            ax.set_ylabel(r'$log_{10} ' + r'\left( \varpi_{' + r'{}'.format(transitions[0]) + r'} \ / \ \varpi_{'
                          + r'{}'.format(transitions[1]) + r'} \right)$', fontsize=fontsize)
        if ylabel == '':
            ax.set_ylabel(r'$log_{10} ' + r'\left( T_\mathrm{B, ' + r'{}'.format(transitions[0]) + r'} \ / \ \varpi_{'
                          + r'{}'.format(transitions[1]) + r'} \right)$', fontsize=fontsize)
        else:
            ax.set_ylabel(ylabel)
    else:
        comp = ''
        data = [survey_ratios[0][0], *survey_ratios[1]]
        ax.violinplot(data, positions=violin_positions, widths=violin_width)
        if ylabel == '' and comp_type == 'integrated':
            ax.set_ylabel(r'$\varpi_{' + r'{}'.format(transitions[0]) + r'} \ / \ \varpi_{'
                          + r'{}'.format(transitions[1]) + r'}$', fontsize=fontsize)
        elif ylabel == '':
            ax.set_ylabel(r'$T_\mathrm{B, ' + r'{}'.format(transitions[0]) + r'} \ / \ \varpi_{'
                          + r'{}'.format(transitions[1]) + r'}$', fontsize=fontsize)
        else:
            ax.set_ylabel(ylabel)
    if len(ylim):
        ax.set_ylim(ylim)
    ax.set_xticks(violin_positions, labels=labels)
    ax.xaxis.set_tick_params(labelrotation=label_rotation, labelsize=label_fontsize)
    #xticknames = ax.get_xticklabels()
    #plt.setp(xticknames, rotation=label_rotation, fontsize=label_fontsize)
    if len(ylim):
        ax.set_ylim(ylim)
    if title == '':
        ax.set_title('{} / {}'.format(*transitions), fontsize=fontsize)
    else:
        ax.set_title(title)
    #plt.tight_layout()
    if save_plot:
        if output_file == '':
            filename = '_'.join(file_format.replace('/', '').split('_')[1:]).replace('{}', '_')
            current_output_file = 'boxplot_ratio_{}-{}_{}-{}_{}'.format(survey_paths[0][0].split('/')[-1].split('.')[0],
                                                                        transitions[0],
                                                                        survey_paths[0][1].split('/')[-1].split('.')[0],
                                                                        transitions[1], filename) + comp
        else:
            current_output_file = output_file
        plt.savefig(model_path + 'fit_results/Plots/{}.{}'.format(current_output_file, output_format),
                    format=output_format)
    else:
        plt.show()
    plt.close()

    return


def violin_comparison(path='/mnt/hpc_backup/yanitski/projects/pdr/observational_data/MilkyWay/', 
                      file_format='',
                      surveys=[], model_param=[[]], log_comp=True, lat=None, comp_type='integrated',
                      ylim=[], ylabel='', title='', violin_width=0.5, violin_spacing=1, showextrema=False, 
                      label_rotation=30, label_fontsize=16, fontsize=20,
                      figsize=None, save_plot=False, output_file='', output_format='png',
                      debug=False, verbose=False, **kwargs):
    '''
    Compare models and observation as violin plots.
    '''

    # Check that the missions are specified properly.
    if surveys == '' or surveys == None or surveys == []:
        surveys = os.listdir(path)
    elif isinstance(surveys, str):
        surveys = [surveys]
    elif isinstance(surveys, list):
        # Use both COBE instruments if simply 'COBE' is specified
        if 'COBE' in surveys:
            survey.remove('COBE')
            if 'COBE-FIRAS' in surveys: surveys.remove('COBE-FIRAS')
            if 'COBE-DIRBE' in surveys: surveys.remove('COBE-DIRBE')
            surveys.append('COBE-FIRAS')
            surveys.append('COBE-DIRBE')
    else:
        print('Please specify a list of missions to compare the models.')
        return

    if file_format == '' or model_param == [[]]:
        print('Please specify both file_format and model_param.')
        return

    if log_comp:
        comp = comp_type + '_logT'
    else:
        comp = comp_type

    model_param_grid = np.meshgrid(*np.asarray(model_param, dtype=object))
    model_params = zip(*[model_param_grid[n].flatten() for n in range(len(model_param_grid))])

    for survey in surveys:

        if survey == 'Plots':
            continue

        obs_directory = path + survey + '/regridded/temp/'

        # if verbose:
        print('\n {}'.format(survey))

        if (survey + '_files') in kwargs.keys():
            if isinstance(kwargs[survey + '_files'], list):
                survey_files = kwargs[survey + '_files']
            elif isinstance(kwargs[survey + '_files'], str):
                survey_files = [kwargs[survey + '_files']]
            else:
                print('Incorrect files specified for survey {}'.format(survey))
                continue
        else:
            survey_files = os.listdir(path + survey + '/regridded/temp/')

        print(survey_files)

        if not os.path.exists(path.replace('observational_data', 'KT3_history') + 'fit_results/' + survey):
            os.mkdir(path.replace('observational_data', 'KT3_history') + 'fit_results/' + survey)

        if not 'Plots' in os.listdir(path.replace('observational_data', 'KT3_history') + 'fit_results/' + survey):
            os.mkdir(path.replace('observational_data', 'KT3_history') + 'fit_results/' + survey + '/Plots/')

        for survey_file in survey_files:

            if 'error' in survey_file or (not '.fits' in survey_file and not '.idl' in survey_file):
                continue

            print('\n{}'.format(survey_file))

            if '.fits' in survey_file or '.FITS' in survey_file:
                obs = fits.open(obs_directory + survey_file)
                obs_error = fits.open(obs_directory + survey_file.replace('.fits', '_error.fits'))[0].data
            elif '.sav' in survey_file or '.idl' in survey_file:
                obs = readsav(obs_directory + survey_file)
            else:
                print(survey_file)
                print('Invalid: Unknown file type in observation folder.')
                continue

            if survey == 'COBE-DIRBE' or survey == 'COBE-FIRAS':
                if survey_file == 'craig.idl':

                    # This file requires a comparison to the galactic plane
                    # lat = 0

                    # these values are hard-coded since the data is not in the file
                    linfrq = np.array([115.3, 230.5, 345.8, 424.8, 461.0, 492.2, 556.9, 576.3, 691.5, 808.1,
                                       1113, 1460, 2226, 1901, 2060, 2311, 2459, 2589, 921.8])
                    transitions = np.array(['CO 1', 'CO 2', 'CO 3', 'CO 4', 'C 3', 'CO 5',
                                            'CO 6', 'CO 7 + C 1', 'C+ 1', 'O 1', 'CO 8'])
                    transition_indeces = np.array([0, 1, 2, 4, 5, 7, 8, 9, 13, 14, 18])

                    obs_data = obs['amplitude'] / (2.9979**3) * (linfrq**3) * 2 * 1.38 / 10 ** 8
                    obs_error = obs['sigma'] / (2.9979**3) * (linfrq**3) * 2 * 1.38 / 10 ** 8
                    lon_obs = obs['long']
                    lat_obs = np.array([0])
                    i_lat_obs_init = np.ones(1, dtype=bool)
                    comp_type = 'integrated'

                else:

                    obs_data = obs[0].data
                    if debug:
                        pprint(obs[0].header)
                    # obs_error = obs[0].header['OBSERR'].split()
                    transitions = obs[0].header['TRANSL'].split(', ')
                    transition_indeces = obs[0].header['TRANSI'].split(', ')
                    lon_obs = np.linspace(obs[0].header['CRVAL1'] - obs[0].header['CDELT1'] * (obs[0].header['CRPIX1'] - 1),
                                          obs[0].header['CRVAL1'] + obs[0].header['CDELT1'] * (
                                                      obs[0].header['NAXIS1'] - obs[0].header['CRPIX1']),
                                          num=obs[0].header['NAXIS1'])
                    lat_obs = np.linspace(obs[0].header['CRVAL2'] - obs[0].header['CDELT2'] * (obs[0].header['CRPIX2'] - 1),
                                          obs[0].header['CRVAL2'] + obs[0].header['CDELT2'] * (
                                                      obs[0].header['NAXIS2'] - obs[0].header['CRPIX2']),
                                          num=obs[0].header['NAXIS2'])
                    i_lat_obs_init = obs_data.any(0).any(1)
                vel = None
                comp_type = 'pv'

            elif survey == 'Planck':
                obs_data = obs[0].data
                if debug:
                    pprint(obs[0].header)
                # obs_error = obs[0].header['OBSERR'].split()
                transitions = obs[0].header['TRANSL'].split(', ')
                transition_indeces = obs[0].header['TRANSI'].split(', ')
                lon_obs = np.linspace(obs[0].header['CRVAL1'] - obs[0].header['CDELT1'] * (obs[0].header['CRPIX1'] - 1),
                                      obs[0].header['CRVAL1'] + obs[0].header['CDELT1'] * (
                                                  obs[0].header['NAXIS1'] - obs[0].header['CRPIX1']),
                                      num=obs[0].header['NAXIS1'])
                lat_obs = np.linspace(obs[0].header['CRVAL2'] - obs[0].header['CDELT2'] * (obs[0].header['CRPIX2'] - 1),
                                      obs[0].header['CRVAL2'] + obs[0].header['CDELT2'] * (
                                                  obs[0].header['NAXIS2'] - obs[0].header['CRPIX2']),
                                      num=obs[0].header['NAXIS2'])
                i_lat_obs_init = obs_data.any(1)
                vel = None
                comp_type = 'pv'

            else:
                if debug:
                    pprint(obs[0].header)
                # obs_error = obs[0].header['OBSERR']
                transitions = obs[0].header['TRANSL'].split(', ')
                transition_indeces = obs[0].header['TRANSI'].split(', ')
                lon_obs = np.linspace(obs[0].header['CRVAL1'] - obs[0].header['CDELT1'] * (obs[0].header['CRPIX1'] - 1),
                                      obs[0].header['CRVAL1'] + obs[0].header['CDELT1'] * (
                                                  obs[0].header['NAXIS1'] - obs[0].header['CRPIX1']),
                                      num=obs[0].header['NAXIS1'])
                lat_obs = np.linspace(obs[0].header['CRVAL2'] - obs[0].header['CDELT2'] * (obs[0].header['CRPIX2'] - 1),
                                      obs[0].header['CRVAL2'] + obs[0].header['CDELT2'] * (
                                                  obs[0].header['NAXIS2'] - obs[0].header['CRPIX2']),
                                      num=obs[0].header['NAXIS2'])
                vel_obs = np.linspace(obs[0].header['CRVAL3'] - obs[0].header['CDELT3'] * (obs[0].header['CRPIX3'] - 1),
                                      obs[0].header['CRVAL3'] + obs[0].header['CDELT3'] * (
                                                  obs[0].header['NAXIS3'] - obs[0].header['CRPIX3']),
                                      num=obs[0].header['NAXIS3'])
                if comp_type == 'pv':
                    obs_data = obs[0].data
                else:
                    obs_data = np.trapz(obs[0].data, vel_obs, axis=0)
                i_lat_obs_init = obs_data.any(1)

            for i_trans_obs, transition in enumerate(transitions):

                print(transition)

                map_list = [deepcopy(obs_data)]
                map_labels = [survey + ' ' + transition]
                i_model = []

                for param in deepcopy(model_params):

                    map_labels.append(file_format.format(*param))
                    model_dir = path.replace('observational_data', 'KT3_history') + map_labels[-1] \
                                + 'synthetic_intensity.fits'

                    model = fits.open(model_dir)

                    # Create arrays for the longitude and velocity axes
                    lat_model = np.linspace(
                        model[1].header['CRVAL3']
                        - model[1].header['CDELT3'] * (model[1].header['CRPIX3'] - 0.5),
                        model[1].header['CRVAL3']
                        + model[1].header['CDELT3'] * (model[1].header['NAXIS3'] - model[1].header['CRPIX3'] - 0.5),
                        num=model[1].header['NAXIS3']) * 180 / np.pi
                    lon_model = np.linspace(
                        model[1].header['CRVAL2']
                        - model[1].header['CDELT2'] * (model[1].header['CRPIX2'] - 0.5),
                        model[1].header['CRVAL2']
                        + model[1].header['CDELT2'] * (model[1].header['NAXIS2']-model[1].header['CRPIX2'] - 0.5),
                        num=model[1].header['NAXIS2']) * 180 / np.pi
                    # if not (mission == 'COGAL'):
                    #     lon_model[lon_model<0] += 360
                    vel_model = np.linspace(
                        model[1].header['CRVAL4']
                        - model[1].header['CDELT4'] * (model[1].header['CRPIX4'] - 0.5),
                        model[1].header['CRVAL4']
                        + model[1].header['CDELT4'] * (model[1].header['CRPIX4'] - 0.5),
                        num=model[1].header['NAXIS4'])

                    if transition == 'Dust':
                        # map_list.append(deepcopy(model[2].data[0, :, :, 0]))
                        map_list.append(deepcopy(model[2].data[:, :, 0]))
                    else:
                        i_transition = np.where(np.asarray(model[1].header['SPECIES'].split(', ')) == transition)[0]
                        if i_transition.size and comp_type == 'pv':
                            map_list.append(model[1].data[:, :, :, i_transition][:, :, :, 0] -
                                            # model[2].data[:, :, :, 0], vel_model, axis=0))
                                            model[2].data[:, :, 0])
                        elif i_transition.size and comp_type == 'integrated':
                            map_list.append(np.trapz(model[1].data[:, :, :, i_transition][:, :, :, 0] -
                                                     # model[2].data[:, :, :, 0], vel_model, axis=0))
                                                     model[2].data[:, :, 0], vel_model, axis=0))
                        else:
                            print('Transition {} not found in file {}'.format(transition, survey_file))
                            break

                    if ((isinstance(lat, int) or isinstance(lat, float)) and comp_type == 'pv' and
                        not (survey == 'COBE-FIRAS' or survey == 'COBE-DIRBE' or survey == 'Planck')):
                        lat_min = lat
                        lat_max = lat
                    else:
                        lat_min = lat_obs[i_lat_obs_init].min()
                        lat_max = lat_obs[i_lat_obs_init].max()

                    i_vel_obs = np.linspace(0, vel_obs.size-1, num=vel_obs.size, dtype=int)
                    i_vel_model = np.linspace(0, vel_model.size-1, num=vel_model.size, dtype=int)
                    i_lon_obs = np.linspace(0, lon_obs.size-1, num=lon_obs.size, dtype=int)
                    i_lon_model = np.linspace(0, lon_model.size-1, num=lon_model.size, dtype=int)
                    if survey == 'COBE-FIRAS' or survey == 'DOBE-DIRBE' or survey == 'Planck':
                        i_lat_model = (lat_model >= lat_min) \
                                      & (lat_model <= lat_max)
                        i_lat_obs = (lat_obs >= lat_model[i_lat_model].min()) \
                                    & (lat_obs <= lat_model[i_lat_model].max())
                    else:
                        i_lat_model = np.where((lat_model >= lat_min)
                                               & (lat_model <= lat_max))[0]
                        i_lat_obs = np.where((lat_obs >= lat_model[i_lat_model].min())
                                             & (lat_obs <= lat_model[i_lat_model].max()))[0]
                    if not survey == 'COBE-FIRAS' and not survey == 'COBE-DIRBE' and comp_type == 'pv':
                        i_obs = np.ix_(i_vel_obs, i_lat_obs, i_lon_obs)
                    elif not survey == 'COBE-FIRAS' and not survey == 'COBE-DIRBE' and comp_type == 'integrated':
                        i_obs = np.ix_(i_lat_obs, i_lon_obs)
                    elif not survey_file == 'craig.idl':
                        i_obs = np.ix_([i_trans_obs], i_lat_obs, i_lon_obs)
                    else:
                        i_obs = np.ix_(i_lon_obs, [i_trans_obs])
                    if comp_type == 'pv' and not survey in ['COBE-FIRAS', 'COBE-DIRBE', 'Planck']:
                        i_model.append(np.ix_(i_vel_model, i_lat_model, i_lon_model))
                    elif comp_type == 'integrated':
                        i_model.append(np.ix_(i_lat_model, i_lon_model))

                if len(map_list) == 1:
                    print('No models used in comparison...')
                    continue

                if debug:
                    print('number of maps and labels: {}, {}'.format(len(map_list), len(map_labels)))
                    print('observation')
                    print('  full map size: {}, indexed map size: {}'.format(map_list[0].size,
                                                                             map_list[0][i_obs].size))
                    if len(map_list) > 1:
                        print('model')
                        print('  full map size: {}, indexed map size: {}'.format(map_list[1].size,
                                                                                 map_list[1][i_model].size))

                fig, ax = plt.subplots(1, 1, figsize=(violin_spacing*len(map_list), 7))
                violin_positions = np.arange(0, violin_spacing*len(map_labels), violin_spacing)
                i_nan_obs = (map_list[0][i_obs] <= 0) | np.isnan(map_list[0][i_obs])
                i_nan_model = [((m[i_model[i]] <= 0) | np.isnan(m[i_model[i]]))
                        for i, m in enumerate(map_list[1:])]
                if log_comp == True:
                    # map_list[0] = np.nan_to_num(map_list[0], nan=0)
                    data = [np.log10(map_list[0][i_obs][~i_nan_obs].flatten()),
                            *[np.log10(m[i_model[i]][~i_nan_model[i]].flatten())
                              for i, m in enumerate(map_list[1:])]]
                    ax.violinplot(data, positions=violin_positions, widths=violin_width, 
                                  showextrema=showextrema)
                    if ylabel == '' and comp_type == 'integrated':
                        ax.set_ylabel(r'$log_{10} (\varpi) \ \left( K \frac{km}{s} \right)$', 
                                      fontsize=fontsize)
                    elif ylabel == '':
                        ax.set_ylabel(r'$log_{10} (T_\mathrm{B}) \ \left( K \right)$',
                                      fontsize=fontsize)
                    else:
                        ax.set_ylabel(ylabel)
                else:
                    data = [map_list[0][i_obs][~i_nan_obs].flatten(), 
                            *[m[i_model[i]][~i_nan_model[i]].flatten()
                              for i, m in enumerate(map_list[1:])]]
                    ax.violinplot(data, positions=violin_positions, widths=violin_width,
                                  showextrema=showextrema)
                    if comp_type == 'integrated':
                        ax.set_ylabel(r'$\varpi \ \left( K \frac{km}{s} \right)$', 
                                      fontsize=fontsize)
                    else:
                        ax.set_ylabel(r'$T_\mathrm{B} \ \left( K \right)$', 
                                      fontsize=fontsize)
                ax.set_xticks(violin_positions, labels=map_labels)
                ax.xaxis.set_tick_params(labelrotation=label_rotation, labelsize=label_fontsize)
                ax.yaxis.set_tick_params(labelsize=label_fontsize)
                #xticknames = ax.get_xticklabels()
                #plt.setp(xticknames, rotation=label_rotation, fontsize=label_fontsize)
                if len(ylim):
                    ax.set_ylim(ylim)
                if title == '':
                    ax.set_title(survey + ' -- ' + transition, fontsize=fontsize)
                elif isinstance(title, str):
                    ax.set_title(title, fontsize=fontsize)
                else:
                    pass
                plt.tight_layout()
                plt.subplots_adjust(bottom=0.15, left=0.15, right=1, top=0.85)
                if save_plot:
                    if output_file == '':
                        current_output_file = path.replace('observational_data', 'KT3_history') + \
                                              f'fit_results/{survey}/Plots/violinplot_{survey_file}-{transition}' + comp
                    else:
                        current_output_file = output_file
                    plt.savefig('{}.{}'.format(current_output_file, output_format),
                                format=output_format)
                    plt.close()
                else:
                    plt.show()

    return


def double_line_plot(obs_path='/mnt/hpc_backup/yanitski/projects/pdr/observational_data/MilkyWay/',
                     model_path='/mnt/hpc_backup/yanitski/projects/pdr/KT3_history/fit_results/MilkyWay/',
                     missions=[], transitions=[], survey_files=[], file_format='', model_param=[[]],
                     log_comp=True, lat=None, comp_type='integrated', bins=50, density=False, 
                     cmax=None, cmin=None, vmax=None, vmin=None, figsize=(), 
                     xlabel='', ylabel='', ylim=(), title='',
                     label_rotation=30, label_fontsize=16, fontsize=20,
                     save_plot=False, output_file='', output_format='png',
                     debug=False, verbose=False, **kwargs):
    '''
    Plot two lines against each other. Not setup well for comparing model grids.
    '''

    survey_paths = [[], []]  # list for paths to plotted data (in order observed, synthetic)
    survey_maps = [[], []]   # list for unprocessed (pre-calculation) data to plot (in order observed, synthetic)
    survey_i_lat_init = []   # list for initial observed lattitude indeces where there is an observation
    survey_i_vel_init = []   # list for initial observed velocity indeces where there is an observation
    survey_ix = [[], []]     # list of indeces used in comparison (in order observed, synthetic)
    survey_ratios = [[], []] # list of processed (post-calculation) data to plot (in order observed, synthetic)
    # list of labels for each ratio (in order observed, synthetic)
    survey_labels = [['/'.join(missions)], []]
    # list of available transitions in the data (in order observed, synthetic; used for debugging)
    survey_transitions = [[], []]
    survey_lons = [[], []]   # list of data longitude (in order observed, synthetic)
    survey_lats = [[], []]   # list of data lattitude (in order observed, synthetic)
    survey_vels = [[], []]   # list of data velocities (in order observed, synthetic; 0 for 2D maps)

    if ('COBE-FIRAS' in missions or 'COBE-DIRBE' in missions or 'Planck' in missions):
        comp_type = 'integrated'

    for i, survey in enumerate(missions):
        print(survey)
        survey_paths[0].append(obs_path + missions[i] + '/regridded/temp/' + survey_files[i])
        if '.fits' in survey_paths[0][i]:
            obs = fits.open(survey_paths[0][i])
            survey_transitions[0].append(np.asarray(obs[0].header['TRANSL'].split(', ')))
            i_trans_obs = survey_transitions[0][i] == transitions[i]
            if i_trans_obs.any() == False:
                print('Invalid transition {} for file {}'.format(transitions[i], survey_paths[i]))
                return
            survey_lons[0].append(np.linspace(obs[0].header['CRVAL1'] - obs[0].header['CDELT1'] * (
                                              obs[0].header['CRPIX1'] - 1),
                                              obs[0].header['CRVAL1'] + obs[0].header['CDELT1'] * (
                                              obs[0].header['NAXIS1'] - obs[0].header['CRPIX1']),
                                              num=obs[0].header['NAXIS1']))
            survey_lats[0].append(np.linspace(obs[0].header['CRVAL2'] - obs[0].header['CDELT2'] * (
                                              obs[0].header['CRPIX2'] - 1),
                                              obs[0].header['CRVAL2'] + obs[0].header['CDELT2'] * (
                                              obs[0].header['NAXIS2'] - obs[0].header['CRPIX2']),
                                              num=obs[0].header['NAXIS2']))
            if survey == 'COBE-FIRAS':
                survey_vels[0].append(np.asarray([0]))
                survey_maps[0].append(obs[0].data[i_trans_obs, :, :][0, :, :])
                survey_i_lat_init.append(survey_maps[0][i].any(1))
            elif survey == 'COBE-DIRBE':
                survey_vels[0].append(np.asarray([0]))
                survey_maps[0].append(obs[0].data[i_trans_obs, :, :][0, :, :])
                survey_i_lat_init.append(survey_maps[0][i].any(1))
            elif survey == 'Planck':
                survey_vels[0].append(np.asarray([0]))
                survey_maps[0].append(obs[0].data[:, :])
                survey_i_lat_init.append(survey_maps[0][i].any(1))
            elif ('COBE-FIRAS' in missions or 'Planck' in missions) or comp_type == 'integrated':
                survey_vels[0].append(np.linspace(obs[0].header['CRVAL3'] - obs[0].header['CDELT3'] * (
                                                  obs[0].header['CRPIX3'] - 1),
                                                  obs[0].header['CRVAL3'] + obs[0].header['CDELT3'] * (
                                                  obs[0].header['NAXIS3'] - obs[0].header['CRPIX3']),
                                                  num=obs[0].header['NAXIS3']))
                survey_maps[0].append(np.trapz(obs[0].data, survey_vels[0][i], axis=0))
                survey_i_lat_init.append(survey_maps[0][i].any(1))
            elif comp_type == 'pv':
                survey_vels[0].append(np.linspace(obs[0].header['CRVAL3'] - obs[0].header['CDELT3'] * (
                                                  obs[0].header['CRPIX3'] - 1),
                                                  obs[0].header['CRVAL3'] + obs[0].header['CDELT3'] * (
                                                  obs[0].header['NAXIS3'] - obs[0].header['CRPIX3']),
                                                  num=obs[0].header['NAXIS3']))
                survey_maps[0].append(obs[0].data)
                survey_i_lat_init.append(survey_maps[0][i].any(0).any(1))
                survey_i_vel_init.append(survey_maps[0][i].any(2).any(1))
            else:
                print('''Enter a valid string for `comp_type`: 'integrated' for a 2D comparison of maps '''
                      + '''and 'pv' for a 3D comparison of spectroscopic data.''')
                return
        elif '.idl' in survey_paths[0][i]:
            obs = readsav(survey_paths[0][i])
            survey_lons[0].append(obs['long'])
            survey_lats[0].append(np.asarray([0]))
            survey_transitions[0].append(transitions[i])
            i_trans_obs = cobe_idl_transitions == transitions[i]
            if i_trans_obs.any() == False:
                print()
                return
            survey_maps[0].append((obs['amplitude'][:, i_trans_obs] / (2.9979**3)
                                   * (cobe_idl_linfrq[i_trans_obs]**3) * 2 * 1.38 / 10 ** 8).reshape(0, -1))
            survey_i_lat_init.append(np.ones(1, dtype=bool))
        else:
            print('Invalid mission path {}.'.format(survey_paths[0][i]))
            return

    print(survey_maps[0][0].shape, survey_maps[0][1].shape)

    i_lon_obs = []
    i_lon_model = []
    i_lat_obs = []
    i_lat_model = []
    if ((isinstance(lat, int) or isinstance(lat, float)) and comp_type == 'pv' and
        not ('COBE-FIRAS' in missions or 'COBE-DIRBE' in missions or 'Planck' in missions)):
        lat_min = lat
        lat_max = lat
    else:
        lat_min = np.amin([survey_lats[0][i][survey_i_lat_init[i]] for i in range(2)])
        lat_max = np.amax([survey_lats[0][i][survey_i_lat_init[i]] for i in range(2)])
    if comp_type == 'pv':
        vel_min = np.amin([survey_vels[0][i][survey_i_vel_init[i]] for i in range(2)])
        vel_max = np.amax([survey_vels[0][i][survey_i_vel_init[i]] for i in range(2)])

    model_param_grid = np.meshgrid(*np.asarray(model_param, dtype=object))
    model_params = zip(*[model_param_grid[n].flatten() for n in range(len(model_param_grid))])

    for param in model_params:

        survey_labels[1].append(file_format.format(*param))
        model_dir = model_path + survey_labels[1][-1] + 'synthetic_intensity.fits'

        model = fits.open(model_dir)

        # Create arrays for the longitude and velocity axes of the synthetic observation if required
        if len(survey_lats[1]) == 0:
            survey_lats[1].append(
                np.linspace(model[1].header['CRVAL3']
                            - model[1].header['CDELT3'] * (model[1].header['CRPIX3'] - 0.5),
                            model[1].header['CRVAL3']
                            + model[1].header['CDELT3'] * (model[1].header['NAXIS3'] - model[1].header['CRPIX3'] - 0.5),
                            num=model[1].header['NAXIS3']) * 180 / np.pi)
            survey_lons[1].append(
                np.linspace(model[1].header['CRVAL2']
                            - model[1].header['CDELT2'] * (model[1].header['CRPIX2'] - 0.5),
                            model[1].header['CRVAL2']
                            + model[1].header['CDELT2'] * (model[1].header['NAXIS2'] - model[1].header['CRPIX2'] - 0.5),
                            num=model[1].header['NAXIS2']) * 180 / np.pi)
            # if not (mission == 'COGAL'):
            #     lon_model[lon_model<0] += 360
            survey_vels[1].append(
                np.linspace(model[1].header['CRVAL4']
                            - model[1].header['CDELT4'] * (model[1].header['CRPIX4'] - 0.5),
                            model[1].header['CRVAL4']
                            + model[1].header['CDELT4'] * (model[1].header['CRPIX4'] - 0.5),
                            num=model[1].header['NAXIS4']))

        survey_maps[1].append([])

        for i, transition in enumerate(transitions):
            if transition == 'Dust':
                survey_transitions[1].append(np.asarray(model[2].header['DUST'].split(', ')))
                # survey_maps[1][-1].append(deepcopy(model[2].data[0, :, :, 0]))
                survey_maps[1][-1].append(deepcopy(model[2].data[:, :, 0]))
            else:
                survey_transitions[1].append(np.asarray(model[1].header['SPECIES'].split(', ')))
                i_transition = np.where(survey_transitions[1][-1] == transition)[0]
                if i_transition.size == 1 and (comp_type == 'integrated'
                                               or not ('COBE-FIRAS' in missions or 'Planck' in missions)):
                    survey_maps[1][-1].append(np.trapz(model[1].data[:, :, :, i_transition][:, :, :, 0]
                                                       # - model[2].data[:, :, :, 0], survey_vels[1][0], axis=0))
                                                       - model[2].data[:, :, 0], survey_vels[1][0], axis=0))
                elif i_transition.size == 1:
                    survey_maps[1][-1].append(model[1].data[:, :, :, i_transition][:, :, :, 0]
                                              # - model[2].data[:, :, :, 0])
                                              - model[2].data[:, :, 0])
                else:
                    print('Transition {} not found in model'.format(transition))
                    return

            if len(survey_ix[0]) < 2:
                i_lon_obs = np.linspace(0, survey_lons[0][i].size-1, num=survey_lons[0][i].size, dtype=int)
                i_lon_model = np.linspace(0, survey_lons[1][0].size-1, num=survey_lons[1][0].size, dtype=int)
                i_lat_model = np.where((survey_lats[1][0] >= lat_min)
                                       & (survey_lats[1][0] <= lat_max))[0]
                i_lat_obs = np.where((survey_lats[0][i] >= survey_lats[1][0][i_lat_model].min())
                                     & (survey_lats[0][i] <= survey_lats[1][0][i_lat_model].max()))[0]
                if comp_type == 'integrated': #compare 2D maps
                    survey_ix[0].append(np.ix_(i_lat_obs, i_lon_obs))
                    survey_ix[1].append(np.ix_(i_lat_model, i_lon_model))
                elif comp_type == 'pv': #compare 3D maps
                    i_vel_model = np.where((survey_vels[1][0] >= vel_min)
                                           & (survey_vels[1][0] <= vel_max))[0]
                    i_vel_obs = np.where((survey_vels[0][i] >= survey_vels[1][0][i_vel_model].min())
                                         & (survey_vels[0][i] <= survey_vels[1][0][i_vel_model].max()))[0]
                    survey_ix[0].append(np.ix_(i_vel_obs, i_lat_obs, i_lon_obs))
                    survey_ix[1].append(np.ix_(i_vel_model, i_lat_model, i_lon_model))

        if len(survey_ratios[0]) == 0:
            survey_lines[0].append((survey_maps[0][0][survey_ix[0][0]].flatten(),
                                    survey_maps[0][1][survey_ix[0][1]].flatten()))
        # for i in range(len(survey_maps[1])):
        survey_lines[1].append((survey_maps[1][-1][0][survey_ix[1][0]].flatten(),
                                survey_maps[1][-1][1][survey_ix[1][1]].flatten()))

    if len(figsize) == 0:
        figsize = (1*(len(survey_ratios[0])+len(survey_ratios[1])), 10)
    #fig, ax = plt.subplots(1, 1, figsize=figsize)
    labels = [*survey_labels[0], *survey_labels[1]]
    if log_comp == True:
        comp = '_logT'
        i_nan_obs = (survey_ratios[0][0] <= 0) | np.isnan(survey_ratios[0][0])
        data = [np.log10(survey_lines[0][0][~i_nan_obs]),
                *[np.log10(survey_lines[1][i]) for i in range(len(survey_ratios[1]))]]
        #ax.boxplot(data, labels=labels, widths=boxwidth, notch=notch)
        if ylabel == '':
            xlabel = r'$log_{10} ' + r'\left( \varpi_{' + r'{}'.format(transitions[0]) \
                     + r'} \right)$'
            ylabel = r'$log_{10} ' + r'\left( \varpi_{' + r'{}'.format(transitions[1]) \
                     + r'} \right)$'
        #else:
        #    ax.set_xlabel(xlabel)
        #    ax.set_ylabel(ylabel)
    else:
        comp = ''
        data = [survey_lines[0][0], *survey_lines[1]]
        #ax.boxplot(data, labels=labels, widths=boxwidth, notch=notch)
        if ylabel == '':
            xlabel = r'$\varpi_{' + r'{}'.format(transitions[0]) + r'}$'
            ylabel = r'%\varpi_{' + r'{}'.format(transitions[1]) + r'}$'
        #else:
        #    ax.set_xlabel(xlabel)
        #    ax.set_ylabel(ylabel)
    for i, lines in enumerate(data):
        fig, ax = plt.subplots(1, 1, figsize=figsize)
        plt.hist2d(lines[0], lines[1], bins=bins, density=density, vmin=vmin, vmax=vmax, cmin=cmin, cmax=cmax)
        ax.set_xlabel(xlabel)
        ax.set_ylabel(ylabel)
        if len(ylim):
            ax.set_ylim(ylim)
        #xticknames = ax.get_xticklabels()
        #plt.setp(xticknames, rotation=label_rotation, fontsize=label_fontsize)
        if title == '':
            ax.set_title('{} / {}'.format(*transitions), fontsize=fontsize)
        else:
            ax.set_title(title)
        plt.tight_layout()
        if save_plot:
            if output_file == '':
                #filename = '_'.join(file_format.replace('/', '').split('_')[1:]).replace('{}', '_')
                output_args = (survey_paths[0][0].split('/')[-1].split('.')[0], transitions[0],
                               survey_paths[0][1].split('/')[-1].split('.')[0], transitions[1], labels[i])
                current_output_file = 'line_comp_{}-{}_{}-{}==>{}'.format(*output_args) + comp
            else:
                current_output_file = output_file + '==>{}'.format(labels[i])
            plt.savefig(model_path + 'fit_results/Plots/{}.{}'.format(current_output_file, output_format),
                        format=output_format)
        else:
            plt.show()
        plt.close()

    return



def plot_comparison(path='/mnt/hpc_backup/yanitski/projects/pdr/KT3_history/MilkyWay/fit_results/', file_format='',
                    missions=[], model_param=[[]], transitions=[], i_x=0, i_y=1, stat_lim=1e100, 
                    comp_type='pv', log_comp=True, likelihood=True, log=False, normalise=False,
                    contour=True, levels=10, fraction=0.1, cb_aspect=20, cmap='viridis',
                    xlabel='', ylabel='', xticks=[], yticks=[], xrot=0, yrot=0, supxlabel='', supylabel='', 
                    clabel='', clabel_xa=0.98, clabel_ha='left',
                    title='', offset_idx=(1.0, 1.0), fontsize=24, labelsize=16, ax_aspect=1.5, 
                    pad=1.0, pad_left=0, pad_right=0.125, pad_bottom=0, pad_top=0.12, wspace=0, hspace=0,
                    figsize=None, save_plot=False, output_file='', prefix='', output_format='png', transparent=False,
                    verbose=False, debug=False, **kwargs):
    '''
    Plot heatmaps of a test statistic to exaimine over model grid.
    '''

    # Check that the missions are specified properly.
    if missions == '' or missions == None or missions == []:
        missions = os.listdir(path)
    elif isinstance(missions, str):
        missions = [missions]
    elif isinstance(missions, list):
        # Use both COBE instruments if simply 'COBE' is specified
        if 'COBE' in missions:
            if 'COBE-FIRAS' in missions: missions.remove('COBE-FIRAS')
            if 'COBE-DIRBE' in missions: missions.remove('COBE-DIRBE')
            missions.append('COBE-FIRAS')
            missions.append('COBE-DIRBE')
    else:
        print('Please specify a list of missions to compare the models.')
        return

    transitions_all = copy(transitions)

    if file_format == '' or model_param == [[]]:
        print('Please specify both file_format and model_param.')
        return

    if log_comp:
        comp = comp_type + '_logT'
    else:
        comp = comp_type

    if not clabel:
        if likelihood:
            clabel = r'$log_{10}(\mathcal{L})$'
        else:
            clabel = r'$\chi^2_\mathrm{min}$'

    dimensions = len(model_param)
    naxis = np.zeros((dimensions), dtype=bool)
    naxis[[i_x, i_y]] = True
    # print(model_param)
    lenaxis = np.asarray([len(arr) for arr in model_param])

    model_param_grid = np.meshgrid(*np.asarray(model_param, dtype=object)[naxis])
    if i_x < i_y:
        x_grid = model_param_grid[0]
        y_grid = model_param_grid[1]
    else:
        x_grid = model_param_grid[1]
        y_grid = model_param_grid[0]

    if naxis.size == 2:
        sub_params = zip([0], [0])
    elif naxis.size == 3:
        sub_params = zip(np.asarray(model_param, dtype=object)[~naxis])
    else:
        sub_params = zip(*tuple(p.flatten() for p in np.meshgrid(*np.asarray(model_param, dtype=object)[~naxis])))

    # Detemine the size of the subplot grid, for now allow a maximum of 2 additional parameters.
    if naxis.size == 2:
        sub_x, sub_y = (1, 1)
    elif naxis.size == 3:
        sub_x = lenaxis[~naxis]
        sub_y = 1
    elif naxis.size == 4:
        sub_x, sub_y = lenaxis[~naxis]
    else:
        print('Too many dimensions in grid.')
        return

    # Initialise the tick labels
    if len(xticks):
        xtick_base = np.linspace(0, lenaxis[i_x]-1, num=len(xticks))
        xtick_labels = xticks
    else:
        xtick_base = np.linspace(0, lenaxis[i_x]-1, num=lenaxis[i_x])
        xtick_labels = [str(par) for par in model_param[i_x]]
    if len(yticks):
        ytick_base = np.linspace(0, lenaxis[i_y]-1, num=len(yticks))
        ytick_labels = yticks
    else:
        ytick_base = np.linspace(0, lenaxis[i_y]-1, num=lenaxis[i_y])
        ytick_labels = [str(par) for par in model_param[i_y]]
    if contour:
        pass
    else:
        xtick_base += 0.5
        ytick_base += 0.5

    # Initialise likelihood and figure for plot analysing all missions
    loglikelihood_overall = np.zeros((*x_grid.shape, sub_y, sub_x))
    overall_dof = 1#0
    fig_overall, axes_overall = plt.subplots(sub_y, sub_x, figsize=figsize)
    if not isinstance(axes_overall, np.ndarray):
        axes_overall = np.asarray([[axes_overall]])
    elif axes_overall.ndim == 1:
        axes_overall.resize(-1, 1)

    for survey in missions:

        if survey == 'Plots':
            continue

        # if verbose:
        print('\n  {}\n'.format(survey))

        if (survey + '_files') in kwargs.keys():
            if isinstance(kwargs[survey + '_files'], list):
                survey_files = kwargs[survey + '_files']
            elif isinstance(kwargs[survey + '_files'], str):
                survey_files = [kwargs[survey + '_files']]
            else:
                print('Incorrect files specified for survey {}'.format(survey))
                continue
        else:
            survey_files = os.listdir(path + survey + '/')

        if 'Plots' in survey_files:
            survey_files.remove('Plots')

        if not 'Plots' in os.listdir(path + survey + '/'):
            os.mkdir(path + survey + '/Plots/')

        file_transitions = []
        file_transition_plots = []
        file_transition_log_likelihood = []

        if debug:
            print('survey files: {}\n'.format(survey_files))

        if not figsize:
            subgrid_aspect = sub_x/sub_y
            figsize = (subgrid_aspect*15, 10)
        fig, axes = plt.subplots(sub_y, sub_x, figsize=figsize)
        if not isinstance(axes, np.ndarray):
            axes = np.asarray([[axes]])
        elif axes.ndim == 1:
            axes.resize(-1, 1)

        for param in deepcopy(sub_params):

            # Calculate likelihood
            # --------------------

            if debug:
                print(param)

            survey_dof = 1#0
            file_dof = [1]#[0] * len(survey_files)

            log_likelihood = np.zeros(x_grid.shape)

            for f,survey_file in enumerate(survey_files):

                if '.png' in survey_file or 'Plots' in survey_file or '.npy' in survey_file \
                        or '.npz' in survey_file:
                    continue
                if not survey_file in os.listdir(path + survey + '/'):
                    print('File {} not available for survey {}'.format(survey_file, survey))
                    continue

                if verbose:
                    print(survey_file)

                if len(transitions_all):
                    transitions = copy(transitions_all)
                    transitions_skipped = []
                else:
                    transitions = os.listdir(path + survey + '/' + survey_file + '/')
                    for _ in copy(transitions):
                        if '.npy' in _ or '.npz' in _:
                            transitions.remove(_)
                    #if f'{output_file}_comparison.npy' in transitions:
                    #    transitions.remove(f'{output_file}_comparison.npy')

                    transitions_skipped = []
                    for t in copy(transitions):
                        if debug:
                            print(t)
                        if (len(os.listdir(path + survey + '/' + survey_file + '/' + t + '/')) == 0):
                            transitions.remove(t)
                            transitions_skipped.append(t)
                        if debug:
                            print('  transitions compared: {}\n  transitions skipped: {}'.format(transitions,
                                                                                                 transitions_skipped))
                file_transitions.append(copy(transitions))
                #survey_dof += len(transitions)
                #file_dof[f] = len(transitions)
                if len(transitions_skipped) and verbose:
                    print('  transitions {} not available.'.format(', '.join(transitions_skipped)))

                transition_plots = []
                transition_log_likelihood = []
                for _ in transitions:
                    t_fig, t_axes = plt.subplots(sub_y, sub_x, figsize=figsize)
                    if not isinstance(t_axes, np.ndarray):
                        t_axes = np.asarray([[t_axes]])
                    elif t_axes.ndim == 1:
                        t_axes.resize(-1, 1)
                    transition_plots.append((copy(t_fig), copy(t_axes)))
                    transition_log_likelihood.append(np.zeros(x_grid.shape))
                    plt.close(t_fig)

                for i in range(log_likelihood.shape[0]):
                    for j in range(log_likelihood.shape[1]):

                        file_format_param = np.zeros(naxis.size, dtype=object)
                        file_format_param[i_x] = x_grid[i, j]
                        file_format_param[i_y] = y_grid[i, j]

                        if verbose:
                            print(file_format.format(*file_format_param) \
                                  + '_{}_chi2.npy'.format(comp))

                        if naxis.size > 2:
                            file_format_param[~naxis] = param
                        if likelihood:
                            filename = file_format.format(*file_format_param) \
                                       + '_{}_loglikelihood.npy'.format(comp_type)
                        else:
                            filename = file_format.format(*file_format_param) \
                                       + '_{}_chi2.npy'.format(comp)
                            # filename = file_format.format(x_grid[i, j], y_grid[i, j]) \
                            #            + '_{}_loglikelihood.npy'.format(comp_type)
                        param_log_likelihood = [np.load(path + survey + '/' + survey_file + '/' +
                                                            t + '/' + filename) for t in transitions]
                        # filter out unnecessary pixels
                        for _ in range(len(param_log_likelihood)):
                            param_log_likelihood[_][param_log_likelihood[_] == 0] = np.nan
                            param_log_likelihood[_][param_log_likelihood[_] >= stat_lim] = np.nan
                        if likelihood:
                            log_likelihood[i, j] = log_likelihood[i, j] + np.nansum(param_log_likelihood)
                        else:
                            log_likelihood[i, j] = log_likelihood[i, j] + np.nanprod(param_log_likelihood)
                        for t in range(len(transitions)):
                            transition_log_likelihood[t][i, j] = transition_log_likelihood[t][i, j] \
                                                                 + np.nansum(param_log_likelihood[t])

                file_transition_plots.append(copy(transition_plots))
                file_transition_log_likelihood.append(deepcopy(transition_log_likelihood))
                np.savez(f'{path}{survey}/{survey_file}/{prefix}{output_file}_{comp}_comparison.npz', trans=file_transitions, chi2=file_transition_log_likelihood)

            if normalise:
                if (log_likelihood<0).all():
                    log_likelihood = log_likelihood.max() / log_likelihood
                else:
                    log_likelihood = log_likelihood / log_likelihood.max()
                for f in range(len(survey_files)):
                    for t in range(len(file_transitions[f])):
                        if (log_likelihood<0).all():
                            file_transition_log_likelihood[f][t] = file_transition_log_likelihood[f][t].max() \
                                                                   / file_transition_log_likelihood[f][t]
                        else:
                            file_transition_log_likelihood[f][t] = file_transition_log_likelihood[f][t] \
                                                                   / file_transition_log_likelihood[f][t].max()
            if log:
                if (log_likelihood<0).all():
                    log_likelihood = np.log10(-log_likelihood)
                else:
                    log_likelihood = np.log10(log_likelihood)
                for f in range(len(survey_files)):
                    for t in range(len(file_transitions[f])):
                        if (file_transition_log_likelihood[f][t]<0).all():
                            file_transition_log_likelihood[f][t] = np.log10(-file_transition_log_likelihood[f][t])
                        else:
                            file_transition_log_likelihood[f][t] = np.log10(file_transition_log_likelihood[f][t])

            # Plot subplots
            # -------------

            if len(param) == 1:
                sub_indeces = (0, np.asarray(model_param, dtype=object)[~naxis].index(param))
                x_param = param[0]
                y_param = ''
            elif axes.size > 1:
                sub_indeces = tuple(arr.index(param[i]) for i,arr in
                                    enumerate(np.asarray(model_param, dtype=object)[~naxis]))[::-1]
                x_param = param[1]
                y_param = param[0]
            else:
                sub_indeces = (0, 0)
                x_param = ''
                y_param = ''

            # Add likelihood in overall array
            loglikelihood_overall[:, :, sub_indeces[0], sub_indeces[1]] \
                = loglikelihood_overall[:, :, sub_indeces[0], sub_indeces[1]] + log_likelihood
            
            # Update overall degrees of freedom
            #overall_dof += np.sum(file_dof)

            # Minimum chi^2
            minval = (log_likelihood/np.sum(file_dof)).min()
            # minval = np.around(minval/(10**np.floor(np.log10(minval))), 3) * 10**np.floor(np.log10(minval))

            if contour:
                cm = axes[sub_indeces].contourf((log_likelihood/np.sum(file_dof)-minval)/minval*100, antialiased=True, levels=levels, cmap=cmap)
            else:
                cm = axes[sub_indeces].imshow(log_likelihood/np.sum(file_dof),
                                              extent=[0, lenaxis[i_x], 0, lenaxis[i_y]], cmap=cmap)
            axes[sub_indeces].set_aspect(lenaxis[i_x]/ax_aspect/lenaxis[i_y])
            cb = fig.colorbar(cm, ax=axes[sub_indeces], fraction=fraction, aspect=cb_aspect)
            # cb.ax.ticklabel_format(useOffset=minval, style='plain', useMathText=True)
            cb.ax.ticklabel_format(useOffset=False, style='plain', useMathText=True)
            #cb.ax.yaxis.offsetText.set_text('')
            #offsetx, offsety = cb.ax.yaxis.offsetText.get_position()
            # axes[sub_indeces].text(1.1, 1.0, f'{cb.ax.yaxis.major.formatter.offset:1.3e}', transform=axes[sub_indeces].transAxes, fontsize=labelsize, ha='left', va='bottom')
            axes[sub_indeces].text(*offset_idx, f'{minval:1.3e} +', transform=axes[sub_indeces].transAxes, fontsize=labelsize, ha='left', va='bottom')
            #cb.ax.yaxis.major.formatter.set_useOffset(False)
            cb.ax.tick_params(labelsize=labelsize)
            #cb.ax.yaxis.offsetText.set_fontsize(labelsize)
            #cb.ax.yaxis.offsetText.set_ha('left')
            #cb.ax.yaxis.offsetText.set_va('bottom')
            for c in cm.collections:
                c.set_edgecolor("face")
                c.set_linewidth(1e-16)

            axes[sub_indeces].set_xticks(xtick_base)
            axes[sub_indeces].set_xticklabels(xtick_labels, rotation=xrot)
            axes[sub_indeces].set_yticks(ytick_base)
            axes[sub_indeces].set_yticklabels(ytick_labels, rotation=yrot)
            # if contour:
                # axes[sub_indeces].set_xticks(np.arange(lenaxis[i_x]))
                # axes[sub_indeces].set_xticklabels([str(par) for par in model_param[i_x]])
                # axes[sub_indeces].set_yticks(np.arange(lenaxis[i_y]))
                # axes[sub_indeces].set_yticklabels([str(par) for par in model_param[i_y]])
            # else:
                # axes[sub_indeces].set_xticks(np.arange(lenaxis[i_x])+0.5)
                # axes[sub_indeces].set_xticklabels([str(par) for par in model_param[i_x]])
                # axes[sub_indeces].set_yticks(np.arange(lenaxis[i_y])+0.5)
                # axes[sub_indeces].set_yticklabels([str(par) for par in model_param[i_y][::-1]])

            axes[sub_indeces].set_xlabel(xlabel, fontsize=fontsize)
            axes[sub_indeces].set_ylabel(ylabel, fontsize=fontsize)
            axes[sub_indeces].tick_params(labelsize=labelsize)
            axes[sub_indeces].tick_params(axis='x', pad=10)
            axes[sub_indeces].tick_params(axis='y', pad=10)
            if axes.size > 1:
                axes[sub_indeces].set_title('{} {}, {} {}'.format(supylabel, y_param, supxlabel, x_param),
                                            fontsize=fontsize-4)
            # cb.ax.set_ylabel(clabel, fontsize=fontsize-4)

            if debug:
                print('File transitions: {}'.format(file_transitions))

            # For parameter comparisons split by transition and survey
            for f in range(len(survey_files)):
                for t in range(len(file_transitions[f])):
                    
                    # Minimum chi^2
                    minval = file_transition_log_likelihood[f][t].min()
                    # minval = np.around(minval/(10**np.floor(np.log10(minval))), 3) * 10**np.floor(np.log10(minval))

                    if contour:
                        cm = file_transition_plots[f][t][1][sub_indeces].contourf((file_transition_log_likelihood[f][t]-minval)/minval*100,
                                                                                  antialiased=True, levels=levels, cmap=cmap)
                    else:
                        cm = file_transition_plots[f][t][1][sub_indeces].imshow(file_transition_log_likelihood[f][t],
                                                                                extent=[0, lenaxis[i_x],
                                                                                        0, lenaxis[i_y]],
                                                                                cmap=cmap)
                    file_transition_plots[f][t][1][sub_indeces].set_aspect(lenaxis[i_x]/ax_aspect/lenaxis[i_y])
                    cb = file_transition_plots[f][t][0].colorbar(cm, ax=file_transition_plots[f][t][1][sub_indeces],
                                                                 fraction=fraction, aspect=cb_aspect)
                    # cb.ax.ticklabel_format(useOffset=minval, style='plain', useMathText=True)
                    cb.ax.ticklabel_format(useOffset=False, style='plain', useMathText=True)
                    #cb.ax.yaxis.offsetText.set_text('')
                    #offsetx, offsety = cb.ax.yaxis.offsetText.get_position()
                    # cb.ax.text(-1.5, 1.0, f'{cb.ax.yaxis.major.formatter.offset:1.3e}', fontsize=labelsize, ha='center', va='bottom')
                    file_transition_plots[f][t][1][sub_indeces].text(*offset_idx, f'{minval:1.3e} +', transform=file_transition_plots[f][t][1][sub_indeces].transAxes, 
                                                                     fontsize=labelsize, ha='left', va='bottom')
                    #cb.ax.yaxis.major.formatter.set_useOffset(False)
                    cb.ax.tick_params(labelsize=labelsize)
                    #cb.ax.yaxis.offsetText.set_fontsize(labelsize)
                    #cb.ax.yaxis.offsetText.set_ha('left')
                    #cb.ax.yaxis.offsetText.set_va('bottom')
                    for c in cm.collections:
                        c.set_edgecolor("face")
                        c.set_linewidth(1e-16)

                    file_transition_plots[f][t][1][sub_indeces].set_xticks(xtick_base)
                    file_transition_plots[f][t][1][sub_indeces].set_xticklabels(xtick_labels, rotation=xrot)
                    file_transition_plots[f][t][1][sub_indeces].set_yticks(ytick_base)
                    file_transition_plots[f][t][1][sub_indeces].set_yticklabels(ytick_labels, rotation=yrot)
                    # if contour:
                    #     file_transition_plots[f][t][1][sub_indeces].set_xticks(np.arange(lenaxis[i_x]))
                    #     file_transition_plots[f][t][1][sub_indeces].set_xticklabels([str(param) for
                    #                                                                  param in model_param[i_x]])
                    #     file_transition_plots[f][t][1][sub_indeces].set_yticks(np.arange(lenaxis[i_y]))
                    #     file_transition_plots[f][t][1][sub_indeces].set_yticklabels([str(param) for
                    #                                                                  param in model_param[i_y]])
                    # else:
                    #     file_transition_plots[f][t][1][sub_indeces].set_xticks(np.arange(lenaxis[i_x])+0.5)
                    #     file_transition_plots[f][t][1][sub_indeces].set_xticklabels([str(param) for
                    #                                                                  param in model_param[i_x]])
                    #     file_transition_plots[f][t][1][sub_indeces].set_yticks(np.arange(lenaxis[i_y])+0.5)
                    #     file_transition_plots[f][t][1][sub_indeces].set_yticklabels([str(param) for
                    #                                                                  param in model_param[i_y][::-1]])

                    file_transition_plots[f][t][1][sub_indeces].set_xlabel(xlabel, fontsize=fontsize)
                    file_transition_plots[f][t][1][sub_indeces].set_ylabel(ylabel, fontsize=fontsize)
                    file_transition_plots[f][t][1][sub_indeces].tick_params(labelsize=labelsize)
                    file_transition_plots[f][t][1][sub_indeces].tick_params(axis='x', pad=10)
                    file_transition_plots[f][t][1][sub_indeces].tick_params(axis='y', pad=10)
                    if file_transition_plots[f][t][1].size > 1:
                        file_transition_plots[f][t][1][sub_indeces].set_title('{} {}, {} {}'.format(supylabel, y_param,
                                                                                                    supxlabel, x_param),
                                                                              fontsize=fontsize-4)

        if title == '':
            if likelihood:
                suptitle = survey + ' likelihood'
            else:
                suptitle = survey + r' $\chi^2$'
            if normalise:
                suptitle += ', normalised'
            if log:
                suptitle += ', logged'
        else:
            suptitle = copy(title)

        if isinstance(suptitle, str):
            fig.suptitle(suptitle, fontsize=fontsize)
        if clabel == '' or clabel == None:
            clabel = 'Value'
        # fig.supylabel(clabel, x=clabel_xa, ha=clabel_ha, fontsize=fontsize)
        fig.text(clabel_xa, 0.5, clabel, rotation=90, va='center', ha=clabel_ha, fontsize=fontsize)

        if suptitle:
            fig_top = 1 - pad*pad_top
        if clabel:
            fig_right = 1 - pad*pad_right
        fig.subplots_adjust(left=pad*pad_left, right=fig_right, bottom=pad*pad_bottom, top=fig_top,
                            wspace=wspace, hspace=hspace)
        # fig.tight_layout(pad=pad)

        plt.figure(fig)
        if save_plot:
            if output_file == None or output_file == '':
                output_file = file_format.replace('{}', '_')
            if likelihood:
                plt.savefig(path + survey + '/Plots/{}_{}_loglikelihood.{}'.format(output_file, comp, output_format),
                            format=output_format, transparent=transparent)
            else:
                plt.savefig(path + survey + '/Plots/{}_{}_chi2.{}'.format(output_file, comp, output_format),
                            format=output_format, transparent=transparent)
        else:
            plt.show()

        plt.close(fig)

        # For parameter comparisons split by transition and survey
        for f,file in enumerate(survey_files):
            for t in range(len(file_transitions[f])):

                if title == '':
                    if likelihood:
                        suptitle = survey + ' ' + file_transitions[f][t] + ' likelihood'
                    else:
                        suptitle = survey + ' ' + file_transitions[f][t] + r' $\chi^2$'
                    if normalise:
                        suptitle += ', normalised'
                    if log:
                        suptitle += ', logged'
                else:
                    suptitle = copy(title)

                if isinstance(suptitle, str):
                    file_transition_plots[f][t][0].suptitle(suptitle, fontsize=fontsize)
                if clabel == '' or clabel == None:
                    clabel = 'Value'
                # file_transition_plots[f][t][0].supylabel(clabel, x=clabel_xa, ha=clabel_ha, fontsize=fontsize)
                file_transition_plots[f][t][0].text(clabel_xa, 0.5, clabel, rotation=90, ha=clabel_ha, va='center', fontsize=fontsize)

                if suptitle:
                    fig_top = 1 - pad*pad_top
                if clabel:
                    fig_right = 1 - pad*pad_right
                # file_transition_plots[f][t][0].tick_params(labelsize=labelsize)
                file_transition_plots[f][t][0].subplots_adjust(left=pad*pad_left, right=fig_right,
                                                               bottom=pad*pad_bottom, top=fig_top,
                                                               wspace=wspace, hspace=hspace)
                # file_transition_plots[f][t][0].tight_layout()

                plt.figure(file_transition_plots[f][t][0])
                if save_plot:
                    if output_file == None or output_file == '':
                        output_file = file_format.replace('{}', '_')
                    if likelihood:
                        plt.savefig(path + survey + '/Plots/{}_{}-{}_{}_loglikelihood.{}'
                                    .format(output_file, file, file_transitions[f][t].replace(' ', '-'),
                                            comp, output_format),
                                    format=output_format, transparent=transparent)
                    else:
                        plt.savefig(path + survey + '/Plots/{}_{}-{}_{}_chi2.{}'
                                    .format(output_file, file, file_transitions[f][t].replace(' ', '-'),
                                            comp, output_format),
                                    format=output_format, transparent=transparent)
                else:
                    plt.show()

                plt.close(file_transition_plots[f][t][0])

    # For parameter comparisons including all transitions/surveys
    for param in deepcopy(sub_params):

        if len(param) == 1:
            sub_indeces = (0, np.asarray(model_param, dtype=object)[~naxis].index(param))
            x_param = param[0]
            y_param = ''
        elif axes.size > 1:
            sub_indeces = tuple(
                arr.index(param[i]) for i, arr in enumerate(np.asarray(model_param, dtype=object)[~naxis]))[::-1]
            x_param = param[1]
            y_param = param[0]
        else:
            sub_indeces = (0, 0)
            x_param = ''
            y_param = ''

        # minimum chi^2
        minval = (loglikelihood_overall[:, :, sub_indeces[0], sub_indeces[1]]/overall_dof).min()
        # minval = np.around(minval/(10**np.floor(np.log10(minval))), 3) * 10**np.floor(np.log10(minval))
        
        if contour:
            cm = axes_overall[sub_indeces].contourf((loglikelihood_overall[:, :, sub_indeces[0], sub_indeces[1]]
                                                     / overall_dof - minval)/minval*100, antialiased=True, levels=levels, cmap=cmap)
        else:
            cm = axes_overall[sub_indeces].imshow(loglikelihood_overall[:, :, sub_indeces[0], sub_indeces[1]]
                                                  / overall_dof, extent=[0, lenaxis[i_x], 0, lenaxis[i_y]], cmap=cmap)
        axes_overall[sub_indeces].set_aspect(lenaxis[i_x]/ax_aspect/lenaxis[i_y])
        cb = fig_overall.colorbar(cm, ax=axes_overall[sub_indeces], fraction=fraction, aspect=cb_aspect)
        cb.ax.ticklabel_format(useOffset=False, style='plain', useMathText=True)
        #cb.ax.yaxis.offsetText.set_text('')
        #offsetx, offsety = cb.ax.yaxis.offsetText.get_position()
        # cb.ax.text(-1.5, 1.0, f'{cb.ax.yaxis.major.formatter.offset:1.3e}', fontsize=labelsize, ha='center', va='bottom')
        axes_overall[sub_indeces].text(*offset_idx, f'{minval:1.3e} +', fontsize=labelsize, ha='center', va='bottom')
        #cb.ax.yaxis.major.formatter.set_useOffset(False)
        cb.ax.tick_params(labelsize=labelsize)
        #cb.ax.yaxis.offsetText.set_fontsize(labelsize)
        #cb.ax.yaxis.offsetText.set_ha('left')
        #cb.ax.yaxis.offsetText.set_va('bottom')
        for c in cm.collections:
            c.set_edgecolor("face")
            c.set_linewidth(1e-16)

        axes_overall[sub_indeces].set_xticks(xtick_base)
        axes_overall[sub_indeces].set_xticklabels(xtick_labels, rotation=xrot)
        axes_overall[sub_indeces].set_yticks(ytick_base)
        axes_overall[sub_indeces].set_yticklabels(ytick_labels, rotation=yrot)
        # if contour:
        #     axes_overall[sub_indeces].set_xticks(np.arange(lenaxis[i_x]))
        #     axes_overall[sub_indeces].set_xticklabels([str(param) for param in model_param[i_x]])
        #     axes_overall[sub_indeces].set_yticks(np.arange(lenaxis[i_y]))
        #     axes_overall[sub_indeces].set_yticklabels([str(param) for param in model_param[i_y]])
        # else:
        #     axes_overall[sub_indeces].set_xticks(np.arange(lenaxis[i_x]) + 0.5)
        #     axes_overall[sub_indeces].set_xticklabels([str(param) for param in model_param[i_x]])
        #     axes_overall[sub_indeces].set_yticks(np.arange(lenaxis[i_y]) + 0.5)
        #     axes_overall[sub_indeces].set_yticklabels([str(param) for param in model_param[i_y][::-1]])

        axes_overall[sub_indeces].set_xlabel(xlabel, fontsize=fontsize)
        axes_overall[sub_indeces].set_ylabel(ylabel, fontsize=fontsize)
        axes_overall[sub_indeces].tick_params(labelsize=labelsize)
        axes_overall[sub_indeces].tick_params(axis='x', pad=10)
        axes_overall[sub_indeces].tick_params(axis='y', pad=10)
        if axes_overall.size > 1:
            axes_overall[sub_indeces].set_title('{} {}, {} {}'.format(supylabel, y_param, supxlabel, x_param),
                                                fontsize=fontsize - 4)

    if title == '':
        if likelihood:
            suptitle = 'Total (all lines/missions) likelihood'
        else:
            suptitle = 'Total (all lines/missions) ' + r'$\chi^2$'
        if normalise:
            suptitle += ', normalised'
        if log:
            suptitle += ', logged'
    else:
        suptitle = copy(title)

    if isinstance(suptitle, str):
        fig_overall.suptitle(suptitle, fontsize=fontsize)
    if clabel == '' or clabel == None:
        clabel = 'Value'
    # fig_overall.supylabel(clabel, x=clabel_xa, ha=clabel_ha, fontsize=fontsize)
    fig_overall.text(clabel_xa, 0.5, clabel, rotation=90, va='center', ha=clabel_ha, fontsize=fontsize)

    if suptitle:
        fig_top = 1 - pad * pad_top
    if clabel:
        fig_right = 1 - pad * pad_right
    fig_overall.subplots_adjust(left=pad * pad_left, right=fig_right, bottom=pad * pad_bottom, top=fig_top,
                                wspace=wspace, hspace=hspace)
    # fig_overall.tight_layout(pad=pad)

    plt.figure(fig_overall)
    if save_plot:
        if not 'Plots' in os.listdir(path):
            os.mkdir(path + '/Plots/')
        if output_file == None or output_file == '':
            output_file = file_format.replace('{}', '_')
        if likelihood:
            plt.savefig(path + '/Plots/{}_{}_loglikelihood.{}'.format(output_file, comp, output_format),
                        format=output_format, transparent=transparent)
        else:
            plt.savefig(path + '/Plots/{}_{}_chi2.{}'.format(output_file, comp, output_format),
                        format=output_format, transparent=transparent)
    else:
        plt.show()

    plt.close('all')

    return

