import cygrid
import numpy as np
import healpy as hp
import astrokit
from astropy.io import fits
from astropy.wcs import WCS
from astropy.convolution import Gaussian1DKernel, Gaussian2DKernel, convolve
from astropy.visualization.wcsaxes.frame import EllipticalFrame
from scipy.interpolate import interp1d, interp2d
from scipy.integrate import trapz
import matplotlib.pyplot as plt
import matplotlib.colors as colors
from functools import lru_cache
from spectres import spectres
from pprint import pprint
from copy import copy, deepcopy
import os


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


def determine_rms(hdul, mission='', file=''):
    if mission == 'COGAL':

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
                elif (attribute == 'WCSAXES'):
                    hdul_rms[0].header['WCSAXES'] = 2

        hdul_rms[0].header['NAXIS']=2

        if mean:
            med = np.nanmean(cube_co_fix[0].data, axis=0)
        #     print('mean', med, 'K')
        else:
            med = np.nanmedian(cube_co_fix[0].data, axis=0)
        #     print('median', med, 'K')
        print('std for data < median:', cube_co_fix[0].data[cube_co_fix[0].data < med.reshape(1, *med.shape)].std(), 'K')

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


def regrid_observations(path='/media/hpc_backup/yanitski/projects/pdr/observational_data/MilkyWay/', mission=None,
                        target_header=None, target_kernel=None, output_file='obs_regridded.fits'):
    '''
    This function will regrid the specified mission data to the specified target header. This prepares the
    observational data for comparison to simulated data. For now, the target header must

    :param path: The directory storing the mission folders. Each mission folder should contain the relevant fits files.
    :param mission: The mission you wish to regrid. The default functionality is to regrid all missions.
    :param target_header: The target header, which must have longitude on axis 1 and latitude on axis 2. An optional
        axis 3 may be used for the velocity. If there is not a third dimension, the data will be regridded with a
        spacing of 1 km/s.
    :param target_kernel: The kernel to use for regridding.
    :return: Float -- 0 for success; 1 for fail
    '''

    if isinstance(mission, str):
        mission = [mission]

    if mission==None or mission[0]=='':
        if path[-1] != '/': path += '/'
        mission = os.listdir(path)
    elif mission[0] in os.listdir(path):
        print('Regridding {} survey'.format(mission[0]))
    else:
        print('Invalid mission... {} observations do not exist or are not downloaded'.format(mission[0]))
        return 1

    if target_header==None:
        print('Error: Please specify target_header to regrid the observations.')
        print('Regrid failed.')
        return 1

    if target_header['NAXIS'] == 3:
        target_vel = np.linspace(target_header['CRVAL3'] - target_header['CDELT3'] * (target_header['CRPIX3'] - 1),
                                 target_header['CRVAL3'] + target_header['CDELT3'] *
                                 (target_header['NAXIS3'] - target_header['CRPIX3']),
                                 num=target_header['NAXIS3'])
    else:
        print('Error: Please enter a target velocity to regrid the spectroscopic observations.')
        return 1
        target_vel = None

    min_vel = 0
    max_vel = 0

    CO1 = False
    CO2 = False
    iCO1 = False
    iCO2 = False
    dust = False
    cobe = False

    for survey in mission:

        if survey == 'COBE':
            survey = 'COBE-FIRAS'

        temp_header = copy(target_header)

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
            co1_gridder_err = cygrid.WcsGrid(target_header)
            co1_gridder_err.set_kernel(*target_kernel)
            CO1 = True
        elif survey == 'Mopra':
            co1_gridder = cygrid.WcsGrid(target_header)
            co1_gridder.set_kernel(*target_kernel)
            co1_gridder_err = cygrid.WcsGrid(target_header)
            co1_gridder_err.set_kernel(*target_kernel)
            CO1 = True
            ico1_gridder = cygrid.WcsGrid(target_header)
            ico1_gridder.set_kernel(*target_kernel)
            ico1_gridder_err = cygrid.WcsGrid(target_header)
            ico1_gridder_err.set_kernel(*target_kernel)
            iCO1 = True
        elif survey == 'ThrUMMS':
            co1_gridder = cygrid.WcsGrid(target_header)
            co1_gridder.set_kernel(*target_kernel)
            co1_gridder_err = cygrid.WcsGrid(target_header)
            co1_gridder_err.set_kernel(*target_kernel)
            CO1 = True
            ico1_gridder = cygrid.WcsGrid(target_header)
            ico1_gridder.set_kernel(*target_kernel)
            ico1_gridder_err = cygrid.WcsGrid(target_header)
            ico1_gridder_err.set_kernel(*target_kernel)
            iCO1 = True
        elif survey == 'SEDIGISM':
            ico2_gridder = cygrid.WcsGrid(target_header)
            ico2_gridder.set_kernel(*target_kernel)
            ico2_gridder_err = cygrid.WcsGrid(target_header)
            ico2_gridder_err.set_kernel(*target_kernel)
            iCO2 = True
        elif survey == 'Planck':
            dust = True
        elif survey == 'COBE-FIRAS':
            cobe = True
        else:
            print('Survey {} not available. Choose another.'.format(survey))
            continue

        for file in files:

            if file == 'regridded' or file == 'temp' or 'RMS' in file or survey == 'HiGAL':
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
                lon = np.linspace(obs[0].header['CRVAL2'] - obs[0].header['CDELT2'] * (obs[0].header['CRPIX2'] - 1),
                                  obs[0].header['CRVAL2'] + obs[0].header['CDELT2'] * (
                                              obs[0].header['NAXIS2'] - obs[0].header['CRPIX2']),
                                  num=obs[0].header['NAXIS2'])
                lat = np.linspace(obs[0].header['CRVAL3'] - obs[0].header['CDELT3'] * (obs[0].header['CRPIX3'] - 1),
                                  obs[0].header['CRVAL3'] + obs[0].header['CDELT3'] * (
                                              obs[0].header['NAXIS3'] - obs[0].header['CRPIX3']),
                                  num=obs[0].header['NAXIS3'])
                vel = np.linspace(obs[0].header['CRVAL1'] - obs[0].header['CDELT1'] * (obs[0].header['CRPIX1'] - 1),
                                  obs[0].header['CRVAL1'] + obs[0].header['CDELT1'] * (
                                              obs[0].header['NAXIS1'] - obs[0].header['CRPIX1']),
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

                obs_data = spectres(target_vel, vel, np.nan_to_num(obs[0].data, nan=0), fill=0)
                obs_data = np.nan_to_num(obs_data.reshape(-1, obs_data.shape[-1]), nan=0)
                obs_error = determine_rms(obs, mission=survey)[0].data.reshape(-1, 1)
                # obs_error = np.swapaxes(obs_error, 0, 2)
                # obs_error = np.swapaxes(obs_error, 0, 1)

                # Grid
                lon_mesh, lat_mesh = np.meshgrid(lon, lat)
                co1_gridder.grid(lon_mesh.flatten(), lat_mesh.flatten(), obs_data)
                co1_gridder_err.grid(lon_mesh.flatten(), lat_mesh.flatten(), obs_error)

                if vel.min() < min_vel:
                    min_vel = vel.min()
                if vel.max() > max_vel:
                    max_vel = vel.max()
            elif survey == 'Mopra' and ('_Vfull.fits' in file) and (('12CO' in file) or ('13CO' in file)):

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
                obs_error = np.nanmean(fits.open(path + survey + '/' + file.replace('_Vfull', '.sigma'))[0].data)

                # Get axes
                lon = np.linspace(obs[0].header['CRVAL1'] - obs[0].header['CDELT1'] * (obs[0].header['CRPIX1'] - 1),
                                  obs[0].header['CRVAL1'] + obs[0].header['CDELT1'] * (
                                              obs[0].header['NAXIS1'] - obs[0].header['CRPIX1']),
                                  num=obs[0].header['NAXIS1'])
                lat = np.linspace(obs[0].header['CRVAL2'] - obs[0].header['CDELT2'] * (obs[0].header['CRPIX2'] - 1),
                                  obs[0].header['CRVAL2'] + obs[0].header['CDELT2'] * (
                                              obs[0].header['NAXIS2'] - obs[0].header['CRPIX2']),
                                  num=obs[0].header['NAXIS2'])
                vel = np.linspace(obs[0].header['CRVAL3'] - obs[0].header['CDELT3'] * (obs[0].header['CRPIX3'] - 1),
                                  obs[0].header['CRVAL3'] + obs[0].header['CDELT3'] * (
                                              obs[0].header['NAXIS3'] - obs[0].header['CRPIX3']),
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
                obs_data = spectres(target_vel, vel, np.nan_to_num(obs[0].data.T, nan=0), fill=0)
                obs_data = obs_data.reshape(-1, obs_data.shape[-1])
                obs_error = fits.open(path + survey + '/' + file.replace('_Vfull', '.sigma'))[0].data.reshape(-1, 1)
                i_nan = np.isnan(obs_error)
                del obs

                # Grid
                lon_mesh, lat_mesh = np.meshgrid(lon, lat)
                if '12CO' in transition:
                    co1_gridder.grid(lon_mesh.flatten(), lat_mesh.flatten(), obs_data)
                    co1_gridder_err.grid(lon_mesh.flatten()[~i_nan.flatten()], lat_mesh.flatten()[~i_nan.flatten()], obs_error[~i_nan])
                elif '13CO' in transition:
                    ico1_gridder.grid(lon_mesh.flatten(), lat_mesh.flatten(), obs_data)
                    ico1_gridder_err.grid(lon_mesh.flatten()[~i_nan.flatten()], lat_mesh.flatten()[~i_nan.flatten()], obs_error[i_nan])

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
                lon = np.linspace(obs[0].header['CRVAL1'] - obs[0].header['CDELT1'] * (obs[0].header['CRPIX1'] - 1),
                                  obs[0].header['CRVAL1'] + obs[0].header['CDELT1'] * (
                                              obs[0].header['NAXIS1'] - obs[0].header['CRPIX1']),
                                  num=obs[0].header['NAXIS1'])
                lat = np.linspace(obs[0].header['CRVAL2'] - obs[0].header['CDELT2'] * (obs[0].header['CRPIX2'] - 1),
                                  obs[0].header['CRVAL2'] + obs[0].header['CDELT2'] * (
                                              obs[0].header['NAXIS2'] - obs[0].header['CRPIX2']),
                                  num=obs[0].header['NAXIS2'])
                vel = np.linspace(obs[0].header['CRVAL3'] - obs[0].header['CDELT3'] * (obs[0].header['CRPIX3'] - 1),
                                  obs[0].header['CRVAL3'] + obs[0].header['CDELT3'] * (
                                              obs[0].header['NAXIS3'] - obs[0].header['CRPIX3']),
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
                obs_data = spectres(target_vel, vel, np.nan_to_num(obs_data, nan=0), fill=0)
                obs_data = obs_data.reshape(-1, obs_data.shape[-1])
                if '12' in transition:
                    obs_error = np.full(obs_data.shape, 1.3)  # from Barnes et al. (2015)
                elif '13' in transition:
                    obs_error = np.full(obs_data.shape, 0.7)  # from Barnes et al. (2015)
                del obs

                # Grid
                lon_mesh, lat_mesh = np.meshgrid(lon, lat)
                if '12' in transition:
                    co1_gridder.grid(lon_mesh.flatten(), lat_mesh.flatten(), obs_data)
                    co1_gridder_err.grid(lon_mesh.flatten(), lat_mesh.flatten(), obs_error)
                elif '13' in transition:
                    ico1_gridder.grid(lon_mesh.flatten(), lat_mesh.flatten(), obs_data)
                    ico1_gridder_err.grid(lon_mesh.flatten(), lat_mesh.flatten(), obs_error)

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
                lon = np.linspace(obs[0].header['CRVAL1'] - obs[0].header['CDELT1'] * (obs[0].header['CRPIX1'] - 1),
                                  obs[0].header['CRVAL1'] + obs[0].header['CDELT1'] * (
                                              obs[0].header['NAXIS1'] - obs[0].header['CRPIX1']),
                                  num=obs[0].header['NAXIS1'])
                lat = np.linspace(obs[0].header['CRVAL2'] - obs[0].header['CDELT2'] * (obs[0].header['CRPIX2'] - 1),
                                  obs[0].header['CRVAL2'] + obs[0].header['CDELT2'] * (
                                              obs[0].header['NAXIS2'] - obs[0].header['CRPIX2']),
                                  num=obs[0].header['NAXIS2'])
                vel = np.linspace(obs[0].header['CRVAL3'] - obs[0].header['CDELT3'] * (obs[0].header['CRPIX3'] - 1),
                                  obs[0].header['CRVAL3'] + obs[0].header['CDELT3'] * (
                                              obs[0].header['NAXIS3'] - obs[0].header['CRPIX3']),
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
                obs_data = spectres(target_vel, vel, np.nan_to_num(obs_data, nan=0), fill=0)
                obs_data = obs_data.reshape(-1, obs_data.shape[-1])
                obs_error = determine_rms(obs, mission=survey, file=file)[0].data.reshape(-1, 1)
                i_nan = np.isnan(obs_error)
                # print(np.nanmean(obs_error), np.mean(obs_error[~i_nan]))
                # np.save('/home/yanitski/obs_error.npy', obs_error)
                # exit()
                del obs

                # Grid
                lon_mesh, lat_mesh = np.meshgrid(lon, lat)
                ico2_gridder.grid(lon_mesh.flatten(), lat_mesh.flatten(), obs_data)
                ico2_gridder_err.grid(lon_mesh.flatten()[~i_nan.flatten()], lat_mesh.flatten()[~i_nan.flatten()],
                                      obs_error[~i_nan])

                if vel.min() < min_vel:
                    min_vel = vel.min()
                if vel.max() > max_vel:
                    max_vel = vel.max()
            elif survey == 'Planck':

                print(file)

                # Specify observation type
                transitions = ['Dust']

                # Open file and get data
                obs = fits.open(path + survey + '/' + file)
                obs_data = obs[1].data['I_ML'] * 1e-6
                obs_error = obs[1].data['I_RMS'] * 1e-6

                nside = 256
                npix = hp.nside2npix(nside)

                # fix header
                if 'CDELT3' in temp_header.keys():
                    temp_header['NAXIS'] = 2
                    del temp_header['NAXIS3']
                    del temp_header['CTYPE3']
                    del temp_header['CDELT3']
                    del temp_header['CRVAL3']
                    del temp_header['CRPIX3']

                # Grid
                lat_mesh, lon_mesh = np.degrees(hp.pix2ang(nside=nside, ipix=np.arange(npix), nest='True'))
                lat_mesh = 90 - lat_mesh
                dust_gridder = cygrid.WcsGrid(temp_header)
                dust_gridder.set_kernel(*target_kernel)
                dust_gridder.grid(lon_mesh, lat_mesh, obs_data)
                dust_gridder_err = cygrid.WcsGrid(temp_header)
                dust_gridder_err.set_kernel(*target_kernel)
                dust_gridder_err.grid(lon_mesh, lat_mesh, obs_error)
            elif survey == 'COBE-FIRAS':

                print(file)

                # Specify transitions
                if 'HIGH' in file.split('.')[0]:
                    transitions = ['CO 6', 'C 2', 'H2O f1113', 'N+ 1', 'H2O 2', 'C+ 1', 'O 1', 'Si', 'N+ 2', 'CH 2']   #'CO 6', 'O 2'
                    transition_indeces = [0, 1, 2, 3, 4, 5, 6, 7, 8, 9]   #0, 6
                elif 'HRES' in file.split('.')[0]:
                    transitions = ['CO 1', 'CO 2', 'CO 3', 'O2 13', 'CO 4', 'C 1', 'H2O f557', 'CO 5']
                    transition_indeces = [0, 1, 2, 4, 5, 7]
                elif 'LOWF' in file.split('.')[0]:
                    transitions = ['CO 1', 'CO 2', 'CO 3', 'CO 4', 'C 1', 'CO 5']
                    transition_indeces = [0, 1, 2, 4, 5, 7]

                # Open data and convert to brightness temperature
                obs = fits.open(path + survey + '/' + file)
                linfrq = np.array([obs[0].header[key] for key in obs[0].header.keys() if 'LINFRQ' in key])
                # obs_data = np.nan_to_num(obs[1].data['LINE_FLU'], nan=0) * (2.9979**3) \
                #            / (linfrq**3) / 2 / 1.38 * 10 ** 8
                # obs_data = np.nan_to_num(obs[1].data['LINE_FLU'], nan=0) / (2.9979**3) \
                #            * (linfrq**3) * 2 * 1.38 / 10 ** 8
                obs_data = obs[1].data['LINE_FLU']
                obs_error = (np.nan_to_num(obs[1].data['LINE_FL2'], nan=0) * (2.9979**3) / (linfrq**3) / 2
                             / 1.38 * 10 ** 8)

                # Get axes
                lon_mesh = obs[1].data['GAL_LON']
                lat_mesh = obs[1].data['GAL_LAT']

                if 'HIGH' in file and False:
                    lon_test = np.linspace(2.5, 357.5, num=72)
                    lat_test = np.zeros(lon_test.shape)
                    i_cii = np.zeros(lon_test.shape)
                    for i in range(i_cii.size):
                        idx = (np.abs(lon_mesh-lon_test[i])<2.5) & (np.abs(lat_mesh-lat_test[i])<1)
                        i_cii[i] = np.average(obs_data[idx, 5])
                    fig, ax = plt.subplots(1, 1)
                    ax.step(lon_test, i_cii)
                    plt.show()

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
                temp_header['TRANSI'] = ', '.join('{}'.format(_) for _ in np.arange(len(transitions)))
                grid_hdu = fits.PrimaryHDU(data=gridder.get_datacube(), header=fits.Header(temp_header))
                grid_hdu_err = fits.PrimaryHDU(data=gridder_err.get_datacube(), header=fits.Header(temp_header))
                grid_hdu.writeto(path + survey + '/regridded/temp/' + file.replace('.FITS', '_regridded.fits'),
                                 overwrite=True, output_verify='ignore')
                grid_hdu_err.writeto(path + survey + '/regridded/temp/' + file.replace('.FITS', '_regridded_error.fits'),
                                     overwrite=True, output_verify='ignore')
            else:
                # print('The specified survey {} is unavailable. Please add it to ´model selection´.'.format(survey))
                continue

        if CO1:
            temp_header['TRANSL'] = 'CO 1'
            temp_header['TRANSI'] = '0'
            grid_hdu = fits.PrimaryHDU(data=co1_gridder.get_datacube(), header=fits.Header(temp_header))
            grid_hdu_err = fits.PrimaryHDU(data=co1_gridder_err.get_datacube(), header=fits.Header(temp_header))
            grid_hdu.writeto(path + survey + '/regridded/temp/' +
                             'co1_test_regridded.fits', overwrite=True, output_verify='ignore')
            grid_hdu_err.writeto(path + survey + '/regridded/temp/' +
                                 'co1_test_regridded_error.fits', overwrite=True, output_verify='ignore')
        if CO2:
            temp_header['TRANSL'] = 'CO 2'
            temp_header['TRANSI'] = '0'
            grid_hdu = fits.PrimaryHDU(data=co2_gridder.get_datacube(), header=fits.Header(temp_header))
            grid_hdu_err = fits.PrimaryHDU(data=co2_gridder_err.get_datacube(), header=fits.Header(temp_header))
            grid_hdu.writeto(path + survey + '/regridded/temp/' +
                             'co2_test_regridded.fits', overwrite=True, output_verify='ignore')
            grid_hdu_err.writeto(path + survey + '/regridded/temp/' +
                                 'co2_test_regridded_error.fits', overwrite=True, output_verify='ignore')
        if iCO1:
            temp_header['TRANSL'] = '13CO 1'
            temp_header['TRANSI'] = '0'
            grid_hdu = fits.PrimaryHDU(data=ico1_gridder.get_datacube(), header=fits.Header(temp_header))
            grid_hdu_err = fits.PrimaryHDU(data=ico1_gridder_err.get_datacube(), header=fits.Header(temp_header))
            grid_hdu.writeto(path + survey + '/regridded/temp/' +
                             '13co1_test_regridded.fits', overwrite=True, output_verify='ignore')
            grid_hdu_err.writeto(path + survey + '/regridded/temp/' +
                                 '13co1_test_regridded_error.fits', overwrite=True, output_verify='ignore')
        if iCO2:
            temp_header['TRANSL'] = '13CO 2'
            temp_header['TRANSI'] = '0'
            grid_hdu = fits.PrimaryHDU(data=ico2_gridder.get_datacube(), header=fits.Header(temp_header))
            grid_hdu_err = fits.PrimaryHDU(data=ico2_gridder_err.get_datacube(), header=fits.Header(temp_header))
            grid_hdu.writeto(path + survey + '/regridded/temp/' +
                             '13co2_test_regridded.fits', overwrite=True, output_verify='ignore')
            grid_hdu_err.writeto(path + survey + '/regridded/temp/' +
                                 '13co2_test_regridded_error.fits', overwrite=True, output_verify='ignore')
        if dust:
            temp_header['TRANSL'] = 'Dust'
            temp_header['TRANSI'] = '0'
            grid_hdu = fits.PrimaryHDU(data=dust_gridder.get_datacube(), header=fits.Header(temp_header))
            grid_hdu_err = fits.PrimaryHDU(data=dust_gridder_err.get_datacube(), header=fits.Header(temp_header))
            grid_hdu.writeto(path + survey + '/regridded/temp/' + 'planck_dust_regridded.fits',
                             overwrite=True, output_verify='ignore')
            grid_hdu_err.writeto(path + survey + '/regridded/temp/' + 'planck_dust_regridded_error.fits',
                                 overwrite=True, output_verify='ignore')

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
    if (mission != 'COBE-FIRAS') and (mission != 'Planck') and not ('error' in filename):
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


def model_selection(path='/mnt/hpc_backup/yanitski/projects/pdr/KT3_history/MilkyWay', mission=None, lat=None,
                    model_dir='', model_param=[[]], cmap='gist_ncar', PLOT=False, PRINT=False, debug=False):
    '''
    This function will compare the Milky Way models to the position-velocity observations depending on the
      dataset. It will utilise the sightlines of the dataset.

      :param path: directory containing the model folders.
      :param model_dir: directory format for each model in the grid.
      :param model_param: a list of lists containing the parameter space of the models.
      :param resolution: the voxel size (in pc) of the models used in the comparison; used in model folder name.
      :param mission: as this is setup for the Milky Way, the mission selection is
             - COBE-FIRAS
             - COGAL
             - Mopra
             - ThrUMMS
             - SEDIGISM
             - Planck
             - Hi-GAL (in development)
      :param cmap: color map for the plots.
      :param PLOT: flag to save plots.
      :param PRINT: flag to enable verbose output.
      :param debug: flag to run in debug mode.

      currently, this will look solely at the galactic plane.
    '''

    if isinstance(mission, str):
        mission = [mission]
    elif isinstance(mission, list):
        pass
    elif mission == '' or mission == None:
        mission = os.listdir(path.replace('KT3_history', 'observational_data'))
    else:
        print('Please specify a list of missions to compare the models.')
        return

    if model_dir == '' or model_param == [[]]:
        print('Please specify both model_dir and model_param.')
        return

    model_params = np.meshgrid(*model_param)
    # model_params = zip(np.transpose([model_params[n].flatten() for n in range(len(model_param))]))
    model_params = np.transpose([model_params[n].flatten() for n in range(len(model_param))])

    obs_data = []
    lon = []
    vel = []

    for survey in mission:

        print('\n\n', survey)

        if survey == 'HiGAL':
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

            if not '.fits' in file or 'error' in file:
                continue

            # print(file)

            if '.fits' in file:
                obs = fits.open(directory + file)
                if survey == 'COGAL' or survey == 'ThrUMMS':
                    obs_error = fits.open(directory + file.replace('.fits', '_error.fits'))[0].data[0, :, :]
                else:
                    obs_error = fits.open(directory + file.replace('.fits', '_error.fits'))[0].data
            elif '.sav' in file or '.idl' in file:
                obs = readsav(directory + file)
            else:
                continue

            if survey == 'COBE' or survey == 'COBE-FIRAS':
                # if 'HIGH' in file.split('.')[0]:
                #     transitions = ['C 2', 'C+ 1']   #'CO 6', 'O 2'
                #     transition_indeces = [1, 5]   #0, 6
                # elif 'HRES' in file.split('.')[0]:
                #     transitions = ['CO 1', 'CO 2', 'CO 3', 'CO 4', 'C 1', 'CO 5']
                #     transition_indeces = [0, 1, 2, 4, 5, 7]
                # elif 'LOWF' in file.split('.')[0]:
                #     transitions = ['CO 1', 'CO 2', 'CO 3', 'CO 4', 'C 1', 'CO 5']
                #     transition_indeces = [0, 1, 2, 4, 5, 7]
                # obs = fits.open(directory + file)
                # linfrq = np.array([obs[0].header[key] for key in obs[0].header.keys() if 'LINFRQ' in key])
                # obs_data = obs[1].data['LINE_FLU'] * (2.9979**3) / (linfrq**3) / 2 / 1.38 * 10 ** 8
                # obs_error = obs[1].data['LINE_FL2'] * (2.9979**3) / (linfrq**3) / 2 / 1.38 * 10 ** 8
                # lon = obs[1].data['GAL_LON']
                # lat = obs[1].data['GAL_LAT']
                # vel = None
                # sightline_indeces = np.abs(lat) <= 2
                if file == 'craig.idl':
                    # these values are manually coded-in since the data is not in the file
                    linfrq = np.array([115.3, 230.5, 345.8, 424.8, 461.0, 492.2, 556.9, 576.3, 691.5, 808.1,
                                       1113, 1460, 2226, 1901, 2060, 2311, 2459, 2589, 921.8])
                    transitions = np.array(['CO 1', 'CO 2', 'CO 3', 'CO 4', 'C 3', 'CO 5', 'CO 6', 'CO 7 + C 1', 'C+ 1', 'O 1', 'CO 8'])
                    transition_indeces = np.array([0, 1, 2, 4, 5, 7, 8, 9, 13, 14, 18])

                    obs_data = obs['amplitude'] / (2.9979**3) * (linfrq**3) * 2 * 1.38 / 10 ** 8
                    obs_error = obs['sigma'] / (2.9979**3) * (linfrq**3) * 2 * 1.38 / 10 ** 8
                    lon_obs = obs['long']
                    lat_obs = np.array([0])
                else:
                    obs_data = obs[0].data
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

                if os.path.isdir(
                        path + 'fit_results/{}/{}/{}/'.format(survey, file.replace('.fits', ''), transition)) == False:
                    os.makedirs(path + 'fit_results/{}/{}/{}/'.format(survey, file.replace('.fits', ''), transition))

                if (survey == 'COBE') or (survey == 'COBE-FIRAS'):
                    vmin = obs_data[int(transition_indeces[i]), :, :].min()
                    vmax = obs_data[int(transition_indeces[i]), :, :].max()
                    print(transition)
                    print(vmin, vmax)
                else:
                    vmin = np.nanmin(obs_data)
                    vmax = np.nanmax(obs_data)

                for param in model_params:

                    if PRINT:
                        print(model_dir.format(*param))

                    # Check existance of model
                    dir_model = model_dir.format(*param) + 'channel_intensity.fits'
                        # '/r{}_cm{}-{}_d{}_uv{}/channel_intensity.fits'.format(resolution, fcl,
                        #                                                               ficl, fden, fuv)
                    if not os.path.isfile(path + dir_model): continue

                    # Open the model
                    params.append(param)
                    model = fits.open(path + dir_model)
                    # pprint(model[1].header)

                    # Locate the species transition (only one dust value is considered; constant background)
                    if survey == 'Planck':
                        i_spec = np.array([0])
                    else:
                        i_spec = np.where(np.asarray(model[1].header['SPECIES'].split(', ')) == transition)[0]
                    if i_spec.size == 0:
                        print('{} not in model.'.format(transition))
                        continue
                    else:
                        i_spec = i_spec[0]

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
                    if isinstance(lat, int) or isinstance(lat, float):
                        lat_min = lat
                        lat_max = lat
                    else:
                        lat_min = lat_obs[i_lat_obs_init].min()
                        lat_max = lat_obs[i_lat_obs_init].max()
                    if survey == 'COBE-FIRAS' or survey =='Planck':
                        i_lat_model = (lat_model >= lat_min) \
                                      & (lat_model <= lat_max)
                        i_lat_obs = (lat_obs >= lat_model[i_lat_model].min()) \
                                    & (lat_obs <= lat_model[i_lat_model].max())
                    else:
                        i_lat_model = (lat_model >= lat_min) \
                                      & (lat_model <= lat_max)
                        i_lat_obs = (lat_obs >= lat_model[i_lat_model].min()) \
                                    & (lat_obs <= lat_model[i_lat_model].max())

                    # Interpolate at the observational axes
                    # print('grid shape:', model[1].data[:, i_lat, :, i_spec].shape)
                    # print('axes shape:', vel_model.size, lon_model.size)
                    if (survey == 'COBE') or (survey == 'COBE-FIRAS'):
                        # print(vel_model.shape,
                        #       model[1].data[:, i_lat_model, :, i_spec].shape,
                        #       model[2].data[:, i_lat_model, :, 0].shape,
                        #       i_lat_obs.size)
                        model_interp = trapz(model[1].data[:, i_lat_model, :, i_spec]
                                             - model[2].data[:, i_lat_model, :, 0],
                                             vel_model, axis=1)
                        # f_model_interp = interp2d(lon_model, lat_model, model_int[:, :, i_spec])
                        # model_interp = []
                        # for _ in range(len(lat)):
                        #     model_interp.append(f_model_interp(lon[_], lat[_]))
                        # model_interp = copy(obs_data[:, int(transition)])
                        #     print(model_interp[-1])
                        # print(model_interp)
                        # input()

                    elif survey == 'Planck':
                        model_interp = model[2].data[0, i_lat_model, :, i_spec]

                    else:
                        # print(vel)
                        # print(lon)
                        # print(vel_model)
                        # print(lon_model)

                        print(vel_model.shape,
                              model[1].data[:, i_lat_model, :, i_spec].shape,
                              model[2].data[:, i_lat_model, :, 0].shape,
                              obs_data[:, i_lat_obs, :].shape,
                              np.swapaxes(obs_data[:, i_lat_obs, :], 0, 1).shape)
                        model_interp = (model[1].data[:, i_lat_model, :, i_spec]
                                        - model[2].data[:, i_lat_model, :, 0])
                        obs_data_temp = np.swapaxes(obs_data[:, i_lat_obs, :], 0, 1)
                        obs_error_temp = np.swapaxes([obs_error[i_lat_obs, :]]*obs_data_temp.shape[1],
                                                     0, 1)

                        # if (mission == 'SEDIGISM') and (lon.mean() < 0.5):
                        #     lon_model[-10:] -= 360
                        # elif (mission == 'SEDIGISM') and (lon.mean() > 359.5):
                        #     lon_model[:10] += 360
                        #
                        # i_sort = np.argsort(lon_model)
                        # lon_model = lon_model[i_sort]
                        # model_data = model_data[:, i_sort]
                        #
                        # #print(vel_model[450], 'km/s', lon_model[80], 'deg', model_data[450, 80])
                        #
                        # f_model_interp = interp2d(lon_model, vel_model, model_data, kind='linear')
                        #
                        # #print(vel_model[450], 'km/s', lon_model[80], 'deg', f_model_interp(lon_model[80], vel_model[450]))
                        # model_data[model_data == 0] = np.nan
                        #
                        # # print(f_model_interp(100, -25))
                        # # print(f_model_interp(100, 25))
                        #
                        # #print(vel[0], 'km/s', lon[0], 'deg')
                        # model_interp = f_model_interp(lon, vel)
                        #
                        # # reverse longitude dimension, which was sorted when interpolated
                        # if lon[0]>lon[-1]:
                        #     model_interp = model_interp[:, ::-1]
                        #
                        # # reverse velocity dimension, which was sorted when interpolated
                        # if vel[0]>vel[-1]:
                        #     model_interp = model_interp[::-1, :]

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

                    # Calculate the chi**2 and overall likelihood
                    if ((survey == 'COBE') or (survey == 'COBE-FIRAS')) and '.fits' in file:
                        chi2 = (obs_data[int(transition_indeces[i]), i_lat_obs, :] - model_interp) ** 2 / \
                                obs_error[int(transition_indeces[i]), i_lat_obs, :] ** 2
                        likelihood = np.exp(-chi2 / 2) / \
                                     np.sqrt(2 * np.pi * obs_error[int(transition_indeces[i]), i_lat_obs, :])
                    elif ((survey == 'COBE') or (survey == 'COBE-FIRAS')) and ('.idl' in file or '.sav' in file):
                        chi2 = (obs_data[int(transition_indeces[i])] - model_interp) ** 2 / \
                                obs_error[int(transition_indeces[i])] ** 2
                        likelihood = np.exp(-chi2 / 2) / \
                                     np.sqrt(2 * np.pi * obs_error[i_lat_obs, :])
                    elif survey == 'Planck':
                        chi2 = (obs_data[i_lat_obs, :] - model_interp) ** 2 / \
                                obs_error[i_lat_obs, :] ** 2
                        likelihood = np.exp(-chi2 / 2) / \
                                     np.sqrt(2 * np.pi * obs_error[i_lat_obs, :])
                    else:
                        chi2 = (obs_data_temp - model_interp) ** 2 / obs_error_temp ** 2
                        likelihood = np.exp(-chi2 / 2) / np.sqrt(2 * np.pi * obs_error_temp)

                    # print(model_interp)
                    # print(obs_data[~np.isnan(obs_data)])
                    # print(chi2)
                    chi2 = np.nan_to_num(chi2, nan=0)
                    likelihood = np.nan_to_num(likelihood, nan=1)
                    likelihood[likelihood == 0] = 1e-100
                    # input()
                    np.save(path + 'fit_results/{}/{}/{}/{}_chi2.npy'.format(survey, file.replace('.fits', ''),
                                                                             transition, dir_model.replace('/', '')
                                                                             .replace('channel_intensity.fits', '')),
                            chi2)
                    np.save(path + 'fit_results/{}/{}/{}/{}_loglikelihood.npy'.format(survey, file.replace('.fits', ''),
                                                                                      transition,
                                                                                      dir_model.replace('/', '')
                                                                                      .replace('channel_intensity.fits',
                                                                                               '')),
                            np.log10(likelihood))

                    chi2_grid.append(chi2.sum())
                    loglikelihood_grid.append(np.log10(likelihood).sum())
                    #                             print('  ', likelihood.min())

                    # print(obs_data.min(), obs_data.max())

                    if PLOT:
                        if (survey == 'COBE') or (survey == 'COBE-FIRAS'):
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
                    del likelihood
                    del lon_model
                    del lat_model
                    del vel_model

                if np.size(loglikelihood_grid) == 0:
                    continue
                # print(loglikelihood_grid)

                i_bestfit2 = np.where(loglikelihood_grid == np.max(loglikelihood_grid))[0][0]
                # print(i_bestfit2)
                print('The best-fitting model for {} with transition {}\n'.format(file.replace('.fits', ''), transition) +
                      '  has parameters ' + model_dir.format(*params[i_bestfit2]))

                np.save(path + 'fit_results/{}/{}/{}/chi2.npy'.format(survey, file.replace('.fits', ''), transition),
                        chi2_grid)
                np.save(path + 'fit_results/{}/{}/{}/loglikelihood.npy'.format(survey, file.replace('.fits', ''), transition),
                        loglikelihood_grid)
                np.save(path + 'fit_results/{}/{}/{}/parameters.npy'.format(survey, file.replace('.fits', ''), transition),
                        params)

            obs.close()
            del obs_data
            del obs_error
            del lon_obs
            del lat_obs
            if not vel is None:
                del vel

        try:
            i_bestfit1 = np.where(chi2_grid == np.min(chi2_grid))[0][0]
            print('The best-fitting model has parameters' +
                  model_dir.format(*params[i_bestfit1]))
        except ValueError:
            pass
        # i_bestfit2 = np.where(loglikelihood_grid == np.max(loglikelihood_grid))[0][0]


    return


def plot_comparison(path='/mnt/hpc_backup/yanitski/projects/pdr/KT3_history/MilkyWay/fit_results/', file_format='',
                    model_param=[[]], missions=[], i_x=0, i_y=1,
                    xlabel='', ylabel='', clabel='', title='', fontsize=16, save_plot=False):
    #

    # Check that the missions are specified properly.
    if isinstance(missions, str):
        mission = [missions]
    elif missions == '' or missions == None or missions == []:
        missions = os.listdir(path)
    elif isinstance(missions, list):
        pass
    else:
        print('Please specify a list of missions to compare the models.')
        return

    if file_format == '' or model_param == [[]]:
        print('Please specify both file_format and model_param.')
        return

    dimensions = len(model_param)
    naxis = np.zeros(dimensions)
    naxis[[i_x, i_y]] = 1
    lenaxis = [len(arr) for arr in model_param]

    model_param_grid = np.meshgrid(*model_param)
    x_grid = model_param_grid[i_x]
    y_grid = model_param_grid[i_y]
    # model_params = zip(np.transpose([model_params[n].flatten() for n in range(len(model_param))]))
    # model_params = np.transpose([model_params[n].flatten() for n in range(len(model_params))])

    for survey in missions:

        survey_files = os.listdir(path + survey + '/')

        log_likelihood = np.zeros(x_grid.shape)

        for survey_file in survey_files:

            transitions = os.listdir(path + survey + '/' + survey_file + '/')

            transitions_skipped = []
            for t in transitions:
                if len(os.listdir(path + survey + '/' + survey_file + '/' + t + '/')) == 0:
                    transitions.remove(t)
                    transitions_skipped.append(t)
            if len(transitions_skipped):
                print('  transitions {} not available.'.format(', '.join(transitions_skipped)))

            for i in range(log_likelihood.shape[0]):
                for j in range(log_likelihood.shape[1]):

                    filename = file_format.format(x_grid[i, j], y_grid[i, j]) + '_loglikelihood.npy'
                    param_log_likelihood = [np.load(path + survey + '/' + survey_file + '/' + t + '/' + filename)
                                            for t in transitions]
                    log_likelihood[i, j] = log_likelihood[i_x, i_y] + np.nansum(param_log_likelihood)

        fig, axes = plt.subplots(1, 1, figsize=(10, 10))
        if not isinstance(axes, np.ndarray):
            axes = np.asarray([axes])

        cm = axes[0].imshow(log_likelihood, extent=[0, lenaxis[i_x], 0, lenaxis[i_y]])
        axes[0].set_aspect(lenaxis[i_x]/lenaxis[i_y])
        cb = fig.colorbar(cm, ax=axes[0], fraction=0.1, aspect=20)

        axes[0].set_xticks(np.arange(lenaxis[i_x])+0.5)
        axes[0].set_xticklabels([str(param) for param in model_param[i_x]])
        axes[0].set_yticks(np.arange(lenaxis[i_y])+0.5)
        axes[0].set_yticklabels([str(param) for param in model_param[i_y][::-1]])

        axes[0].set_xlabel(xlabel, fontsize=fontsize)
        axes[0].set_ylabel(ylabel, fontsize=fontsize)
        cb.ax.set_ylabel(clabel, fontsize=fontsize)
        if title == '':
            axes[0].set_title(survey, fontsize=fontsize)

        plt.show()

    return
