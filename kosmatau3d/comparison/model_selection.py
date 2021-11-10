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
from scipy.io import readsav
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


def model_selection(path='/mnt/hpc_backup/yanitski/projects/pdr/KT3_history/MilkyWay', missions=None, lat=None,
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

    # Check that the missions are specified properly.
    if missions == '' or missions == None or missions == []:
        missions = os.listdir(path.replace('KT3_history', 'observational_data'))
    elif isinstance(missions, str):
        missions = [missions]
    elif isinstance(missions, list):
        pass
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

            if not ('.fits' in file or '.idl' in file) or 'error' in file:
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

                if os.path.isdir(
                        path + 'fit_results/{}/{}/{}/'.format(survey, file.replace('.fits', ''), transition)) == False:
                    os.makedirs(path + 'fit_results/{}/{}/{}/'.format(survey, file.replace('.fits', ''), transition))

                if (survey == 'COBE') or (survey == 'COBE-FIRAS'):
                    if '.idl' in file:
                        vmin = obs_data[:, int(transition_indeces[i])].min()
                        vmax = obs_data[:, int(transition_indeces[i])].max()
                    else:
                        vmin = obs_data[int(transition_indeces[i]), :, :].min()
                        vmax = obs_data[int(transition_indeces[i]), :, :].max()
                    # print(transition)
                    # print(vmin, vmax)
                else:
                    vmin = np.nanmin(obs_data)
                    vmax = np.nanmax(obs_data)

                for param in model_params:

                    if PRINT:
                        print('  ' + model_dir.format(*param) + '    \r', end='')

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
                        print('  {} not in model.            '.format(transition))
                        break
                    # else:
                    #     i_spec = i_spec[0]

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
                        not (survey == 'COBE-FIRAS' or survey == 'Planck')):
                        lat_min = lat
                        lat_max = lat
                    else:
                        lat_min = lat_obs[i_lat_obs_init].min()
                        lat_max = lat_obs[i_lat_obs_init].max()

                    if survey == 'COBE-FIRAS' or survey == 'Planck':
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
                    if (survey == 'COBE') or (survey == 'COBE-FIRAS'):
                        idx_t = np.ix_(np.arange(vel_model.size), i_lat_model, i_lon_model, i_spec)
                        idx_d = np.ix_(np.arange(vel_model.size), i_lat_model, i_lon_model, [0])
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
                            obs_data_final = obs_data[idx_obs].reshape(-1)#[:, int(transition_indeces[i])]
                            obs_error_final = obs_error[idx_obs].reshape(-1)#[:, int(transition_indeces[i])]
                            model_interp = model_interp.reshape(-1)
                        else:
                            idx_obs = np.ix_([int(transition_indeces[i])], i_lat_obs, np.arange(obs_data.shape[2]))
                            obs_data_final = obs_data[idx_obs][0, :, :]#[int(transition_indeces[i]), i_lat_obs, :]
                            obs_error_final = obs_error[idx_obs][0, :, :]#[int(transition_indeces[i]), i_lat_obs, :]
                            model_interp = model_interp[:, :, 0]
                        # print(model_interp.shape, obs_data_final.shape)

                    elif survey == 'Planck':
                        idx_d = np.ix_([0], i_lat_model, np.arange(lon_model.size), [0])
                        model_interp = model[2].data[idx_d][0, :, :, 0]#[0, i_lat_model, :, i_spec]
                        idx_obs = np.ix_(i_lat_obs, np.arange(obs_data.shape[1]))
                        obs_data_final = obs_data[idx_obs]#[i_lat_obs, :]
                        obs_error_final = obs_error[idx_obs]#[i_lat_obs, :]
                        # print(model_interp.shape, obs_data_final.shape)

                    else:
                        # print(vel_model.shape,
                        #       model[1].data[:, i_lat_model, :, i_spec].shape,
                        #       model[2].data[:, i_lat_model, :, 0].shape,
                        #       obs_data[:, i_lat_obs, :].shape,
                        #       np.swapaxes(obs_data[:, i_lat_obs, :], 0, 1).shape)
                        # input()
                        model_data = (model[1].data[:, i_lat_model, :, i_spec]
                                        - model[2].data[:, i_lat_model, :, 0])
                        obs_data_temp = np.swapaxes(obs_data[:, i_lat_obs, :], 0, 1)
                        if comp_type == 'integrated':
                            model_interp = np.trapz(model_data, vel_model, axis=model_data.shape.index(vel_model.size))
                            obs_data_final = np.trapz(obs_data_temp, vel_obs, axis=obs_data_temp.shape.index(vel_obs.size))
                            obs_error_final = obs_error[i_lat_obs, :]
                        elif comp_type == 'pv':
                            model_interp = np.swapaxes(model_data, 0, 1)
                            obs_data_final = copy(obs_data[:, i_lat_obs, :])
                            obs_error_final = np.asarray([obs_error[i_lat_obs, :]]*obs_data_temp.shape[1])
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
                        i_signal = np.where(obs_data_final >= obs_data_final.min())
                    else:
                        i_signal = np.where(obs_data_final > obs_error_final)
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
                    chi2[i_signal] = (obs_data_final[i_signal] - model_interp[i_signal]) ** 2 / \
                                     obs_error_final[i_signal] ** 2 / (obs_data_final[i_signal].size - len(param))
                    loglikelihood[i_signal] = -chi2[i_signal] / 2 \
                                              - 0.5 * np.log(np.sqrt(2 * np.pi) * obs_error_final[i_signal])
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


def plot_comparison(path='/mnt/hpc_backup/yanitski/projects/pdr/KT3_history/MilkyWay/fit_results/', file_format='',
                    missions=[], model_param=[[]], i_x=0, i_y=1,
                    comp_type='pv', log_comp=True, likelihood=True, log=False, normalise=False,
                    contour=True, levels=10, fraction=0.1, aspect=20, cmap='viridis',
                    xlabel='', ylabel='', supxlabel='', supylabel='', clabel='', clabel_xa=0.98, clabel_ha='left',
                    title='', fontsize=20,
                    pad=1.0, pad_left=0, pad_right=0.125, pad_bottom=0, pad_top=0.12, wspace=0, hspace=0,
                    figsize=None, save_plot=False, output_file='', output_format='png', transparent=False,
                    verbose=False, debug=False):
    #

    # Check that the missions are specified properly.
    if missions == '' or missions == None or missions == []:
        missions = os.listdir(path)
    elif isinstance(missions, str):
        missions = [missions]
    elif isinstance(missions, list):
        pass
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

    dimensions = len(model_param)
    naxis = np.zeros((dimensions), dtype=bool)
    naxis[[i_x, i_y]] = True
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
        sub_y =1
    elif naxis.size == 4:
        sub_x, sub_y = lenaxis[~naxis]
    else:
        print('Too many dimensions in grid.')
        return

    # Initialise likelihood and figure for plot analysing all missions
    loglikelihood_overall = np.zeros((*x_grid.shape, sub_x, sub_y))
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

        survey_files = os.listdir(path + survey + '/')
        if 'Plots' in survey_files:
            survey_files.remove('Plots')
        else:
            os.mkdir(path + survey + '/Plots/')

        file_transitions = []
        file_transition_plots = []
        file_transition_log_likelihood = []

        if debug:
            print('survey files: {}\n'.format(survey_files))

        if not figsize:
            subgrid_aspect = sub_x/sub_y
            figsize = (subgrid_aspect*10, 10)
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

            log_likelihood = np.zeros(x_grid.shape)

            for survey_file in survey_files:

                if '.png' in survey_file or 'Plots' in survey_file:
                    continue

                if verbose:
                    print(survey_file)

                transitions = os.listdir(path + survey + '/' + survey_file + '/')

                transitions_skipped = []
                for t in copy(transitions):
                    if debug:
                        print(t)
                    if (len(os.listdir(path + survey + '/' + survey_file + '/' + t + '/')) < 369):
                        transitions.remove(t)
                        transitions_skipped.append(t)
                    if debug:
                        print(transitions, transitions_skipped)
                file_transitions.append(copy(transitions))
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
                        if naxis.size > 2:
                            file_format_param[~naxis] = param
                        if likelihood:
                            filename = file_format.format(*file_format_param) \
                                       + '_{}_loglikelihood.npy'.format(comp_type)
                        else:
                            filename = file_format.format(*file_format_param) \
                                       + '_{}_chi2.npy'.format(comp_type)
                            # filename = file_format.format(x_grid[i, j], y_grid[i, j]) \
                            #            + '_{}_loglikelihood.npy'.format(comp_type)
                        param_log_likelihood = [np.load(path + survey + '/' + survey_file + '/' +
                                                            t + '/' + filename) for t in transitions]
                        log_likelihood[i, j] = log_likelihood[i, j] + np.nansum(param_log_likelihood)
                        for t in range(len(transitions)):
                            transition_log_likelihood[t][i, j] = transition_log_likelihood[t][i, j] \
                                                                 + np.nansum(param_log_likelihood[t])

                file_transition_plots.append(copy(transition_plots))
                file_transition_log_likelihood.append(deepcopy(transition_log_likelihood))

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

            if contour:
                cm = axes[sub_indeces].contourf(log_likelihood, levels=levels, cmap=cmap)
            else:
                cm = axes[sub_indeces].imshow(log_likelihood, extent=[0, lenaxis[i_x], 0, lenaxis[i_y]], cmap=cmap)
            axes[sub_indeces].set_aspect(lenaxis[i_x]/lenaxis[i_y])
            cb = fig.colorbar(cm, ax=axes[sub_indeces], fraction=fraction, aspect=aspect)

            if contour:
                axes[sub_indeces].set_xticks(np.arange(lenaxis[i_x]))
                axes[sub_indeces].set_xticklabels([str(param) for param in model_param[i_x]])
                axes[sub_indeces].set_yticks(np.arange(lenaxis[i_y]))
                axes[sub_indeces].set_yticklabels([str(param) for param in model_param[i_y]])
            else:
                axes[sub_indeces].set_xticks(np.arange(lenaxis[i_x])+0.5)
                axes[sub_indeces].set_xticklabels([str(param) for param in model_param[i_x]])
                axes[sub_indeces].set_yticks(np.arange(lenaxis[i_y])+0.5)
                axes[sub_indeces].set_yticklabels([str(param) for param in model_param[i_y][::-1]])

            axes[sub_indeces].set_xlabel(xlabel, fontsize=fontsize-4)
            axes[sub_indeces].set_ylabel(ylabel, fontsize=fontsize-4)
            if axes.size > 1:
                axes[sub_indeces].set_title('{} {}, {} {}'.format(supylabel, y_param, supxlabel, x_param),
                                            fontsize=fontsize-4)
            # cb.ax.set_ylabel(clabel, fontsize=fontsize-4)

            for f in range(len(survey_files)):
                for t in range(len(file_transitions[f])):

                    if contour:
                        cm = file_transition_plots[f][t][1][sub_indeces].contourf(file_transition_log_likelihood[f][t],
                                                                                  levels=levels, cmap=cmap)
                    else:
                        cm = file_transition_plots[f][t][1][sub_indeces].imshow(file_transition_log_likelihood[f][t],
                                                                                extent=[0, lenaxis[i_x],
                                                                                        0, lenaxis[i_y]],
                                                                                cmap=cmap)
                    file_transition_plots[f][t][1][sub_indeces].set_aspect(lenaxis[i_x]/lenaxis[i_y])
                    cb = file_transition_plots[f][t][0].colorbar(cm, ax=file_transition_plots[f][t][1][sub_indeces],
                                                         fraction=fraction, aspect=aspect)

                    if contour:
                        file_transition_plots[f][t][1][sub_indeces].set_xticks(np.arange(lenaxis[i_x]))
                        file_transition_plots[f][t][1][sub_indeces].set_xticklabels([str(param) for
                                                                                     param in model_param[i_x]])
                        file_transition_plots[f][t][1][sub_indeces].set_yticks(np.arange(lenaxis[i_y]))
                        file_transition_plots[f][t][1][sub_indeces].set_yticklabels([str(param) for
                                                                                     param in model_param[i_y]])
                    else:
                        file_transition_plots[f][t][1][sub_indeces].set_xticks(np.arange(lenaxis[i_x])+0.5)
                        file_transition_plots[f][t][1][sub_indeces].set_xticklabels([str(param) for
                                                                                     param in model_param[i_x]])
                        file_transition_plots[f][t][1][sub_indeces].set_yticks(np.arange(lenaxis[i_y])+0.5)
                        file_transition_plots[f][t][1][sub_indeces].set_yticklabels([str(param) for
                                                                                     param in model_param[i_y][::-1]])

                    file_transition_plots[f][t][1][sub_indeces].set_xlabel(xlabel, fontsize=fontsize-4)
                    file_transition_plots[f][t][1][sub_indeces].set_ylabel(ylabel, fontsize=fontsize-4)
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
        fig.supylabel(clabel, x=clabel_xa, ha=clabel_ha, fontsize=fontsize)

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
                output_file = file_format
            if likelihood:
                plt.savefig(path + survey + '/Plots/{}_{}_loglikelihood.{}'.format(output_file, comp, output_format),
                            format=output_format, transparent=transparent)
            else:
                plt.savefig(path + survey + '/Plots/{}_{}_chi2.{}'.format(output_file, comp, output_format),
                            format=output_format, transparent=transparent)
        else:
            plt.show()

        plt.close(fig)

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
                file_transition_plots[f][t][0].supylabel(clabel, x=clabel_xa, ha=clabel_ha, fontsize=fontsize)

                if suptitle:
                    fig_top = 1 - pad*pad_top
                if clabel:
                    fig_right = 1 - pad*pad_right
                file_transition_plots[f][t][0].subplots_adjust(left=pad*pad_left, right=fig_right,
                                                               bottom=pad*pad_bottom, top=fig_top,
                                                               wspace=wspace, hspace=hspace)
                # transition_plots[t][0].tight_layout(pad=pad)

                plt.figure(file_transition_plots[f][t][0])
                if save_plot:
                    if output_file == None or output_file == '':
                        output_file = file_format
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

        if contour:
            cm = axes_overall[sub_indeces].contourf(loglikelihood_overall[:, :, sub_indeces[0], sub_indeces[1]],
                                                    levels=levels, cmap=cmap)
        else:
            cm = axes_overall[sub_indeces].imshow(loglikelihood_overall[:, :, sub_indeces[0], sub_indeces[1]],
                                                  extent=[0, lenaxis[i_x], 0, lenaxis[i_y]], cmap=cmap)
        axes_overall[sub_indeces].set_aspect(lenaxis[i_x] / lenaxis[i_y])
        cb = fig_overall.colorbar(cm, ax=axes_overall[sub_indeces], fraction=fraction, aspect=aspect)

        if contour:
            axes_overall[sub_indeces].set_xticks(np.arange(lenaxis[i_x]))
            axes_overall[sub_indeces].set_xticklabels([str(param) for param in model_param[i_x]])
            axes_overall[sub_indeces].set_yticks(np.arange(lenaxis[i_y]))
            axes_overall[sub_indeces].set_yticklabels([str(param) for param in model_param[i_y]])
        else:
            axes_overall[sub_indeces].set_xticks(np.arange(lenaxis[i_x]) + 0.5)
            axes_overall[sub_indeces].set_xticklabels([str(param) for param in model_param[i_x]])
            axes_overall[sub_indeces].set_yticks(np.arange(lenaxis[i_y]) + 0.5)
            axes_overall[sub_indeces].set_yticklabels([str(param) for param in model_param[i_y][::-1]])

        axes_overall[sub_indeces].set_xlabel(xlabel, fontsize=fontsize - 4)
        axes_overall[sub_indeces].set_ylabel(ylabel, fontsize=fontsize - 4)
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
    fig_overall.supylabel(clabel, x=clabel_xa, ha=clabel_ha, fontsize=fontsize)

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
            output_file = file_format
        if likelihood:
            plt.savefig(path + '/Plots/{}_{}_loglikelihood.{}'.format(output_file, comp, output_format),
                        format=output_format, transparent=transparent)
        else:
            plt.savefig(path + '/Plots/{}_{}_chi2.{}'.format(output_file, comp, output_format),
                        format=output_format, transparent=transparent)
    else:
        plt.show()

    plt.close(fig_overall)

    return
