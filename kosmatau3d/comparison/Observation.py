import numpy as np
import os

from astropy.io import fits
from scipy.io import readsav

class Observation(object):
    '''
    This is an object to load individual observations. This is merely for
    the convenience of opening all of the observations in a consistent manner.
    There is an optional argument when initialising to set a base directory, which
    makes it easier to load multiple models in succession.
    '''

    cobe_idl_linfrq = np.array([115.3, 230.5, 345.8, 424.8, 461.0, 492.2, 556.9, 576.3, 691.5, 808.1,
                                1113, 1460, 2226, 1901, 2060, 2311, 2459, 2589, 921.8])
    cobe_idl_transitions = np.array(['CO 1', 'CO 2', 'CO 3', 'CO 4', 'C 3', 'CO 5',
                                     'CO 6', 'CO 7 + C 1', 'C+ 1', 'O 1', 'CO 8'])
  
    def __init__(self, base_dir='', regridded_dir=''):
        '''
        This initialises the object along with the base directory.
        The owned objects of `base_dir` and `files` are created.
        `files` can be modified again when loading a model, but for now it
        has the default filenames created with `kosmatau3d`.
  
        :param base_dir: the base directory to use when loading observations. Default: `''`.
  
  
        '''
  
        self.base_dir = base_dir
        self.regridded_dir = regridded_dir
        self.obs = []
        self.obs_header = []
        self.obs_data = []
        self.obs_lon = []
        self.obs_lat = []
        self.obs_vel = []
        self.obs_iid = []
        self.obs_wavelength = []
        self.obs_frequency = []
        # self.files = {'intensity': 'synthetic_intensity',
        #               'optical_depth': 'synthetic_optical_depth',
        #               'dust_absorption': 'dust_absorption',
        #               'dust_emissivity': 'dust_emissivity',
        #               'species_absorption': 'species_absorption',
        #               'species_emissivity': 'species_emissivity',
        #               'density': 'voxel_density',
        #               'ensemble_dispersion': 'voxel_ensemble_dispersion', 
        #               'ensemble_mass': 'voxel_ensemble_mass', 
        #               'fuv_absorption': 'voxel_FUVabsorption', 
        #               'fuv': 'voxel_fuv', 
        #               'position': 'voxel_position', 
        #               'velocity': 'voxel_velocity', 
        #               'los_count': 'sightlines', 
        #               'log': 'log', }
          
        return

    def load_survey(self, directory=None, survey=None):
        '''
        '''

        if survey is None:
            print('ERROR: Please specify an observational survey')

        self.survey = survey

        if os.path.exists(self.base_dir + directory + self.survey):
            pass
        else:
            print(f'ERROR: Survey {self.base_dir + directory + survey} does not exist. Ensure '
                  + 'that base_dir was initialised correctly and that you have the correct survey')
            return

        self.directory = directory

        if os.path.exists(self.base_dir + directory + survey + self.regridded_dir):
            full_path = self.base_dir + self.directory + self.survey + self.regridded_dir
        else:
            print('ERROR: Either the survey has not been regridded or the directory differs from '
                  + f'{self.regridded_dir}.')
            return

        self.files = list(f.name for f in os.scandir(full_path) 
                          if (f.is_file() and not '_error' in f.name))
        for f in self.files:
            if '.idl' in f:
                self.obs.append(readsav(full_path + f))
                self.obs_header.append[None]
                self.obs_iid.append(deepcopy(cobe_idl_transitions)
                self.obs_frequency.append(np.array(cobe_idl_linfrq)*1e9)
                self.obs_wavelength.append(299792458/self.obs_frequency[-1])
                self.obs_data.append(obs['amplitude']/(2.9979**3)*(self.obs_frequency**3)*2*1.38/10**-1)
                self.obs_lons.append(obs['long'])
                self.obs_lats.append(np.asarray([0]))
            else:
                self.obs.append(fits.open(full_path + f))
                self.obs_header.append[self.obs[0].header]
                self.obs_data.append(self.obs[-1][0].data)
                self.obs_lon.append(np.linspace(obs[0].header['CRVAL1']
                                                - obs[0].header['CDELT1']*(obs[0].header['CRPIX1']-1),
                                                obs[0].header['CRVAL1'] 
                                                + obs[0].header['CDELT1']*(obs[0].header['NAXIS1']
                                                                           -obs[0].header['CRPIX1']),
                                                num=obs[0].header['NAXIS1']))
                self.obs_lat.append(np.linspace(obs[0].header['CRVAL2']
                                                - obs[0].header['CDELT2']*(obs[0].header['CRPIX2']-1),
                                                obs[0].header['CRVAL2'] 
                                                + obs[0].header['CDELT2']*(obs[0].header['NAXIS2']
                                                                           -obs[0].header['CRPIX2']),
                                                num=obs[0].header['NAXIS2']))
                if self.obs_header[-1]['NAXIS'] == 3 and not self.survey == 'COBE-FIRAS':
                    self.obs_vel.append(np.linspace((self.obs_header[-1]['CRVAL3']
                                                     - self.obs_header[-1]['CDELT3']
                                                     *(self.obs_header[-1]['CRPIX3']-1)),
                                                    (self.obs_header[-1]['CRVAL3'] 
                                                     + self.obs_header[-1]['CDELT3']
                                                     *(self.obs_header[-1]['NAXIS3']
                                                       -self.obs_header[-1]['CRPIX3'])),
                                                    num=self.obs_header[-1]['NAXIS3']))
                self.obs_iid.append(self.obs_header[-1])
