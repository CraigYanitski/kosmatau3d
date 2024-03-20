'''
This is a subpackage containing a class to compare methods of interpolating the
KOSMA-:math:`\\tau` grid.
'''

import inspect
import matplotlib.pyplot as plt
import numpy as np
import os

from copy import copy
from scipy.interpolate import LinearNDInterpolator
from sklearn.ensemble import ExtraTreesRegressor


class CompareInterpolation():
    '''
    This is a class to facilitate the comparison of interpolation methods.
    '''

    def __init__(self, tmb_file='clump_Tmb_LineCenter.dat', tau_file='clump_tau_LineCenter.dat',
                 n_param=4):
        '''
        Initialise attributes and open grid files.
        '''

        filename = inspect.getframeinfo(inspect.currentframe()).filename
        package_path = os.path.abspath(os.path.dirname(filename)+'/../../')  # kosmatau3d dir
        self.path = package_path + '/grid/'  # grid dir
        self.tmb_file = tmb_file
        self.tau_file = tau_file
        self.open_files(n_param=n_param)
        self.tmb_orig = []
        self.tmb_interp_lin = []
        self.tmb_interp_ml = []
        self.tau_orig = []
        self.tau_interp_lin = []
        self.tau_interp_ml = []

        return

    def open_files(self, n_param=4):
        '''
        Parse header of grid files and load data into memory.
        '''

        header = []
        with open(self.path+self.tmb_file) as tmb:
            header.append(tmb.readline())
            header.append(tmb.readline())
        molecules = header[1].split(': ')[1]
        self.species = []
        for molecule in molecules.split(', '):
            for transition in np.arange(1, int(molecule.split(' ')[1])+1):
                self.species.append('{} {}'.format(molecule.split(' ')[0], transition))
        tmb = np.genfromtxt(self.path+self.tmb_file)
        self.tmb_data = (tmb[:, :n_param], tmb[:, n_param:])
        tau = np.genfromtxt(self.path+self.tau_file)
        self.tau_data = (tau[:, :n_param], tau[:, n_param:])

        return

    def interpolate(self, transition=None, full_grid=False, all_species=False, savedir='', verbose=False):
        '''
        Perform interpolation of model grid at each point in multiple ways.
        Currently this performs a linear interpolation as well as predicting
        the interpolated point using extremely-randomised trees.
        '''

        if not transition:
            transition = self.species
        elif isinstance(transition, (str, np.ndarray, tuple)):
            transition = list(transition)

        self.tmb_orig = []
        self.tmb_interp_lin = []
        self.tmb_interp_ml = []
        self.tau_orig = []
        self.tau_interp_lin = []
        self.tau_interp_ml = []
        if full_grid:
            print('full grid')
            tmb_data = (self.tmb_data[0], self.tmb_data[1])
            tau_data = (self.tau_data[0], self.tau_data[1])
        else:
            print('partial grid')

        for i, params in enumerate(self.tmb_data[0]):

            if verbose:
                print(params)#, end='\r')

            tmb_orig = []
            tmb_interp_lin = []
            tmb_interp_ml = []
            tau_orig = []
            tau_interp_lin = []
            tau_interp_ml = []
            idx = np.all([self.tmb_data[0][:, p]==params[p]
                          for p in range(len(params))], axis=0)

            if not full_grid:
                tmb_data = (copy(self.tmb_data[0][~idx]), copy(self.tmb_data[1][~idx]))
                tau_data = (copy(self.tau_data[0][~idx]), copy(self.tau_data[1][~idx]))

            if all_species:
                i_species = np.arange(len(self.species))
                # print(tmb_data[1][:, i_species].shape)
                # print(tau_data[1][:, i_species].shape)
                tmb_orig = self.tmb_data[1][idx]
                tau_orig = self.tau_data[1][idx]
                lin_interp = LinearNDInterpolator(tmb_data[0],
                                                  tmb_data[1][:, i_species])
                tmb_interp_lin.append(lin_interp(params))
                ml_interp = ExtraTreesRegressor(random_state=0)
                ml_interp.fit(tmb_data[0], tmb_data[1][:, i_species])
                tmb_interp_ml.append(ml_interp.predict(np.asarray(params).reshape(1, -1)))
                lin_interp = LinearNDInterpolator(tau_data[0],
                                                  tau_data[1][:, i_species])
                tau_interp_lin.append(lin_interp(params))
                ml_interp = ExtraTreesRegressor(random_state=0)
                ml_interp.fit(tau_data[0], tau_data[1][:, i_species])
                tau_interp_ml.append(ml_interp.predict(np.asarray(params).reshape(1, -1)))

            for species in transition:

                if all_species:
                    continue

                if verbose:
                    print(species + '\t\t', end='\r')

                i_species = self.species.index(species)

                # Interpolate intensity
                tmb_orig.append(self.tmb_data[1][idx, i_species])
                # - linear
                lin_interp = LinearNDInterpolator(tmb_data[0],
                                                  tmb_data[1][:, i_species])
                tmb_interp_lin.append(lin_interp(params))
                # - ML
                ml_interp = ExtraTreesRegressor(random_state=0)
                ml_interp.fit(tmb_data[0], tmb_data[1][:, i_species])
                tmb_interp_ml.append(ml_interp.predict(np.asarray(params).reshape(1, -1)))

                # Interpolate optical depth
                tau_orig.append(self.tau_data[1][idx, i_species])
                # - linear
                lin_interp = LinearNDInterpolator(tau_data[0],
                                                  tau_data[1][:, i_species])
                tau_interp_lin.append(lin_interp(params))
                # - ML
                ml_interp = ExtraTreesRegressor(random_state=0)
                ml_interp.fit(tau_data[0], tau_data[1][:, i_species])
                tau_interp_ml.append(ml_interp.predict(np.asarray(params).reshape(1, -1)))

            if verbose:
                print('\r')

            self.tmb_orig.append(tmb_orig)
            self.tmb_interp_lin.append(tmb_interp_lin)
            self.tmb_interp_ml.append(tmb_interp_ml)
            self.tau_orig.append(tau_orig)
            self.tau_interp_lin.append(tau_interp_lin)
            self.tau_interp_ml.append(tau_interp_ml)

        self.tmb_orig = np.asarray(self.tmb_orig)
        self.tmb_interp_lin = np.asarray(self.tmb_interp_lin)
        self.tmb_interp_ml = np.asarray(self.tmb_interp_ml)
        self.tau_orig = np.asarray(self.tau_orig)
        self.tau_interp_lin = np.asarray(self.tau_interp_lin)
        self.tau_interp_ml = np.asarray(self.tau_interp_ml)

        self.save_results(directory=savedir, full_grid=full_grid)
        print(f'\nsaved comparison to {savedir}')

        return

    def save_results(self, directory='', full_grid=True):
        '''
        Save original and interpolated data in separate numpy binary files.
        '''
        if directory:
            savedir = (directory+'/') if not (directory[-1]=='/') else directory
        else:
            os.makedirs(self.path + 'comparison/')
            savedir = self.path + 'comparison/'
        if full_grid:
            suffix = 'full'
        else:
            suffix = 'partial'
        np.save(f"{savedir+self.tmb_file.split('.')[0]}_tmb_orig_{suffix}.npy", self.tmb_orig)
        np.save(f"{savedir+self.tmb_file.split('.')[0]}_tmb_interp_lin_{suffix}.npy", self.tmb_interp_lin)
        np.save(f"{savedir+self.tmb_file.split('.')[0]}_tmb_interp_ml_{suffix}.npy", self.tmb_interp_ml)
        np.save(f"{savedir+self.tau_file.split('.')[0]}_tau_orig_{suffix}.npy", self.tau_orig)
        np.save(f"{savedir+self.tau_file.split('.')[0]}_tau_interp_lin_{suffix}.npy", self.tau_interp_lin)
        np.save(f"{savedir+self.tau_file.split('.')[0]}_tau_interp_ml_{suffix}.npy", self.tau_interp_ml)

        return

    def plot_result(self, ax=None, transition=None, **kwargs):
        '''
        Plot ratio of the interpolated and original data as a function of parameter.
        '''
        if not transition:
            warnings.error('Please specify a transition')
            exit()
        if not ax:
            fig, ax = plt.subplots(figsize=(10, 10))
        i_species = self.species.index(transition)
        idx = np.where((self.tmb_orig[:, i_species] >= 1e-16) | (self.tau_orig[:, i_species] >= 1e-16))[0]
        ix = np.ix_(idx, i_species)
        ax.semilogy(np.arange(self.tmb_data[0].shape[0]), (self.tmb_interp_lin/self.tmb_orig)[ix], ls='-', color='xkcd:maroon', **kwargs)
        ax.semilogy(np.arange(self.tmb_data[0].shape[0]), (self.tmb_interp_ml/self.tmb_orig)[ix], ls='--', color='xkcd:maroon', **kwargs)
        ax.semilogy(np.arange(self.tau_data[0].shape[0]), (self.tau_interp_lin/self.tau_orig)[ix], ls='-', color='xkcd:sapphire', **kwargs)
        ax.semilogy(np.arange(self.tau_data[0].shape[0]), (self.tau_interp_ml/self.tau_orig)[ix], ls='--', color='xkcd:sapphire', **kwargs)
        ax.set_xlabel('Parameter', fontsize=24)
        ax.set_ylabel('Relative error', fontsize=24)
        return ax
