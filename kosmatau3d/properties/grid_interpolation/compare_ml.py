import inspect
import matplotlib.pyplot as plt
import numpy as np
import os

from scipy.interpolate import LinearNDInterpolator
from sklearn.tree import DecisionTreeRegressor, plot_tree
from sklearn.ensemble import RandomForestRegressor, ExtraTreesRegressor


class CompareInterpolation():

    def __init__(self, tmb_file='clump_Tmb_LineCenter.dat', tau_file='clump_tau_LineCenter.dat', 
                 n_param=4):
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

    def interpolate(self, transition, full_grid=False, all_species=False):
        
        if all_species:
            i_species = np.arange(self.species.size)
        else:
            i_species = self.species.index(transition)

        self.tmb_orig = []
        self.tmb_interp_lin = []
        self.tmb_interp_mb = []
        self.tau_orig = []
        self.tau_interp_lin = []
        self.tau_interp_ml = []

        for i, params in enumerate(self.tmb_data[0]):

            if full_grid:
                idx = np.ones(self.tmb_data[0].shape[0], dtype=bool)
            else:
                idx = np.all([self.tmb_data[0][:, p]==params[p] 
                              for p in range(len(params))], axis=1)
            
            # Interpolate intensity
            self.tmb_orig.append(self.tmb_data[1][idx, i_species])
            # - linear
            lin_interp = LinearNDInterpolator(self.tmb_data[0][~idx, :], 
                                              self.tmb_data[1][:, i_species][~idx])
            self.tmb_interp_lin.append(lin_interp(params))
            # - ML
            ml_interp = ExtraTreesRegressor(random_state=0)
            ml_interp.fit(self.tmb_data[0][~idx, :], self.tmb_data[1][:, i_species][~idx])
            self.tmb_interp_ml.append(ml_interp.predict(np.asarray(params).reshape(1, -1)))
            
            # Interpolate optical depth
            self.tau_orig.append(self.tau_data[1][idx, i_species])
            # - linear
            lin_interp = LinearNDInterpolator(self.tau_data[0][~idx, :], 
                                              self.tau_data[1][:, i_species][~idx])
            self.tau_interp_lin.append(lin_interp(params))
            # - ML
            ml_interp = ExtraTreesRegressor(random_state=0)
            ml_interp.fit(self.tau_data[0][~idx, :], self.tau_data[1][:, i_species][~idx])
            self.tau_interp_ml.append(ml_interp.predict(np.asarray(params).reshape(1, -1)))
        
        return

    def plot_result(self, ax=None):
        if not ax:
            fig, ax = plt.subplots(figsize=(10, 10))
        ax.semilogy(np.arange(self.tmb_data[0].shape[0]), np.asarray(self.tmb_interp_lin)/np.asarray((self.tmb_orig)), ls='-', color='xkcd:maroon')
        ax.semilogy(np.arange(self.tmb_data[0].shape[0]), np.asarray(self.tmb_interp_ml)/np.asarray((self.tmb_orig)), ls='--', color='xkcd:maroon')
        ax.semilogy(np.arange(self.tau_data[0].shape[0]), np.asarray(self.tau_interp_lin)/np.asarray((self.tau_orig)), ls='-', color='xkcd:sapphire')
        ax.semilogy(np.arange(self.tau_data[0].shape[0]), np.asarray(self.tau_interp_ml)/np.asarray((self.tau_orig)), ls='--', color='xkcd:sapphire')
        ax.set_xlabel('Parameter', fontsize=24)
        ax.set_ylabel('Relative error', fontsize=24)
        return ax
