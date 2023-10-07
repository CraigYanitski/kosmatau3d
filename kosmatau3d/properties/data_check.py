'''
A module to check the structure of the KOSMA-tau grid files.
'''

import numpy as np
import pandas as pd


def grid_data_check(file=''):
    '''
    Parse and return parameters used in KOSMA-tau grid from header.
    '''
    with open(file) as gd:
        header = gd.readline()
    params = header.split(': ')[1].split(', molecules')[0].split(', ')
    gd = np.genfromtxt(file)
    param_vals = {}
    for i, p in enumerate(params):
        param_vals[p] = np.unique(gd[:, i])
    return param_vals

def grid_data_limits(file=''):
    '''
    Return limits of KOSMA-tau grid parameters.
    '''
    params = grid_data_check(file)
    for p in params.keys():
        print(p + ' range: 10^{:.2f} -> 10^{:.2f}'.format(params[p].min()/10, params[p].max()/10))
    return

