'''
The :code:`comparison` subpackage is useful both to load the observational data used for comparison
and to perform a comparison to :code:`kosmatau3d` models.
It is currently setup to compare to the Milky Way models from Yanitski (2023), but it will
be extended soon to arbitrary model setup.
'''

from .model_selection import *
from .observation import *


def help():
      print('This submodule is the location of the model comparison methods I wrote to '
            'compare the kosmatau3d models to observational data. It regrids the '
            'observations to a common resolution, calculates the observation error, '
            'and performs multiple goodness-of-fit plots. These can be overall '
            'model grid likelihood, chi-squared, or line ratios.')
      return
