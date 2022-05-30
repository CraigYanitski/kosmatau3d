from .model_selection import *
from .Observation import *


def help():
      print('This submodule is the location of the model comparison methods I wrote to '
            'compare the kosmatau3d models to observational data. It regrids the '
            'observations to a common resolution, calculates the observation error, '
            'and performs multiple goodness-of-fit plots. These can be overall '
            'model grid likelihood, chi-squared, or line ratios.')
      return
