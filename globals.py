# Identify the working directories needed
import os
import inspect
filename = inspect.getframeinfo(inspect.currentframe()).filename
KOSMAPATH = os.path.dirname(filename)+'/'
COMPOUNDPATH = KOSMAPATH+'compounds/disk/'
ORIONPATH = KOSMAPATH+'compounds/Orion_Bar/'
FITSPATH = KOSMAPATH+'fits_file/data/'
TEMPPATH = KOSMAPATH+'temp/'
INPUTPATH = KOSMAPATH+'input_files/'
# Initialise a global dictionary to be used throughout the code.
_globals = {}
velocities = []
plotFlags = {           'Total mass':  False,   # 1 shows plot
                        'Clump mass':  False,
                   'Interclump mass':  False,   # if mass is there
                   'Total intensity':  False,
                   'Clump intensity':  False,
              'Interclump intensity':  False,
                           'Density':  False,
                        'Z velocity':  False,
                    'Total velocity':  False,
                        'Dispersion':  False,
                  'Total dispersion':  False,
                  'FUV distribution':  False}
