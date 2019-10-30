import os
from shutil import copyfile
import numpy as np

directory_name = str(raw_input('name of new directory: '))

path = '../results_paper2_version2/' + directory_name
try: os.mkdir(path)
except OSError, e:
    if e.errno == 17 and os.path.isdir(path): pass
    # errno == 17 : path exists
    # os.path.isdir returns true if path is en exisiting directory

path = '../results_paper2_version2/' + directory_name + '/input/'
try: os.mkdir(path)
except OSError, e:
    if e.errno == 17 and os.path.isdir(path): pass
copyfile('./compounds/Orion_Bar/model.py', path + 'model.py')
copyfile('./compounds/Orion_Bar/model_FUV.py', path + 'model_FUV.py')

path = '../results_paper2_version2/' + directory_name + '/maps/'
try: os.mkdir(path)
except OSError, e:
    if e.errno == 17 and os.path.isdir(path): pass
os.system('cp ./maps/*.dat' + ' ' + path)

path = '../results_paper2_version2/' + directory_name + '/spectra/'
try: os.mkdir(path)
except OSError, e:
    if e.errno == 17 and os.path.isdir(path): pass
os.system('cp ./spectra/*.dat' + ' ' + path)

path = '../results_paper2_version2/' + directory_name + '/temp/'
try: os.mkdir(path)
except OSError, e:
    if e.errno == 17 and os.path.isdir(path): pass
os.system('cp ./temp/*.dat' + ' ' + path)

path = '../results_paper2_version2/' + directory_name + '/pics/'
try: os.mkdir(path)
except OSError, e:
    if e.errno == 17 and os.path.isdir(path): pass
os.system('cp ./pics/*' + ' ' + path)




print 'done'

