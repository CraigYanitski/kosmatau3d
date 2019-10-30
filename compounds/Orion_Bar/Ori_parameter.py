
######### line and observation specific parameters ########
import numpy as np
from PDR import _globals
### Numbering and naming for species and transition as defined in
# SetUpOrionBarModelEnvironment.v.1.1.nb ########################
# _globals['compound']['species'] = ['C+', 'CO']
_globals['compound']['species'] = ['C+', 'CO', 'CO', 'CO', 'CO', 'CO', '13CO', '13CO', '13CO', '13CO', 'HCO+', 'HCO+']

# List of Species (string). possible choices are: 
# C+, C, O, CO, 13CO, 13C+, 13C, HCO+, H13CO+, H3O+, (CH+)
# In all files and tabels species are arranged in this order
# CH+ does not yet work (due to missing tau-information)

_globals['compound']['transition'] = [1, 2, 3, 6, 10, 16, 3, 5, 6, 10, 3, 6]
# _globals['compound']['transition'] = [1, 2]
# provide transitions (for example 10 indicats transition (10->9))
# lines per species: 1, 3, 3, 49, 49, 1, 3, 15, 30, 17, 10
# List needs to have the same length as "_globals['compound']['species']"
_globals['compound']['dbeam'] = [11.2, 30.5, 21.9, 10.6, 18.4, 11.5, 21.9, 38.5, 10.6, 19.3, 30.5, 39.6] 
# _globals['compound']['dbeam'] = [11.2, 30.5]
# List of Beamsizes. same Length again
_globals['compound']['names'] = ['[CII]', 'CO 2-1', 'CO 3-2', 'CO 6-5', 'CO 10-9', 'CO 16-15', '$^{13}$CO 3-2', '$^{13}$CO 5-4', '$^{13}$CO 6-5', '$^{13}$CO 10-9', 'HCO$^+$ 3-2', 'HCO$^+$ 6-5']
# species +  transitions dencently set for plot axes ...only needed for nice 
# pics, can contain any string 
_globals['compound']['g_sigma'] = 3 
# g_sigma * FWHM(beam) will be accounted for when creating the 
# gauss kernel for the convolution of the map

_globals['statistic'] = {}
_globals['statistic']['Nmin_clumps'] = 1
_globals['statistic']['Nmin_inter'] = 100
# minimum number of large clumps for clump and interclump medium.
# only relevant for binomial option


##### sanity check##########################################################

nspe1 = len(_globals['compound']['species'])
nspe2 = len(_globals['compound']['transition'])
nspe3 = len(_globals['compound']['dbeam'])
if nspe1 == nspe2 and  nspe1 == nspe3: pass
else: 
    import sys    
    sys.exit('Lists of species, transition numbers and beamsizes do not have the same length ...exiting...')
_globals['compound']['nspe'] = nspe1

for i in np.arange(_globals['compound']['npoints']):
    _globals['compound']['ens'][i].epsilon = np.zeros(nspe1)
    _globals['compound']['ens'][i].epsilon_tot = np.zeros(nspe1)
    _globals['compound']['ens'][i].kappa = np.zeros(nspe1)
    _globals['compound']['ens'][i].kasppa_tot = np.zeros(nspe1)

print 'choosen species: ', _globals['compound']['species']
print 'choosen transitions: ', _globals['compound']['transition']
print 'beamsizes: ', _globals['compound']['dbeam']
pause = raw_input('...ok?...')

#########################################################
# C+, 13C+, O, C, 13C, CO, 13CO, HCO+, H13CO+, H3O+ ????
# 1, 1, 3, 3, 3, 49, 49, 15, 30, 17 ????
# if map in Kkm/s is needed, the frequency is necessary for conversion
# READ IN ONLY TRAINSITION would of course be more elegant!!!!
#_globals['compound']['frequency'] = 550.926 # 13CO(5-4) in GHz
#_globals['compound']['frequency'] = 1101.35 # 13CO(10-9) in GHz
#_globals['compound']['frequency'] = 1900.537 # CII in GHz
