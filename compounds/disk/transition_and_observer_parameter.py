######### line and observation specific parameters ########
import numpy as np
#from PDR import gbl._globals
# edit: Craig, 14.10.2019
import globals as gbl
#end edit
### Numbering and naming for species and transition:

#names = ['C+_1-0', 'C_1-0', 'C_2-1', 'CO_1-0', 'CO_2-1', 'CO_3-2', 'CO_4-3', 'CO_5-4', 'CO_6-5'] #firas and Dust1-22
#selected_dust = [1,6,10,12,16,21] #Dust interpol
#gbl._globals['compound']['species'] = ['13CO','C+', 'CO']
gbl._globals['compound']['species'] = ['13CO','C+','CO','13CO','O','Dust','Dust','Dust','Dust','Dust']
#gbl._globals['compound']['species'] = ['C+', 'C', 'C', 'CO', 'CO', 'CO', 'CO', 'CO', 'CO', 'Dust', 'Dust', 'Dust', 'Dust', 'Dust', 'Dust']
#gbl._globals['compound']['species'] = ['Dust', 'Dust', 'Dust', 'Dust', 'Dust', 'Dust', 'Dust', 'Dust', 'Dust', 'Dust', 'Dust', 'Dust', 'Dust', 'Dust', 'Dust', 'Dust', 'Dust', 'Dust', 'Dust', 'Dust', 'Dust', 'Dust']
#gbl._globals['compound']['species'] = ['Dust']
# List of Species (string). possible choices are:
# C+, C, O, CO, 13CO, 13C+, 13C, HCO+, H13CO+, H3O+, (CH+) and Dust
# In all files and tabels species are arranged in this order
# CH+ does not yet work (due to missing tau-information)

#gbl._globals['compound']['transition'] = [10, 1, 1]
gbl._globals['compound']['transition'] = [10, 1, 1, 1, 2, 1, 2, 3, 4, 5]
#gbl._globals['compound']['transition'] = [1, 1, 2, 1, 2, 3, 4, 5, 6, 1, 6, 10, 12, 16, 21]
#gbl._globals['compound']['transition'] = [1,2]
#gbl._globals['compound']['transition'] = range(1,23) #22 dust transitions#[1, 1]#
#gbl._globals['compound']['transition'] = [1]

# provide transitions (for example 10 indicats transition (10->9))
#gbl._globals['compound']['transition'] = [1,2,3,4,5,6,7,8,1,10] ###CO-LINES
# lines per species: 1, 3, 3, 49, 49, 1, 3, 15, 30, 17, 10
# List needs to have the same length as "gbl._globals['compound']['species']"

gbl._globals['compound']['g_sigma'] = 3
# g_sigma * FWHM(beam) will be accounted for when creating the
# gauss kernel for the convolution of the map

gbl._globals['statistic'] = {}
gbl._globals['statistic']['Nmin_clumps'] = 1
gbl._globals['statistic']['Nmin_inter'] = 100
# minimum number of large clumps for clump and interclump medium.
# only relevant for binomial option


##### sanity check##########################################################

nspe1 = len(gbl._globals['compound']['species'])
nspe2 = len(gbl._globals['compound']['transition'])
#nspe3 = len(gbl._globals['compound']['dbeam'])
if nspe1 == nspe2 : pass
else:
    import sys
    sys.exit('Lists of species and transition numbers do not have the same length ...exiting...')
gbl._globals['compound']['nspe'] = nspe1


print('choosen species: ', gbl._globals['compound']['species'])
print('choosen transitions: ', gbl._globals['compound']['transition'])
#print 'beamsizes: ', gbl._globals['compound']['dbeam']
#pause = input('...ok?...')
