# KOSMA-tau_3D
# INPUT: model.py, model_FUV.py, Ori_parameter.py
#
# COMMENTS / THINGS TO BE IMPROVED
# MASER LINES:
#   At the moment negative kappa (i.e. maser lines) are set to zero.
#   The radiative transfer part is capable of treating WEAK (tau<0.1) maser
#   lines (MASER lines are treated in linear approximation). 
#   however, the interpolation on the KOSMA-tau grid does not work for
#   negative values as it is done in "log-space"

# RE-GRID
#   currently the code is not regridding. 
#   regrid will not work for the new version of the code with
#   velocity dependent intensities and optical depth
#   hence the reference line of sight has to go through 0/0.
#   (as translation works parallel to z should be sufficient. TRY OUT)
#   if rotation of grid is ever (re-)included it needs to go before 
#   line_absorption_and_emission*.py, because the velocity component 
#   in z-direction is needed there
#   regrid need to be performend AFTER line_absorption_and_emission*.py

######################################################################
# iclude Ori_parametes into model.py?

# edit: Craig, 14.10.2019
import globals as gbl
import sys
sys.path.append(gbl.KOSMAPATH)
sys.path.append(gbl.COMPOUNDPATH)
#end edit
#sys.path.insert(1,'./compounds/disk')
# Do you want an exhaustive output of each menial calculation??
gbl._globals['verbose'] = False
# add directorys containing prstoject files to sys.path variable
gbl._globals['Files'] = {}
gbl._globals['Files']['compound_filename'] = 'model' 
gbl._globals['Files']['parameter_filename'] = 'transition_and_observer_parameter'

#reuse_FUV = 1 # (=1): use FUV field for each pixel which has been calculated 
# in the previous simulation (read in temp/FUV.dat). can only be used of geometry has not changed
#reuse_kappa_epsilon = 0
# enter name of parameter file
gbl._globals['Files']['fuv_filename'] = 'model_FUV' #
# enter name of file which calculates the FUV intensity for the ensembles.
# not nessesary if reuse_FUV != 1

# give the element from the velocity interval which gives 
#   the velocity at which the line intensities are plotted. i.e. for 21 elements
#   (counting from 0 to 20) 10 is the central element. not needed for FUV plot

# Set to True to use existing .fits files
gbl._globals['Flags'] = {}
gbl._globals['Flags']['reuse_kappa_epsilon'] = False
gbl._globals['Flags']['reuse_mass'] = False
gbl._globals['Flags']['reuse_clump'] = False
gbl._globals['Flags']['reuse_inter'] = False
gbl._globals['Flags']['reuse_vz'] = False
gbl._globals['Flags']['reuse_density'] = False
gbl._globals['Flags']['reuse_FUV'] = False
gbl._globals['Flags']['reuse_disp'] = False

# read in fits files
from astropy.io import fits
if gbl._globals['Flags']['reuse_mass']:
  a=1
  #hdulist = fits.open('./fits_file/data/mass.fits')
  #data = hdulist[0].data

# IMPORTS
from tkinter import *
from copy import copy
import datetime
import time
import os
import numpy as np
norm = np.linalg.norm

import PDRfunctions as pf  # collection of self-written routines
from astropy.io import fits  # for fit files


#from intersection import intersection
#from signum import signum
#from mod import mod
#from arrayops import array_equal
#from makeList import makeList
  
# Modify global dictionary
gbl._globals['plots'] = {}  # parameters for plots

gbl._globals['runtimes'] = {}  #for runtime control

# Define Constants
gbl._globals['constants'] = {}
gbl._globals['constants']['pc'] = 3.0856776 * 10 **(18)  #pc in cm

gbl._globals['constants']['kb'] = 1.3806 * 10 **(-23)  #Boltzmann constant (SI units)

gbl._globals['constants']['c'] =  2.998*10**8  #speed of light (m/s)

gbl._globals['constants']['alpha'] = None  #exponent for mass spectrum. will be overwrittern 
                                           #by alpha_cl and alpha_inter respectively

gbl._globals['constants']['alpha_cl'] = 1.84  #exponent for mass spectrum 

gbl._globals['constants']['alpha_intercl'] = 1.84  #exponent for mass spectrum, for interclump medium. only necessary if 
                                                   #more than one type of clumps is used

gbl._globals['constants']['gamma'] = 2.31  #exponent for mass-size relation
                                           #exponents: from Heithausen et al.

gbl._globals['constants']['M_sol'] = 1.98892*10**33 #solar mass in g

gbl._globals['constants']['M_H'] = 1.008 * 1.6605 * 10**(-24)  #hydrogen mass in g (1.008 u)

gbl._globals['constants']['nsigma'] = 3  #the n*sigma interval is accounted for during the opacity calculation
                                         #it can be increased for higher precision, decreased for computing speed

gbl._globals['threshold'] = {}

gbl._globals['threshold']['epsgrid'] = 10  #numerical threshold

gbl._globals['npar'] = 4  # number of parameter which have to be interpolated using regrid. 

gbl._globals['Etilde'] = {}  #directory for tables f1 and f2

gbl._globals['Etilde']['switch'] = True  #tables Ereal and Eimag still need to be read in of needed

# Define the flags
gbl._globals['check_negative'] = True  #(=0): check for negative optical depth when reading in KOSMA-\tau grid.
                                       #the current code functionality cannot treat negative optical depth
                                       #because the interpolations are done in log-space.

# Initialize the directory pathes
gbl._globals['directories'] = {}

gbl._globals['directories'] = pf.directoryManager(gbl._globals['directories'], \
                              'temp', gbl.KOSMAPATH+'temp/') 

# create temp directory if not existent and add to dictionary
gbl._globals['directories'] = pf.directoryManager(gbl._globals['directories'], \
                              'maps', gbl.KOSMAPATH+'maps/') 

# Run the main program
print('\n                                                                   _____      ______  ' + '\33[36m' + \
      '\n |   /   -------   -----  |\    /|      /\                        /     \     |     \ ' + '\33[35m' + \
      '\n |  /    |     |  /       | \  / |     /  \                              \    |      \' + '\33[34m' + \
      '\n | /     |     |  \____   |  \/  |    /    \         -----            ___/    |      |' + '\33[33m' + \
      '\n |/      |     |       \  |      |   /------\   ---    |                 \    |      |' + '\33[32m' + \
      '\n |  \    |     |       /  |      |  /        \         |                 /    |      /' + '\33[31m' + \
      '\n |   \   -------  -----   |      | /          \        |          \_____/     |_____/ ' + '\33[37m' + \
      '\n\n\n')
pf.runKOSMAt()