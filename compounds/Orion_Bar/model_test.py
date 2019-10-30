# OriBar_small: geometry of the Orion bar PDR

# set up project: 
#     - define geometry of the object

########## initializations ##########################################
import numpy as np
from PDRfunctions import Ensemble
from PDR import _globals

#x_length = 55 # edge length of map (projection)
x_length = 3 
y_length = 3 # edge length of map (projection)
z_length = 3 #66 + 14 - 20 - 10
x_offset = 0
y_offset = 5
z_offset = 0
_globals['compound'] = {} 
# directory entries for compound

coords = [] # liste where coordinates [x,y,z] will be appended
coords_index = np.zeros([x_length, y_length, z_length])
coords_index[:] = -1
#####################################################################
################# parameters to be changed ##########################
_globals['compound']['name'] = 'OriBar_hogerheijde; 1 pixel = 0.01 pc' 
# name of compound
_globals['compound']['mapSize'] = {}
_globals['compound']['mapSize']['x'] = x_length
_globals['compound']['mapSize']['y'] = y_length
_globals['compound']['mapSize']['z'] = z_length
_globals['compound']['mapPixels'] = x_length * y_length
_globals['compound']['offsets'] = {}
_globals['compound']['offsets']['x'] = x_offset
_globals['compound']['offsets']['y'] = y_offset
_globals['compound']['offsets']['z'] = z_offset
# map Pixels = number of lines of sight
########### geometry parameters #####################################
thick = 100.   # PDR thickness 0.05 oder 0.06 pc. choose: 3pixels
theta = 3. # 15 # angle between los and edge-on part of the pdr. in degree
# IF at y=0

# choose coorinates [x,y,z] at which ensembles shall be put. append to coords
# for x in np.arange(10) :
number = 0
for x in np.arange(x_length) - x_offset:
    for y in np.arange(y_length) - y_offset:
        for z in np.arange(z_length) - z_offset:
            if z < 0: 
                if y >= np.tan(theta * np.pi/180.) * abs(z) - thick: 
                    coords.append([x, y, z])
                    coords_index[x + x_offset][y + y_offset][z + z_offset] = number
                    number = number + 1
            elif z >= 30 - z_offset: 
                if y <= -np.tan(theta * np.pi/180.) * z: 
                    coords.append([x, y, z])
                    coords_index[x + x_offset][y + y_offset][z + z_offset] = number
                    number = number + 1
            else: 
                if y >= -np.tan(theta * np.pi/180.) * z - thick and \
                   y <= -np.tan(theta * np.pi/180.) * z:
                       coords.append([x, y, z])
                       coords_index[x + x_offset][y + y_offset][z + z_offset] = number
                       number = number + 1
"""
            else: 
                if y <= -np.tan(theta * np.pi/180.) * z:
                       coords.append([x, y, z])
                       coords_index[x + x_offset][y + 32][z + 6] = number
                       number = number + 1

"""
_globals['compound']['coordinates'] = np.array(coords)
_globals['compound']['npoints'] = _globals['compound']['coordinates'].size/3
_globals['compound']['coords_index'] = coords_index

print 'number of pixels', _globals['compound']['coordinates'].size/3
print 'doublecheck...', number
for i in np.arange(_globals['compound']['npoints']):
    if np.alltrue(np.equal([0,0,0], coords[i])):
        print 'i', i
# test = coords_index[0 + x_offset][0 + y_offset][0 + z_offset]
"""
test0_m2_0 = coords_index[0 + x_offset][-2 + y_offset][50 + z_offset] 
test0_m5_0 = coords_index[0 + x_offset][-5 + y_offset][50 + z_offset] 
test0_m7_0 = coords_index[0 + x_offset][-7 + y_offset][50 + z_offset] 
test0_m10_0 = coords_index[0 + x_offset][-10 + y_offset][50 + z_offset] 
test0_m15_18 = coords_index[0 + x_offset][-15 + y_offset][50 + z_offset] 
#test2_m15_18 = coords_index[2 + x_offset][-15 + y_offset][18 + z_offset] 
test0_m16_18 = coords_index[0 + x_offset][-16 + y_offset][50 + z_offset] 
#test2_m16_18 = coords_index[2 + x_offset][-16 + y_offset][18 + z_offset] 
test0_m20_18 = coords_index[0 + x_offset][-20 + y_offset][50 + z_offset] 
#test2_m20_18 = coords_index[2 + x_offset][-20 + y_offset][18 + z_offset] 
test0_m32_18 = coords_index[0 + x_offset][-32 + y_offset][50 + z_offset] 
#test2_m32_18 = coords_index[2 + x_offset][-32 + y_offset][18 + z_offset]
 
print 'test0_m2_0', test0_m2_0
print 'test0_m5_0', test0_m5_0
print 'test0_m7_0', test0_m7_0
print 'test0_m10_0', test0_m10_0
print 'test0_m15_18', test0_m15_18
#print 'test2_m15_18', test2_m15_18 
print 'test0_m16_18', test0_m16_18
#print 'test2_m16_18', test2_m16_18 
print 'test0_m20_18', test0_m20_18
#print 'test2_m20_18', test2_m20_18
print 'test0_m32_18', test0_m32_18
#print 'test2_m32_18', test2_m32_18 
print 'COORDS[TEST]' , coords[int(test)]

pause = raw_input('...')
"""
#############################################################

_globals['compound']['npoints'] = _globals['compound']['coordinates'].size/3
_globals['compound']['pixelsize'] = 0.01 # pixel size in pc

# print 'number of ensembles', _globals['compound']['npoints']
# pause = raw_input('...')

_globals['binomial_poisson'] = 1
# are the clump statistics in the source best described by a 
# binomial or a poisson distribution? 1 means binomial, other
# input means poisson 

_globals['compound']['los_start'] = np.array([0, 0, -22], int)
_globals['compound']['los_end'] =  np.array([0, 0, 62], int)
# define start and endpoint of line of sight.
# currently rotation of compound is not working so line
# of sight has to be parallel to z axis (in positive direction)
# TEST: needs to go through 0/0 ?

_globals['compound']['losoffsets'] = np.array([[0, -5]])
#_globals['compound']['losoffsets'] = np.array([[0, 1],[0, 0],[0, -1],[0, -2],[0, -3],[0, -4], [0, -5],[0, -6],[0, -7],[0, -8],[0, -9], [0, -10], [0, -15]])
# give x-y-offsets of position(s) at which you want to see spectra.
# offsets are relative to line of sight AFTER translation and rotation
# if not given the default is [0,0].
# inputs needs to have the form np.array([[1,1],[2,2],..])

_globals['compound']['I_bg'] = 0 # Background Intensity 
# distance to source in pc. needed for conversion of beam size into pixels when creating maps
_globals['compound']['distance'] = 414
# side length of pixels. value used for line of sight integration (dint = ds). in cm!!!!
_globals['compound']['dint'] = _globals['compound']['pixelsize'] *\
     _globals['constants']['pc']

# parameters for plots
# _globals['plots'] = {}
# set plotrange for geomety plot: 
# [[x_min, x_max], [y_min, y_max], [z_min, z_max]]
_globals['plots']['plotrange'] = np.array([[-30, 30], [-35, 25], [-20, 40]])

# information on sampling velocities: range: vel +- d_vel. 
# and number of steps: nstep 
_globals['compound']['vel'] = float(11.3)
_globals['compound']['d_vel'] = float(5)
_globals['compound']['nstep'] = float(21) 
# 
_globals['sigma'] = {}
_globals['sigma']['sigma_cl_j'] = np.array([1/2.3548, 1/2.3548, 1/2.3548, 1/2.3548], float) 
# intrinsic clump linewidth (standard deviation)
_globals['sigma']['sigma_inter_j'] = np.array([1/2.3548], float)
# intrinsic clump linewidth for clumps representing interclump medium
_globals['sigma']['sigma_ens_cl_j'] = np.array([2/2.3548, 2/2.3548, 2/2.3548, 2/2.3548], float)
# ensemble "inter-clump" velocity spread (standard deviation)
_globals['sigma']['sigma_ens_inter_j'] = np.array([4/2.3548], float)
# ensemble "inter-clump" velocity spread for interclump medium (sigma???)
################## (FWHM = 2.3548 sigma)

############################################################################
_globals['plots']['unit'] = 0
# set to 1 to get advanced map in erg/(cm^2 s sr). any other value
# to leave map in Kkm/s

#############################################################################
# Define ensemble (clump and interclump) parameters for each pixel
ens = [] 
# initialize array for clumpy ensembles
# CHOOSE THE FOLLOWING PARAMETERS
# rho_ens_clumps, rho_ens_inter: ensemble-averaged density in cm^-3
# for clump and interclump medium  
# (surface density = rho_ens/1.91)
# ......
# Mu, Ml  are the discrete highest and lowest mass. not the
# ensemble interval limits. can be 10^-3. 10^-2,.... 10^3 with Mu>Ml
# Mu and Ml can be set to the same value, if only one type of clumps
# shall be used 
# B: for logarithmic scaling B=M_(i+1)/M_i. At the moment only B=10 has
# been tried out
# velocity is the systematic velocity of the pixel/ensemble. If the 
#          pixel is not rotated only the z-component is used.

for i in np.arange(_globals['compound']['npoints']):
    ens.append(Ensemble(i, \
               rho_ens_clumps  = 4. * 10**6, \
               rho_ens_inter   = 4. * 10**4, \
               Mens_clumps     = 0., \
               Mens_inter      = 0.00173 * 0.45, \
               Ml_cl           = 0.001, \
               Mu_cl           = 1., \
               Ml_inter        = 0.01, \
               Mu_inter        = 0.01, \
               B               = 10 ,\
               velocity        = np.array([0, 0, 11.3], float)\
               ))
################## position dependent clump/interclump ratio ##############

f = None
#fi = None # interclump
Mens_clumps = []
for i in np.arange(_globals['compound']['npoints']):  
    # print i
    y = coords[i][1]
    z = coords[i][2]
    if z >= 0 and y <= 0: 
        if y <= -np.tan(theta * np.pi/180.) * z \
            and y > -np.tan(theta * np.pi/180.) * z - 3:  
                f = 0.0
                # fi = 5./10.
            # IF
        #elif y <= -np.tan(theta * np.pi/180.) * z - 3: f = 3
        elif y <= -np.tan(theta * np.pi/180.) * z - 3\
            and y > -np.tan(theta * np.pi/180.) * z - 4:  
                f = 5./10.
                # fi = 5./10.
        elif y <= -np.tan(theta * np.pi/180.) * z - 4\
            and y > -np.tan(theta * np.pi/180.) * z - 5:  
                f = 5./10.
                # fi = 5./10.
        elif y <= -np.tan(theta * np.pi/180.) * z - 5: 
            f = 5./10. 
            # fi = 10./10.
        # else: f = 1./6.       
        # pause = raw_input()
    elif z < 0 and y > 0: 
        if z >= -3: f = 0.0
        elif z < -3 and z >= -4: 
            f = 5./10.
            # fi = 5./10.
        elif z < -4 and z >= -5: 
            f = 5./10.
            # fi = 5./10.
        # elif z < -3 and z >= -4: f = 2./10
        elif z < -5: 
            f = 5./10. 
            # fi = 10./10.
    elif z <= 0 and y <= 0: 
        if (z >= -3 and y > -3):  
            f = 0.0
            # fi = 5./10.
        #else: f = 1   
        #if z == 0 and y == 0: f = 0
        elif (z == -4 and y >= -3) or (y == -3 and z >= -4):  
            f = 5./10.
            # fi = 5./10.
        elif (z == -5 and y >= -4) or (y == -4 and z >= -5):  
            f = 5./10. 
            # fi = 5./10.
        else: 
            f = 5./10. 
            # fi = 10./10.
    # print 'f', f
    ens[i].Mens_clumps = 0.00173 * 4 * f
    # ens[i].Mens_inter = 0.00173 * fi 
    # print 'ens[i].Mens_clumps', ens[i].Mensigma)
    # print 'ens[i].Mens_inter', ens[i].Mens_inter
    Mens_clumps.append(ens[i].Mens_clumps)

"""
import PDRfunctions as pf 
import pylab as p
fig2 = pf.plot_geo(_globals['compound']['coordinates'], \
           wei = Mens_clumps, limits = _globals['plots']['plotrange'], \
           title = 'mass distribution, clumps', cbarLabel = "mass")
p.show()
"""
#################### interclump #########################################
"""
fi = None # interclump
Mens_inter = []
for i in np.arange(_globals['compound']['npoints']):  
    # print i
    y = coords[i][1]
    z = coords[i][2]
    if z >= 0 and y <= 0: 
        if y <= -np.tan(theta * np.pi/180.) * z \
            and y > -np.tan(theta * np.pi/180.) * z - 5:  
                fi = 1.
        elif y <= -np.tan(theta * np.pi/180.) * z - 5: 
            fi = 2. 
    elif z < 0 and y > 0: 
        if z >= -5: fi = 1.
        elif z < -5: 
            fi = 2. 
    elif z <= 0 and y <= 0: 
        if (z >= -5 and y > -5):  
            fi = 1.
        else: 
            fi = 2. 
    ens[i].Mens_inter = 0.00173 * 0.5 * fi 
    Mens_inter.append(0.00173 * 0.5 * fi)

fig2 = pf.plot_geo(_globals['compound']['coordinates'], \
           wei = Mens_inter, limits = _globals['plots']['plotrange'], \
           title = 'mass distribution, interclump', cbarLabel = "mass")
p.show()
"""
#####################################################################

_globals['compound']['ens'] = ens
# save ensemble set-up
# pause = raw_input()

######################### END OF DEFINITIONS ###############################
#############################################################################
# set axis labels for 3D geometry plots. Do not need to be changed.
_globals['plots']['xlabel'] = 'X offset [' + str(_globals['compound']['pixelsize']) + 'pc]'
_globals['plots']['ylabel'] = 'Y offset [' + str(_globals['compound']['pixelsize']) + 'pc]'
_globals['plots']['zlabel'] = 'to observer [' + str(_globals['compound']['pixelsize']) +  'pc]'
#############################################################################
_globals['compound']['dgrid'] = 1 
# factor defining streching of grid in dgrid
# currently not used
# DO NOT CHANGE
