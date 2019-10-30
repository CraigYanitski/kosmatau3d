# set up project: 
#     - define geometry of the object
########## initializations ##########################################
import random as rand
import numpy as np
import time
from PDRfunctions import Ensemble
from PDR import _globals
from astropy.io import fits
from scipy.interpolate import griddata
_globals['compound'] = {}  #directory entries for compound


r = 1.5 * 1000#18    #disk radius in pc

scale = _globals['compound']['pixelsize'] = 333.3333333333333333333 #pixel size in pc
_globals['compound']['pixelsize'] = scale


x_length = 3#37  number of elements along x axis,r*2+1 # odd-number
y_length = 3#37 # edge length of map (projection)
z_length = 3#37 #
x_offset = (x_length-1)/2 #  to zentralize around zero
y_offset = (y_length-1)/2 #  (length-1)/2
z_offset = (z_length-1)/2

#x_axis = np.zeros(x_length)
x_axis = np.linspace( -(x_length-1)/2.*scale , (x_length-1)/2.*scale , num=x_length)
y_axis = np.linspace( -(y_length-1)/2.*scale , (y_length-1)/2.*scale , num=y_length)
z_axis = np.linspace( -(z_length-1)/2.*scale , (z_length-1)/2.*scale , num=z_length)
_globals['compound']['x_axis'] = x_axis
_globals['compound']['y_axis'] = y_axis
_globals['compound']['z_axis'] = z_axis

#print 'xAxis', x_axis
#print 'yAxis', y_axis
#print 'zAxis', z_axis

r_disk = r # radius disk
z_disk = 0.25 *1000 # thickness /  scaleheight of disk
inc_disk = 0.   /180.*np.pi # inclination of disk [rad] relative to y-axis
#phi_disk = 0 # ancle of disk in x,y-plane


#################################### Plotting options
show_mass_plot = 0
show_clump_mass_plot = 0 #1 shows plot
show_inter_mass_plot = 0
show_density_plot = 0
show_vz_plot = 0
show_disp_plot = 0
_globals['namesuffix'] = 'FUV=70_m1' #'FUV=70complement'

_globals['compound']['r_disk'] = r_disk
_globals['compound']['z_disk'] = z_disk
_globals['compound']['inc_disk'] = inc_disk
#_globals['compound']['phi_disk'] = phi_disk

import math as math
def calc_h(x,y,z,inc) : #calculates height h from koordinates relative to disk-plain
    h = (y - np.tan(inc)*z) * np.cos(inc)
    return h

def calc_r (x,y,z,inc) : #calculates the radius in the diskplain 
    h = (y - np.tan(inc)*z) * np.cos(inc)
    Px = x ; Py = y - np.cos(inc)*h ; Pz = z + np.tan(inc)*np.cos(inc)*h
    r = (Px**2 + Py**2 + Pz**2)**0.5
    return r
#print "calc_r", calc_r(0,0,2,45./180*np.pi)

#def calc_phi (x,y)
#    phi = np.arctan2 ( (x,(y**2+z**2)**0.5) )
#    return phi

coords = [] # liste where coordinates [x,y,z] will be appended
#coords_index = np.zeros([x_length, y_length, z_length])
#coords_index[:] = -1
#####################################################################
################# parameters to be changed ##########################
_globals['compound']['name'] = 'Galaxy' 
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
#thick = 100.   # PDR thickness 0.05 oder 0.06 pc. choose: 3pixels

# choose coorinates [x,y,z] at which ensembles shall be put. append to coords if inside disk

number = 0
for x in x_axis: #np.arange(x_length) - x_offset:
  #print "x", x
  for y in y_axis: #np.arange(y_length) - y_offset:
        for z in z_axis: #np.arange(z_length) - z_offset:
            if calc_r (x,y,z,inc_disk) < r_disk and abs(calc_h(x,y,z,inc_disk)) < z_disk: # and z < 8.5*1000:
                coords.append([x, y, z])
                #coords_index[x + x_offset][y + y_offset][z + z_offset] = number
                 
                #print 'cacl_r' , calc_r (x,y,z,inc_disk)
                #print 'h' , calc_h(x,y,z,inc_disk)
                #print x,y,z
                #print 'number', number

print "coords", coords
#print "number\n", number
_globals['compound']['coordinates'] = np.array(coords)
_globals['compound']['npoints'] = _globals['compound']['coordinates'].size/3
#_globals['compound']['coords_index'] = coords_index


print 'number of voxels', _globals['compound']['coordinates'].size/3
#print 'doublecheck...', number
#for i in np.arange(_globals['compound']['npoints']):
#    if np.alltrue(np.equal([0,0,0], coords[i])):
#        print 'i', i

#############################################################
# print 'number of ensembles', _globals['compound']['npoints']
# pause = raw_input('...')

_globals['binomial_poisson'] = 1 # 1
# are the clump statistics in the source best described by a 
# binomial or a poisson distribution? 1 means binomial, other
# input means poisson 

_globals['compound']['los_start'] =  np.array([0, 0, -z_length], int)
_globals['compound']['los_end']   =  np.array([0, 0, z_length], int)
# define start and endpoint of line of sight.
# line of sight has to be parallel to z axis (in positive direction)

#_globals['compound']['losoffsets'] = np.array([[-1, -1]]) #obsolete, fit cubes contains all possible LoS
# if not given the default is [0,0]. inputs needs to have the form np.array([[1,1],[2,2],..])

_globals['compound']['I_bg'] = 0 # Background Intensity 

# distance to source in pc. needed for conversion of beam size into pixels when creating maps
_globals['compound']['distance'] = 2.5*10**6 # 414    1kpc

# side length of pixels. value used for line of sight integration (dint = ds). in cm!!!!
_globals['compound']['dint'] = _globals['compound']['pixelsize'] *\
     _globals['constants']['pc']


# parameters for plots
# _globals['plots'] = {}
# set plotrange for geomety plot: 
# [[x_min, x_max], [y_min, y_max], [z_min, z_max]]
#_globals['plots']['plotrange'] = np.array([[-25, 25], [-25, 25], [-25, 25]])
_globals['plots']['plotrange'] = np.array([[-x_offset*scale, x_offset*scale], \
[-y_offset*scale,y_offset*scale ], [-z_offset*scale,z_offset*scale]])

_globals['sigma'] = {}
_globals['sigma']['sigma_cl_j'] = np.array([1.67/2.3548, 1.67/2.3548, 1.67/2.3548, 1.67/2.3548], float) 
# intrinsic clump linewidth (standard deviation)
# im Moment muessen die im code fuer alle massenpunkte identisch sein, da nur der erste Eintrag zur Umrechnung
# in peak werte genutzt wird. AENDERN!!!
_globals['sigma']['sigma_inter_j'] = np.array([1.67/2.3548, 1.67/2.3548, 1.67/2.3548, 1.67/2.3548], float)
# intrinsic clump linewidth for clumps representing interclump medium
_globals['sigma']['sigma_ens_cl_j'] = np.array([1.1/2.3548, 1.1/2.3548, 1.1/2.3548, 1.1/2.3548], float)
_globals['sigma']['sigma_ens_cl_j'] = np.array([10.1/2.3548, 10.1/2.3548, 10.1/2.3548, 10.1/2.3548], float) #broader change back
# ensemble "inter-clump" velocity spread (standard deviation)
_globals['sigma']['sigma_ens_inter_j'] = np.array([3.635/2.3548, 3.635/2.3548, 3.635/2.3548, 3.635/2.3548], float)
_globals['sigma']['sigma_ens_inter_j'] = np.array([10.635/2.3548, 10.635/2.3548, 10.635/2.3548, 10.635/2.3548], float)
# ensemble "inter-clump" velocity spread for interclump medium (sigma???)
################## (FWHM = 2.3548 sigma)


# information on sampling velocities: range: vel +- d_vel. 
# and number of steps: nstep # km s^-1
_globals['compound']['vel'] = float(0.) # central v
_globals['compound']['d_vel'] = float(330) *np.cos(inc_disk) +11.# max v
_globals['compound']['nstep'] = float(200) #number of v-baskets should be uneven
#
print 'min Velocity Sampling [km/s]:',-_globals['compound']['d_vel'] + _globals['compound']['vel']
print 'max Velocity Sampling [km/s]:',_globals['compound']['d_vel'] + _globals['compound']['vel'],'\n'
Dv = 2*_globals['compound']['d_vel'] / (_globals['compound']['nstep']-1)
print 'bin-distance km/s' , Dv

##################Abtastrate der Geschwindigkeitsbins 4 Koerbe pro Sigma
#bin size *2 sollte  kleiner totaler halbwertsbreite sein
sigma_cl_tot = ( np.array(_globals['sigma']['sigma_cl_j'])**2 + np.array(_globals['sigma']['sigma_ens_cl_j'])**2)**0.5
_globals['sigma']['sigma_ens_cl_j_total'] = sigma_cl_tot
sigma_inter_tot = (np.array(_globals['sigma']['sigma_inter_j'])**2 + np.array(_globals['sigma']['sigma_ens_inter_j'])**2)**0.5
_globals['sigma']['sigma_inter_j_total'] = sigma_inter_tot
print 'Sigma cl tot', sigma_cl_tot[0], 'Sigma_inter', sigma_inter_tot[0], 'km/s'

print 'adviced bin-number:',round(_globals['compound']['d_vel']*2 / min(sigma_cl_tot)/2.3548 *6.,0)


############################################################################
_globals['plots']['unit'] = 0
# set to 1 to get advanced map in erg/(cm^2 s sr). any other value
# to leave map in Kkm/s

#############################################################################
# Define ensemble (clump and interclump) parameters for each voxel
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
# FUV-strenght could be set directly
# velocity is the systematic velocity of the pixel/ensemble. If the 
#          pixel is not rotated only the z-component is used.
# v_dispersion is the velocity dispersion in z, created by different 
#          systematic velocity inside a voxel which is not resolved


for i in np.arange(_globals['compound']['npoints']):
    ens.append(Ensemble(i, \
               rho_ens_clumps  = 0. ,\
               rho_ens_inter   = 0. ,\
               Mens_clumps     = 0. ,\
               Mens_inter      = 0. ,\
               Ml_cl           = 10 ,\
               Mu_cl           = 10 ,\
               Ml_inter        = 1 ,\
               Mu_inter        = 1 ,\
               B               = 10 ,\
               FUV             = 0 ,\
               velocity        = np.array([0, 0, 0], float),\
               v_dispersion    = 0. \
               ))

##################################### read in rotation curve
convert = 1000. #if r_data are in kpc or in pc
rot_list_r =[]
rot_list_v =[]
with open ('compounds/disk/rot_curve.dat','r') as rot_file:
    next(rot_file) #skip header
    for line in rot_file:
        rot_list_r.append(line.split()[0])
        rot_list_v.append(line.split()[1])

rot_curve = np.zeros(  (len(rot_list_r),2)) #creates rot array

for i in np.arange(len(rot_list_r)):
    rot_curve[i,0] = float(rot_list_r[i])*convert #convert from kPC to PC
    rot_curve[i,1] = float(rot_list_v[i])
#print rot_curve
#pause = raw_input('rot_curve')


##################################### read in mass profile clumps
mass_list_r =[]
mass_list_v =[]
with open ('compounds/disk/mass_profile.dat','r') as mass_file:
  next(mass_file) #skip header
  for line in mass_file:
    mass_list_r.append(line.split()[0])
    mass_list_v.append(line.split()[1])

mass_profile = np.zeros(  (len(mass_list_r),2)) #creates mass array

for i in np.arange(len(mass_list_r)):
  mass_profile[i,0] = float(mass_list_r[i])*convert
  mass_profile[i,1] = float(mass_list_v[i])
#print mass_profile
#pause = raw_input('mass_profile')


##################################### read in density profile of clumps
dens_list_r =[]
dens_list_v =[]
with open ('compounds/disk/densities_clouds.dat','r') as dens_file:
  next(dens_file) #skip header
  for line in dens_file:
    dens_list_r.append(line.split()[0])
    dens_list_v.append(line.split()[1])

dens_profile = np.zeros(  (len(dens_list_r),2)) #creates dens array

for i in np.arange(len(dens_list_r)):
  dens_profile[i,0] = float(dens_list_r[i])*convert #kpc -> PC
  dens_profile[i,1] = float(dens_list_v[i])
#print dens_profile
#pause = raw_input('dens_profile')

##################################### read in mass profile interclumpmedium
mass_list_r =[]
mass_list_v =[]
with open ('compounds/disk/mass_profile_inter.dat','r') as mass_file:
  next(mass_file) #skip header
  for line in mass_file:
    mass_list_r.append(line.split()[0])
    mass_list_v.append(line.split()[1])

mass_profile_inter = np.zeros(  (len(mass_list_r),2)) #creates mass array

for i in np.arange(len(mass_list_r)):
  mass_profile_inter[i,0] = float(mass_list_r[i])*convert #kpc->PC
  mass_profile_inter[i,1] = float(mass_list_v[i])
#print mass_profile_inter
#pause = raw_input('mass_profile_inter')

#####################################################################
##################################### define speed and masses for npoints
total_clump_mass = 0
total_inter_mass = 0
for i in np.arange(_globals['compound']['npoints']):  
    ################### define all  coordinates
    x = coords[i][0]
    y = coords[i][1]
    z = coords[i][2]
    #print x,y,z
    #raw_input('sad')
    r = calc_r (x,y,z,inc_disk)
    h = calc_h (x,y,z,inc_disk)
    #phi = np.arctan2 (x,y)
    phi = np.arctan2 (x,(y**2+z**2)**0.5)
    #print 'x', x, 'y', y, 'z', z
    #print 'umlaufwinkel', phi/3.14*180
    ###################################### Setting Vz
    ####interpolate from rot_curve for current r
    v = griddata(rot_curve[:,0], rot_curve[:,1], r, method='linear')
    if x == 0 and y == 0 and z == 0:
      v = 233 # center speed
    v_z = v * np.cos(inc_disk) * np.sin(phi) # calc v_z-component
    ens[i].velocity[2] = v_z
    #print "v_z [km/s]", ens[i].velocity[2]
    

    ####################################### voxel mass clump (H2)
    # Gesamt-Dichteprofil in R-Richtung fuer massen der voxel
    mass_surface = griddata(mass_profile[:,0], mass_profile[:,1], r, method='linear') #s-density of voxel (not clumps)
    mass_density = mass_surface / _globals['compound']['pixelsize'] # -> in SM / PC^3
    
    ########interpolate surrounding mass of current r for more precise density and calculate v_dispersion
 
    bins = 5
    end_dens = 0.
    end_disp = 0.
    for loc_x in np.linspace(-_globals['compound']['pixelsize']/2,+_globals['compound']['pixelsize']/2, bins):
        for loc_y in np.linspace(-_globals['compound']['pixelsize']/2,+_globals['compound']['pixelsize']/2, bins):
            #print 'x',x, 'y',y, 'r', r, 'r|'   #environment r_s around r central r of voxel 
            alpha = np.arctan2 (x,y) #absolute ancle 
            loc_alpha = np.arctan2 (loc_x,loc_y) #local ancle
            absolute_x =  np.cos(alpha)*r * _globals['compound']['pixelsize'] + np.cos(loc_alpha) * (loc_x**2+loc_y**2)**0.5
            absolute_y =  np.sin(alpha)*r * _globals['compound']['pixelsize'] + np.sin(loc_alpha) * (loc_x**2+loc_y**2)**0.5
            #in pc not kpc!
            #print 'absolute x', absolute_x, 'loc_x', loc_x
            #print 'absolute y', absolute_y, 'loc_y', loc_y
            absolute_r = (absolute_x**2 + absolute_y**2)**0.5 / _globals['compound']['pixelsize']
            #print absolute_r, 'rrrrrrrrrrrrrrrr'
            loc_dens = griddata(mass_profile[:,0], mass_profile[:,1], absolute_r, method='linear') #actual surface density
            loc_v = griddata(rot_curve[:,0], rot_curve[:,1], absolute_r, method='linear')
            #print 'loc_dens', loc_dens
            end_dens = end_dens + loc_dens
            end_disp = end_disp + ( (loc_v - v)**2)**0.5
    end_dens = end_dens / bins**2 / _globals['compound']['pixelsize']
    end_disp = abs( end_disp / bins**2  * np.cos(inc_disk) * np.sin(phi) )# normalize to binnumber and make calc v_z-component
    #print 'x,y,z', x,y,z
    #print 'end_disp', end_disp
    #print 'central v', v_z 
    #a = raw_input ('vvvv')
    ens[i].v_dispersion = end_disp # dispersion by systematic speed which can't be resolved
    #print mass_density, 'mass_density'

    #ens[i].Mens_clumps = mass_density * _globals['compound']['pixelsize']**3 # -> total voxel mass
    ens[i].Mens_clumps = end_dens * _globals['compound']['pixelsize']**3 # -> total voxel mass with additional interpol
    #print "Voxel_mass:", ens[i].Mens_clumps ,"[Sm]" ,"at r", r      
    #SM/pc->N/cm^3
    #density_cm = ens[i].Mens_clumps * _globals['constants']['M_sol'] / _globals['constants']['M_H'] /_globals['constants']['pc']**3
    
   
   ######################################### Setting clump H2 surface density n(r)
    clump_rho = griddata(dens_profile[:,0], dens_profile[:,1], r, method='linear')
    if x == 0 and y == 0 and z == 0:
      clump_rho = 15000. # set centermass
    ens[i].rho_ens_clumps = clump_rho # N/cm^3

    if np.isnan(clump_rho) == True:
      print 'voxel number', i
      print 'r', r, 'h' , h
      pause = raw_input('NaN in clumpdensity found')
    #print "Clump density in N/cm^3", ens[i].rho_ens_clumps , "at r in kpc", r , "h in kpc", h 


    ####################################### voxel mass interclump (H1)
    # Gesamt-Dichteprofil in R-Richtung fuer massen der voxel
    mass_surface_inter = griddata(mass_profile_inter[:,0], mass_profile_inter[:,1], r, method='linear')
    mass_density_inter = mass_surface_inter / _globals['compound']['pixelsize'] # -> in SM / PC^3
    ens[i].Mens_inter = mass_density_inter * _globals['compound']['pixelsize']**3
    #print  ens[i].Mens_inter, ''

    ####################################### density inter (H1)
    ens[i].rho_ens_inter = ens[i].rho_ens_clumps  ####now same, but with different input grid should asume inter medium is less dense

    total_clump_mass = total_clump_mass + ens[i].Mens_clumps
    total_inter_mass = total_inter_mass + ens[i].Mens_inter
print "Total clump-mass:", total_clump_mass/10**9, 'in G SunM'
print "Total inter-mass:", total_inter_mass/10**9, 'in G SunM'

_globals['compound']['ens'] = ens #save ensamble to globals
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

###################################################### Mass-Plot
import PDRfunctions as pf 
import pylab as p
if show_mass_plot == 1:
  Mplot = []
  for i in np.arange(_globals['compound']['npoints']):  
    #Mplot.append(ens[i].Mens_clumps + ens[i].Mens_inter)
    Mplot.append(np.log10(ens[i].Mens_clumps + ens[i].Mens_inter))
    #print "total Voxel mass", ens[i].Mens_clumps , ens[i].Mens_inter
  #print Mplot
  fig2 = pf.plot_geo(_globals['compound']['coordinates'], \
           wei = Mplot, limits = _globals['plots']['plotrange'], \
           title = 'total mass distribution', cbarLabel = "log S_mass per voxel", \
           labelx ='x [pc]', labely ='y [pc]', labelz='z [pc]')
  p.show()

###################################################### clump_Mass-Plot
if show_clump_mass_plot == 1:
  Mplot = []
  for i in np.arange(_globals['compound']['npoints']):  
    Mplot.append(np.log10(ens[i].Mens_clumps))
  fig2 = pf.plot_geo(_globals['compound']['coordinates'], \
           wei = Mplot, limits = _globals['plots']['plotrange'], \
           title = 'inter mass distribution', cbarLabel = "log Smass of clumps per voxel", \
           labelx ='x [pc]', labely ='y [pc]', labelz='z [pc]')
  p.show()

###################################################### inter_Mass-Plot
if show_inter_mass_plot == 1:
  Mplot = []
  for i in np.arange(_globals['compound']['npoints']):  
    Mplot.append(np.log10(ens[i].Mens_inter))
  fig2 = pf.plot_geo(_globals['compound']['coordinates'], \
           wei = Mplot, limits = _globals['plots']['plotrange'], \
           title = 'inter mass distribution', cbarLabel = "log S_mass intermedium per voxel", \
           labelx ='x [kpc]', labely ='y [kpc]', labelz='z [kpc]')
  p.show()


##################################################### Density-Plot
#show_density_plot = 0 #1 shows plot
if show_density_plot == 1:
  Mplot = []
  for i in np.arange(_globals['compound']['npoints']):  
    Mplot.append(np.log10(ens[i].rho_ens_clumps + ens[i].rho_ens_inter))
  #print Mplot
  fig4 = pf.plot_geo(_globals['compound']['coordinates'], \
           wei = Mplot, limits = _globals['plots']['plotrange'], \
           title = 'density distribution', cbarLabel = "log density [/cm^3]", \
           labelx ='x [pc]', labely ='y [pc]', labelz='z [pc]')
  p.show()

#rho_ens_clumps

##################################################### Vz-Plot
#show_vz_plot = 0   # 1 shows plot
if show_vz_plot == 1:
  Mplot = []
  for i in np.arange(_globals['compound']['npoints']):
      Mplot.append(ens[i].velocity[2])
  fig3 = pf.plot_geo(_globals['compound']['coordinates'], \
           wei = Mplot, limits = _globals['plots']['plotrange'], \
           title = 'speed distribution', cbarLabel = "speed of z-component [km/s]", \
           labelx ='x [pc]', labely ='y [pc]', labelz='z [pc]')
  p.show()


##################################################### v_dispersion-Plot
if show_disp_plot == 1:
  Mplot = []
  for i in np.arange(_globals['compound']['npoints']):
      Mplot.append(ens[i].v_dispersion)
  fig3 = pf.plot_geo(_globals['compound']['coordinates'], \
           wei = Mplot, limits = _globals['plots']['plotrange'], \
           title = 'systematic v dispersion', cbarLabel = "dispersion in z-component [km/s]", \
           labelx ='x [pc]', labely ='y [pc]', labelz='z [pc]')
  p.show()


################################################### Setting FUV-field
calc_fuv = 1
if calc_fuv == 1:
  print "calulation disk-FUV"
  fuv_file = open ('temp/FUV.dat','w')
  #fuv_file = open ('temp/FUV_chris.dat','w')
  
  for i in np.arange(_globals['compound']['npoints']):
    # print i
    x = coords[i][0]
    y = coords[i][1]
    z = coords[i][2]
    r = calc_r (x,y,z,inc_disk)
    h = calc_h (x,y,z,inc_disk)
    #phi = np.arctan2 (x,(y**2+z**2)**0.5)
    #print r
    r_OB = 1/6**0.5 *1000 # distance on average of OBs in pc
    d_MC = rand.random()*30+20  #random for distances of the clouds to next OB assosiation in pc (20-50) 
    UV = (r_OB**2 / d_MC**2)/3
    UV = 70    
    #print "lokal UV-field",UV

    #write out like in Silkes method: (should become mostly obsolete)        
    fuv_file.write(str(UV) + str(" ") + str(  _globals['compound']['coordinates'][i] ) +"\n"  )

    #fill directly into ensamble, and create fitsfile later with:
    ens[i].FUV = UV

  fuv_file.close()



################################################### Writing fits files

write_fits = 1
if write_fits == 1:
  print "creating fits-file"
  
  mass_data = np.zeros( (z_length,y_length,x_length) , dtype=np.float32)  #initialize fit grids with zeros
  mass_clump_data = np.zeros( (z_length,y_length,x_length) , dtype=np.float32)
  mass_inter_data = np.zeros( (z_length,y_length,x_length) , dtype=np.float32)
  vz_data = np.zeros( (z_length,y_length,x_length), dtype=np.float32 )
  density_data = np.zeros( (z_length,y_length,x_length), dtype=np.float32 )
  FUV_data = np.zeros( (z_length,y_length,x_length), dtype=np.float32 )
  disp_data = np.zeros( (z_length,y_length,x_length), dtype=np.float32 )

  for i in np.arange(_globals['compound']['npoints']): # fill grids with values
    x = coords[i][0]#+x_offset #transfrom real position to array coordinates
    y = coords[i][1]#+y_offset
    z = coords[i][2]#+z_offset
    #print 'x,y,z in PC' , x,y,z
    x_i = coords[i][0]/scale + (_globals['compound']['mapSize']['x']-1)/2. #calculates correct array positions
    y_i = coords[i][1]/scale + (_globals['compound']['mapSize']['y']-1)/2.
    z_i = coords[i][2]/scale + (_globals['compound']['mapSize']['z']-1)/2.
    #print 'x_i, y_i, z_i', x_i, y_i, z_i
    #raw_input('Koords')
    mass_data[z_i][y_i][x_i] = ens[i].Mens_clumps + ens[i].Mens_inter # write mass data to array coordiantes
    mass_clump_data[z_i][y_i][x_i] = ens[i].Mens_clumps
    mass_inter_data[z_i][y_i][x_i] = ens[i].Mens_inter
    vz_data[z_i][y_i][x_i] = ens[i].velocity[2] # velocitydata in z
    density_data[z_i][y_i][x_i] = ens[i].rho_ens_clumps # densitydata
    FUV_data[z_i][y_i][x_i] = ens[i].FUV #FUV data
    disp_data[z_i][y_i][x_i] = ens[i].v_dispersion

  #print 'FUV_data', FUV_data
  #a = raw_input ('pause')

  hdu = fits.PrimaryHDU(mass_data) # sets data as hdu-list
  #################some Header content
  hdu.header['OBSERVER'] = "Just a Simulation"
  hdu.header['NAME'] = "galactic Disk"
  hdu.header['COMMENT'] = " easter egg Fr. 13.5 "
  #hdu.header['BZERO'] = 0
  #hdu.header['BSCALE'] = 1

  hdu.header['CRPIX1'] = 1   #ref Pixel
  hdu.header['CTYPE1'] = 'x-dim [pc]' #Name 
  hdu.header['CRVAL1'] = ((-x_length-1)/2.)*scale #Reference Value in kpc
  hdu.header['CDELT1'] = 1. #stepsize per counting index in kpc

  hdu.header['CRPIX2'] = 1
  hdu.header['CTYPE2'] = "y-dim [pc]"
  hdu.header['CRVAL2'] = ((-y_length-1)/2.)*scale 
  hdu.header['CDELT2'] = 1. 
 
  hdu.header['CRPIX3'] = 1
  hdu.header['CTYPE3'] = "z-dim [pc]"
  hdu.header['CRVAL3'] = ((-z_length-1)/2.)*scale
  hdu.header['CDELT3'] = 1. 

  #### clearing old files
  open('fits_file/data/mass.fits', 'w').close()
  open('fits_file/data/mass_clump.fits', 'w').close()
  open('fits_file/data/mass_inter.fits', 'w').close()
  open('fits_file/data/vz.fits', 'w').close()
  open('fits_file/data/density.fits', 'w').close()
  open('fits_file/data/FUV.fits', 'w').close()
  open('fits_file/data/disp.fits', 'w').close()

##### first writing
  hdu.header['BTYPE'] = 'all mass'
  hdu.header['BUNIT'] = " SM/voxel"
  hdu.writeto('fits_file/data/mass.fits', clobber = True)

  hdu.data = mass_clump_data
  hdu.header['BTYPE'] = 'clump mass'
  hdu.writeto('fits_file/data/mass_clump.fits', clobber = True)

  hdu.data = mass_inter_data
  hdu.header['BTYPE'] = 'inter mass'
  hdu.writeto('fits_file/data/mass_inter.fits', clobber = True)

###velocity
  hdu.header['BTYPE'] = 'velocity'
  hdu.header['BUNIT'] = " km/s"
  hdu.data = vz_data
  hdu.writeto('fits_file/data/vz.fits', clobber = True)
### density
  hdu.header['BTYPE'] = 'density of clouds'
  hdu.header['BUNIT'] = " particle / cm^3"
  hdu.data = density_data
  hdu.writeto('fits_file/data/density.fits', clobber = True)
### FUV
  hdu.header['BTYPE'] = 'FUV-field'
  hdu.header['BUNIT'] = " drain-field"
  hdu.data = FUV_data
  hdu.writeto('fits_file/data/FUV.fits', clobber = True)
### v_dispersion
  hdu.header['BTYPE'] = 'unresolved systematic vz_dispersion'
  hdu.header['BUNIT'] = " km/s"
  hdu.data = disp_data
  hdu.writeto('fits_file/data/disp.fits', clobber = True)

_globals['runtimes']['end_model'] = time.time()

################ write out voxel sanity check
write_voxel = 1
if write_voxel == 1:
  print "writing Voxel"
  voxel_file = open ('temp/voxel.dat','w')
    
  for i in np.arange(_globals['compound']['npoints']):
    # print i
    x = coords[i][0]
    y = coords[i][1]
    z = coords[i][2]
    r = calc_r (x,y,z,inc_disk)
    h = calc_h (x,y,z,inc_disk)
    phi = np.arctan2 (x,(y**2+z**2)**0.5)
    #print r
    k=1 #kommastellen
    voxel_file.write( 'x,y,z' + str(  _globals['compound']['coordinates'][i] ) +\
    ' r '+str(round(r,k)) + ' phi ' + str(round(phi/np.pi*180,k)) + ' vz' + str(round(ens[i].velocity[2],k)) + \
    ' vz_disp ' + str(round(ens[i].v_dispersion,k)) + \
    ' rho_clumps [/cm^3] ' + str(round(ens[i].rho_ens_clumps,k)) + ' mens_clump [SM]' + str(round(ens[i].Mens_clumps,k))  + \
 
    "\n"  )

  voxel_file.close()

print "model.py durchgelaufen!"
pause = raw_input('OK')
