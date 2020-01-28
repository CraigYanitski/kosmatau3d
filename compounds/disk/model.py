# set up project: 
#     - define geometry of the object
# Initializations
import random as rand
import numpy as np
import time
import PDRfunctions as pf
#from PDR import gbl._globals
from astropy.io import fits
from scipy.interpolate import griddata
# edit: Craig, 14.10.2019
import globals as gbl
#end edit
gbl._globals['compound'] = {}  #directory entries for compound
r = 18*1000#18    #disk radius in pc
#r = 18.1*1000
factor = 1 #larger makes finer voxelgrid
scale = gbl._globals['compound']['pixelsize'] = 1000.0 / factor #voxel size in pc
#gbl._globals['compound']['pixelsize'] = scale
#print 'cells in one dimension', (36*factor+1)
x_length = 36*factor+1#37  number of elements/voxels along x axis,r*2+1 # odd-number
y_length = 36*factor+1#37 # edge length of map (projection) can be lower for think galactic disk in 90 degree ancle
y_length = 0*factor+1#37 # edge length of map (projection - increases last parts speed radiative transfer
z_length = 36*factor+1#37 #
x_offset = (x_length-1)/2 #  to zentralize around zero
y_offset = (y_length-1)/2 #  (length-1)/2
z_offset = (z_length-1)/2
r_disk = r # radius disk
z_disk = 1*1000 # thickness /  scaleheight of disk
inc_disk = 0.   /180.*np.pi # inclination of disk [rad] relative to y-axis
#phi_disk = 0 # ancle of disk in x,y-plane
# Naming and setting for testcases
global_UV = 10 #in drain-field - max UV
clump_mass_l = 10.#0.1 # used clump masses in SM, lower
clump_mass_u = 100# upper boundary
global_densityfactor = 2. #multiplies surface density  1 bis 3
global_massfactor = 1. #multiplies mass   1 bis 4
# with intermedium parameters
clump_there = 'yes'   # 'no' puts not mass to exclude medium for testing 
inter_there = 'yes'
massratio_clump_to_inter = 1.0
inter_mass_l = 0.01
inter_mass_u = 0.01
UV_interpol = 'yes' #yes = natale sonst = const
cut_switch = 'total' #'(inner vs outer    or total)
#cut_switch = 'outer'
#cut_switch = 'total'
gbl._globals['compound']['interclump_FUV'] = 1 # one drainfield
gbl._globals['namesuffix'] = '_' + str( int(gbl._globals['compound']['pixelsize']) )+ '_' \
                              + 'FUV' + str(global_UV)+ '_mass' + str(clump_mass_l) +'_'+ str(clump_mass_u) + '_densfac' + str(global_densityfactor) \
                              + '_mass_factor'+ str(global_massfactor) +'_cut_'+ str(cut_switch) + 'dust-UVnatale-yesinter-oneclump'      ######UV2max1_clump0.01' #'maxFUV10withinter' #'inner/outer complement'
print(gbl._globals['namesuffix'])
#a = input('pause')
filename = gbl.COMPOUNDPATH+'rot_milki2018_14.dat'
#filename = 'compounds/disk/alternative_rot/rot_curve_max.dat' #alt
#filename = 'compounds/disk/rot_curve.dat'
# Plotting options
# - located in globals.py -
gbl._globals['compound']['r_disk'] = r_disk
gbl._globals['compound']['z_disk'] = z_disk
gbl._globals['compound']['inc_disk'] = inc_disk
#gbl._globals['compound']['phi_disk'] = phi_disk
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
#def calc_phi (x,y,z)
#    phi = np.arctan2 ( (x,(y**2+z**2)**0.5) )
#    return phi
coords = [] # liste where coordinates [x,y,z] will be appended
abs_coords = [] #list where absolute coordinates will be appended
real_coords = [] # list where real_coordinates [x,y,z] will be appended
# Parameters to be changed 
gbl._globals['compound']['name'] = 'Galaxy'  # name of compound
gbl._globals['compound']['mapSize'] = {}
gbl._globals['compound']['mapSize']['x'] = x_length
gbl._globals['compound']['mapSize']['y'] = y_length
gbl._globals['compound']['mapSize']['z'] = z_length
gbl._globals['compound']['mapPixels'] = x_length * y_length
gbl._globals['compound']['offsets'] = {}
gbl._globals['compound']['offsets']['x'] = x_offset
gbl._globals['compound']['offsets']['y'] = y_offset
gbl._globals['compound']['offsets']['z'] = z_offset
# map Pixels = number of lines of sight
# geometry parameters 
# choose coorinates [x,y,z] at which ensembles shall be put. append to coords if inside disk
# this way you don't need to plan in advance how many points will be inside list
R0 = 8500 #sun radial position in Parsec
for x in np.arange(x_length) - x_offset:
  for y in np.arange(y_length) - y_offset:
    for z in np.arange(z_length) - z_offset:
      x_pos = (x + x_offset) * scale  - (x_length-1)/2*scale #real coordiantes
      y_pos = (y + y_offset) * scale  - (y_length-1)/2*scale
      z_pos = (z + z_offset) * scale  - (z_length-1)/2*scale
      

      if cut_switch == 'inner':
        if calc_r (x_pos,y_pos,z_pos,inc_disk) < r_disk and abs(calc_h(x_pos,y_pos,z_pos,inc_disk)) < z_disk and (R0**2-x_pos**2) >= 0.  and z_pos <= (R0**2-x_pos**2)**0.5: #cutdisk inner
          coords.append([x, y, z])
          abs_coords.append([x_pos,y_pos,z_pos])
      if cut_switch == 'outer':
        if calc_r (x_pos,y_pos,z_pos,inc_disk) < r_disk and abs(calc_h(x_pos,y_pos,z_pos,inc_disk)) < z_disk and (R0**2-x_pos**2) >= 0. and -z_pos >= (R0**2-x_pos**2)**0.5: #cut disk outer down correct-for radiative transfer  
          coords.append([x, y, z])
          abs_coords.append([x_pos,y_pos,z_pos])
      if cut_switch == 'total':
        if calc_r (x_pos,y_pos,z_pos,inc_disk) < r_disk and abs(calc_h(x_pos,y_pos,z_pos,inc_disk)) < z_disk: #and z_pos < R0: #total
          coords.append([x, y, z])
          abs_coords.append([x_pos,y_pos,z_pos])
#print "coords", coords
#print "number\n", number
gbl._globals['compound']['coordinates'] = np.array(coords) #list/tuples to matrix
gbl._globals['compound']['abs_coordinates'] = np.array(abs_coords)
#print 'array', gbl._globals['compound']['coordinates']
gbl._globals['compound']['npoints'] = int(gbl._globals['compound']['coordinates'].size/3)
print('Number of voxels: ', int(gbl._globals['compound']['coordinates'].size/3))
npoints = range(gbl._globals['compound']['npoints'])
#for i in range(gbl._globals['compound']['npoints']):
    #print i
    #if np.alltrue(np.equal([0,0,0], coords[i])):
        #print 'i', i
#a = input('pause')
#gbl._globals['compound']['npoints'] = gbl._globals['compound']['coordinates'].size/3
# print 'number of ensembles', gbl._globals['compound']['npoints']
gbl._globals['binomial_poisson'] = True # 1
# are the clump statistics in the source best described by a 
# binomial or a poisson distribution? 1 means binomial, other
# input means poisson 
gbl._globals['compound']['los_start'] =  np.array([0, 0, -y_length], int)
gbl._globals['compound']['los_end']   =  np.array([0, 0, y_length], int)
# define start and endpoint of line of sight.
# line of sight has to be parallel to z axis (in positive direction)
#gbl._globals['compound']['losoffsets'] = np.array([[-1, -1]]) #obsolete, fit cubes contains all possible LoS
# if not given the default is [0,0]. inputs needs to have the form np.array([[1,1],[2,2],..])
gbl._globals['compound']['I_bg'] = 0 # Background Intensity 
# distance to source in pc. needed for conversion of beam size into pixels when creating maps
gbl._globals['compound']['distance'] = 2.5*10**6 # 414    1kpc
# side length of pixels. value used for line of sight integration (dint = ds). in cm!!!!
gbl._globals['compound']['dint'] = gbl._globals['compound']['pixelsize'] *\
     gbl._globals['constants']['pc']
# parameters for plots
# gbl._globals['plots'] = {}
# set plotrange for geomety plot: 
# [[x_min, x_max], [y_min, y_max], [z_min, z_max]]
gbl._globals['plots']['plotrange'] = np.array([[-x_offset*scale, x_offset*scale], [-y_offset*scale, y_offset*scale ], [-z_offset*scale, z_offset*scale]])
#gbl._globals['plots']['plotrange'] = np.array([[-x_offset*scale, x_offset*scale], [-y_offset*scale, y_offset*scale ], [-z_offset*scale/5., z_offset*scale/5.]])
gbl._globals['sigma'] = {}
gbl._globals['sigma']['sigma_cl_j'] = np.array([1.67/2.3548, 1.67/2.3548, 1.67/2.3548, 1.67/2.3548, 1.67/2.3548], dtype=np.float32) 
# intrinsic clump linewidth (standard deviation)
# im Moment muessen die im code fuer alle massenpunkte identisch sein, da nur der erste Eintrag zur Umrechnung
# in peak werte genutzt wird. AENDERN!!!
gbl._globals['sigma']['sigma_inter_j'] = np.array([1.67/2.3548, 1.67/2.3548, 1.67/2.3548, 1.67/2.3548, 1.67/2.3548], dtype=np.float32)
# intrinsic clump linewidth for clumps representing interclump medium
# internal sigma of clumps used in "line_absorption...py"
#gbl._globals['sigma']['sigma_ens_cl_j'] = np.array([1.1/2.3548, 1.1/2.3548, 1.1/2.3548, 1.1/2.3548], float)
gbl._globals['sigma']['sigma_ens_cl_j'] = np.array([10./2.3548, 10./2.3548, 10./2.3548, 10./2.3548, 10./2.3548], dtype=np.float32)
# ensemble "inter-clump" velocity spread (standard deviation)
#gbl._globals['sigma']['sigma_ens_inter_j'] = np.array([3.635/2.3548, 3.635/2.3548, 3.635/2.3548, 3.635/2.3548], float)
gbl._globals['sigma']['sigma_ens_inter_j'] = np.array([10./2.3548, 10./2.3548, 10./2.3548, 10./2.3548, 10./2.3548], dtype=np.float32)
# ensemble "inter-clump" velocity spread for interclump medium (sigma???)
# (FWHM = 2.3548 sigma)
# information on sampling velocities: range: vel +- d_vel. 
# and number of steps: nstep # km s^-1
gbl._globals['compound']['vel'] = float(0.) # central v
gbl._globals['compound']['d_vel'] = float(350) *np.cos(inc_disk) +10.# max v
gbl._globals['compound']['nstep'] = int(101) #101 #number of v-baskets should be uneven
#gbl._globals['compound']['nstep'] = int(101) #101 #number of v-baskets should be uneven
print('min Velocity Sampling [km/s]: ', -gbl._globals['compound']['d_vel'] + gbl._globals['compound']['vel'])
print('max Velocity Sampling [km/s]: ', gbl._globals['compound']['d_vel'] + gbl._globals['compound']['vel'], '\n')
Dv = 2*gbl._globals['compound']['d_vel'] / (gbl._globals['compound']['nstep']-1)
print('bin-distance km/s: ' , Dv)
# Abtastrate der Geschwindigkeitsbins 4 Koerbe pro Sigma
#bin size *2 sollte  kleiner totaler halbwertsbreite sein
minv= 30 # in km/s
sigma_cl_tot = ( np.array(gbl._globals['sigma']['sigma_cl_j'])**2 + np.array(gbl._globals['sigma']['sigma_ens_cl_j'])**2 + minv**2 )**0.5
gbl._globals['sigma']['sigma_ens_cl_j_total'] = sigma_cl_tot
sigma_inter_tot = (np.array(gbl._globals['sigma']['sigma_inter_j'])**2 + np.array(gbl._globals['sigma']['sigma_ens_inter_j'])**2 + minv**2)**0.5
gbl._globals['sigma']['sigma_inter_j_total'] = sigma_inter_tot
print('Sigma cl tot', sigma_cl_tot[0], 'Sigma_inter', sigma_inter_tot[0], 'km/s')
print('adviced bin-number:',round(gbl._globals['compound']['d_vel']*2 / min(sigma_cl_tot)/2.3548 *6.,0))
gbl._globals['plots']['unit'] = 0
#set to 1 to get advanced map in erg/(cm^2 s sr). any other value
#to leave map in Kkm/s
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
for i in npoints:
  ens.append(pf.Ensemble(i, \
              rho_ens_clumps  = 0. ,\
              rho_ens_inter   = 0. ,\
              Mens_clumps     = 0. ,\
              Mens_inter      = 0. ,\
              Ml_cl           = clump_mass_l ,\
              Mu_cl           = clump_mass_u ,\
              Ml_inter        = inter_mass_l ,\
              Mu_inter        = inter_mass_u ,\
              B               = 1 ,\
              FUV             = 0 ,\
              velocity        = np.array([0, 0, 0], float),\
              v_dispersion    = 0. \
              ))
# Read in rotation curve
convert = 1000 #if r_data are in kpc (*1000) or in pc (*1)
rot_list_r =[]
rot_list_v =[]
#filename = 'compounds/disk/rot_curve.dat'
#with open ('compounds/disk/rot_curve.dat','r') as rot_file:
with open(filename,'r') as rot_file:
  next(rot_file) #skip header
  for line in rot_file:
    rot_list_r.append(line.split()[0])
    rot_list_v.append(line.split()[1])
rot_curve = np.zeros(  (len(rot_list_r),2)) #creates rot array
for i in range(len(rot_list_r)):
  rot_curve[i,0] = float(rot_list_r[i])*convert
  rot_curve[i,1] = float(rot_list_v[i])
#print rot_curve
#pause = input('rot_curve')
# Read in mass profile clumps
mass_list_r =[]
mass_list_v =[]
with open(gbl.COMPOUNDPATH+'mass_profile.dat','r') as mass_file:
  next(mass_file) #skip header
  for line in mass_file:
    mass_list_r.append(line.split()[0])
    mass_list_v.append(line.split()[1])
mass_profile = np.zeros(  (len(mass_list_r),2)) #creates mass array
for i in range(len(mass_list_r)):
  mass_profile[i,0] = float(mass_list_r[i])*convert
  mass_profile[i,1] = float(mass_list_v[i])
#print mass_profile
#pause = input('mass_profile')
# Read in density profile of clumps
dens_list_r =[]
dens_list_v =[]
with open(gbl.COMPOUNDPATH+'densities_clouds.dat','r') as dens_file:
  next(dens_file) #skip header
  for line in dens_file:
    dens_list_r.append(line.split()[0])
    dens_list_v.append(line.split()[1])
dens_profile = np.zeros(  (len(dens_list_r),2)) #creates dens array
for i in range(len(dens_list_r)):
  dens_profile[i,0] = float(dens_list_r[i])*convert
  dens_profile[i,1] = float(dens_list_v[i])
#print dens_profile
#pause = input('dens_profile')
# Read in mass profile interclumpmedium
mass_list_r =[]
mass_list_v =[]
with open(gbl.COMPOUNDPATH+'mass_profile_inter.dat','r') as mass_file:
  next(mass_file) #skip header
  for line in mass_file:
    mass_list_r.append(line.split()[0])
    mass_list_v.append(line.split()[1])
mass_profile_inter = np.zeros(  (len(mass_list_r),2)) #creates mass array
for i in range(len(mass_list_r)):
  mass_profile_inter[i,0] = float(mass_list_r[i])*convert
  mass_profile_inter[i,1] = float(mass_list_v[i])
  #print mass_profile_inter
#pause = input('mass_profile_inter')
# Define speed and masses for npoints
total_clump_mass = 0
total_inter_mass = 0
minimum_roh = 10**10 #just a to high start value
for i in npoints:  
  # define all  coordinates
  x = abs_coords[int(i)][0]
  y = abs_coords[int(i)][1]
  z = abs_coords[int(i)][2]
  r = calc_r (x,y,z,inc_disk)
  h = calc_h (x,y,z,inc_disk)
  phi = np.arctan2 (x,(y**2+z**2)**0.5)
  #print 'x', x, 'y', y, 'z', z
  # Setting Vz
  #interpolate from rot_curve for current r
  v = griddata(rot_curve[:,0], rot_curve[:,1], r, method='linear')
  v_z = v * np.cos(inc_disk) * np.sin(phi) # calc v_z-component
  v_x = v * np.cos(inc_disk) * np.cos(phi) # calc v_x-component
  ens[int(i)].velocity[0] = v_x
  ens[int(i)].velocity[2] = v_z
  if x == 0 and y == 0 and z == 0:
    v = 0 # central velocity just in case it should not be zero
    #print 'v_center', v_z
    #pause = input ( 'v_center' )
  #print "v_z [km/s]", ens[i].velocity[2]
  # voxel mass clump (H2)
  # Gesamt-Dichteprofil in R-Richtung fuer massen der voxel
  mass_surface = griddata(mass_profile[:,0], mass_profile[:,1], r, method='linear') #s-density of voxel (not clumps)
  mass_density = mass_surface / gbl._globals['compound']['pixelsize'] # -> in SM / PC^3
  # interpolate surrounding mass of current r for more precise density and calculate v_dispersion
  bins = 2 #number of interpolation points on axis for precise mass and v, dispersion
  end_dens = 0.
  disp = 0.
  grid_points = ((np.arange(bins) + 0.5)/bins - 0.5 ) * gbl._globals['compound']['pixelsize'] #for interpol
  #print 'interpol', grid_points
  #pause = input('tets')
  for loc_x in grid_points:
    for loc_y in grid_points: #range(1): #grid_points:
      for loc_z in grid_points:
        #print 'x',x, 'y',y, 'z', z, 'r', r, 'r|\n' 
        #print 'locx', loc_x, 'locy', loc_y, 'locz', loc_z
        total_x = x + loc_x
        total_y = y + loc_y
        total_z = z + loc_z
        disk_r = calc_r (total_x,total_y,total_z,inc_disk)
        disk_phi = np.arctan2 (total_x,(total_y**2+total_z**2)**0.5)
        #absolute_r = (absolute_x**2 + absolute_y**2 + absolute_z**2)**0.5
        loc_v =  griddata(rot_curve[:,0], rot_curve[:,1], abs(disk_r), method='linear')
        loc_dens = griddata(mass_profile[:,0], mass_profile[:,1], disk_r, method='linear') #actual surface density
        loc_v_z = loc_v * np.cos(inc_disk) * np.sin(disk_phi) # calc v_z-component                
        #print 'absolute_x', total_x
        #print 'absolute_y', total_y
        #print 'absolute_z', total_z
        #print 'disk_phi', disk_phi * 180. /np.pi
        #print 'disk_R', disk_r
        #print 'loc_v', loc_v
        #print 'loc_dens', loc_dens
        #print 'v_z',v_z # of voxel center
        #print 'loc_v_z', loc_v_z # of interpolated point
        #print 'dispersions beitrag', ((loc_v_z - v_z)**2)**0.5
        disp = disp + ( (loc_v_z - v_z)**2) #add dispersion of interpol point relative to center
        end_dens = end_dens + loc_dens
  end_disp =  (disp / bins**3)**0.5  #divided by number of interpolations
  #print 'dispersion of voxel', end_disp
  end_dens = end_dens / bins**3 / gbl._globals['compound']['pixelsize'] #normalizing to binnumber and pix3D density
  #pause = input('pause')
  ens[int(i)].v_dispersion = end_disp # dispersion by systematic speed which can't be resolved  
  ens[int(i)].Mens_clumps = end_dens * gbl._globals['compound']['pixelsize']**3 # -> total voxel mass with additional interpol    
  ens[int(i)].Mens_clumps = ens[int(i)].Mens_clumps * global_massfactor #* 2
  #print "Voxel_mass:", ens[i].Mens_clumps ,"[Sm]" ,"at r", r      
  #SM/pc->N/cm^3
  #density_cm = ens[i].Mens_clumps * gbl._globals['constants']['M_sol'] / gbl._globals['constants']['M_H'] /gbl._globals['constants']['pc']**3
  # Setting clump H2 surface density n(r)
  clump_rho = griddata(dens_profile[:,0], dens_profile[:,1], r, method='linear')
  if x == 0 and y == 0 and z == 0:
    clump_rho = 15000. # set centermass
  ens[int(i)].rho_ens_clumps = clump_rho # N/cm^3
  ens[int(i)].rho_ens_clumps = ens[int(i)].rho_ens_clumps * global_densityfactor  #densities for different test cases
  if ens[int(i)].rho_ens_clumps < minimum_roh: minimum_roh = ens[int(i)].rho_ens_clumps #memorises lowest roh
  if np.isnan(clump_rho) == True:
    print('\nVoxel number: ', i)
    print('r, ', r, '\nh: ', h)
    pause = input('NaN in clumpdensity found')
  #print "Clump density in N/cm^3", ens[i].rho_ens_clumps , "at r in kpc", r , "h in kpc", h 
  # interclump (H1)
  # Gesamt-Dichteprofil in R-Richtung fuer massen der voxel
  mass_surface_inter = griddata(mass_profile_inter[:,0], mass_profile_inter[:,1], r, method='linear')
  mass_density_inter = mass_surface_inter / gbl._globals['compound']['pixelsize'] # -> in SM / PC^3
  ens[int(i)].Mens_inter = mass_density_inter * gbl._globals['compound']['pixelsize']**3
  #print  ens[i].Mens_inter, ''
  ens[int(i)].Mens_inter = mass_density_inter * gbl._globals['compound']['pixelsize']**3
  if inter_there == 'yes':
    #totalmass = ens[i].Mens_clumps #has to distribute mass in clump and inter mass
    #ens[i].Mens_clumps = totalmass * massratio_clump_to_inter # mass clumps
    #ens[i].Mens_inter  = totalmass * (1-massratio_clump_to_inter) # mass intermedium
    ens[int(i)].rho_ens_inter = 1911#const minimum density=1000  ens[i].rho_ens_clumps  ####density inter same as clump
  # check if mass excluded
  if clump_there == 'no' : ens[int(i)].Mens_clumps = 0
  if inter_there == 'no' : ens[int(i)].Mens_inter = 0
  total_clump_mass = total_clump_mass + ens[int(i)].Mens_clumps
  total_inter_mass = total_inter_mass + ens[int(i)].Mens_inter
print('Total clump-mass: ', total_clump_mass/10**9, ' in G M_sol')
print('Total inter-mass: ', total_inter_mass/10**9, ' in G M_sol')
print('Minimum ens-rho: ', minimum_roh)
gbl._globals['compound']['ens'] = ens #save ensamble to globals
# END OF DEFINITIONS
# set axis labels for 3D geometry plots. Do not need to be changed.
gbl._globals['plots']['xlabel'] = 'X offset [' + str(gbl._globals['compound']['pixelsize']) + 'pc]'
gbl._globals['plots']['ylabel'] = 'Y offset [' + str(gbl._globals['compound']['pixelsize']) + 'pc]'
gbl._globals['plots']['zlabel'] = 'to observer [' + str(gbl._globals['compound']['pixelsize']) +  'pc]'

gbl._globals['compound']['dgrid'] = 1 
# factor defining streching of grid in dgrid
# currently not used
# DO NOT CHANGE
## Setting FUV-field
# read in mass profile interclumpmedium
uv_list_r =[]
uv_list_v =[]
with open(gbl.COMPOUNDPATH+'galactic_FUV.dat','r') as file:
  next(file) #skip header
  for line in file:
    uv_list_r.append(float(line.split()[0]))
    uv_list_v.append(float(line.split()[1]))
uv_r = np.array (uv_list_r)
uv_v = np.array (uv_list_v)
#print uv_r
#print uv_v
#print type(uv_r[0])
#pause = input('UV OK?')
norm_UV = 1.06*1e39 * 2.7305
calc_fuv = 1
if calc_fuv == 1:
  print('Calulation disk-FUV')
  fuv_file = open (gbl.KOSMAPATH+'temp/FUV.dat','w')
  for i in npoints:
    # print i
    x = abs_coords[int(i)][0]
    y = abs_coords[int(i)][1]
    z = abs_coords[int(i)][2]
    r = calc_r (x,y,z,inc_disk)
    h = calc_h (x,y,z,inc_disk)
    #phi = np.arctan2 (x,(y**2+z**2)**0.5)
    #print r
    r_OB = 1/6**0.5 *1000 # distance on average of OBs in pc
    d_MC = rand.random()*30+20  #random for distances of the clouds to next OB assosiation in pc (20-50) 
    UV = (r_OB**2 / d_MC**2)/3
    UV = global_UV #if constant and not random method
    if UV_interpol == 'yes':
      temp = griddata( uv_r, uv_v, r, method='linear') #UV interpol for actual position
      #print 'temp', temp
      UV = temp/norm_UV*global_UV
      #print 'r',r, 'FUV', UV
      if UV < 1.:
        UV = 1.
      #pause = input('UV OK?')
    #print "lokal UV-field",UV
    #write out like in Silkes method: (should become mostly obsolete)        
    fuv_file.write(str(UV) + str(" ") + str(  gbl._globals['compound']['coordinates'][int(i)] ) +"\n"  )
    #fill directly into ensamble, and create fitsfile later with:
    ens[int(i)].FUV = UV
  fuv_file.close()
# Writing fits files
write_fits = 1
if write_fits == 1:
  print('Creating .fits files . . .')
  mass_data = np.zeros( (z_length,y_length,x_length) , dtype=np.float32)  #initialize fit grids with zeros
  mass_clump_data = np.zeros( (z_length,y_length,x_length) , dtype=np.float32)
  mass_inter_data = np.zeros( (z_length,y_length,x_length) , dtype=np.float32)
  vz_data = np.zeros( (z_length,y_length,x_length), dtype=np.float32 )
  density_data = np.zeros( (z_length,y_length,x_length), dtype=np.float32 )
  FUV_data = np.zeros( (z_length,y_length,x_length), dtype=np.float32 )
  # TEST
  disp_data = np.zeros( (z_length,y_length,x_length), dtype=np.float32 )
  for i in npoints: # fill grids with values
    x = int(coords[int(i)][0]+x_offset) #transfrom to array coordinates
    y = int(coords[int(i)][1]+y_offset)
    z = int(coords[int(i)][2]+z_offset)
    #x_i = abs_coords[i][0]/scale + (gbl._globals['compound']['mapSize']['x']-1)/2. #calculates correct array positions
    #y_i = abs_coords[i][1]/scale + (gbl._globals['compound']['mapSize']['y']-1)/2.
    #z_i = abs_coords[i][2]/scale + (gbl._globals['compound']['mapSize']['z']-1)/2.
    #print "z", z - z_offset
    mass_data[z][y][x] = ens[int(i)].Mens_clumps + ens[int(i)].Mens_inter # write mass data to array coordiantes
    mass_clump_data[z][y][x] = ens[int(i)].Mens_clumps
    mass_inter_data[z][y][x] = ens[int(i)].Mens_inter
    vz_data[z][y][x] = ens[int(i)].velocity[2] # velocitydata in z
    density_data[z][y][x] = ens[int(i)].rho_ens_clumps # densitydata
    FUV_data[z][y][x] = ens[int(i)].FUV #FUV data
    disp_data[z][y][x] = ens[int(i)].v_dispersion
  hdu = fits.PrimaryHDU(mass_data) # sets data as hdu-list
  # Header content
  hdu.header['OBSERVER'] = "Just a Simulation"
  hdu.header['NAME'] = "galactic Disk"
  hdu.header['COMMENT'] = " easter egg Fr. 13.5 "
  #hdu.header['BZERO'] = 0
  #hdu.header['BSCALE'] = 1

  hdu.header['CRPIX1'] = 1   #ref Pixel
  hdu.header['CTYPE1'] = 'x-dim [pc]' #Name 
  hdu.header['CRVAL1'] = ((-x_length+1)/2.)  *scale #Reference Value in kpc
  hdu.header['CDELT1'] = 1. * scale #stepsize per counting index in pc

  hdu.header['CRPIX2'] = 1
  hdu.header['CTYPE2'] = "y-dim [pc]"
  hdu.header['CRVAL2'] = ((-y_length+1)/2.)  *scale
  hdu.header['CDELT2'] = 1. * scale
 
  hdu.header['CRPIX3'] = 1
  hdu.header['CTYPE3'] = "z-dim [pc]"
  hdu.header['CRVAL3'] = ((-z_length+1)/2.)  *scale
  hdu.header['CDELT3'] = 1. * scale 

  # Clearing old files
  open(gbl.FITSPATH+'mass.fits', 'w').close()
  open(gbl.FITSPATH+'mass_clump.fits', 'w').close()
  open(gbl.FITSPATH+'mass_inter.fits', 'w').close()
  open(gbl.FITSPATH+'vz.fits', 'w').close()
  open(gbl.FITSPATH+'density.fits', 'w').close()
  open(gbl.FITSPATH+'FUV.fits', 'w').close()
  open(gbl.FITSPATH+'disp.fits', 'w').close()

# First writing
  hdu.header['BTYPE'] = 'all mass'
  hdu.header['BUNIT'] = " SM/voxel"
  hdu.writeto(gbl.FITSPATH+'mass.fits', overwrite=True)

  hdu.data = mass_clump_data
  hdu.header['BTYPE'] = 'clump mass'
  hdu.writeto(gbl.FITSPATH+'mass_clump.fits', overwrite=True)

  hdu.data = mass_inter_data
  hdu.header['BTYPE'] = 'inter mass'
  hdu.writeto(gbl.FITSPATH+'mass_inter.fits', overwrite=True)

  # Velocity
  hdu.header['BTYPE'] = 'velocity'
  hdu.header['BUNIT'] = " km/s"
  hdu.data = vz_data
  hdu.writeto(gbl.FITSPATH+'vz.fits', overwrite=True)
  # Density
  hdu.header['BTYPE'] = 'density of clouds'
  hdu.header['BUNIT'] = " particle / cm^3"
  hdu.data = density_data
  hdu.writeto(gbl.FITSPATH+'density.fits', overwrite=True)
  # FUV
  hdu.header['BTYPE'] = 'FUV-field'
  hdu.header['BUNIT'] = " drain-field"
  hdu.data = FUV_data
  hdu.writeto(gbl.FITSPATH+'FUV.fits', overwrite=True)
  # V_dispersion
  hdu.header['BTYPE'] = 'sigma voxel unresolved systematic vz_disp'
  hdu.header['BUNIT'] = " km/s"
  hdu.data = disp_data
  hdu.writeto(gbl.FITSPATH+'disp.fits', overwrite=True)

gbl._globals['runtimes']['end_model'] = time.time()
# Write out voxel sanity check
write_voxel = 1
if write_voxel == 1:
  print('Writing voxel...')
  voxel_file = open(gbl.TEMPPATH+'voxel.dat','w')
  for i in npoints:
    # print i
    x_i = coords[int(i)][0]
    y_i = coords[int(i)][1]
    z_i = coords[int(i)][2]
    x = abs_coords[int(i)][0]
    y = abs_coords[int(i)][1]
    z = abs_coords[int(i)][2]
    r = calc_r (x,y,z,inc_disk)
    h = calc_h (x,y,z,inc_disk)
    phi = np.arctan2 (x,(y**2+z**2)**0.5)
    #print r
    k=1 #kommastellen
    voxel_file.write( 'x,y,z' + str(  gbl._globals['compound']['abs_coordinates'][int(i)] ) +\
    ' r '+str(round(r,k)) + ' phi ' + str(round(phi/np.pi*180,k)) + ' vz' + str(round(ens[int(i)].velocity[2],k)) + \
    ' vz_disp ' + str(round(ens[int(i)].v_dispersion,k)) + \
    ' rho_clumps [/cm^3] ' + str(round(ens[int(i)].rho_ens_clumps,k)) + ' mens_clump [SM]' + str(round(ens[int(i)].Mens_clumps,k))  + \
    "\n"  )
  voxel_file.close()
print('model.py durchgelaufen!')
#pause = input('OK?')
