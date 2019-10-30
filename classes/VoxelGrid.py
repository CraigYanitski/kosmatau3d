class VoxelGrid(Model):
  '''
  This is a class to handle all of the voxels in KOSMA-tau^3. It contains a
  specified arrangement of voxels, and must coordinate with the Dimensions class
  to make the Shape class functional.
  '''
  # PRIVATE
  def __init__(self, indeces=1):
    for i in indeces: self.__voxels.append(Voxel(i))
    self.__map = {}       #dictionary object to map the voxel indeces to the correct location
    return

  # PUBLIC
  #def createGrid(self, indeces):
  #  for i in indeces: self.__voxels.append(Voxel(i))
  #  return
  def initialiseVoxels(self):
    for voxel in self.__voxels:
      voxel.createClump()
      voxel.createInterclump()
      voxel.
    return
  def allVoxels(self):
    return self.__voxels






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
    #UV_interpol = 'yes'
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