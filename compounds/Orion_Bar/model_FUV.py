# model_FUV
#########################################################################
# calculate incident FUV field strength for each voxel (in Units of Draine
# field) based on the Luminosity of the illuminating source.
# The luminosity is first calculated basen on the r^-2 - law and then
# corrected for absorption by dust
# Currently FUV field strength < 1 Draine Field are set to 1 because
# input for clumps with fields < 1 is not provided (grid needs to be
# changed
#########################################################################
# THIS PART (model_FUV.py) OF THE CODE NEEDS TO BE CHANGED IF:
# - the position of the illuminating source needs to be changed -> pos_trap
# - the strength of the illuminating source needs to be changed -> luminosity
# - the size of the sub-voxel grid needs to be changed -> subpixel_step
# - the distance between the sample points on the connecting line between
#   FUV source and needs to be changed -> step
#########################################################################
#########################################################################
import numpy as np
from numpy.linalg import norm
import globals as gbl

ds = gbl._globals['compound']['pixelsize'] # pixel size in pc
coords = gbl._globals['compound']['coordinates'] #coordinates
coords_index = gbl._globals['compound']['coords_index']
x_offset = gbl._globals['compound']['offsets']['x']
y_offset = gbl._globals['compound']['offsets']['y']
z_offset = gbl._globals['compound']['offsets']['z']
#########################################################################
################### parameters to be adjusted ###########################

step = 1/20.
# distance between sample points measured in pixels
subpixel_step = 3 
# finer gitter for FUV calculation

# position of Theta 1 Ori C. z position unsure. possibilities
pos_trap = np.array([0, 22.3, 40.]) 
# z-position from pellegrini 
# pos_trap = np.array([0, 11.15, 9.1]) 
# z-position from van der werf
# pos_trap = np.array([0, 11.15, 30.]) 
# directly in front of bar (upper edge)

#  I found : 1.4 * 10**4
for i in np.arange(gbl._globals['compound']['npoints']):
    gbl._globals['compound']['ens'][i].FUV = max(5.5 * 10**4 * ((0.223/ds)/\
        (np.linalg.norm(pos_trap - \
         np.array(gbl._globals['compound']['coordinates'][i]))))**2., 1)
# FUV intensity (before attenuation by dust) ; units of Draine field
# 0.223 is the distance to the source in pc
# point of reference: 4.4*10**4 draine fields flux for a distance of
# 0.223pc. all other distances need to be scaled.

#########################################################################


#######################################################################
 
# correct for fuv absorption
#### read in extinction A(lambda)
# from scipy.interpolate import interp1d
# lamb = []
# Al = []
# with open("input_files/AlAV.dat","r") as A:
#     i = 0
#     for line in A:
#         i = i+1
#         if i > 3:
#             l=line.split()
#             lamb.append(float(l[0]))
#             Al.append(float(l[1]))
# print Al[1]
# print lamb[0]
# lamb.reverse()
# Al.reverse()
# Alint = interp1d(lamb, Al)
"""
print 'coords_index..'
print coords_index[x_offset-3][32-1][10]
print coords_index[x_offset-3][32-1][9]
print coords_index[x_offset-3][32-1][8]
print coords_index[x_offset-3][32-1][7]
print coords_index[x_offset-3][32-1][6]
print coords_index[x_offset-3][32-1][5]
print coords_index[x_offset-3][32-1][10]
print coords_index[x_offset-4][32-1][9]
print coords_index[x_offset-4][32-1][8]
print coords_index[x_offset-4][32-1][7]
print coords_index[x_offset-4][32-1][6]
print coords_index[x_offset-4][32-1][5]
# pause = raw_input('..')
#print coords[coords_index[23][3+32][6]]
print coords[np.round(coords_index[x_offset-3][32-1][10])]
print coords[np.round(coords_index[x_offset-3][32-1][9])]
print coords[np.round(coords_index[x_offset-3][32-1][8])]
print coords[np.round(coords_index[x_offset-3][32-1][7])]
"""
# pause = raw_input('..')
#######################################################################
# calculate H-column density between star and different pixels ########
for i in np.arange(gbl._globals['compound']['npoints']):
#for i in [11134, 11133, 11132, 11131,10600, 10599,10598, 10597]:
#for i in [11134, 10600]:
#for i in np.arange(2):#+12760:
    print('pixel number', i)
    print("coords[i]", coords[i])
    # pause=raw_input('..')
    A = 0 # hydrogen column density
    subcoords = []
    if subpixel_step <= 1 : subcoords.append(coords[i])
    else:
        for subx in np.arange(subpixel_step):
            for suby in np.arange(subpixel_step):
                for subz in np.arange(subpixel_step):
                    subcoords.append([coords[i][0]-0.5+subx/(subpixel_step-1.),\
                                      coords[i][1]-0.5+suby/(subpixel_step-1.),\
                                      coords[i][2]-0.5+subz/(subpixel_step-1.)])
    print('subcoords', subcoords)
    # pause= raw_input('...')
    Apixel = 0
    for subc in subcoords:
        found = 0 
        # first pixel not yet found
        direc = np.array([subc[0], subc[1], subc[2]]) - \
            np.array([pos_trap[0], pos_trap[1], pos_trap[2]])
        distance = norm(direc)
        # distance between star and pixel
        direc = direc/norm(direc)
        # normalised direction vector pointing from star to pixel 
        # print 'distance', distance
        # print 'direc', direc
        # print 'step', step
        # pause = raw_input('...')
        for j in np.arange(round(distance/step)) + 1:
            position = np.round(pos_trap + direc * j * step)
            # position on connecting line
            # print 'j, position', j, position
            try: pos_coords = coords_index[position[0] + x_offset]\
                              [position[1] + y_offset][position[2] + z_offset]
            except IndexError: pos_coords = -1
            # print 'pos_coords', pos_coords
            if pos_coords != -1: 
                if found == 0: 
                    found = 1
                    # if np.alltrue(np.equal(position, coords[i])): break
                    # break loop if first position is final pixel position
                    # (otherwise first pixel might be found twice)
                elif found == 1: 
                    # pause = raw_input('..A..')
                    # print 'pos_coords', pos_coords
                    Apixel = Apixel + \
                      gbl._globals['compound']['ens'][int(pos_coords)].Afuv_tot*step
    A = Apixel/float(subpixel_step**3)
    print('A_final', A)
    # pause = raw_input('...')
    print('FUV without correction', gbl._globals['compound']['ens'][i].FUV)
    # pause = raw_input('...')
    gbl._globals['compound']['ens'][i].FUV = \
            max(gbl._globals['compound']['ens'][i].FUV*10**(-A/2.5), 1)
    print('FUV corrected', gbl._globals['compound']['ens'][i].FUV)
    # pause = raw_input('...')

########### correct single pixel fuv intensity ######################
"""
theta = 15 # angle between los and edge-on part of the pdr. in degree
coords = gbl._globals['compound']['coordinates']
for i in np.arange(gbl._globals['compound']['npoints']):
#for i in np.arange(10):
    # figure out how deep the pixel lies inside the pdr
    print "coords[i]", coords[i]
    if pos_trap[2] == coords[i][2]:
        # star and pixel have the same z position
        intersec = -np.tan(theta * np.pi/180.) * coords[i][2]
        distance = abs(intersec - coords[i][1]) * ds
        print 'distance in pc', distance
    elif coords[i][2] > pos_trap[2]:
        direc = np.array([coords[i][0], coords[i][1], coords[i][2]]) - \
                np.array([pos_trap[0], pos_trap[1], pos_trap[2]])
        direc = direc/norm(direc)
        print 'direc', direc
        la = (-pos_trap[1] - pos_trap[2] * np.tan(theta * np.pi/180.))\
            /(np.tan(theta * np.pi/180.) * direc[2] + direc[1])
        print 'la', la
        intersec = np.array([pos_trap[0] + la * direc[0], \
                     pos_trap[1] + la * direc[1], pos_trap[2] + la*direc[2]])
        print 'intersec', intersec
        distance = norm(intersec - coords[i]) * ds
        print 'distance in pc', distance
    else:
        direc = np.array([coords[i][0], coords[i][1], coords[i][2]]) - \
                np.array([pos_trap[0], pos_trap[1], pos_trap[2]])
        direc = direc/norm(direc)
        yzero = pos_trap[1] - (pos_trap[2]/float(direc[2])) *\
                float(direc[1])
        if yzero > 0:
            print 'yzero > 0'
            xzero = pos_trap[0] - (pos_trap[2]/float(direc[2])) *\
                    float(direc[0])
            intersec = np.array([xzero, yzero, 0])
            print 'intersec', intersec
            distance = norm(intersec - coords[i]) * ds
            print 'distance in pc', distance
        else:
            print 'yzero < 0'
            la = (-pos_trap[1] - pos_trap[2] * np.tan(theta * np.pi/180.))\
                /(np.tan(theta * np.pi/180.) * direc[2] + direc[1])
            print 'la', la
            intersec = np.array([pos_trap[0] + la * direc[0], \
                     pos_trap[1] + la * direc[1], pos_trap[2] + la*direc[2]])
            print 'intersec', intersec
            distance = norm(intersec - coords[i]) * ds
            print 'distance in pc', distance
    
    # A = AlAV * 5.3 * 10**(-22) * 10**4 * distance * gbl._globals['constants']['pc']
    A = AlAV * 9.92 * 10**(-22) * 10**5 * distance * gbl._globals['constants']['pc']
    print 'FUV', gbl._globals['compound']['ens'][i].FUV
    # gbl._globals['compound']['ens'][gbl._globals['compound']['npoints']-i-1].FUV = max(gbl._globals['compound']['ens'][i].FUV*10**(-A/2.5), 1)
    gbl._globals['compound']['ens'][i].FUV = max(gbl._globals['compound']['ens'][i].FUV*10**(-A/2.5), 1)
    # gbl._globals['compound']['ens'][i].FUV = max(gbl._globals['compound']['ens'][i].FUV*10**(-1000/2.5), 1)
    # gbl._globals['compound']['ens'][i].FUV =10**4
    print 'FUV corrected', gbl._globals['compound']['ens'][i].FUV
"""
