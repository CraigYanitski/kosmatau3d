import os
import sys
import inspect
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
pathname = inspect.getframeinfo(inspect.currentframe()).filename
print(os.path.dirname(pathname))
sys.path.append(os.path.abspath(os.path.dirname(pathname))+'/classes/')
from Model import *
species = ['13CO 10', 'C+ 1', 'CO 1']#, '13CO 1', 'O 2', 'Dust 1', 'Dust 2', 'Dust 3', 'Dust 4', 'Dust 5']
kosma = Model(36, 36, 3,modelType='disk',resolution=500)
kosma.addSpecies(species)
kosma.initialiseModel()
kosma.calculateEmission()
#for voxel in kosma.getGrid().allVoxels():
#  position = voxel.getPosition()
#  print('Voxel {} at position x={}, y={}, z={}.'.format(voxel.getIndex(), position[0], position[1], position[2]))
#kosma.setLOS()
kosma.calculateObservation()
#kosma.printIntensityMap()
positions, intensity = kosma.getIntensityMap()
positions = positions.T
positions = kosma.getGrid().getVoxelPositions()
intensity = kosma.getGrid().totalEmission()[0]
fig = plt.figure()
ax = fig.add_subplot(111, projection='3d')
ax.scatter(positions[0], positions[1], positions[2], c=(intensity.max(2)).sum(1)/intensity.max(), cmap=plt.cm.jet)
plt.show()
#print('at position {}, the intensity is:\n{}\n'.format(position[50],intensity[50]))
#kosma.getGrid().printVoxels()
##
#for i in range(kosma.getGrid().getVoxelNumber()): print(i, ':', kosma.getGrid().allVoxels()[i].getPosition())
##
#for i in range(kosma.getGrid().getVoxelNumber()): print(kosma.getGrid().allVoxels()[i])