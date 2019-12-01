import os
import sys
import inspect
pathname = inspect.getframeinfo(inspect.currentframe()).filename
print(os.path.dirname(pathname))
sys.path.append(os.path.dirname(pathname)+'/classes/')
from Model import *
species = ['13CO 10', 'C+ 1', 'CO 1', '13CO 1', 'O 2', 'Dust 1', 'Dust 2', 'Dust 3', 'Dust 4', 'Dust 5']
kosma = Model(36.2, 36.2, 0,modelType='disk',resolution=1000)
kosma.addSpecies(species)
kosma.initialiseModel()
kosma.calculateEmission()
kosma.setLOS()
kosma.calculateObservation()
##
for i in range(kosma.getGrid().getVoxelNumber()): print(i, ':', kosma.getGrid().allVoxels()[i].getPosition())
##
for i in range(kosma.getGrid().getVoxelNumber()): print(kosma.getGrid().allVoxels()[i])