from VoxelGrid import *
from Dimensions import *
import numpy as np
d=Dimensions(10,10,10)
v=VoxelGrid(d)
##
from Model import *
species = ['13CO 10', 'C+ 1', 'CO 1']#, 'C+ 1', '13CO 1', 'O 2', 'Dust 1', 'Dust 2', 'Dust 3', 'Dust 4', 'Dust 5']
kosma = Model(12, 12, 0,modelType='disk',resolution=1000)
kosma.addSpecies(species)
kosma.initialiseModel()
kosma.calculateEmission()
kosma.setLOS()
kosma.calculateObservation()
##
for i in range(kosma.getGrid().getVoxelNumber()): print(i, ':', kosma.getGrid().allVoxels()[i].getPosition())
##
for i in range(kosma.getGrid().getVoxelNumber()): print(kosma.getGrid().allVoxels()[i])