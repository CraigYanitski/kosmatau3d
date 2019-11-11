from VoxelGrid import *
from Dimensions import *
import numpy as np
d=Dimensions(10,10,10)
v=VoxelGrid(d)
##
#import os
#os.chdir('/home/craig/Desktop/Köln/kosma-tau^3-develop/kosma-tau-3d/classes')
from Model import *
species = ['C 1', 'CO 1', 'C+ 1', '13CO 1', 'O 2']
kosma = Model(10,10,2,'disk')
kosma.addSpecies(species)
kosma.initialiseModel()
kosma.calculateEmission()
##
for i in range(kosma.getGrid().getDimensions().voxelNumber()): print(i, ':', kosma.getGrid().allVoxels()[i].getPosition())