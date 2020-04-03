import os
import sys
import inspect
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D

pathname = inspect.getframeinfo(inspect.currentframe()).filename
sys.path.append(os.path.abspath(os.path.dirname(pathname))+'/classes/')

from Model import *
import constants

complete = False

modelFlag = False
rtFlag = True
cyplotFlag = False

x = 6
y = 6
z = 2

shape = 'disk'

resolution = 400

constants.changeDirectory('MilkyWay')

modelFolder = 'r400_n38094/'

# Factors
constants.clumpMassFactor = 1
constants.interclumpMassFactor = 1
constants.FUVFactor = 1
constants.DensityFactor = 1

# Constant
constants.interclumpLogFUV = 1

# Model masses
constants.clumpLogMassNum = 4
constants.clumpLogMassRange = [-1, 2]
constants.interclumpLogMassNum = 2
constants.interclumpLogMassRange = [-3, -2]

print('KOSMA-tau^3')

species = ['13CO 10', 'C+ 1', 'CO 1', 'CO 2', 'CO 3', 'CO 4', 'CO 5', 'CO 6', 'CO 7', 'CO 8', 'CO 9', 'CO 10', '13CO 1', '13CO 2', '13CO 3', '13CO 4', '13CO 5', '13CO 6', '13CO 7', '13CO 8', '13CO 9', '13CO 10', 'O 2']
kosma = Model(x, y, z, modelType=shape, resolution=resolution)
kosma.addSpecies(species)

if modelFlag:

  kosma.calculateModel()
  kosma.writeEmission()
  
  modelFolder = constants.history

# Calculate integrated intensity maps

if rtFlag:

  import radiativeTransfer

  # radiativeTransfer.plotModel(plot='velocity', directory='r1000.0_n3015/')

  radiativeTransfer.calculateObservation(directory=modelFolder, dim='spherical')

# Create cygrid images

if cyplotFlag:

  import cyplot

  images,wcs = cyplot.convertMap(modelFolder, input_coord='spherical')
  print(images.shape)

  np.nan_to_num(images)

  for i in range(images[:,0,0,0].size):
    plt.imshow(images[i,333,:,:])
    plt.show()

# if False:
#   while(selection==0):

#     print('\nSelect a process')
#     print('1 -> Select molecules <INCOMPLETE>')
#     print('2 -> Calculate intrinsic emission')
#     print('3 -> Write intrinsic emission')
#     print('4 -> Calculate observed emission')
#     print('5 -> Show plots')
#     print('6 -> Exit')

#     selection = input()

#     if selection=='1':
#       selection = 0
#       if complete:
#         while(selection==0):
#           print('Select molecules (separate by commas)')
#           print('1 -> C+ 1->0')
#           print('2 -> C 1->0')
#           print('3 -> C 2->1')
#           print('4 -> C 3->2')
#           print('5 -> O 1->0')
#           print('6 -> O 2->1')
#           print('7 -> O 3->2')
#           print('8 -> CO 1->0')
#           print('9 -> CO 2->1')
#           print('10 -> CO 3->2')
#           print('11 -> CO 4->3')
#           print('12 -> CO 5->4')
#           print('13 -> CO 6->5')
#           print('14 -> CO 7->6')
#           print('15 -> CO 8->7')
#           print('16 -> CO 9->8')
#           print('17 -> CO 10->9')
#           print('18 -> CO 11->10')
#           print('19 -> CO 12->11')
#           print('20 -> CO 13->12')
#           print('21 -> CO 14->13')
#           print('22 -> CO 15->14')
#           print('23 -> CO 16->15')
#           print('24 -> CO 17->16')
#           print('25 -> CO 18->17')
#           print('26 -> CO 19->18')
#           print('27 -> CO 20->19')
#           print('28 -> CO 21->20')
#           print('29 -> CO 22->21')
#           print('30 -> CO 23->22')
#           print('31 -> CO 24->23')
#           print('32 -> CO 25->24')
#           print('33 -> CO 26->25')
#           print('34 -> CO 27->26')
#           print('35 -> CO 28->27')
#           print('36 -> CO 29->28')
#           print('37 -> CO 30->29')
#           print('38 -> CO 31->30')
#           print('39 -> CO 32->31')
#           print('40 -> CO 33->32')
#           print('41 -> CO 34->33')
#           print('42 -> CO 35->34')
#           print('43 -> CO 36->35')
#           print('44 -> CO 37->36')
#           print('45 -> CO 38->37')
#           print('46 -> CO 39->38')
#           print('47 -> CO 40->39')
#           print('48 -> CO 41->40')
#           print('49 -> CO 42->41')
#           print('50 -> CO 43->42')
#           print('51 -> CO 44->43')
#           print('52 -> CO 45->44')
#           print('53 -> CO 46->45')
#           print('54 -> CO 47->46')
#           print('55 -> CO 48->47')
#           print('56 -> CO 49->48')
#           print('57 -> 13CO 1->0')
#           print('58 -> 13CO 2->1')
#           print('59 -> 13CO 3->2')
#           print('60 -> 13CO 4->3')
#           print('61 -> 13CO 5->4')
#           print('62 -> 13CO 6->5')
#           print('63 -> 13CO 7->6')
#           print('64 -> 13CO 8->7')
#           print('65 -> 13CO 9->8')
#           print('66 -> 13CO 10->9')
#           print('67 -> 13CO 11->10')
#           print('68 -> 13CO 12->11')
#           print('69 -> 13CO 13->12')
#           print('70 -> 13CO 14->13')
#           print('71 -> 13CO 15->14')
#           print('72 -> 13CO 16->15')
#           print('73 -> 13CO 17->16')
#           print('74 -> 13CO 18->17')
#           print('75 -> 13CO 19->18')
#           print('76 -> 13CO 20->19')
#           print('77 -> 13CO 21->20')
#           print('78 -> 13CO 22->21')
#           print('79 -> 13CO 23->22')
#           print('80 -> 13CO 24->23')
#           print('81 -> 13CO 25->24')
#           print('82 -> 13CO 26->25')
#           print('83 -> 13CO 27->26')
#           print('84 -> 13CO 28->27')
#           print('85 -> 13CO 29->28')
#           print('86 -> 13CO 30->29')
#           print('87 -> 13CO 31->30')
#           print('88 -> 13CO 32->31')
#           print('89 -> 13CO 33->32')
#           print('90 -> 13CO 34->33')
#           print('91 -> 13CO 35->34')
#           print('92 -> 13CO 36->35')
#           print('93 -> 13CO 37->36')
#           print('94 -> 13CO 38->37')
#           print('95 -> 13CO 39->38')
#           print('96 -> 13CO 40->39')
#           print('97 -> 13CO 41->40')
#           print('98 -> 13CO 42->41')
#           print('99 -> 13CO 43->42')
#           print('100 -> 13CO 44->43')
#           print('101 -> 13CO 45->44')
#           print('102 -> 13CO 46->45')
#           print('103 -> 13CO 47->46')
#           print('104 -> 13CO 48->47')
#           print('105 -> 13CO 49->48')
#           print('106 -> 13C+ 1->0')
#           print('107 -> 13C 1->0')
#           print('108 -> 13C 2->1')
#           print('109 -> 13C 3->2')
#           print('110 -> HCO+ 1->0')
#           print('111 -> HCO+ 2->1')
#           print('112 -> HCO+ 3->2')
#           print('113 -> HCO+ 4->3')
#           print('114 -> HCO+ 5->4')
#           print('115 -> HCO+ 6->5')
#           print('116 -> HCO+ 7->6')
#           print('117 -> HCO+ 8->7')
#           print('118 -> HCO+ 9->8')
#           print('119 -> HCO+ 10->9')
#           print('120 -> HCO+ 11->10')
#           print('121 -> HCO+ 12->11')
#           print('122 -> HCO+ 13->12')
#           print('123 -> HCO+ 14->13')
#           print('124 -> HCO+ 15->14')
#           print('125 -> H13CO+ 1->0')
#           print('126 -> H13CO+ 2->1')
#           print('127 -> H13CO+ 3->2')
#           print('128 -> H13CO+ 4->3')
#           print('129 -> H13CO+ 5->4')
#           print('130 -> H13CO+ 6->5')
#           print('131 -> H13CO+ 7->6')
#           print('132 -> H13CO+ 8->7')
#           print('133 -> H13CO+ 9->8')
#           print('134 -> H13CO+ 10->9')
#           print('135 -> H13CO+ 11->10')
#           print('136 -> H13CO+ 12->11')
#           print('137 -> H13CO+ 13->12')
#           print('138 -> H13CO+ 14->13')
#           print('139 -> H13CO+ 15->14')
#           print('140 -> H13CO+ 16->15')
#           print('141 -> H13CO+ 17->16')
#           print('142 -> H13CO+ 18->17')
#           print('143 -> H13CO+ 19->18')
#           print('144 -> H13CO+ 20->19')
#           print('145 -> H13CO+ 21->20')
#           print('146 -> H13CO+ 22->21')
#           print('147 -> H13CO+ 23->22')
#           print('148 -> H13CO+ 24->23')
#           print('149 -> H13CO+ 25->24')
#           print('150 -> H13CO+ 26->25')
#           print('151 -> H13CO+ 27->26')
#           print('152 -> H13CO+ 28->27')
#           print('153 -> H13CO+ 29->28')
#           print('154 -> H13CO+ 30->29')
#           print('155 -> H3O+ 1->0')
#           print('156 -> H3O+ 2->1')
#           print('157 -> H3O+ 3->2')
#           print('158 -> H3O+ 4->3')
#           print('159 -> H3O+ 5->4')
#           print('160 -> H3O+ 6->5')
#           print('161 -> H3O+ 7->6')
#           print('162 -> H3O+ 8->7')
#           print('163 -> H3O+ 9->8')
#           print('164 -> H3O+ 10->9')
#           print('165 -> H3O+ 11->10')
#           print('166 -> H3O+ 12->11')
#           print('167 -> H3O+ 13->12')
#           print('168 -> H3O+ 14->13')
#           print('169 -> H3O+ 15->14')
#           print('170 -> H3O+ 16->15')
#           print('171 -> H3O+ 17->16')
#           print('172 -> C18O 1->0')
#           print('173 -> C18O 2->1')
#           print('174 -> C18O 3->2')
#           print('175 -> C18O 4->3')
#           print('176 -> C18O 5->4')
#           print('177 -> C18O 6->5')
#           print('178 -> C18O 7->6')
#           print('179 -> C18O 8->7')
#           print('180 -> C18O 9->8')
#           print('181 -> C18O 10->9')
#           print('182 -> C18O 11->10')
#           print('183 -> C18O 12->11')
#           print('184 -> C18O 13->12')
#           print('185 -> C18O 14->13')
#           print('186 -> C18O 15->14')
#           print('187 -> C18O 16->15')
#           print('188 -> C18O 17->16')
#           print('189 -> C18O 18->17')
#           print('190 -> C18O 19->18')
#           print('191 -> C18O 20->19')
#           selection = input()

#       species = ['13CO 10', 'C+ 1', 'CO 1', 'CO 2', 'CO 3', 'CO 4', 'CO 5', 'CO 6', 'CO 7', 'CO 8', 'CO 9', 'CO 10', '13CO 1', '13CO 2', '13CO 3', '13CO 4', '13CO 5', '13CO 6', '13CO 7', '13CO 8', '13CO 9', '13CO 10', 'O 2']
#       kosma = Model(x, y, z, modelName='MilkyWay', modelType=shape, resolution=resolution)
#       kosma.addSpecies(species)

#     elif selection=='2':
#       selection = 0
#       kosma.calculateModel()

#     elif selection=='3':
#       selection = 0
#       kosma.writeEmission()

#     elif selection=='4':
#       selection = 0
#       kosma.calculateObservation()

#     elif selection=='5':
#       plot = 0
#       while plot==0:
#         print('1 -> FUV')
#         print('2 -> Clump mass')
#         print('3 -> Interclump mass')
#         print('4 -> Velocity')
#         print('5 -> Clump intensity')
#         print('6 -> Interclump intensity')
#         print('7 -> Total intensity')
#         print('8 -> Clump optical depth')
#         print('9 -> Interclump optical depth')
#         print('10 -> Total optical depth')
#         plot = input()
#         if plot==1: kosma.radiativeTransfer.plotModel(plot='FUV')
#         elif plot==2: kosma.radiativeTransfer.plotModel(plot='Afuv')
#         elif plot==3: kosma.radiativeTransfer.plotModel(plot='velocity')
#         elif plot==4: kosma.radiativeTransfer.plotModel(plot='clump intensity')
#         elif plot==5: kosma.radiativeTransfer.plotModel(plot='interclump intensity')
#         elif plot==6:
#           molecule = input('Which molecule?')
#           kosma.radiativeTransfer.plotModel(plot='species total intensity', species=molecule)
#         elif plot==7: kosma.radiativeTransfer.plotModel(plot='total intensity')
#         elif plot==8: kosma.radiativeTransfer.plotModel(plot='clump optical depth')
#         elif plot==9: kosma.radiativeTransfer.plotModel(plot='interclump optical depth')
#         elif plot==10: kosma.radiativeTransfer.plotModel(plot='total optical depth')
#         elif plot==11:
#           molecule = input('Which molecule?')
#           kosma.radiativeTransfer.plotModel(plot='species total optical depth', species=molecule)

#     elif selection=='6': break

#     selection = 0
