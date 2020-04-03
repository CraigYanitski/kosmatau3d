import matplotlib.pyplot as plt

import constants
import observations
import interpolations

from Voxel import *

constants.directory = 'MilkyWay/'
constants.resolution = 100
observations.methods.initialise()
interpolations.initialise()

FUV0 = 10**0
v = Voxel(0)
v.setPosition(x=1000, y=0, z=0, r=1000, phi=0)
v.setProperties(values='automatic', clumpMass=250, interclumpMass=250, mass=500, velocity=0, density=15000, FUV=FUV0)
v.calculateEmission()

FUV1 = 10**2
w = Voxel(1)
w.setPosition(x=1000, y=0, z=0, r=1000, phi=0)
w.setProperties(values='manual', clumpMass=250, interclumpMass=250, mass=500, velocity=0, density=15000, FUV=FUV1)
w.calculateEmission()

FUV2 = 10**4
x = Voxel(2)
x.setPosition(x=1000, y=0, z=0, r=1000, phi=0)
x.setProperties(values='manual', clumpMass=250, interclumpMass=250, mass=500, velocity=0, density=15000, FUV=FUV2)
x.calculateEmission()

FUV3 = 10**6
y = Voxel(3)
y.setPosition(x=1000, y=0, z=0, r=1000, phi=0)
y.setProperties(values='manual', clumpMass=250, interclumpMass=250, mass=500, velocity=0, density=15000, FUV=FUV3)
y.calculateEmission()

(cI0,iI0),depth = v.getEmission()
(cI1,iI1),depth = w.getEmission()
(cI2,iI2),depth = x.getEmission()
(cI3,iI3),depth = y.getEmission()

plt.loglog(constants.wavelengths, cI0[0,:], marker='.', ms=2, ls='', lw=0.5, label=r'$r_{gal}=1$ kpc')
plt.loglog(constants.wavelengths, cI1[0,:], marker='', ms=0.5, ls='-', lw=0.5, label=r'$\chi$={:.0e}'.format(FUV1))
plt.loglog(constants.wavelengths, cI2[0,:], marker='', ms=0.5, ls='-', lw=0.5, label=r'$\chi$={:.0e}'.format(FUV2))
plt.loglog(constants.wavelengths, cI3[0,:], marker='', ms=0.5, ls='-', lw=0.5, label=r'$\chi$={:.0e}'.format(FUV3))

plt.legend()

plt.show()