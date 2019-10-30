# quicklook for specs
#######################################################
#from matplotlib.figure import Figure
import matplotlib.cm as cm # colormap
#import pylab as p
from mpl_toolkits.axes_grid1 import make_axes_locatable
import numpy as np
import matplotlib.pyplot as plt

class gauss:
    def __init__(self, v0, sigma, peak):
        self.peak = float(peak) # area below the curve (integrated curve) 
        self.v0 = float(v0) # peak velocity
        self.sigma = float(sigma) # standard derivation
    def gaussfunc(self, v):
        if self.sigma == 0:
            return 0
        else:
            return self.peak * np.exp(-.5 * ( (v - self.v0) / self.sigma)**2 )


#files = ['spec[[ 0 -3]]_13CO_10-9.dat']
#files = ['spec[[0 0]]_C+_1-0.dat','spec[[ 0 -3]]_CO_2-1.dat','spec[[ 0 -3]]_CO_3-2.dat', 'spec[[ 0 -5]]_CO_6-5.dat', 'spec[[0 0]]_CO_10-9.dat', 'spec[[0 0]]_CO_16-15.dat','spec[[ 0 -6]]_13CO_3-2.dat','spec[[ 0 -5]]_13CO_5-4.dat','spec[[ 0 -5]]_13CO_6-5.dat', 'spec[[0 0]]_13CO_10-9.dat', 'spec[[0 0]]_HCO+_3-2.dat','spec[[0 0]]_HCO+_6-5.dat']

files = ['spec[[0 0]]_C+_1-0.dat','spec[[ 0 -4]]_CO_2-1.dat','spec[[ 0 -4]]_CO_3-2.dat', 'spec[[ 0 -4]]_CO_6-5.dat', 'spec[[ 0 -4]]_CO_10-9.dat', 'spec[[ 0 -4]]_CO_16-15.dat','spec[[ 0 -4]]_13CO_3-2.dat','spec[[ 0 -4]]_13CO_5-4.dat','spec[[ 0 -4]]_13CO_6-5.dat', 'spec[[ 0 -4]]_13CO_10-9.dat', 'spec[[ 0 -4]]_HCO+_3-2.dat','spec[[ 0 -4]]_HCO+_6-5.dat']

names = ['[CII], y-offset: 0','CO 2-1, y-offset: -4','CO 3-2 y-offset: -4', 'CO 6-5, y-offset: -4', 'CO 10-9, y-offset: -4', 'CO 16-15, y-offset: -4','$^{13}$CO_3-2, y-offset: -4','$^{13}$CO 5-4, y-offset: -4','$^{13}$CO 6-5, y-offset: -4', '$^{13}$CO 10-9, y-offset: -4', 'HCO$^+$ 3-2, y-offset: -4','HCO$^+$ 6-5, y-offset: -4']

#names = files


l = len(files)
spec = [None] * l

for sp in np.arange(l):
    spec[sp] = []
    with open(files[sp], 'r') as m:
        for line in m:
            spec[sp].append(line.split())
            # read in spectra
    spec[sp] = np.array(spec[sp], float)
    
print spec

# color = ('b', 'g', 'r', 'c', 'm', 'y', 'k')

f, axarr = plt.subplots(l, sharex = True)

f.set_size_inches([10,15])

for sp in np.arange(l):
    xval = np.linspace(11.3-5, 11.3+5, 100)
    g1 = gauss(11.3, 0.85, max(spec[sp][:,1])) # FWHM=2
    ga1 = g1.gaussfunc(xval)
    g2 = gauss(11.3, 1.7, max(spec[sp][:,1])) # FWHM=4
    ga2 = g2.gaussfunc(xval)
    axarr[sp].plot(spec[sp][:,0], spec[sp][:,1], color = 'k')
    axarr[sp].plot(xval, ga1, ls='--')
    axarr[sp].plot(xval, ga2, ls='--')
    axarr[sp].set_title(names[sp]) 
    # if sp == 0:
    axarr[int(l/2.-1)].set_ylabel('intensity [K]', labelpad = 15)   
    axarr[l-1].set_xlabel('velocity [km/s]')
    # axarr[l-1].xlabel.set_position((0.2, 0.1))
    axarr[sp].set_yticks(np.around(np.linspace(0, max(spec[sp][:,1]), 3), 1))
    # axarr[sp].set_yticks(np.array(np.linspace(0, max(spec[sp][:,1]), 3), int))
plt.subplots_adjust(hspace=.5)
#gcf().set_inches_size(4,8)

plt.savefig('test2.eps', format='eps' , figsize = (3,8) , dpi=600, bbox_inches='tight') 

plt.show()

