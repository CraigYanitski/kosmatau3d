import numpy as np
import matplotlib.pyplot as plt
from scipy.interpolate import griddata





"""
##################################### read in mass profile clumps
mass_list_r =[]
mass_list_m =[]
with open ('mass_profile.dat','r') as mass_file:
  next(mass_file) #skip header
  for line in mass_file:
    mass_list_r.append(line.split()[0])
    mass_list_m.append(line.split()[1])
"""
"""
##################################### plotting
plt.figure(1)
plt.plot(mass_list_r, mass_list_m, 'bo')
plt.ylabel(" H2 gas  Smass/pc^2")
plt.xlabel("r [kpc]")
plt.yscale('log')

#plt.show()
plt.savefig('mass_r.png', format='png', dpi=180, bbox_inches='tight')
plt.close()
"""
