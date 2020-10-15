import matplotlib.cm as cm # colormap
from mpl_toolkits.axes_grid1 import make_axes_locatable
import numpy as np
import matplotlib.pyplot as plt
import matplotlib as mpl

from matplotlib.backends.backend_agg import FigureCanvasAgg as FigureCanvas
from matplotlib.figure import Figure
from scipy.interpolate import griddata
import os, sys

# Open a file
path = "data"          # verzeichnis mit den dateien
files = os.listdir( path ) # erstellt liste mit dateinamen

lambdas = [912, 1350, 1500, 1650, 2000, 2200, 2500, 2800, 3650] #lambda [A]
files.sort() #files are now in same order than the lambdas

list = [None]*(len(lambdas)+0) #with files and radiation input#
for i, l in enumerate(list):
  list[i]=[] #leere Liste vorbereiten

for i, filename in enumerate(files):
  with open ('data/' + filename) as file:
    #print type(file)
    for j, line in enumerate(file):
      if j > 81: #81: # start at datapoints    
        temp = line.split() #r[pc], z[pc], urad[erg/pc3/A], eabs[erg/s/pc3/A]
        if float(temp[1]) == 0:
          arr = [temp[0],temp[2],temp[3]]
          list[i].append(arr)
          #list[i].append(str(temp[0] + ' ' + temp[2] + ' ' + temp[3]))
          #print temp[0], temp[2]

#print list[0]

################ prepare plots writing
for i, files in enumerate(list):
    r = []
    radiation = []
    for j, leng in enumerate(list[i]):
        #print j, leng[0], leng[1]
        r.append(float(leng[0]))
        radiation.append(float(leng[1]))
##################################### einzel plots    
    plt.figure(0)
    plt.subplot(111)     #number of plots in x, in y, plot ID number
    plt.ylabel('radiation-energy-dens [erg/pc3/A]')
    plt.xlabel("radius [PC]")
    text = str(round(lambdas[i],1))   #Angst
    plt.title(text) # + file[:-5])
    plt.plot(r,radiation, 'bo')
    plt.savefig('./writen/'+str(lambdas[i])+'A_plot'+'.jpg', format='jpg', dpi=180, bbox_inches='tight')
    plt.close()


#################################### multiplot
################ prepare plots writing
colors = iter(cm.rainbow(np.linspace(0, 1, len(lambdas))))
for i, files in enumerate(list):
    r = []
    radiation = []
    for j, leng in enumerate(list[i]):
        #print j, leng[0], leng[1]
        r.append(float(leng[0]))
        radiation.append(float(leng[1]))
##################################### einzel plots    
    plt.figure(0)
    plt.subplot(111)     #number of plots in x, in y, plot ID number
    plt.ylabel('radiation-engery-dens [erg/pc3/A]')
    plt.xlabel("radius [PC]")
    text = str(round(lambdas[i],1))   #Angst
    plt.title(text)

    legend = [None]* (len(lambdas))
    #legend = np.zeros(len(lamdas))
    for i, id in enumerate(lambdas):
        legend[i] = (str(lambdas[i]) + 'A')
    #print legend
    plt.legend(legend, loc='upper right')


    plt.plot(r,radiation, 'bo', color=next(colors))
    plt.savefig('./writen/'+str('total')+'plot_galactic-UV'+'.jpg', format='jpg', dpi=180, bbox_inches='tight')
plt.close()


meta_r = []
meta_radiation = []

################ prepare plots writing
for i, files in enumerate(list):
    r = []
    radiation = []
    for j, leng in enumerate(list[i]):
        #print j, leng[0], leng[1]
        r.append(float(leng[0]))
        radiation.append(float(leng[1]))
    meta_r.append(r)
    meta_radiation.append(radiation)

#print meta_r
with open('galactic_FUV', 'w') as file:
    file.write('#Header: r[pc] then galactic radiation-field energy-density distributions [erg/pc^3/A] avg and for different wavelenghts [A]')
    file.write( str(lambdas))
    
    for i, df in enumerate(r):
        
        avg = 0
        avg = meta_radiation[0][i] + meta_radiation[1][i] + meta_radiation[2][i] + meta_radiation[3][i] + \
              meta_radiation[4][i] + meta_radiation[5][i] + meta_radiation[6][i] + meta_radiation[7][i] + \
              meta_radiation[8][i]
        avg = avg /9. 
        #print 'avg', (avg)

        file.write('\n' +  str(r[i]) +'\t' +str(avg) +'\t'  + str(meta_radiation[0][i]) +'\t' +str(meta_radiation[1][i]) \
        +'\t' +str(meta_radiation[2][i]) +'\t' +str(meta_radiation[3][i]) +'\t' +str(meta_radiation[4][i])\
        +'\t' +str(meta_radiation[5][i]) +'\t' +str(meta_radiation[6][i]) +'\t' +str(meta_radiation[7][i])\
        +'\t' +str(meta_radiation[8][i])   )




















#################################### multiplot 2
################ prepare plots writing
colors = iter(cm.rainbow(np.linspace(0, 1, len(lambdas))))
for i, files in enumerate(list):
    r = []
    radiation = []
    for j, leng in enumerate(list[i]):
        #print j, leng[0], leng[1]
        r.append(float(leng[0]))
        radiation.append(float(leng[1]))
##################################### einzel plots    
    plt.figure(0)
    plt.subplot(111)     #number of plots in x, in y, plot ID number
    plt.ylabel('radiation-engery-dens [erg/pc3/A]')
    plt.xlabel("radius [PC]")
    text = str(round(lambdas[i],1))   #Angst
    plt.title(text)

    legend = [None]* (len(lambdas))
    #legend = np.zeros(len(lamdas))
    for i, id in enumerate(lambdas):
        legend[i] = (str(lambdas[i]) + 'A')
    #print legend
    plt.legend(legend, loc='upper right')
    plt.plot(r,radiation, 'bo', color=next(colors))
    #plt.plot(r,radiation, 'v')


habi = []
drai = []
#################################### create theory curves
parsec = 3.08567758149137*10**18 #parsec in cm
#lambdas = [912, 1350, 1500, 1650, 2000, 2200, 2500, 2800, 3650] #
# in ergs/cm^3
colors = iter(cm.rainbow(np.linspace(0, 1, len(lambdas))))
for i, lamb in enumerate(lambdas):
  l3 = lamb/1000. # in angststrom / 1000
  habing = 10**-14*(-25./6*l3**3 + 25./2*l3**2- 13./3*l3)
  chi = 1.71
  draine = 4e-14*chi*l3**-5 * (31.016*l3**2 - 49.913*l3 + 19.897)

  draine = draine * parsec**3 / lamb # erg/cm^3 -> erg/pc^3/Angst
  habing = habing * parsec**3 / lamb #
  #print 'draine', draine
  #print 'habing', habing
  habi.append(habing)
  drai.append(draine)
  if i < 10:
    plt.plot([0, 2.5e4], [habing, habing], 'k-', color = next(colors))
    #plt.plot([0, 2.5e4], [draine, draine], 'k-', color = next(colors))
#hab = []
#for x in r:
#  hab.append(


plt.savefig('./writen/'+str('total')+'plot_galactic-UV'+'.jpg', format='jpg', dpi=180, bbox_inches='tight')
plt.close()
