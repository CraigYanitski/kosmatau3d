from tkinter import *
import PDRfunctions as pf
import pylab as p
import numpy as np
from copy import copy
from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg, \
                                              NavigationToolbar2Tk
from maps import maps
from plots import plots
from matplotlib.figure import Figure
import matplotlib.cm as cm # colormap
# edit: Craig
import globals as gbl
#end edit

#import matplotlib.pyplot as plt
class PDR_GUI(): #++maps
  def __init__(self, root): 
    #from PDR import _globals # import global variables
    self.root = root
    self.root.title('3D PDR Simulation: ' + gbl._globals['compound']['name']) 
    # name of project
    self.done = 0

    fig = pf.plot_geo(gbl._globals['compound']['coordinates'], \
          wei = gbl._globals['compound']['weights'][0,:],\
          labelx = gbl._globals['plots']['xlabel'], \
          labely = gbl._globals['plots']['ylabel'], \
          labelz = gbl._globals['plots']['zlabel'],\
          limits = gbl._globals['plots']['plotrange'], \
          title = 'Geometry',\
          los_start = gbl._globals['compound']['los_start'], \
          los_end = gbl._globals['compound']['los_end'],\
          losoffset = gbl._globals['compound']['losoffsets'])
    geoCanvas = FigureCanvasTkAgg(fig, master = self.root)
    geoCanvas.get_tk_widget().grid(row=0, column=0, columnspan=9, \
                                    rowspan=9, sticky=W+E+N+S)  
    mapLabel = Label(self.root, text = "MAPS",bg="slate gray" )
    mapLabel.grid(row=0, column=9, columnspan=4, sticky = N+S+E+W)
    specLabel1 = Label(self.root, text = "SPECTRA",bg="slate gray" )
    specLabel1.grid(row=4, column=9,columnspan=4, sticky = N+S+E+W)
    specLabel1 = Label(self.root, text = "1")
    specLabel1.grid(row=6, column=9)
    specLabel2 = Label(self.root, text = "2")
    specLabel2.grid(row=7, column=9)
#     self.progress = Entry(self.root)
#     self.progress.grid(row = 8, column = 9, columnspan = 3, sticky = N+S+E+W)
#     self.progress.insert(END, "3D PDR simulation ready...")

    vlabel = Label(self.root, text = "velocity [km/s]:")
    vlabel.grid(row=1, column=10, sticky=W)
    dvlabel = Label(self.root, text = "delta velocity [km/s]:")
    dvlabel.grid(row=2, column=10, sticky=W)
    dbeamlabel = Label(self.root, text = "beamsize [arcsec]:")
    dbeamlabel.grid(row=3, column=10, sticky=W)
    offxLabel = Label(self.root, text = "x Offsets")
    offxLabel.grid(row=5, column=10)
    self.x1 = Entry(self.root)
    self.x1.grid(row=6, column=10)
    self.x1.insert(END, 0)
    self.x2 = Entry(self.root)
    self.x2.grid(row=7, column=10)
    self.x2.insert(END, 0)

    self.vEntry = Entry(self.root)
    self.vEntry.grid(row=1, column=11) 
    self.vEntry.insert(END, gbl._globals['compound']['vel'])
    self.dvEntry = Entry(self.root)
    self.dvEntry.grid(row=2, column=11)
    self.dvEntry.insert(END, gbl._globals['compound']['d_vel'])
    self.beamEntry = Entry(self.root)
    self.beamEntry.grid(row=3, column=11)
    self.beamEntry.insert(END, gbl._globals['compound']['dbeam'][0])
    offyLabel = Label(self.root, text = "y Offsets")
    offyLabel.grid(row=5, column=11)
    self.y1 = Entry(self.root)
    self.y1.grid(row=6, column=11)
    self.y1.insert(END, 0)
    self.y2 = Entry(self.root)
    self.y2.grid(row=7, column=11)
    self.y2.insert(END, 0)

    mapButton = Button(self.root, text="refresh map", \
                        command = self.refresh_map_int)
    mapButton.grid(row=1, column=12, rowspan=2, sticky= N+S+E+W)
    convButton = Button(self.root, text="convolve", \
                        command = self.refresh_map_conv)
    convButton.grid(row=3, column=12, sticky= N+S+E+W)
    spec1Button = Button(self.root, text="refresh spec1", \
                          command = self.refresh_spec1)
    spec1Button.grid(row=6, column=12, sticky= N+S+E+W)
    spec2Button = Button(self.root, text="refresh spec2", \
                          command = self.refresh_spec2)
    spec2Button.grid(row=7, column=12, sticky= N+S+E+W)
    quitButton = Button(self.root, text="QUIT", bg="red",\
                        command = self.quit)
    quitButton.grid(row=8, column=12, sticky = N+S+E+W)

############# functions ##################### 

  def quit(self):
    self.root.quit()
    self.root.destroy()

########## refresh non-convolved map #################

    def refresh_map_int(self): # channel map not convolved with beam
 #       self.progress.delete(0, END)
 #       self.progress.insert(END, "creating maps...")

        try: vel = float(self.vEntry.get())
        # line peak velocity
        except: 
            vel = 0
            self.vEntry.delete(0, END) 
            self.vEntry.insert(END, '0')
            print('no peak velocity given for integration. setting vel = 0')
        try: d_vel = float(self.dvEntry.get())
        # velocity integration range
        except: 
            d_vel = 3
            self.dvEntry.delete(0, END)                  
            self.dvEntry.insert(END, '3')
            print('no velocity range given for integration. setting d_vel = 3')
        print('Variable vel: ', vel)
        print('Variable dvel: ', d_vel)
        m = maps()
        m.create_map(vel, d_vel)
        self.done = 1
        # import complete map from maps

        top = Toplevel()
        top.title("Non-convolved map")
        #from PDR import gbl._globals
        map_int = gbl._globals['map_int']
        # plot map   
        mapFrame = Frame(top)
        mapFrame.pack()

        fig = Figure(dpi = 100) 
        ax = fig.add_subplot(111)
        cax = ax.imshow(map_int, cmap=cm.hot, interpolation='nearest',\
                        extent = [gbl._globals['map_limits']['minx'], \
                        gbl._globals['map_limits']['maxx'], \
                        gbl._globals['map_limits']['miny'], \
                        gbl._globals['map_limits']['maxy']])
        maxintense = max(np.array(map_int).flat)
        print('Variable maxintense: ', maxintense)
        ax.set_title('non-convolved...')
        # Add colorbar, make sure to specify tick locations 
        # to match desired ticklabels
        cbar = fig.colorbar(cax)
        # cbar.ax.set_yticklabels(['0', str(maxintense)])# vertically oriented colorbar
     
        mapCanvas = FigureCanvasTkAgg(fig, master = mapFrame)
        mapCanvas.show()
        mapCanvas.get_tk_widget().pack(side = "top", fill = "both", expand = 1)
        # side="top", fill="both", expand=1

        toolbar = NavigationToolbar2TkAgg(mapCanvas, top)
        toolbar.update()
        mapCanvas._tkcanvas.pack(side="top", fill="both", expand=1)   

#######################convolve map############################################################

  def refresh_map_conv(self): 
    # channel map convolved with beam
    if self.done == 0: self.refresh_map_int() 
    # create non-convolved map if not yet existent

#      self.progress.delete(0, END)
#      self.progress.insert(END, "convolving maps...")

    try: dbeam = float(self.beamEntry.get())
    except: 
      dbeam = 20
      self.beamEntry.delete(0, END) 
      self.beamEntry.insert(END, '20')
    print('dbeam: ', dbeam)

#       self.progress.delete(0, END)
#       self.progress.insert(END, "3D PDR simulation ready...")

    # import convolved map
    m = maps()
    m.create_convolved_map(dbeam)
    #from PDR import gbl._globals
    # print gbl._globals['map'] 
    map_conv = gbl._globals['map_conv']   
    # plot convolved map 

    fig = Figure(figsize = (6,6), dpi=100) 
    ax = fig.add_subplot(111)
    ax.set_xlabel(gbl._globals['plots']['xlabel']) # labels cause error? 
    ax.set_ylabel(gbl._globals['plots']['ylabel'])
    cax = ax.imshow(map_conv, cmap = cm.hot, \
                    extent = [gbl._globals['map_limits']['minx'] + gbl._globals['threeSigma'] , \
                    gbl._globals['map_limits']['maxx'] - gbl._globals['threeSigma'], \
                    gbl._globals['map_limits']['miny'] + gbl._globals['threeSigma'], \
                    gbl._globals['map_limits']['maxy'] - gbl._globals['threeSigma']])
    ax.set_title('... and convolved map')
    # Add colorbar, make sure to specify tick locations to match desired ticklabels
    cbar = fig.colorbar(cax)
  #  cbar.ax.set_yticklabels(['0', '50', '100 erg...'])# vertically oriented colorbar

    convCanvas = FigureCanvasTkAgg(fig, master = self.root)
    convCanvas.get_tk_widget().grid(row=9, column=11, columnspan=11, rowspan=12, sticky=N+S+E+W)

    self.convAdvButton = Button(self.root, bg="slate gray", text="Advanced Map", command = self.map_convAdv)
    self.convAdvButton.grid(row = 8, column = 11, sticky = N+S+E+W)

 
###########################################################
#####################spectra################

  def refresh_spec1(self):
    offI0 = True # switch
    try: 
      offx1 = int(float(self.x1.get()))
      self.x1.delete(0, END)
      self.x1.insert(END, offx1)   
    except: 
      iffI0 = False
      self.x1.delete(0, END)

    try: 
      offy1 = int(float(self.y1.get()))
      self.y1.delete(0, END)
      self.y1.insert(END, offy1)   
    except: 
      iffI0 = False
      self.y1.delete(0, END)

    print('Variable offI0: ', offI0)

    if offI0:
      p = plots()
      p.plot_spec_simple(offsets = [offx1, offy1])
      #from PDR import gbl._globals
      xval = gbl._globals['xval']
      yval = gbl._globals['yval']
      f_off = Figure(figsize=(8,3)) 
      a_off = f_off.add_subplot(111)
      a_off.plot(xval, yval)

      offCanvas1 = FigureCanvasTkAgg(f_off, master = self.root)
      offCanvas1.get_tk_widget().grid(row = 9, column = 0, columnspan=11, sticky = N+W)

      self.specAdvButton = Button(self.root ,bg="slate gray", text="Advanced Specs", command = self.spec_Adv)
      self.specAdvButton.grid(row = 8, column = 10, sticky = N+S+E+W)

  def refresh_spec2(self):

    offI0 = True # switch

    try: 
      offx2 = int(float(self.x2.get()))
      self.x2.delete(0, END)
      self.x2.insert(END, offx2)   
    except: 
      iffI0 = False
      self.x2.delete(0, END)

    try: 
      offy2 = int(float(self.y2.get()))
      self.y2.delete(0, END)
      self.y2.insert(END, offy2)   
    except: 
      iffI0 = False
      self.y2.delete(0, END) 

    print('Variable offI0: ', offI0)

    if offI0:
      p2 = plots()
      p2.plot_spec_simple(offsets = [offx2, offy2])
      #from PDR import gbl._globals
      xval = gbl._globals['xval']
      yval = gbl._globals['yval']
      f_off2 = Figure(figsize=(8,3)) 
      a_off2 = f_off2.add_subplot(111)
      a_off2.plot(xval, yval)
      offCanvas2 = FigureCanvasTkAgg(f_off2, master = self.root)
      offCanvas2.get_tk_widget().grid(row = 12, column = 0, columnspan=11)
    #    toolbar = NavigationToolbar2TkAgg(offCanvas2, self.root)
    #    toolbar.update()
      #   offCanvas2._tkcanvas.grid()

####################################Toplevel#################################

  def map_convAdv(self): 
    # open convolved map in a new window which enables saving it
      top = Toplevel()
      top.title("Convolved Map")
      #from PDR import gbl._globals
      map_conv = np.array(gbl._globals['map_conv']) # convolved map  

      # print map_conv
      # pasue = input()
      # l = np.array(map_conv).size
      map_conv_erg = copy(map_conv)

      if gbl._globals['plots']['unit'] == 1:  
      # case conversion to Kkm/s
        for i in np.arange(map_conv.shape[0]):
          for j in np.arange(map_conv.shape[1]):
            map_conv_erg[int(i),int(j)] = map_conv[int(i),int(j)] * 2 * \
            gbl._globals['constants']['kb'] * \
            (gbl._globals['compound']['frequency']*10**9)**3 * 10**6 \
            /(gbl._globals['constants']['c']**3)           
      map_conv = map_conv_erg
      mapFrame = Frame(top)
      mapFrame.pack()

      fig = Figure(dpi=100) 
      ax = fig.add_subplot(111)
      ax.set_xlabel(gbl._globals['plots']['xlabel']) # labels cause error? 
      ax.set_ylabel(gbl._globals['plots']['ylabel'])

      # from scipy import ndimage
      # map_conv = ndimage.rotate(map_conv, -37, reshape = False)
      # rotate image

      # cm.hot (yellow-red)
      # or: cm.afmhot (brown...)
      cax = ax.imshow(map_conv, cmap=cm.hot, \
            extent = [gbl._globals['map_limits']['minx']+gbl._globals['threeSigma'],\
                    gbl._globals['map_limits']['maxx'] - gbl._globals['threeSigma'],\
                    gbl._globals['map_limits']['miny'] + gbl._globals['threeSigma'],\
                    gbl._globals['map_limits']['maxy'] - gbl._globals['threeSigma']])
      # Add colorbar, make sure to specify tick locations to match desired ticklabels
      cbar = fig.colorbar(cax)
  #     ax.autoscale_view()

    # create label
      if gbl._globals['plots']['unit'] == 1:  # case conversion to Kkm/s 
        cbar.set_label('Integrated Intensity [K km/s]')
      else: cbar.set_label('Integrated Intensity [erg/s/cm^2/sr]')

      mapCanvas = FigureCanvasTkAgg(fig, master = mapFrame)
      mapCanvas.show()
      mapCanvas.get_tk_widget().pack(side = "top", fill = "both", expand = 1)#side="top", fill="both", expand=1

      toolbar = NavigationToolbar2Tk(mapCanvas, top)
      toolbar.update()
      mapCanvas._tkcanvas.pack(side="top", fill="both", expand=1)

###############################################################################################

  def spec_Adv(self): # Advanced spectra plotting
    top = Toplevel()
    top.title("Spectra advanced")
    #from PDR import gbl._globals
    Entry1 = Entry(top)
    Entry1.grid(row=1, column=1) 
    Entry2 = Entry(top)
    Entry2.grid(row=2, column=1) 
    Entry3 = Entry(top)
    Entry3.grid(row=3, column=1) 
    Entry4 = Entry(top)
    Entry4.grid(row=4, column=1) 
    Entry5 = Entry(top)
    Entry5.grid(row=5, column=1) 


###########################################

if __name__ == '__main__':
#def main():
  root = Tk()
  g = PDR_GUI(root)
  root.mainloop()
