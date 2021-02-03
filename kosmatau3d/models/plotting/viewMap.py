import sys
import matplotlib
import matplotlib.pyplot as plt
from astropy.io import fits
import numpy as np
from PyQt5.QtCore import Qt
from PyQt5.QtWidgets import (
  QApplication,
  QMainWindow,
  QVBoxLayout,
  QWidget,
  QLabel,
  QPushButton,
  QScrollArea,
  QSlider,
)
from matplotlib.backends.backend_qt5agg import (
  FigureCanvasQTAgg as FigCanvas,
  NavigationToolbar2QT as NabToolbar,
)

from .. import constants

# Make sure that we are using QT5
matplotlib.use('Qt5Agg')

def main():
  directory = constants.HISTORYPATH
  
  app = QApplication(sys.argv)
  
  window = viewMap(directory=directory)
  sys.exit(app.exec_())

class viewMap(QWidget):
  fig = plt.figure(figsize=(9, 6))
  def __init__(self, directory, posXY=(700, 40), windowSize=(925, 750)):
    super().__init__()
    self.title = 'Subplot Switcher'
    self.posXY = posXY
    self.windowSize = windowSize
    self.directory = directory
    self.file = fits.open(directory+'/channel_intensity.fits')
    velocity = np.linspace(self.file[1].header['CRVAL4'] - self.file[1].header['CRPIX4'] * self.file[1].header['CDELT4'],
                           self.file[1].header['CRVAL4'] + self.file[1].header['CRPIX4'] * self.file[1].header['CDELT4'],
                           num=self.file[1].header['NAXIS4'])
    self.species_intensity_integrated = np.trapz(self.file[1].data, velocity, axis=0)
    self.dust_intensity_integrated = np.trapz(self.file[2].data, velocity, axis=0)
    self.species = self.file[1].header['SPECIES'].split(', ')
    self.species_length = len(self.species)
    self.dust = self.file[2].header['Dust'].split(', ')
    self.fig = fig
    self.initUI()
  
  def initUI(self):
    # Initialise the user interface, including widgets.
    
    QMainWindow().setCentralWidget(QWidget())
    
    self.setLayout(QVBoxLayout())
    self.layout().setContentsMargins(0, 0, 0, 0)
    self.layout().setSpacing(0)
    
    self.canvas = FigCanvas(self.fig)
    self.canvas.draw()
    
    scroll = QScrollArea(self)
    scroll.setWidget(self.canvas)
    
    self.VELslider = QSlider(Qt.Horizontal, self)
    self.VELslider.setRange(0, self.file[1].shape[0] - 1)
    self.VELslider.valueChanged[int].connect(self.changeVelocity)
    
    self.VEL = QLabel('Velocity: {}'.format(0))
    
    self.SPEslider = QSlider(Qt.Horizontal, self)
    self.SPEslider.setRange(0, len(self.species) + len(self.dust) - 1)
    self.SPEslider.valueChanged[int].connect(self.changeSpecies)
    
    self.SPE = QLabel('Species: {}'.format(self.species[0]))
    
    button = QPushButton('Integrate')
    button.clicked.connect(self.plotIntegrated)
    
    nav = NabToolbar(self.canvas, self)
    
    self.layout().addWidget(nav)
    self.layout().addWidget(scroll)
    self.layout().addWidget(self.VELslider)
    self.layout().addWidget(self.VEL)
    self.layout().addWidget(self.SPEslider)
    self.layout().addWidget(self.SPE)
    self.layout().addWidget(button)
    self.updatePlot(self.file[1].data[0,:,:,0], '{}'.format(self.species[0]))
    
    self.show_basic()
    
    return
  
  def changeVelocity(self, value):
    # Update the subplot value.
    
    spe = self.SPEslider.value()
    self.VEL.setText('Velocity: {}'.format(value))
    
    if spe<self.species_length:
      self.updatePlot(self.file[1].data[value,:,:,spe], '{}'.format(self.species[spe]))
    else:
      dust = spe-self.species_length
      self.updatePlot(self.file[2].data[value,:,:,dust], '{}'.format(self.dust[dust]))
    
    return
  
  def changeSpecies(self, value):
    # Update the subplot value.
    
    vel = self.VELslider.value()
    
    if value<self.species_length:
      self.SPE.setText('Species: {}'.format(self.species[value]))
      self.updatePlot(self.file[1].data[vel,:,:,value], '{}'.format(self.species[value]))
    else:
      dust = value-self.species_length
      self.SPE.setText('Dust: {}'.format(self.dust[dust]))
      self.updatePlot(self.file[2].data[vel,:,:,dust], '{}'.format(self.dust[dust]))
    
    return
  
  def plotIntegrated(self):
    # Plot the corresponding integrated intensity.
    
    spe = self.SPEslider.value()
    
    if spe<self.species_length:
      self.updatePlot(self.species_intensity_integrated[:,:,spe], 'integrated {}'.format(self.species[spe]), integrated=True)
    else:
      dust = spe-self.species_length
      self.updatePlot(self.dust_intensity_integrated[:,:,dust], 'integrated {}'.format(self.dust[dust]), integrated=True)
  
  def updatePlot(self, data, title, integrated=False):
    # Remove the current Axes object and create a new one of the desired subplot.
    
    axTemp = self.fig.axes
    for ax in axTemp:
      self.fig.delaxes(ax)
    
    value = np.full(data.shape, 1e-100)
    value[data>0] = data[data>0]
    value = np.log10(value)
    # value = np.nan_to_num(value, nan=-100)
    vmax = value.max()
    
    ax = self.fig.add_axes((0.1, 0.1, 0.8, 0.9))
    cm = ax.imshow(value, extent=(-np.pi, np.pi, -np.pi/2, np.pi/2), vmin=-2, vmax=max(0,vmax), cmap='cubehelix')
    cb = fig.colorbar(cm, extend='both', ax=ax, fraction=0.02)
    if integrated:
      cb.ax.set_ylabel('Intensity log(K km/s)', fontsize=20, labelpad=16, rotation=270)
    else:
      cb.ax.set_ylabel('Intensity log(K)', fontsize=20, labelpad=16, rotation=270)
    ax.set_xlabel('Longitude (rad)', fontsize=16)
    ax.set_ylabel('Latitude (rad)', fontsize=16)
    ax.set_title('Galactic '+title, fontsize=24)
    self.fig.tight_layout()
    
    self.canvas.draw()
    
    return
  
  def show_basic(self):
    # Set window geometry and show on screen.
    
    self.setWindowTitle(self.title)
    self.setGeometry(*self.posXY, *self.windowSize)
    self.show()
    
    return


if __name__ == '__main__':
  # Run the example if the script is run on its own.
  
  main()