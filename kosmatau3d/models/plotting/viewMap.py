import sys
import matplotlib
import matplotlib.pyplot as plt
import matplotlib.colors as colors
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
  
  # app = QApplication()
  
  window = viewMap(directory=directory)
  
  sys.exit(window.exec_())

class viewMap(QApplication):
  
  def __init__(self, directory='', posXY=(700, 40), windowSize=(925, 750)):
    if not directory:
      print('please specify the directory of the model history files.')
    super().__init__(sys.argv)
    self.window = QWidget()
    self.title = 'Subplot Switcher'
    self.posXY = posXY
    self.windowSize = windowSize
    self.directory = directory
    self.file = fits.open(directory+'/channel_intensity.fits')
    print(self.file[1].header)
    # lon = np.linspace(-np.pi, np.pi, self.file[1].shape[2])
    # lat = np.linspace(-np.pi/2, np.pi/2, self.file[1].shape[1])
    lon = np.linspace(self.file[1].header['CRVAL2'] - self.file[1].header['CRPIX2'] * self.file[1].header['CDELT2'],
                           self.file[1].header['CRVAL2'] + self.file[1].header['CRPIX2'] * self.file[1].header['CDELT2'],
                           num=self.file[1].header['NAXIS2'])
    lat = np.linspace(self.file[1].header['CRVAL3'] - self.file[1].header['CRPIX3'] * self.file[1].header['CDELT3'],
                           self.file[1].header['CRVAL3'] + self.file[1].header['CRPIX3'] * self.file[1].header['CDELT3'],
                           num=self.file[1].header['NAXIS3'])
    self.lon,self.lat = np.meshgrid(lon,lat)
    self.velocity = np.linspace(self.file[1].header['CRVAL4'] - self.file[1].header['CRPIX4'] * self.file[1].header['CDELT4'],
                           self.file[1].header['CRVAL4'] + self.file[1].header['CRPIX4'] * self.file[1].header['CDELT4'],
                           num=self.file[1].header['NAXIS4'])
    self.species_intensity_integrated = np.trapz(self.file[1].data, self.velocity, axis=0)
    self.dust_intensity_integrated = np.trapz(self.file[2].data, self.velocity, axis=0)
    self.species = self.file[1].header['SPECIES'].split(', ')
    self.species_length = len(self.species)
    self.dust = self.file[2].header['Dust'].split(', ')
    self.fig = plt.figure(figsize=(9, 6))
    self.initUI()
  
  def initUI(self):
    # Initialise the user interface, including widgets.
    
    QMainWindow().setCentralWidget(QWidget())
    
    self.window.setLayout(QVBoxLayout())
    self.window.layout().setContentsMargins(0, 0, 0, 0)
    self.window.layout().setSpacing(0)
    
    self.canvas = FigCanvas(self.fig)
    self.canvas.draw()
    
    scroll = QScrollArea(self.window)
    scroll.setWidget(self.canvas)
    
    self.VELslider = QSlider(Qt.Horizontal, self.window)
    self.VELslider.setRange(0, self.file[1].shape[0] - 1)
    self.VELslider.valueChanged[int].connect(self.changeVelocity)
    
    self.VEL = QLabel('Velocity: {:.2f}'.format(self.velocity[0]))
    
    self.SPEslider = QSlider(Qt.Horizontal, self.window)
    self.SPEslider.setRange(0, len(self.species) + len(self.dust) - 1)
    self.SPEslider.valueChanged[int].connect(self.changeSpecies)
    
    self.SPE = QLabel('Species: {}'.format(self.species[0]))
    
    button = QPushButton('Integrate')
    button.clicked.connect(self.plotIntegrated)
    
    nav = NabToolbar(self.canvas, self.window)
    
    self.window.layout().addWidget(nav)
    self.window.layout().addWidget(scroll)
    self.window.layout().addWidget(self.VELslider)
    self.window.layout().addWidget(self.VEL)
    self.window.layout().addWidget(self.SPEslider)
    self.window.layout().addWidget(self.SPE)
    self.window.layout().addWidget(button)
    self.updatePlot(self.file[1].data[0,:,:,0], '{}'.format(self.species[0]))
    
    self.show_basic()
    
    return
  
  def changeVelocity(self, value):
    # Update the subplot value.
    
    spe = self.SPEslider.value()
    self.VEL.setText('Velocity: {:.2f}'.format(self.velocity[value]))
    
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
    
    vmax = data.max()
    vmin = data.min()
    
    ax = self.fig.add_axes((0.1, 0.1, 0.8, 0.9), projection='mollweide')
    # cm = ax.imshow(data, extent=(-np.pi, np.pi, -np.pi/2, np.pi/2), norm=colors.SymLogNorm(linthresh=0.1, vmin=vmin, vmax=vmax), cmap='cubehelix')
    cm = ax.pcolormesh(self.lon, self.lat, data, norm=colors.SymLogNorm(linthresh=0.1, vmin=vmin, vmax=vmax), cmap='cubehelix')
    cb = self.fig.colorbar(cm, extend='both', ax=ax, fraction=0.02)
    if integrated:
      cb.ax.set_ylabel(r'Intensity $log_{10} \left( K \frac{km}{s} \right)$', fontsize=20, labelpad=16, rotation=270)
    else:
      cb.ax.set_ylabel(r'Intensity $log_{10} \left( K \right)$', fontsize=20, labelpad=16, rotation=270)
    ax.set_xlabel(r'Longitude $\left( rad \right)$', fontsize=16)
    ax.set_ylabel(r'Latitude $\left( rad \right)$', fontsize=16)
    ax.set_title('Galactic '+title, fontsize=24)
    self.fig.tight_layout()
    
    self.canvas.draw()
    
    return
  
  def show_basic(self):
    # Set window geometry and show on screen.
    
    self.window.setWindowTitle(self.title)
    self.window.setGeometry(*self.posXY, *self.windowSize)
    self.window.show()
    
    return


if __name__ == '__main__':
  # Run the example if the script is run on its own.
  
  main()