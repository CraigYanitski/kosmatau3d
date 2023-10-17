import sys
from functools import partial
from copy import copy
from pprint import pprint

import numpy as np
from astropy.io import fits
import matplotlib
import matplotlib.pyplot as plt
import matplotlib.colors as colors
# import QtQuick.Controls.Universal 2.2
# from PyQt5.QtCore import Qt
from PyQt5.QtWidgets import (
  QApplication,
  QMainWindow,
  QHBoxLayout,
  QVBoxLayout,
  QGridLayout,
  QWidget,
  QLabel,
  QPushButton,
  QRadioButton,
  QCheckBox,
  QScrollArea,
  QSlider,
  QDoubleSpinBox,
)
from matplotlib.backends.backend_qt5agg import (
  FigureCanvasQTAgg as FigCanvas,
  NavigationToolbar2QT as NabToolbar,
)

from kosmatau3d.models import constants


# Make sure that we are using QT5
matplotlib.use('Qt5Agg')


def main():
    directory = constants.HISTORYPATH
    
    # app = QApplication()
    
    window = ViewMap(directory=directory)
    
    sys.exit(window.exec_())


class ViewMap(QApplication):
  
    def __init__(self, directory='', file='/synthetic_intensity.fits', posXY=(700, 40), windowSize=(925, 850),
                 cmap='cubehelix', cmappv='jet'):
        if not directory:
            print('please specify the directory of the model history files.')
        super().__init__(sys.argv)
        self.window = QWidget()
        self.title = 'Channel Map Viewer'
        self.posXY = posXY
        self.windowSize = windowSize
        self.directory = directory
        self.file = fits.open(directory+file)
        # pprint(self.file[1].header)
        self.longitude = np.linspace(
          self.file[1].header['CRVAL2'] - (self.file[1].header['CRPIX2']-0.5) * self.file[1].header['CDELT2'],
          self.file[1].header['CRVAL2'] + (self.file[1].header['CRPIX2']-0.5) * self.file[1].header['CDELT2'],
          num=self.file[1].header['NAXIS2'])
        self.latitude = np.linspace(
          self.file[1].header['CRVAL3'] - (self.file[1].header['CRPIX3']-0.5) * self.file[1].header['CDELT3'],
          self.file[1].header['CRVAL3'] + (self.file[1].header['CRPIX3']-0.5) * self.file[1].header['CDELT3'],
          num=self.file[1].header['NAXIS3'])
        self.velocity = np.linspace(
          self.file[1].header['CRVAL4'] - self.file[1].header['CRPIX4'] * self.file[1].header['CDELT4'],
          self.file[1].header['CRVAL4'] + self.file[1].header['CRPIX4'] * self.file[1].header['CDELT4'],
          num=self.file[1].header['NAXIS4'])
        self.lonGrid, self.latGrid = np.meshgrid(self.longitude, self.latitude)
        self.velPV, self.lonPV = np.meshgrid(self.velocity, self.longitude)
        print(self.file[1].shape, self.file[2].shape)
        self.species_intensity_integrated = np.trapz([self.file[1].data[:, :, :, i] 
                                                      - self.file[1].data[:, :, :, i].min(0)
                                                      for i in range(self.file[1].shape[3])],
                                                     self.velocity, axis=1)
        self.dust_intensity_integrated = self.file[2].data
        self.species = self.file[1].header['SPECIES'].split(', ')
        self.species_length = len(self.species)
        self.dust = self.file[2].header['Dust'].split(', ')
        self.fig = plt.figure(figsize=(9, 6))
        self.cmap = cmap
        self.cmappv = cmappv
        self.initUI()
  
    def initUI(self):
        # Initialise the user interface, including widgets.
        
        QMainWindow().setCentralWidget(QWidget())
        # print(self.window.styleSheet())
        
        # Setup a scrollable canvas for the plot
        self.window.setLayout(QVBoxLayout())
        self.window.layout().setContentsMargins(0, 0, 0, 0)
        self.window.layout().setSpacing(0)
        
        self.canvas = FigCanvas(self.fig)
        self.canvas.draw()
        
        scroll = QScrollArea(self.window)
        scroll.setWidget(self.canvas)
        
        # Min/max boxes for the intensity colour map
        self.minIlabel = QLabel('minimum I:')
        self.minI = QDoubleSpinBox()
        self.minI.setMaximum(100.0001)
        self.minI.setMinimum(-9.9999)
        self.minI.setSingleStep(0.0001)
        self.minI.setDecimals(4)
        self.minI.setValue(0.01)
        # self.minI.valueChanged.connect(partial(self.plotDifference, minmax=True))
        self.maxIlabel = QLabel('maximum I:')
        self.maxI = QDoubleSpinBox()
        self.maxI.setMaximum(9999.9999)
        self.maxI.setMinimum(0.0001)
        self.maxI.setSingleStep(0.0001)
        self.maxI.setDecimals(4)
        self.maxI.setValue(10)
        # self.maxI.valueChanged.connect(partial(self.plotDifference, minmax=True))
        
        # Button to rescale plot
        self.rescaleButton = QCheckBox('rescale')
        self.rescaleButton.clicked.connect(partial(self.plotDifference, minmax=True))
        
        # Velocity slider (only for channel maps)
        self.VELslider = QSlider(Qt.Horizontal, self.window)
        self.VELslider.setRange(0, self.velocity.size - 1)
        self.VELslider.valueChanged[int].connect(partial(self.plotDifference, velocity=True))
        self.VEL = QLabel('Velocity: {:.2f} km/s'.format(self.velocity[0]))
        
        # Species/dust wavelength slider
        self.SPEslider = QSlider(Qt.Horizontal, self.window)
        self.SPEslider.setRange(0, len(self.species) + len(self.dust) - 1)
        self.SPEslider.valueChanged[int].connect(partial(self.plotDifference, species=True))
        self.SPE = QLabel('Species: {}'.format(self.species[0]))
        
        # Latitude slider (only for position-velocity plot)
        self.LATslider = QSlider(Qt.Horizontal, self.window)
        self.LATslider.setRange(0, self.latitude.size - 1)
        self.LATslider.setValue(int(self.latitude.size/2))
        self.LATslider.valueChanged[int].connect(partial(self.plotDifference, latitude=True))
        self.LATslider.setEnabled(False)
        self.LAT = QLabel('Latitude: {:.1f} degrees'.format(self.latitude[self.LATslider.value()]*180/np.pi))
        
        # Button to plot channel intensity (enabled by default)
        self.channelButton = QRadioButton('Channel intensity')
        self.channelButton.toggle()
        self.channelButton.clicked.connect(partial(self.plotDifference))
        
        # Button to plot integrated intensity
        self.integrateButton = QRadioButton('Integrated intensity')
        self.integrateButton.clicked.connect(partial(self.plotDifference))
        
        # Button to plot position-velocity diagram
        self.PVButton = QRadioButton('Position-velocity Diagram')
        self.PVButton.clicked.connect(partial(self.plotDifference))
        
        # Button to plot position-velocity diagram
        self.tempButton = QPushButton('Print window size')
        self.tempButton.clicked.connect(self.printSize)
        
        nav = NabToolbar(self.canvas, self.window)
        self.window.layout().addWidget(nav)
        self.window.layout().addWidget(scroll)
        
        grid = QGridLayout()
        
        layout = QGridLayout()
        layout.addWidget(self.minIlabel, 0, 0)
        layout.addWidget(self.minI, 0, 1)
        layout.addWidget(self.maxIlabel, 1, 0)
        layout.addWidget(self.maxI, 1, 1)
        layout.addWidget(self.rescaleButton, 0, 2, 2, 2)
        grid.addLayout(layout, 0, 0, 2, 2)
        
        grid.addWidget(self.VEL, 2, 0)
        grid.addWidget(self.VELslider, 2, 1, 1, 6)
        grid.addWidget(self.SPE, 3, 0)
        grid.addWidget(self.SPEslider, 3, 1, 1, 6)
        grid.addWidget(self.LAT, 4, 0)
        grid.addWidget(self.LATslider, 4, 1, 1, 6)
        
        layout = QHBoxLayout()
        layout.addWidget(self.channelButton)
        layout.addWidget(self.integrateButton)
        layout.addWidget(self.PVButton)
        # layout.addWidget(self.tempButton)
        grid.addLayout(layout, 5, 0, 1, 7)
        
        self.window.layout().addLayout(grid)
        
        self.updatePlot(self.file[1].data[0, :, :, 0], '{}'.format(self.species[0]))
        
        self.show_basic()
        
        return
  
    def printSize(self, value):
        print(value)
        print(self.window.sizeHint())
        print('height:', self.window.height())
        print('width:', self.window.width())
        return
  
    def plotDifference(self, **kwargs):
      
        if self.channelButton.isChecked():
            self.VELslider.setEnabled(True)
            self.SPEslider.setEnabled(True)
            self.LATslider.setEnabled(False)
            self.plotChannel(channel=True, **kwargs)
        elif self.integrateButton.isChecked():
            self.VELslider.setEnabled(False)
            self.SPEslider.setEnabled(True)
            self.LATslider.setEnabled(False)
            self.plotIntegrated(integrated=True, **kwargs)
        elif self.PVButton.isChecked():
            self.VELslider.setEnabled(False)
            self.SPEslider.setEnabled(True)
            self.LATslider.setEnabled(True)
            self.plotPV(PV=True, **kwargs)
        
        return
  
    def plotChannel(self, **kwargs):
        # Update the subplot value.
    
        spe = self.SPEslider.value()
        vel = self.VELslider.value()
        self.VEL.setText('Velocity: {:.2f} km/s'.format(self.velocity[vel]))
        
        if spe < self.species_length:
            self.SPE.setText('Species: {}'.format(self.species[spe]))
            self.updatePlot(self.file[1].data[vel, :, :, spe]
                            - self.file[1].data[:, :, :, spe].min(0), 
                            '{}'.format(self.species[spe]), **kwargs)
        else:
            dust = spe-self.species_length
            self.SPE.setText('Dust: {}'.format(self.dust[dust]))
            self.updatePlot(self.file[2].data[vel, :, :, dust], '{}'.format(self.dust[dust]), **kwargs)
        
        return
  
    def changeSpecies(self, **kwargs):
        # Update the subplot value.
    
        if self.PVButton.isChecked() is True:
            self.plotPV(True)
            return
        elif self.integrateButton.isChecked() is True:
            self.plotIntegrated()
            return
        else:
            self.VELslider.setEnabled(True)
            self.LATslider.setEnabled(False)
        
        if value is True:
            value = self.SPEslider.value()
        
        vel = self.VELslider.value()
        
        if value < self.species_length:
            self.SPE.setText('Species: {}'.format(self.species[value]))
            self.updatePlot(self.file[1].data[vel, :, :, value], '{}'.format(self.species[value]))
        else:
            dust = value-self.species_length
            self.SPE.setText('Dust: {}'.format(self.dust[dust]))
            self.updatePlot(self.file[2].data[:, :, dust], '{}'.format(self.dust[dust]))
        
        return
  
    def plotIntegrated(self, **kwargs):
        # Plot the corresponding integrated intensity.
        
        spe = self.SPEslider.value()
        
        if spe < self.species_length:
            self.SPE.setText('Species: {}'.format(self.species[spe]))
            self.updatePlot(self.species_intensity_integrated[spe, :, :], 
                            'integrated {}'.format(self.species[spe]),
                            **kwargs)
        else:
            dust = spe-self.species_length
            self.SPE.setText('Dust: {}'.format(self.dust[dust]))
            self.updatePlot(self.dust_intensity_integrated[:, :, dust], 
                            'integrated {}'.format(self.dust[dust]),
                            **kwargs)
        
        return
  
    def plotPV(self, **kwargs):
        # Plot the corresponding integrated intensity.
        
        spe = self.SPEslider.value()
        lat = self.LATslider.value()
        self.LAT.setText('Latitude: {:.1f} degrees'.format(self.latitude[lat]*180/np.pi))
        
        if spe < self.species_length:
            self.SPE.setText('Species: {}'.format(self.species[spe]))
            i_bg = np.where(self.file[1].data[:, lat, :, spe] == self.file[2].data[lat, :, 0])
            data = (self.file[1].data[:, lat, :, spe] 
                    - self.file[1].data[:, lat, :, spe].min(0))
                    # - np.asarray([self.file[2].data[lat, :, 0]]*self.velocity.size))
            data[i_bg] = np.nan
            self.updatePlot(data.T, 'PV {}'.format(self.species[spe]), **kwargs)
        else:
            dust = spe-self.species_length
            self.SPE.setText('Dust: {}'.format(self.dust[dust]))
            self.updatePlot(np.asarray([self.file[2].data[lat, :, dust]]*self.velocity.size).T, 
                            'PV {}'.format(self.dust[dust]), **kwargs)
        
        return
  
    def updatePlot(self, data, title, integrated=False, PV=False, minmax=False, **kwargs):
        # Remove the current Axes object and create a new one of the desired subplot.
        
        axTemp = self.fig.axes
        for ax in axTemp:
            self.fig.delaxes(ax)
        
        if self.rescaleButton.isChecked():
            vmin = self.minI.value()
            vmax = self.maxI.value()
        else:
            vmax = np.max(data)
            vmin = np.min(data)
            self.minI.setValue(vmin)
            self.maxI.setValue(vmax)
        
        if PV:
            data = copy(data)
            data[data == 0] = np.nan
            if not self.rescaleButton.isChecked():
                vmax = np.nanmax(data)
                vmin = np.nanmin(data)
                self.minI.setValue(vmin)
                self.maxI.setValue(vmax)
            ax = self.fig.add_axes((0.1, 0.1, 0.8, 0.9))
            cm = ax.pcolormesh(self.lonPV*180/np.pi, self.velPV, data, shading='auto',
                               norm=colors.SymLogNorm(linthresh=0.1, vmin=vmin, vmax=vmax, base=10), 
                               cmap=self.cmappv)
        else:
            ax = self.fig.add_axes((0.1, 0.1, 0.8, 0.9), projection='mollweide')
            cm = ax.pcolormesh(self.lonGrid, self.latGrid, data, shading='auto',
                               norm=colors.SymLogNorm(linthresh=0.1, vmin=vmin, vmax=vmax, base=10), 
                               cmap=self.cmap)
        # cm = ax.imshow(data, extent=(-np.pi, np.pi, -np.pi/2, np.pi/2),
        #                norm=colors.SymLogNorm(linthresh=0.1, vmin=vmin, vmax=vmax), cmap='cubehelix')
        cb = self.fig.colorbar(cm, ax=ax, fraction=0.02)
        if integrated:
            self.minIlabel.setText('minimum I km/s:')
            self.maxIlabel.setText('maximum I km/s:')
            cb.ax.set_ylabel(r'Intensity $log_{10} \left( K \frac{km}{s} \right)$',
                             fontsize=20, labelpad=16, rotation=270)
        else:
            self.minIlabel.setText('minimum I:')
            self.maxIlabel.setText('maximum I:')
            cb.ax.set_ylabel(r'Intensity $log_{10} \left( K \right)$',
                             fontsize=20, labelpad=16, rotation=270)
        ax.set_xlabel(r'Longitude $\left( rad \right)$', fontsize=16)
        if PV:
            ax.set_ylabel(r'Velocity $\left( \frac{km}{s} \right)$', fontsize=16)
            ax.invert_xaxis()
            # ax.set_aspect('equal')
            self.fig.tight_layout()
        else:
            ax.set_ylabel(r'Latitude $\left( rad \right)$', fontsize=16)
        ax.set_title('Galactic '+title, fontsize=24)
        
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
