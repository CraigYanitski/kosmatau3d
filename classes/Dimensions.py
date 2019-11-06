import numpy as np
class Dimensions(object):
  '''
  This is a class to contain the dimensions of VoxelGrid(). It will
  coordinate with class VoxelGrid to arrange the voxels. the input dimensions must be
  given in kpc.

  At the moment this will just work in the plane of a disk.
  '''
  # PRIVATE
  def __init__(self, x, y, z, i=0, resolution=1000):
    self.__scale = resolution
    self.__i = i
    self.__x = x*1000
    self.__y = y*1000
    self.__z = z*1000
    self.__xOffset = self.__x/2.
    self.__yOffset = self.__y/2.
    self.__zOffset = self.__z/2.
    self.__xRange = np.arange(0, self.__x+self.__scale, self.__scale)
    self.__yRange = np.arange(0, self.__y+self.__scale, self.__scale)
    self.__zRange = np.arange(0, self.__z+self.__scale, self.__scale)
    grid = np.meshgrid(self.__xRange, self.__yRange, self.__zRange)
    self.__xPositions = grid[0].flatten()
    self.__yPositions = grid[1].flatten()
    self.__zPositions = grid[2].flatten()
    # for ix in range(len(self.__x)):
    #   for iy in range(len(self.__y)):
    #     for iz in range(len(self.__z)):
    #       # Center the coordinates about the origin
    #       self.__xPosition[ix] = self.__x[ix]*self.__scale  - (self.__x[-1]-1)/2*self.__scale
    #       self.__yPosition[iy] = self.__y[iy]*self.__scale  - (self.__y[-1]-1)/2*self.__scale
    #       self.__zPosition[iz] = self.__z[iz]*self.__scale  - (self.__z[-1]-1)/2*self.__scale
    '''
    if cut_switch == 'inner':
      if self.calc_r (x_pos,y_pos,z_pos,inc_disk) < r_disk and abs(calc_h(x_pos,y_pos,z_pos,inc_disk)) < z_disk and (R0**2-x_pos**2) >= 0.  and z_pos <= (R0**2-x_pos**2)**0.5: #cutdisk inner
        coords.append([x, y, z])
        abs_coords.append([x_pos,y_pos,z_pos])
    if cut_switch == 'outer':
      if self.calc_r (x_pos,y_pos,z_pos,inc_disk) < r_disk and abs(calc_h(x_pos,y_pos,z_pos,inc_disk)) < z_disk and (R0**2-x_pos**2) >= 0. and -z_pos >= (R0**2-x_pos**2)**0.5: #cut disk outer down correct-for radiative transfer  
        coords.append([x, y, z])
        abs_coords.append([x_pos,y_pos,z_pos])
    if cut_switch == 'total':
      if self.calc_r (x_pos,y_pos,z_pos,inc_disk) < r_disk and abs(calc_h(x_pos,y_pos,z_pos,inc_disk)) < z_disk: #and z_pos < R0: #total
        coords.append([x, y, z])
        abs_coords.append([x_pos,y_pos,z_pos])
    '''
    #x = abs_coords[int(i)][0]
    #y = abs_coords[int(i)][1]
    #z = abs_coords[int(i)][2]
    self.__r, self.__phi = self.__convertToPolar()
    #self.__h = self.hCalc(self.__x,self.__y,self.__z,self.__i)
    return
  def __convertToPolar(self):
    # This is a function that calculates the polar coordinates of each voxel
    r = ((self.__xPositions-self.__xOffset)**2 + (self.__yPositions-self.__yOffset)**2)**0.5
    phi = np.arcsin(self.__yPositions/self.__xPositions)
    return r, phi
  def __str__(self):
    return 'The dimensions of the model are {}x{}x{}'.format(self.__x[-1], self.__y[-1], self.__z[-1])

  # PUBLIC
  def voxelNumber(self):
    # Return the number of voxels required
    return len(self.__xPositions)
  def voxelCartesianPosition(self):
    # Return the Cartesian coordinates of each voxel
    return (self.__xPositions, self.__yPositions, self.__zPositions, self.__scale)
  def voxelPolarPosition(self):
    # Return the polar coordinates of each voxel
    return (self.__r, self.__phi)
  def hCalc(self):
    #calculates height h from koordinates relative to disk-plain

    # This will soon be rewritten.
    h = (self.__yPosition - np.tan(self.__i)*self.__zPosition) * np.cos(self.__i)
    return h
