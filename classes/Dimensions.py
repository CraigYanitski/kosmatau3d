class Dimensions():
  '''
  This is a class to contain the dimensions of super-class Shape. It will
  coordinate with class VoxelGrid to arrange the voxels. the input dimensions must be
  given in kpc.
  '''
  # PRIVATE
  def __init__(self, x, y, z, i=0, resolution=1000):
    self.__scale = resolution
    self.__i = i
    self.__x = np.arange(x)
    self.__y = np.arange(y)
    self.__z = np.arange(z)
    self.__xOffset = self.__x/2.
    self.__yOffset = self.__y/2.
    self.__zOffset = self.__z/2.
    self.__xPosition = np.zeros(len(self.__x))
    self.__yPosition = np.zeros(len(self.__y))
    self.__zPosition = np.zeros(len(self.__z))
    for ix in range(len(self.__x)):
      for iy in range(len(self.__y)):
        for iz in range(len(self.__z)):
          self.__xPosition[ix] = self.__x[ix]*self.__scale  - (self.__x-1)/2*self.__scale #real coordiantes
          self.__yPosition[iy] = self.__y[iy]*self.__scale  - (self.__y-1)/2*self.__scale
          self.__zPosition[iz] = self.__z[iz]*self.__scale  - (self.__z-1)/2*self.__scale
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
    self.__r = calc_r (self.__x,self.__y,self.__z,self.__i)
    self.__h = calc_h (self.__x,self.__y,self.__z,self.__i)
    self.__phi = np.arctan2 (x,(y**2+z**2)**0.5)
    return

  # PUBLIC
  def indeces(self):
    return len(self.__x)+len(self.__y)+len(self.__z)
  def calc_h(x,y,z,inc) : #calculates height h from koordinates relative to disk-plain
    h = (y - np.tan(inc)*z) * np.cos(inc)
    return h
  def calc_r (x,y,z,inc) : #calculates the radius in the diskplain 
    h = (y - np.tan(inc)*z) * np.cos(inc)
    Px = x ; Py = y - np.cos(inc)*h ; Pz = z + np.tan(inc)*np.cos(inc)*h
    r = (Px**2 + Py**2 + Pz**2)**0.5
    return r
