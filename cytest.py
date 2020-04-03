#!/usr/bin/env python
# coding: utf-8

# # Test of using cygrid with the data
# 
# The capability of cygrid when interfacing with the integrated intensity maps is tested here. If it goes well, a self-contained module will be soon created.

# In[164]:


import cygrid
from astropy.io import fits
import numpy as np


# In[166]:


hdul = fits.open('history/MilkyWay/r1000_n3015/integrated_intensity.fits')


# In[167]:


for hdu in hdul:
    print()
    print(hdu.header['VELOCITY'])
    print(hdu.header['DIREC'])
    print(hdu.data.shape)


# In[168]:


from cygrid.tests import test_cygrid as tc
test1 = tc.TestWcsGrid()
test1.setup()
test2 = tc.TestSlGrid()
test2.setup()


# In[169]:


test1.test_gridding()


# In[170]:


test2.test_gridding()


# In[171]:


from cygrid.tests import test_healpix as th
test = th.TestHealpix()
test.setup()


# In[172]:


test.test_pix2ring()


# In[173]:


header = {
   'NAXIS': 3,
   'NAXIS1': 101,
   'NAXIS2': 101,
   'NAXIS3': 1024,
   'CTYPE1': 'GLON-SFL',
   'CTYPE2': 'GLAT-SFL',
   'CDELT1': -0.1,
   'CDELT2': 0.1,
   'CRPIX1': 51,
   'CRPIX2': 51,
   'CRVAL1': 12.345,
   'CRVAL2': 3.14,
   }

kernelsize_sigma = 0.02
kernel_type = 'gauss1d'
kernel_params = (kernelsize_sigma,)
kernel_support = 3*kernelsize_sigma
hpx_maxres = kernelsize_sigma/2

mygridder = cygrid.WcsGrid(header)
mygridder.set_kernel(
  kernel_type,
  kernel_params,
  kernel_support,
  hpx_maxres
  )

glon = np.random.uniform(170, 180, 101)
glat = np.random.uniform(80, 90, 101)
# glon,glat = np.meshgrid(glon,glat)
sig = np.random.rand(101,1024)

mygridder.grid(glon.flatten(),glat.flatten(),sig)


# In[174]:


data = mygridder.get_datacube()

