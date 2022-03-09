# **kosmatau3d**

This is the current working version of `kosmatau3d`. It uses a series of sub-modules to perform most of the calculations. The reason for doing this is to reduce the memory footprint of the code to increase the overall efficiency. This also increases the maximum number of voxels that can be evaluated, since each voxel no longer owns as much memory.

## Installation
In the directory of this repository, install this package in bash with,

```
pip install .
```

A `Voxel` instance can be initialised using,

```
from kosmatau3d import models
vox = models.Voxel()
```

There are many parameters that must be specified in order to initialise and simulate the clumpy ensemble. For a detailed explanation of the properties that can be defined/accessed with a voxel instance, see the `jupyter` notebook in `./tests/voxel.ipynb`.

## Functionality

#### Single-voxel models

This is the basic component of `kosmatau3d`. It is made available as a self-sufficient object for use in other subgridding models. Given a volume-filling factor, mass, and FUV field, the single voxel object calculates the wavelength-dependant intensity, optical depth, absorption, and emissivity (assuming no background intensity) for a clumpy PDR ensemble. The explanation of this model is thoroughly-explained in `./tests/voxel.ipynb`.

#### 3D models

The full subgrid model to simulate entire 3-dimensional structures. Voxel data can be streamed into fits files and the radiative transfer portion is a self-contained module to save on computational time.

It is currently setup for the Milky Way model initially developed by Christoph Bruckmann as that will be its first application. This galactic model can also be used in a more generalised application for distant galaxies.

The objects that will modelled with this are:

  - Milky Way
    - an approximate description compared to COBE-FIRAS, Planck, CfA, Mopra, ThrUMMS, and SEDIGISM data
  - IC 1396
    - first application of directly comparing single voxels with an observational map

## Code Corrections

The major changes to the `kosmatau3d` model are described in the document treatise.pdf, and the major changes to the Milky Way model will also appear in the upcoming Yanitski et al. (2022) paper. There will be other documents to state the other various corrections and developments made.

## Developmental Progress

* [x] Correct voxel-averaged intensity calculation
* [x] Ensure 3D model saves relevant data
  * [ ] implement multiprocessing
* [x] Ensure correct calculation of sythetic datacube
  * [x] implement multiprocessing
  * [ ] optimise
* [x] Modify single-voxel instance to calculate a column when $f_V > 1$
  * [ ] Allow for arbitrary $f_V$ ($<1$)
* [ ] Fixed observational regridding error maps
* [ ] Implement FUV absorption calculation in full 3D model
* [ ] Modify the `radiativeTransfer` module to work for arbitrary geometry
* [ ] Implement `numba` more fully to optimise the computation time
  * [ ] use this to parallelise the code
* [ ] Create a GUI to make it easier to setup/alter the model
