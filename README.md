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

For a detailed explanation of the properties that can be accessed with a voxel instance, see the `jupyter` notebook in `./tests/voxel.ipynb`.

## Functionality

#### Single-voxel models

This is the basic component of `kosmatau3d`. It is made available as a self-sufficient object for use in other subgridding models. Given a volume-filling factor, mass, and FUV field, the single voxel object calculates its emissivity, absorption, and intensity (assuming no background intensity) for a clumpy PDR ensemble. The explanation of this model is thoroughly-explained in `./tests/voxel.ipynb`.

#### `kosmatau3d` model

The full subgrid model is almost functional. Voxel data can be streamed into fits files and the radiative transfer portion is a self-contained module. What still remains is a few fixes the the radiative transfer module and to correctly map the integrated intensity using ``cygrid``. The current version fixes most of the errors inherent in the old code. It is setup for the galaxy model developed by Christoph Bruckmann as that will be its first application.

The objects that will modelled with this are:

  - Milky Way: an approximate description compared to COBE-FIRAS, Planck, CfA, Mopra, ThrUMMS, and SEDIGISM data
  - Another galaxy...

## Code Corrections

The major changes to the `kosmatau3d` model are described in the document treatise.pdf, and the major changes to the Milky Way model will also appear in the upcoming Yanitski et al. (2022) paper. There will be other documents to state the other various corrections and developments made.

## Other Features

It was necessary to create a version to initialise a single voxel. This is now the default version of initialising a `Voxel` instance, and it requires passing all of the needed arguments as kwargs. This feature can hopefully be used in other models that are unable to see the effect of the inhomogeneity of PDRs.

## Ongoing Development

* Fix observational regridding error maps
* Implement FUV absorption calculation in full 3D model
* Modify the `radiativeTransfer` module to work for arbitrary geometry
* Implement `numba` more fully to optimise the computation time
  * use this to parallelise the code
* Create a GUI to make it easier to setup/alter the model
