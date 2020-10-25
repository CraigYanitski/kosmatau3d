# $KOSMA-tau^3$

This is the current working version of `KOSMA-tau^3` where I have converted some of the classes from the ***class-based*** branch into modules. The reason for doing this is to reduce the memory footprint of the code to increase the overall efficiency. This also increases the maximum number of voxels that can be evaluated, since each voxel no longer owns as much memory.

## Installation
In the directory containing `KOSMA-tau^3`, install this package in bash with,

```
python -m pip install ./KOSMA-tau^3
```

**NOTE:** pip 20.2 has certain depreciations that make this install very slow. This is currently being investigated, so use pip 20.1 for now.

A voxel instance can be initialised using,

```
from kosmatau3d import models
vox = models.Voxel()
```

For a detailed explanation of the properties that can be accessed with a voxel instance, see the test notebook in `KOSMA-tau^3/tests/voxel.ipynb`.

## Functionality

#### Single-voxel models

This is the basic component of KOSMA-tau^3. It is made available as a self-sufficient object for use in other subgridding models. Given a mass, density, and FUV field, the single voxel object calculates its emissivity, absorption, and intensity (assuming no background intensity).

#### KOSMA-tau^3 model

It is almost functional. Voxel data can be streamed into fits files and the radiative transfer portion is a self-contained module. What still remains is a few fixes the the radiative transfer module and to correctly map the integrated intensity using ``cygrid``. The current version fixes most of the errors inherent in the old code. It is setup for the galaxy model developed by Christof Bruckmann as that will be its first application.

## Code Corrections

The major changes to the KOSMA-tau 3D model are described in the document treatise.pdf, and the major changes to the Milky Way model will also appear in the upcoming Bruckmann et al. (2020) paper. There will be other documents to state the other various corrections and developments made.

## Branches

There are a few branches used for the development of this code. The ***master*** branch contains the current working version that is continuously developed and improved. The old code from Christoph Bruckmann and Silke Andree-Labsch (updated to ``python 3`` and slightly optimised) is in the ***old-version*** branch. An intermediate version of the code is stored in the ***class-based*** branch. It contains the first version of the current code that used classes rather than the modules present in the ***master*** branch. It ran into memory issues due to the large number of unnecessary objects. Finally the ***develop*** branch is used to test new features that will be merged to the master branch.

## Other Features

It was necessary to create a version to initialise a single voxel. This is now the default version of initialising a `Voxel` instance, and it requires passing all of the needed arguments as kwargs. This feature can hopefully be used in the subgrid model of Daniel Seifried.

## Ongoing Development

* Correct the ``radiativeTransfer`` module
* Create a module to utilise ``cygrid`` to analyse the integrated intensity fits file
  * integrate the module into the radiative transfer calculation
* Implement ``numba`` more fully to optimise the computation time
  * use this to parallelise the code
* Create a GUI to make it easier to alter the model
