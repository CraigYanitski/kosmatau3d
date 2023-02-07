[![pipeline status](https://git.ph1.uni-koeln.de/yanitski/kosma-tau-3d/badges/master/pipeline.svg)](https://git.ph1.uni-koeln.de/yanitski/kosma-tau-3d/-/commits/master)
[![coverage report](https://git.ph1.uni-koeln.de/yanitski/kosma-tau-3d/badges/master/coverage.svg)](https://git.ph1.uni-koeln.de/yanitski/kosma-tau-3d/-/commits/master)
[![Latest Release](https://git.ph1.uni-koeln.de/yanitski/kosma-tau-3d/-/badges/release.svg)](https://git.ph1.uni-koeln.de/yanitski/kosma-tau-3d/-/releases)

# **kosmatau3d**

This is the current working version of `kosmatau3d`.
It uses a series of sub-modules to perform most of the calculations.
The reason for doing this is to reduce the memory footprint of the code to
increase the overall efficiency.
This also increases the maximum number of voxels that can be evaluated, since
each voxel no longer owns as much memory.

## Installation

### using `pip`

Directly install from git using,

```
pip install git+https://github.com/CraigYanitski/kosmatau3d.git
```

### manually

Download the repository.

```bash
git clone https://github.com/CraigYanitski/kosmatau3d
cd kosmatau3d
```

Now that you are in the root directory of this repository, install this package in bash with,

```bash
pip install .
```

## Creating a voxel

A `Voxel` instance can be initialised using,

```python
>>> from kosmatau3d import models
>>> vox = models.Voxel()
```

There are many parameters that must be specified in order to initialise and
simulate the clumpy ensemble.
For a detailed explanation of the properties that can be defined/accessed with
a voxel instance, see the `jupyter` notebook in `./notebooks/voxel.ipynb`.
If you wish to use the pre-defined properties, you can simply run,

```python
>>> vox.set_properties()
>>> vox.calculate_emission()
```

One can then easily plot different views of the voxel emission using
built-in plotting methods.

```python
>>> vox.plot_spectrum()
```

## Documentation (in prep.)

At the moment there is no CI/CD to automatically compile the documentation files.
In order to properly view the current documentation, from the root directory run,

```bash
cd doc
make html
cd build/html
browse index.html
```

You should now see the homepage of the documentation in its current form.
This will periodically be improved over time.

## Functionality

### Single-voxel models

This is the basic component of `kosmatau3d`.
It is made available as a self-sufficient object for use in other subgridding models.
Given a volume-filling factor, mass, and FUV field, the single voxel object
calculates the wavelength-dependant intensity, optical depth, absorption, and
emissivity (assuming no background intensity) for a clumpy PDR ensemble.
The explanation of this model is thoroughly-explained in `./notebooks/voxel.ipynb`.

The objects that will modelled with this method are:

- IC 1396
  - (Okada et al. in prep) first application of directly comparing single
  voxels with an observational map

### 3D models

The full subgrid model to simulate entire 3-dimensional structures.
Voxel data can be streamed into fits files and the radiative transfer portion
is a self-contained module to save on computational time.

It is currently setup for the Milky Way model initially developed by Christoph
Bruckmann as that will be its first application.
This galactic model can also be used in a more generalised application for
distant galaxies.

The objects that will modelled with this method are:

- Milky Way
  - (Yanitski et al. in prep) an approximate description compared to
  COBE-DIRBE, Planck, COBE-FIRAS, GOT C+, CfA, Mopra, ThrUMMS, and SEDIGISM data

## Code Corrections

The major changes to the functionality of the `kosmatau3d` model over the KOSMA-tau-3D 
model of Silke et al. (2017) are described in the document `./docs/treatise.pdf`, 
and the major changes to the Milky Way model will also
appear in the upcoming Yanitski et al. (2023) paper.
There will be other documents to state the other various corrections and
developments made.

## Developmental Progress

### base development

- 3D model
  - [x] Correct voxel-averaged intensity calculation
  - [x] Ensure 3D model saves relevant data
    - [ ] implement multiprocessing
  - [x] Modify single-voxel instance to calculate a column when $f_V > 1$
    - [ ] allow for arbitrary $f_V$ $(<1)$
- Radiative transfer
  - [x] Ensure correct calculation of sythetic datacube
    - [x] implement multiprocessing
    - [x] optimise
  - [x] Implement opacity map in radiative transfer calculation
  - [ ] Implement FUV absorption calculation in full 3D model
  - [ ] Modify the `radiativeTransfer` module to work for arbitrary geometry
- Miscellaneous
  - [x] Include submodule to explore KOSMA-$\tau$ and `kosmatau3d` properties
  - [x] Include submodule to compare observational data with synthetic datacubes
    - [x] fix regridding of observational error maps
  - [ ] Implement function to mimic two-layer model (*aka* **two-voxel** model)
  - [ ] Implement function to create voxel grid (functions by Yoko Okada)
    - [ ] have function fit observational data
  - [ ] Include submodule covering the Mathematica routines developed by Markus Röllig
  - [ ] Develop the code testing suite
  - [ ] Fully document the code using `sphinx`
  - [ ] Implement CI/CD
  
### future development

- [x] Allow pickling of interpolation functions for faster debugging of the single-voxel model
- [ ] Utilise `cython` to improve code efficiency
- [ ] Implement `numba` more fully to optimise the computation time
  - [ ] use this to parallelise the code
- [ ] Create a GUI to make it easier to setup/alter the model
- [ ] Implement recursion in radiative transfer calculation
