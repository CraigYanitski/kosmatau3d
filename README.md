<!-- PROJECT SHIELDS -->
<!--
*** I'm using markdown "reference style" links for readability.
*** Reference links are enclosed in brackets [ ] instead of parentheses ( ).
*** See the bottom of this document for the declaration of the reference variables
*** for contributors-url, forks-url, etc. This is an optional, concise syntax you may use.
*** https://www.markdownguide.org/basic-syntax/#reference-style-links
-->
[![Forks][forks-shield]][forks-url]
[![Stargazers][stars-shield]][stars-url]
<!-- [![License: GPL v3][(https://img.shields.io/badge/License-GPLv3-blue.svg)]][license-url] -->
[![pylint][pylint-shield]][pylint-url]
[![Contributors][contributors-shield]][contributors-url]
[![Issues][issues-shield]][issues-url]
[![GNU GPLv3 License][license-shield]][license-url]


<!-- [![pipeline status](https://git.ph1.uni-koeln.de/yanitski/kosma-tau-3d/badges/master/pipeline.svg)](https://git.ph1.uni-koeln.de/yanitski/kosma-tau-3d/-/commits/master)
[![coverage report](https://git.ph1.uni-koeln.de/yanitski/kosma-tau-3d/badges/master/coverage.svg)](https://git.ph1.uni-koeln.de/yanitski/kosma-tau-3d/-/commits/master)
[![Latest Release](https://git.ph1.uni-koeln.de/yanitski/kosma-tau-3d/-/badges/release.svg)](https://git.ph1.uni-koeln.de/yanitski/kosma-tau-3d/-/releases) -->

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


<!-- MARKDOWN LINKS & IMAGES -->
<!-- https://www.markdownguide.org/basic-syntax/#reference-style-links -->
[contributors-shield]: https://img.shields.io/github/contributors/CraigYanitski/kosmatau3d.svg
[contributors-url]: https://github.com/CraigYanitski/kosmatau3d/graphs/contributors
[forks-shield]: https://img.shields.io/github/forks/CraigYanitski/kosmatau3d.svg
[forks-url]: https://github.com/CraigYanitski/kosmatau3d/network/members
[stars-shield]: https://img.shields.io/github/stars/CraigYanitski/kosmatau3d.svg
[stars-url]: https://github.com/CraigYanitski/kosmatau3d/stargazers
[issues-shield]: https://img.shields.io/github/issues/CraigYanitski/kosmatau3d.svg
[issues-url]: https://github.com/CraigYanitski/kosmatau3d/issues
[license-shield]: https://img.shields.io/badge/License-GPLv3-blue.svg
<!-- https://img.shields.io/github/license/CraigYanitski/kosmatau3d.svg?style=for-the-badge -->
[license-url]: https://www.gnu.org/licenses/gpl-3.0
<!-- https://github.com/CraigYanitski/kosmatau3d/blob/master/LICENSE.txt -->
[pylint-shield]: https://github.com/CraigYanitski/kosmatau3d/actions/workflows/test-lint.yml/badge.svg
[pylint-url]: https://github.com/CraigYanitski/kosmatau3d/actions/workflows/test-lint.yml
[product-screenshot]: images/screenshot.png
[Next.js]: https://img.shields.io/badge/next.js-000000?style=for-the-badge&logo=nextdotjs&logoColor=white
[Next-url]: https://nextjs.org/
[React.js]: https://img.shields.io/badge/React-20232A?style=for-the-badge&logo=react&logoColor=61DAFB
[React-url]: https://reactjs.org/
[Vue.js]: https://img.shields.io/badge/Vue.js-35495E?style=for-the-badge&logo=vuedotjs&logoColor=4FC08D
[Vue-url]: https://vuejs.org/
[Angular.io]: https://img.shields.io/badge/Angular-DD0031?style=for-the-badge&logo=angular&logoColor=white
[Angular-url]: https://angular.io/
[Svelte.dev]: https://img.shields.io/badge/Svelte-4A4A55?style=for-the-badge&logo=svelte&logoColor=FF3E00
[Svelte-url]: https://svelte.dev/
[Laravel.com]: https://img.shields.io/badge/Laravel-FF2D20?style=for-the-badge&logo=laravel&logoColor=white
[Laravel-url]: https://laravel.com
[Bootstrap.com]: https://img.shields.io/badge/Bootstrap-563D7C?style=for-the-badge&logo=bootstrap&logoColor=white
[Bootstrap-url]: https://getbootstrap.com
[JQuery.com]: https://img.shields.io/badge/jQuery-0769AD?style=for-the-badge&logo=jquery&logoColor=white
[JQuery-url]: https://jquery.com 