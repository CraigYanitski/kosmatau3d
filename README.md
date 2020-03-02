# KOSMA-tau^3

This is a development branch for KOSMA-tau^3 where I have converted some of the classes from the *develop* branch into modules. The reason for doing this is to reduce the memory footprint of the code to increase the overall efficiency. This also increases the maximum number of voxels that can be evaluated, since each voxel no longer owns as much memory.

## Functionality

It is mostly functional, allowing data to be streamed into fits files and separating the radiative transfer portion to a self-contained module. What still remains is to correctly map the integrated intensity using ``healpy``. The master branch currently is functional, and fixes the errors inherent in the old code. The most recent version of KOSMA-tau^3 is setup for the galaxy model developed by Christof Bruckmann. It assumes the galactic disk lies in the x,y plane, and the radiative transfer is calculated from the perspective of Earth.

## Code Corrections

The major changes to the KOSMA-tau model are described in the document treatise.pdf, and the major changes to the Milky Way model will also appear in the upcoming Bruckmann et al. (2020) paper.

## Branches

There are a few branches used for the development of this code. The *master* branch contains the current working version that is continuously developed and improved. The old code from Christoph Bruckmann and Silke Andree-Labsch (updated to ``python 3`` and slightly optimised) is in the *old-version* branch. An intermediate version of the code is stored in the *class-based* branch. It contains the first version of the current code that used classes rather than the modules present in the *master* branch. It ran into memory issues due to the large number of unnecessary objects. Finally the *develop* branch is used to test new features that will be merged to the master branch.

## Ongoing Development

* Create a module to utilise ``healpy`` to analyse the integrated intensity fits file
  * integrate the module into the radiative transfer calculation
* Fix the probability calculations to use ``scipy``
* Implement ``numba`` more fully to optimise the computation time
  * use this to parallelise the code
* Create a GUI to make it easier to alter the model
