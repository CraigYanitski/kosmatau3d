# KOSMA-tau^3

This is a development branch for KOSMA-tau^3 where I have converted some of the classes from the *develop* branch into modules. The reason for doing this is to reduce the memory footprint of the code to increase the overall efficiency. This also increases the maximum number of voxels that can be evaluated, since each voxel no longer owns as much memory.

## Functionality

It is mostly functional, allowing data to be streamed into fits files and separating the radiative transfer portion to a self-contained module. What still remains is to correctly compute the observed velocity of the clump and properly map the integrated intensity using healpy. It has not been converted to the master branch since I am still debugging and ensuring it can reproduce the same result as the old code. The most recent version of KOSMA-tau 3D is the galaxy model developed by Christof Bruckmann. When this branch can create the same model, I will convert it to the master branch. Comparison of the two branches is rather difficult since I am also debugging the old code. The major changes to the KOSMA-tau model are described in the document treatise.pdf, and the major changes to the Milky Way model will also appear in the upcoming Bruckmann et al. (2020) paper.

Among other things, I believe I have successfully taken out the incredibly time-consuming use of the scipy function griddata() for each species for each masspoint in each voxel. I have replaced it with an initial call to either the scipy function LinearNDInterpolation() or Rbf() to create a list of functions that will interpolate the grid of data. This has greatly reduced the execution time of large, high-resolution models. Besides this, I have correctly implemented matrix calculations using numpy arrays.
